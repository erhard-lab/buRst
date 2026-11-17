# run_kinetics_batch.R
# Batch fitting for one or two conditions.

# Internal utilities -----------------------------------------------------------

.require_batch_namespaces <- function(require_h5 = FALSE, use_timeout = FALSE) {
  packages <- c("grandR", "parallel")
  if (isTRUE(require_h5)) packages <- c(packages, "rhdf5")
  if (isTRUE(use_timeout)) packages <- c(packages, "R.utils")

  missing <- packages[!vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )]

  if (length(missing) > 0L) {
    stop(
      "Missing required package(s): ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
}

.is_nonempty_scalar <- function(x) {
  length(x) == 1L && !is.null(x) && !is.na(x) && nzchar(as.character(x))
}

.has_second_condition <- function(cond1, cond2) {
  .is_nonempty_scalar(cond2) && !identical(as.character(cond2), as.character(cond1))
}

.rbind_fill <- function(dfs) {
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  if (length(dfs) == 0L) return(data.frame())

  all_columns <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs <- lapply(dfs, function(x) {
    missing_columns <- setdiff(all_columns, names(x))
    for (column in missing_columns) x[[column]] <- NA
    x[all_columns]
  })

  rownames_result <- NULL
  out <- do.call(rbind, dfs)
  rownames(out) <- rownames_result
  out
}

.function_formals <- function(fun) {
  tryCatch(names(formals(fun)), error = function(e) character())
}

.accepts_argument <- function(fun, argument) {
  fml <- .function_formals(fun)
  argument %in% fml || "..." %in% fml
}

.has_explicit_argument <- function(fun, argument) {
  argument %in% .function_formals(fun)
}

.filter_call_arguments <- function(fun, args) {
  fml <- .function_formals(fun)
  if ("..." %in% fml) return(args)
  args[intersect(names(args), fml)]
}

.resolve_reference_condition <- function(steady, ref_cond, other_cond = NULL) {
  if (isTRUE(steady)) return(NULL)
  if (.is_nonempty_scalar(ref_cond)) return(as.character(ref_cond))
  if (.is_nonempty_scalar(other_cond)) return(as.character(other_cond))
  NULL
}

.h5_gene_ids <- function(h5_path) {
  if (!.is_nonempty_scalar(h5_path)) {
    stop("h5_path is required when require_h5 = TRUE.", call. = FALSE)
  }

  listing <- rhdf5::h5ls(h5_path, recursive = TRUE)
  ids <- listing$name[listing$group == "/genes"]
  unique(as.character(ids))
}

.prepare_batch_gene_map <- function(
    data,
    genes,
    h5_path,
    require_h5,
    keep_only_unique_symbols
) {
  if (length(genes) == 0L) {
    stop("genes must contain at least one gene symbol or ENSG identifier.", call. = FALSE)
  }

  genes <- as.character(genes)
  gi <- grandR::GeneInfo(data)
  needed <- c("Gene", "Symbol")
  if (!all(needed %in% names(gi))) {
    stop("GeneInfo(data) must contain 'Gene' and 'Symbol' columns.", call. = FALSE)
  }

  gi_gene <- as.character(gi$Gene)
  gi_symbol <- as.character(gi$Symbol)

  map_one <- function(input) {
    is_ensg <- grepl("^ENSG\\d+", input)

    if (is_ensg) {
      idx <- which(gi_gene == input)
      ensg <- if (length(idx) > 0L) input else NA_character_
      symbol <- if (length(idx) > 0L) gi_symbol[idx[1L]] else NA_character_
      ambiguous <- FALSE
    } else {
      idx <- which(gi_symbol == input)
      matched_ensg <- unique(gi_gene[idx])
      matched_ensg <- matched_ensg[!is.na(matched_ensg) & nzchar(matched_ensg)]
      ambiguous <- length(matched_ensg) > 1L
      ensg <- if (length(matched_ensg) > 0L) matched_ensg[1L] else NA_character_
      symbol <- input
    }

    data.frame(
      input = input,
      input_is_ensg = is_ensg,
      Symbol = symbol,
      ENSG = ensg,
      ambiguous_symbol = ambiguous,
      ok_map = !is.na(ensg) && nzchar(ensg),
      stringsAsFactors = FALSE
    )
  }

  gene_map <- do.call(rbind, lapply(genes, map_one))
  rownames(gene_map) <- NULL

  if (isTRUE(require_h5)) {
    h5_ids <- .h5_gene_ids(h5_path)
    gene_map$ok_h5 <- gene_map$ENSG %in% h5_ids
  } else {
    gene_map$ok_h5 <- NA
  }

  gene_map$ok <- gene_map$ok_map
  if (isTRUE(require_h5)) gene_map$ok <- gene_map$ok & gene_map$ok_h5
  if (isTRUE(keep_only_unique_symbols)) {
    gene_map$ok <- gene_map$ok & !gene_map$ambiguous_symbol
  }

  gene_map$note <- NA_character_
  gene_map$note[!gene_map$ok_map] <- "unmapped_gene"
  if (isTRUE(require_h5)) {
    gene_map$note[gene_map$ok_map & !gene_map$ok_h5] <- "not_in_h5"
  }
  if (isTRUE(keep_only_unique_symbols)) {
    gene_map$note[gene_map$ambiguous_symbol] <- "dropped_ambiguous_symbol"
  } else {
    gene_map$note[gene_map$ambiguous_symbol & gene_map$ok] <-
      "ambiguous_symbol_using_first_ENSG"
  }

  gene_map$gene_run <- ifelse(
    gene_map$ok,
    ifelse(
      gene_map$input_is_ensg,
      gene_map$ENSG,
      ifelse(!is.na(gene_map$Symbol) & nzchar(gene_map$Symbol),
             gene_map$Symbol,
             gene_map$ENSG)
    ),
    NA_character_
  )

  genes_ok <- unique(stats::na.omit(gene_map$gene_run))
  if (length(genes_ok) == 0L) {
    stop("No valid genes remain after mapping and HDF5 validation.", call. = FALSE)
  }

  list(gene_map = gene_map, genes_ok = genes_ok)
}

.prepare_fit_specs <- function(
    cond1,
    cond2,
    fit1_steady,
    fit2_steady,
    fit1_ref_cond,
    fit2_ref_cond,
    truth_suffix1,
    truth_suffix2
) {
  if (!.is_nonempty_scalar(cond1)) {
    stop("cond1 must be a non-empty condition name.", call. = FALSE)
  }

  two_conditions <- .has_second_condition(cond1, cond2)

  spec1 <- list(
    condition = as.character(cond1),
    steady = isTRUE(fit1_steady),
    ref_cond = .resolve_reference_condition(
      steady = fit1_steady,
      ref_cond = fit1_ref_cond,
      other_cond = if (two_conditions) cond2 else NULL
    ),
    truth_suffix = truth_suffix1,
    fit_number = 1L
  )

  specs <- list(spec1)

  if (two_conditions) {
    specs[[2L]] <- list(
      condition = as.character(cond2),
      steady = isTRUE(fit2_steady),
      ref_cond = .resolve_reference_condition(
        steady = fit2_steady,
        ref_cond = fit2_ref_cond,
        other_cond = cond1
      ),
      truth_suffix = truth_suffix2,
      fit_number = 2L
    )
  }

  for (spec in specs) {
    if (!spec$steady && is.null(spec$ref_cond)) {
      stop(
        "A non-steady-state fit for condition '", spec$condition,
        "' requires a reference condition.",
        call. = FALSE
      )
    }
  }

  list(specs = specs, two_conditions = two_conditions)
}

.prepare_shared_reference <- function(fit_fun, h5_path, reference) {
  if (!is.null(reference)) return(reference)

  # Automatic loading is intentionally restricted to an explicit formal.
  # A custom wrapper using only `...` can still receive a user-supplied object.
  if (!.has_explicit_argument(fit_fun, "reference")) return(NULL)

  if (!.is_nonempty_scalar(h5_path)) {
    stop(
      "fit_fun accepts a reference but h5_path is missing, so the reference ",
      "cannot be loaded automatically.",
      call. = FALSE
    )
  }

  load_10x_reference(dirname(h5_path))
}

.call_fit_fun <- function(
    fit_fun,
    data,
    gene,
    spec,
    t.label,
    t.diff,
    h5_path,
    reference,
    ndraw,
    seed
) {
  args <- list(
    data = data,
    gene = gene,
    fit.cond = spec$condition,
    t.label = t.label,
    steady = spec$steady,
    ref.cond = spec$ref_cond,
    t.diff = if (spec$steady) NULL else t.diff,
    h5_path = h5_path,
    ndraw = ndraw,
    seed = seed,
    truth_suffix = spec$truth_suffix,
    add_pit = TRUE
  )

  if (!is.null(reference) && .accepts_argument(fit_fun, "reference")) {
    args$reference <- reference
  }

  args <- .filter_call_arguments(fit_fun, args)
  do.call(fit_fun, args)
}

.run_with_timeout <- function(fun, timeout_sec, on_timeout) {
  use_timeout <- !is.null(timeout_sec) && is.finite(timeout_sec) && timeout_sec > 0
  if (!use_timeout) return(fun())

  tryCatch(
    R.utils::withTimeout(
      fun(),
      timeout = timeout_sec,
      onTimeout = "error"
    ),
    error = function(e) {
      if (!inherits(e, "TimeoutException")) stop(e)

      message_text <- sprintf("TIMEOUT after %.1f seconds", timeout_sec)
      if (identical(on_timeout, "warn")) warning(message_text, call. = FALSE)
      if (identical(on_timeout, "silent")) message_text <- "TIMEOUT"
      stop(message_text, call. = FALSE)
    }
  )
}

.safe_fit <- function(
    fit_fun,
    data,
    gene,
    spec,
    t.label,
    t.diff,
    h5_path,
    reference,
    ndraw,
    seed,
    timeout_sec,
    on_timeout
) {
  tryCatch(
    .run_with_timeout(
      fun = function() {
        .call_fit_fun(
          fit_fun = fit_fun,
          data = data,
          gene = gene,
          spec = spec,
          t.label = t.label,
          t.diff = t.diff,
          h5_path = h5_path,
          reference = reference,
          ndraw = ndraw,
          seed = seed
        )
      },
      timeout_sec = timeout_sec,
      on_timeout = on_timeout
    ),
    error = function(e) e
  )
}

.is_fit_error <- function(x) inherits(x, "error")

.fit_meta_row <- function(gene, spec, result, runtime_sec = NA_real_) {
  failed <- .is_fit_error(result)
  data.frame(
    gene = gene,
    ENSG = if (!failed && "ENSG" %in% names(result)) result$ENSG[1L] else NA_character_,
    fit_cond = spec$condition,
    steady = spec$steady,
    ref_cond = if (is.null(spec$ref_cond)) NA_character_ else spec$ref_cond,
    ok = !failed,
    error = if (failed) conditionMessage(result) else NA_character_,
    runtime_sec = runtime_sec,
    runtime_min = runtime_sec / 60,
    stringsAsFactors = FALSE
  )
}

.valid_fit_objects <- function(results) {
  results[!vapply(results, .is_fit_error, logical(1))]
}

.first_fit_ensg <- function(results) {
  valid <- .valid_fit_objects(results)
  if (length(valid) == 0L) return(NA_character_)
  as.character(valid[[1L]]$ENSG[1L])
}

.comparison_base_row <- function(gene, specs, results) {
  two_conditions <- length(specs) == 2L
  data.frame(
    gene = gene,
    ENSG = .first_fit_ensg(results),
    cond1 = specs[[1L]]$condition,
    cond2 = if (two_conditions) specs[[2L]]$condition else NA_character_,
    ok = !any(vapply(results, .is_fit_error, logical(1))),
    error_fit1 = if (.is_fit_error(results[[1L]])) {
      conditionMessage(results[[1L]])
    } else {
      NA_character_
    },
    error_fit2 = if (two_conditions && .is_fit_error(results[[2L]])) {
      conditionMessage(results[[2L]])
    } else {
      NA_character_
    },
    error_cmp = NA_character_,
    interpretation = NA_character_,
    plot_path = NA_character_,
    stringsAsFactors = FALSE
  )
}

.flatten_comparison_summary <- function(summary_table) {

  if (
    is.null(summary_table) ||
    !is.data.frame(summary_table) ||
    nrow(summary_table) == 0L
  ) {
    return(NULL)
  }

  # New km/B gate-star output:
  # summary is already a one-row wide table.

  if (!"parameter" %in% names(summary_table)) {

    if (nrow(summary_table) != 1L) {
      warning(
        paste0(
          "comparison$summary has no 'parameter' column and ",
          "contains more than one row; only the first row will be used."
        ),
        call. = FALSE
      )
    }

    out <- summary_table[
      1L,
      ,
      drop = FALSE
    ]

    rownames(out) <- NULL
    return(out)
  }


  wanted <- c(
    "median_cond1",
    "median_cond2",
    "median_log2FC",
    "log2FC_lo",
    "log2FC_hi",
    "p_up",
    "p_down",
    "p_change",
    "direction"
  )

  make_parameter_columns <- function(parameter) {

    row <- summary_table[
      summary_table$parameter == parameter,
      ,
      drop = FALSE
    ]

    if (nrow(row) == 0L) {
      return(NULL)
    }

    row <- row[
      1L,
      ,
      drop = FALSE
    ]

    values <- lapply(
      wanted,
      function(column) {
        if (column %in% names(row)) {
          row[[column]][1L]
        } else {
          NA
        }
      }
    )

    names(values) <- paste0(
      wanted,
      "_",
      parameter
    )

    as.data.frame(
      values,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  pieces <- lapply(
    c("lambda", "d", "km", "B"),
    make_parameter_columns
  )

  pieces <- pieces[
    !vapply(
      pieces,
      is.null,
      logical(1)
    )
  ]

  if (length(pieces) == 0L) {
    return(NULL)
  }

  do.call(
    cbind,
    pieces
  )
}

.safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unnamed")
}

.save_comparison_plot <- function(comparison, plot_dir, gene, cond1, cond2) {
  if (is.null(plot_dir) || is.null(comparison$plot)) return(NA_character_)
  if (!inherits(comparison$plot, "ggplot")) return(NA_character_)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("A comparison plot was returned, but ggplot2 is unavailable.", call. = FALSE)
    return(NA_character_)
  }

  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(
    plot_dir,
    paste0(
      .safe_filename(gene), "__",
      .safe_filename(cond1), "__vs__",
      .safe_filename(cond2), ".pdf"
    )
  )
  ggplot2::ggsave(filename = path, plot = comparison$plot)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.compare_fit_results <- function(
    gene,
    specs,
    results,
    compare_fun,
    compare_args,
    do_compare,
    plot_dir
) {
  row <- .comparison_base_row(gene, specs, results)
  if (length(specs) != 2L || !isTRUE(do_compare) || !isTRUE(row$ok)) return(row)

  comparison <- tryCatch(
    do.call(compare_fun, c(list(results[[1L]], results[[2L]]), compare_args)),
    error = function(e) e
  )

  if (.is_fit_error(comparison)) {
    row$ok <- FALSE
    row$error_cmp <- conditionMessage(comparison)
    return(row)
  }
  if (!is.list(comparison)) {
    row$ok <- FALSE
    row$error_cmp <- "compare_fun must return a list."
    return(row)
  }

  if (!is.null(comparison$interpretation)) {
    row$interpretation <- paste(as.character(comparison$interpretation), collapse = "; ")
  }

  flattened <- .flatten_comparison_summary(comparison$summary)
  if (!is.null(flattened)) row <- cbind(row, flattened)

  row$plot_path <- .save_comparison_plot(
    comparison = comparison,
    plot_dir = plot_dir,
    gene = gene,
    cond1 = specs[[1L]]$condition,
    cond2 = specs[[2L]]$condition
  )
  row
}

.save_fit_rds <- function(fit, directory, gene, condition) {
  if (is.null(directory) || .is_fit_error(fit)) return(NA_character_)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(
    directory,
    paste0(.safe_filename(gene), "__", .safe_filename(condition), ".rds")
  )
  saveRDS(fit, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.lapply_batch <- function(X, FUN, parallel, ncores) {
  if (!isTRUE(parallel) || .Platform$OS.type == "windows" || ncores <= 1L) {
    return(lapply(X, FUN))
  }

  parallel::mclapply(
    X,
    FUN,
    mc.cores = as.integer(ncores),
    mc.preschedule = TRUE
  )
}

.run_kinetics_engine <- function(
    data,
    genes,
    cond1,
    cond2,
    t.label,
    t.diff,
    fit1_steady,
    fit2_steady,
    fit1_ref_cond,
    fit2_ref_cond,
    h5_path,
    reference,
    ndraw,
    seed,
    fit_fun,
    truth_suffix1,
    truth_suffix2,
    verbose_every,
    require_h5,
    keep_only_unique_symbols,
    parallel,
    ncores,
    timeout_sec,
    on_timeout,
    do_compare = FALSE,
    compare_fun = NULL,
    compare_args = list(),
    plot_dir = NULL,
    collect_runtime = FALSE,
    retain_fit_objects = FALSE,
    save_rds_dir = NULL
) {
  on_timeout <- match.arg(on_timeout, c("error", "warn", "silent"))
  use_timeout <- !is.null(timeout_sec) && is.finite(timeout_sec) && timeout_sec > 0
  .require_batch_namespaces(require_h5 = require_h5, use_timeout = use_timeout)

  if (!is.function(fit_fun)) stop("fit_fun must be a function.", call. = FALSE)
  if (isTRUE(do_compare) && !is.function(compare_fun)) {
    stop("compare_fun must be a function when do_compare = TRUE.", call. = FALSE)
  }
  if (!is.numeric(t.label) || length(t.label) != 1L ||
      !is.finite(t.label) || t.label <= 0) {
    stop("t.label must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(ndraw) || length(ndraw) != 1L || !is.finite(ndraw) || ndraw < 1) {
    stop("ndraw must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be one finite number.", call. = FALSE)
  }
  if (!is.list(compare_args)) {
    stop("compare_args must be a list.", call. = FALSE)
  }
  if (!is.numeric(ncores) || length(ncores) != 1L || ncores < 1) {
    stop("ncores must be a positive integer.", call. = FALSE)
  }
  if (length(verbose_every) != 1L || is.na(verbose_every) || verbose_every < 0) {
    stop("verbose_every must be a non-negative integer.", call. = FALSE)
  }

  ndraw <- as.integer(ndraw)
  ncores <- as.integer(ncores)
  verbose_every <- as.integer(verbose_every)
  set.seed(seed)

  fit_setup <- .prepare_fit_specs(
    cond1 = cond1,
    cond2 = cond2,
    fit1_steady = fit1_steady,
    fit2_steady = fit2_steady,
    fit1_ref_cond = fit1_ref_cond,
    fit2_ref_cond = fit2_ref_cond,
    truth_suffix1 = truth_suffix1,
    truth_suffix2 = truth_suffix2
  )
  specs <- fit_setup$specs
  if (!fit_setup$two_conditions) do_compare <- FALSE

  any_nonsteady <- any(!vapply(specs, `[[`, logical(1), "steady"))
  if (any_nonsteady && (is.null(t.diff) || length(t.diff) != 1L ||
                        !is.numeric(t.diff) || !is.finite(t.diff) || t.diff < 0)) {
    stop("t.diff must be one non-negative finite number for non-steady fits.",
         call. = FALSE)
  }

  if (any_nonsteady) {
    .validate_labeling_interval(t.label, steady = FALSE, t.diff = t.diff)
  }

  mapped <- .prepare_batch_gene_map(
    data = data,
    genes = genes,
    h5_path = h5_path,
    require_h5 = require_h5,
    keep_only_unique_symbols = keep_only_unique_symbols
  )

  shared_reference <- .prepare_shared_reference(
    fit_fun = fit_fun,
    h5_path = h5_path,
    reference = reference
  )

  run_one_gene <- function(index) {
    gene <- mapped$genes_ok[index]
    gene_start <- proc.time()[["elapsed"]]

    if (verbose_every > 0L && (index == 1L || index %% verbose_every == 0L)) {
      message(sprintf("[%d/%d] %s", index, length(mapped$genes_ok), gene))
    }

    executions <- lapply(specs, function(spec) {
      fit_start <- proc.time()[["elapsed"]]
      result <- .safe_fit(
        fit_fun = fit_fun,
        data = data,
        gene = gene,
        spec = spec,
        t.label = t.label,
        t.diff = t.diff,
        h5_path = h5_path,
        reference = shared_reference,
        ndraw = ndraw,
        seed = seed,
        timeout_sec = timeout_sec,
        on_timeout = on_timeout
      )
      list(
        result = result,
        runtime_sec = proc.time()[["elapsed"]] - fit_start
      )
    })

    results <- lapply(executions, `[[`, "result")
    fit_runtimes <- vapply(executions, `[[`, numeric(1), "runtime_sec")
    names(results) <- vapply(specs, `[[`, character(1), "condition")

    if (!is.null(save_rds_dir)) {
      for (i in seq_along(results)) {
        .save_fit_rds(results[[i]], save_rds_dir, gene, specs[[i]]$condition)
      }
    }

    fit_meta <- do.call(rbind, lapply(seq_along(results), function(i) {
      .fit_meta_row(
        gene = gene,
        spec = specs[[i]],
        result = results[[i]],
        runtime_sec = fit_runtimes[i]
      )
    }))
    rownames(fit_meta) <- NULL

    comparison <- .compare_fit_results(
      gene = gene,
      specs = specs,
      results = results,
      compare_fun = compare_fun,
      compare_args = compare_args,
      do_compare = do_compare,
      plot_dir = plot_dir
    )

    gene_runtime <- proc.time()[["elapsed"]] - gene_start
    if (isTRUE(collect_runtime)) {
      comparison$runtime_sec <- gene_runtime
      comparison$runtime_min <- gene_runtime / 60
    }

    list(
      gene = gene,
      results = results,
      fit_meta = fit_meta,
      comparison = comparison
    )
  }

  raw_results <- .lapply_batch(
    seq_along(mapped$genes_ok),
    run_one_gene,
    parallel = parallel,
    ncores = ncores
  )

  fit_rows <- unlist(lapply(raw_results, function(x) {
    .valid_fit_objects(x$results)
  }), recursive = FALSE)
  fits <- .rbind_fill(fit_rows)
  fit_meta <- .rbind_fill(lapply(raw_results, `[[`, "fit_meta"))
  comparisons <- .rbind_fill(lapply(raw_results, `[[`, "comparison"))

  fit_objects <- NULL
  if (isTRUE(retain_fit_objects)) {
    fit_objects <- setNames(vector("list", length(raw_results)), mapped$genes_ok)
    for (i in seq_along(raw_results)) {
      valid <- .valid_fit_objects(raw_results[[i]]$results)
      fit_objects[[i]] <- valid
    }
  }

  list(
    gene_map = mapped$gene_map,
    genes_ok = mapped$genes_ok,
    fits = fits,
    fit_meta = fit_meta,
    fit_objects = fit_objects,
    comparisons = comparisons
  )
}


# Public wrappers --------------------------------------------------------------

#' Run bursting inference for many genes
#'
#' Fits one or two conditions per gene and optionally compares posterior draws.
#' Non-steady fits use the explicitly supplied reference condition, or the other
#' fitted condition when two conditions are present.
#'
#' @param data A grandR object.
#' @param genes Character vector of gene symbols or ENSG identifiers.
#' @param cond1 First condition.
#' @param cond2 Optional second condition.
#' @param t.label Labeling duration.
#' @param t.diff Modeled reference-to-observation interval for non-steady fits.
#'   Must equal `t.label`; pre-labeling synthesis is not modeled.
#' @param fit1_steady,fit2_steady Whether each fit uses the steady-state model.
#' @param fit1_ref_cond,fit2_ref_cond Optional reference conditions.
#' @param h5_path HDF5 gene-profile path.
#' @param reference Optional preloaded 10x reference. It is loaded once when
#'   omitted and `fit_fun` explicitly accepts a `reference` argument.
#' @param ndraw Number of Laplace draws per fit.
#' @param seed Seed passed to each fit for backward-compatible reproducibility.
#' @param compare_fun Posterior-comparison function.
#' @param fit_fun Single-gene fitting function.
#' @param do_compare Whether to compare two successful fits.
#' @param truth_suffix1,truth_suffix2 Optional truth-column suffixes.
#' @param plot_dir Optional directory for a `comparison$plot` ggplot.
#' @param verbose_every Print progress every this many genes; zero disables it.
#' @param require_h5 Drop genes absent from the HDF5 file.
#' @param keep_only_unique_symbols Drop symbols mapping to multiple ENSG IDs.
#' @param parallel Use `parallel::mclapply()` on non-Windows systems.
#' @param ncores Number of forked workers.
#' @param timeout_sec Optional timeout per fit.
#' @param on_timeout One of `"error"`, `"warn"`, or `"silent"`.
#' @param compare_args Named arguments passed to `compare_fun` after the two fits.
#'
#' @return A list containing `gene_map`, `genes_ok`, `fits`, and `comparisons`.
#' @export
run_kinetics_batch <- function(
    data,
    genes,
    cond1,
    cond2 = NULL,
    t.label = 2,
    t.diff = 2,
    fit1_steady = TRUE,
    fit2_steady = TRUE,
    fit1_ref_cond = NULL,
    fit2_ref_cond = NULL,
    h5_path,
    ndraw = 20000,
    seed = 1,
    compare_fun = compare_posteriors_gate_star_kmB,
    fit_fun = fit_gene,
    do_compare = TRUE,
    truth_suffix1 = NULL,
    truth_suffix2 = NULL,
    plot_dir = NULL,
    verbose_every = 50,
    require_h5 = TRUE,
    keep_only_unique_symbols = TRUE,
    parallel = FALSE,
    ncores = 2L,
    timeout_sec = NULL,
    on_timeout = c("error", "warn", "silent"),
    compare_args = list(),
    reference = NULL
) {
  result <- .run_kinetics_engine(
    data = data,
    genes = genes,
    cond1 = cond1,
    cond2 = cond2,
    t.label = t.label,
    t.diff = t.diff,
    fit1_steady = fit1_steady,
    fit2_steady = fit2_steady,
    fit1_ref_cond = fit1_ref_cond,
    fit2_ref_cond = fit2_ref_cond,
    h5_path = h5_path,
    reference = reference,
    ndraw = ndraw,
    seed = seed,
    fit_fun = fit_fun,
    truth_suffix1 = truth_suffix1,
    truth_suffix2 = truth_suffix2,
    verbose_every = verbose_every,
    require_h5 = require_h5,
    keep_only_unique_symbols = keep_only_unique_symbols,
    parallel = parallel,
    ncores = ncores,
    timeout_sec = timeout_sec,
    on_timeout = on_timeout,
    do_compare = do_compare,
    compare_fun = compare_fun,
    compare_args = compare_args,
    plot_dir = plot_dir,
    collect_runtime = FALSE,
    retain_fit_objects = FALSE,
    save_rds_dir = NULL
  )

  result[c("gene_map", "genes_ok", "fits", "comparisons")]
}

#' Run bursting fits and retain complete fit objects
#'
#' @inheritParams run_kinetics_batch
#' @param save_rds_dir Optional directory in which each successful fit is saved.
#'
#' @return A list containing the flat summaries, per-condition metadata, and a
#'   nested `fit_objects[[gene]][[condition]]` store.
#' @export
run_kinetics_fits <- function(
    data,
    genes,
    cond1,
    cond2 = NULL,
    t.label = 2,
    t.diff = 2,
    fit1_steady = TRUE,
    fit2_steady = TRUE,
    fit1_ref_cond = NULL,
    fit2_ref_cond = NULL,
    h5_path,
    ndraw = 20000,
    seed = 1,
    fit_fun = fit_gene,
    truth_suffix1 = NULL,
    truth_suffix2 = NULL,
    verbose_every = 50,
    require_h5 = TRUE,
    keep_only_unique_symbols = TRUE,
    parallel = FALSE,
    ncores = 2L,
    timeout_sec = NULL,
    on_timeout = c("error", "warn", "silent"),
    save_rds_dir = NULL,
    reference = NULL
) {
  result <- .run_kinetics_engine(
    data = data,
    genes = genes,
    cond1 = cond1,
    cond2 = cond2,
    t.label = t.label,
    t.diff = t.diff,
    fit1_steady = fit1_steady,
    fit2_steady = fit2_steady,
    fit1_ref_cond = fit1_ref_cond,
    fit2_ref_cond = fit2_ref_cond,
    h5_path = h5_path,
    reference = reference,
    ndraw = ndraw,
    seed = seed,
    fit_fun = fit_fun,
    truth_suffix1 = truth_suffix1,
    truth_suffix2 = truth_suffix2,
    verbose_every = verbose_every,
    require_h5 = require_h5,
    keep_only_unique_symbols = keep_only_unique_symbols,
    parallel = parallel,
    ncores = ncores,
    timeout_sec = timeout_sec,
    on_timeout = on_timeout,
    do_compare = FALSE,
    compare_fun = NULL,
    plot_dir = NULL,
    collect_runtime = FALSE,
    retain_fit_objects = TRUE,
    save_rds_dir = save_rds_dir
  )

  result[c("gene_map", "genes_ok", "fits", "fit_meta", "fit_objects")]
}

#' Run bursting inference with per-gene runtime statistics
#'
#' @inheritParams run_kinetics_batch
#'
#' @return The same structure as [run_kinetics_batch()], with `runtime_sec` and
#'   `runtime_min` columns in `comparisons`.
#' @export
run_kinetics_batch_stat <- function(
    data,
    genes,
    cond1,
    cond2 = NULL,
    t.label = 2,
    t.diff = 2,
    fit1_steady = TRUE,
    fit2_steady = TRUE,
    fit1_ref_cond = NULL,
    fit2_ref_cond = NULL,
    h5_path,
    ndraw = 20000,
    seed = 1,
    compare_fun = compare_posteriors_gate_star_kmB,
    fit_fun = fit_gene,
    do_compare = TRUE,
    truth_suffix1 = NULL,
    truth_suffix2 = NULL,
    plot_dir = NULL,
    verbose_every = 50,
    require_h5 = TRUE,
    keep_only_unique_symbols = TRUE,
    parallel = FALSE,
    ncores = 2L,
    timeout_sec = NULL,
    on_timeout = c("error", "warn", "silent"),
    compare_args = list(),
    reference = NULL
) {
  result <- .run_kinetics_engine(
    data = data,
    genes = genes,
    cond1 = cond1,
    cond2 = cond2,
    t.label = t.label,
    t.diff = t.diff,
    fit1_steady = fit1_steady,
    fit2_steady = fit2_steady,
    fit1_ref_cond = fit1_ref_cond,
    fit2_ref_cond = fit2_ref_cond,
    h5_path = h5_path,
    reference = reference,
    ndraw = ndraw,
    seed = seed,
    fit_fun = fit_fun,
    truth_suffix1 = truth_suffix1,
    truth_suffix2 = truth_suffix2,
    verbose_every = verbose_every,
    require_h5 = require_h5,
    keep_only_unique_symbols = keep_only_unique_symbols,
    parallel = parallel,
    ncores = ncores,
    timeout_sec = timeout_sec,
    on_timeout = on_timeout,
    do_compare = do_compare,
    compare_fun = compare_fun,
    compare_args = compare_args,
    plot_dir = plot_dir,
    collect_runtime = TRUE,
    retain_fit_objects = FALSE,
    save_rds_dir = NULL
  )

  result[c("gene_map", "genes_ok", "fits", "comparisons")]
}
