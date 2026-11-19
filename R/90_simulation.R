# simulation.R
# -------------------------------------------------------------------
# Simulation utilities for buRst:
#   1) Gillespie simulation of new/old/total RNA
#   2) conversion of simulated matrices into a grandR object
#   3) GRAND-SLAM-like read-evidence simulation
#   4) export of a self-contained <prefix>.targets directory containing:
#        - 4sU.Binom.discrete.h5
#        - matrix.mtx.gz
#        - barcodes.tsv.gz
#        - features.tsv.gz
# -------------------------------------------------------------------

#' @useDynLib buRst, .registration = TRUE
NULL


# -------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------

.sim_or <- function(x, fallback) {
  if (is.null(x) || length(x) == 0L || is.na(x[1L]) || !nzchar(as.character(x[1L]))) {
    fallback
  } else {
    x
  }
}


.resolve_ncores <- function(ncores = NULL) {
  if (is.null(ncores)) {
    ncores <- max(1L, parallel::detectCores() - 1L)
  }

  ncores <- as.integer(ncores[1L])
  if (!is.finite(ncores) || ncores < 1L) {
    stop("ncores must be a positive integer.")
  }

  ncores
}


.sim_result_field <- function(result, field, gene_index = NA_integer_) {
  value <- tryCatch(
    result[[field]],
    error = function(e) NULL
  )

  if (is.null(value)) {
    stop(
      "Malformed SimulateReads result",
      if (is.finite(gene_index)) paste0(" for gene index ", gene_index) else "",
      ": missing field '", field, "'. Result class: ",
      paste(class(result), collapse = "/"),
      call. = FALSE
    )
  }

  value
}


.extract_ntr_beta_fit <- function(x, gene = NULL, barcode = NULL) {
  context <- paste0(
    if (!is.null(gene)) paste0(" for gene '", gene, "'") else "",
    if (!is.null(barcode)) paste0(" and cell '", barcode, "'") else ""
  )

  get_named_value <- function(obj, candidates) {
    nm <- names(obj)
    if (is.null(nm)) return(NULL)

    nm_lower <- tolower(nm)
    idx <- match(tolower(candidates), nm_lower, nomatch = 0L)
    idx <- idx[idx > 0L]

    if (length(idx) == 0L) return(NULL)
    obj[[idx[1L]]]
  }

  ntr <- get_named_value(
    x,
    c("ntr", "ntr_hat", "ntr.hat", "estimate", "mean", "map")
  )
  alpha <- get_named_value(
    x,
    c("alpha", "shape1", "beta_alpha", "beta.alpha")
  )
  beta <- get_named_value(
    x,
    c("beta", "shape2", "beta_beta", "beta.beta")
  )

  # Some versions may store the beta approximation in an attribute.
  if (is.null(alpha) || is.null(beta)) {
    beta_attr <- attr(x, "beta")
    if (!is.null(beta_attr)) {
      if (is.null(alpha)) {
        alpha <- get_named_value(beta_attr, c("alpha", "shape1"))
      }
      if (is.null(beta)) {
        beta <- get_named_value(beta_attr, c("beta", "shape2"))
      }
    }
  }

  # Last-resort positional fallback for an unnamed length-3 atomic return.
  if (
    (is.atomic(x) || is.list(x)) &&
    length(x) >= 3L &&
    is.null(ntr) &&
    is.null(alpha) &&
    is.null(beta)
  ) {
    vals <- suppressWarnings(as.numeric(unlist(x[seq_len(3L)])))
    if (length(vals) == 3L && all(is.finite(vals))) {
      ntr <- vals[1L]
      alpha <- vals[2L]
      beta <- vals[3L]
    }
  }

  ntr <- suppressWarnings(as.numeric(ntr)[1L])
  alpha <- suppressWarnings(as.numeric(alpha)[1L])
  beta <- suppressWarnings(as.numeric(beta)[1L])

  if (!is.finite(ntr)) {
    stop(
      "Could not extract an NTR estimate from fit.ntr(beta.approx=TRUE)",
      context,
      ". Return class: ", paste(class(x), collapse = "/"),
      "; names: ", paste(names(x), collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.finite(alpha) || alpha <= 0 || !is.finite(beta) || beta <= 0) {
    stop(
      "Could not extract valid positive alpha and beta values from ",
      "fit.ntr(beta.approx=TRUE)", context,
      ". Return class: ", paste(class(x), collapse = "/"),
      "; names: ", paste(names(x), collapse = ", "),
      call. = FALSE
    )
  }

  list(
    ntr = ntr,
    alpha = alpha,
    beta = beta
  )
}


.simulate_reads_stage <- function(expr, stage) {
  tryCatch(
    force(expr),
    error = function(e) {
      original_call <- conditionCall(e)
      call_text <- if (is.null(original_call)) {
        ""
      } else {
        paste0(" [original call: ", paste(deparse(original_call), collapse = " "), "]")
      }

      stop(
        "SimulateReads failed during '", stage, "': ",
        conditionMessage(e),
        call_text,
        call. = FALSE
      )
    }
  )
}


.normalize_gene_info <- function(gene.info, n_genes) {
  if (!is.data.frame(gene.info)) {
    gene.info <- data.frame(
      Gene = as.character(gene.info),
      Symbol = as.character(gene.info),
      stringsAsFactors = FALSE
    )
  } else {
    gene.info <- as.data.frame(gene.info, stringsAsFactors = FALSE)
  }

  if (!all(c("Gene", "Symbol") %in% names(gene.info))) {
    stop("gene.info must contain columns 'Gene' and 'Symbol'.")
  }

  if (nrow(gene.info) != n_genes) {
    stop(
      "gene.info must contain exactly one row per simulated gene: expected ",
      n_genes, ", got ", nrow(gene.info), "."
    )
  }

  gene.info$Gene <- as.character(gene.info$Gene)
  gene.info$Symbol <- as.character(gene.info$Symbol)

  if (anyNA(gene.info$Gene) || any(!nzchar(gene.info$Gene))) {
    stop("gene.info$Gene must not contain missing or empty identifiers.")
  }

  if (anyDuplicated(gene.info$Gene)) {
    stop("gene.info$Gene must contain unique identifiers.")
  }

  missing_symbol <- is.na(gene.info$Symbol) | !nzchar(gene.info$Symbol)
  gene.info$Symbol[missing_symbol] <- gene.info$Gene[missing_symbol]

  gene.info
}


.make_simulation_barcodes <- function(n) {
  alphabet <- c("A", "C", "G", "T")
  out <- character(n)
  seen <- new.env(hash = TRUE, parent = emptyenv())

  for (i in seq_len(n)) {
    repeat {
      barcode <- paste0(sample(alphabet, 20L, replace = TRUE), collapse = "")
      if (!exists(barcode, envir = seen, inherits = FALSE)) {
        assign(barcode, TRUE, envir = seen)
        out[i] <- barcode
        break
      }
    }
  }

  out
}


.validate_count_matrices <- function(new, old) {
  if (length(dim(new)) != 2L || length(dim(old)) != 2L) {
    stop("new and old must both be two-dimensional matrices.")
  }

  if (!identical(dim(new), dim(old))) {
    stop("new and old must have identical dimensions.")
  }

  if (any(!is.finite(new)) || any(!is.finite(old))) {
    stop("new and old must contain only finite values.")
  }

  if (any(new < 0) || any(old < 0)) {
    stop("new and old counts must be non-negative.")
  }

  invisible(TRUE)
}


.resolve_targets_paths <- function(
    data,
    h5_path,
    output_prefix = NULL,
    targets_dir = NULL,
    h5_filename = "4sU.Binom.discrete.h5",
    write_10x = FALSE
) {
  use_bundle <- isTRUE(write_10x) || !is.null(output_prefix) || !is.null(targets_dir)

  if (!use_bundle) {
    h5_path <- .sim_or(h5_path, "./sim_profiles.h5")

    return(list(
      bundled = FALSE,
      targets_dir = NULL,
      h5_path = h5_path,
      matrix_path = NULL,
      barcodes_path = NULL,
      features_path = NULL
    ))
  }

  if (is.null(targets_dir)) {
    if (is.null(output_prefix)) {
      prefix <- .sim_or(data$prefix, "Simulated")
      output_prefix <- prefix
    }

    output_prefix <- sub("[/\\\\]+$", "", as.character(output_prefix)[1L])
    targets_dir <- if (grepl("\\.targets$", output_prefix)) {
      output_prefix
    } else {
      paste0(output_prefix, ".targets")
    }
  }

  targets_dir <- sub("[/\\\\]+$", "", as.character(targets_dir)[1L])

  if (!nzchar(targets_dir)) {
    stop("targets_dir resolved to an empty path.")
  }

  h5_name <- .sim_or(h5_filename, "4sU.Binom.discrete.h5")

  # In bundled mode, h5_filename controls the file inside <prefix>.targets.
  # h5_path is retained only for legacy, non-bundled output.

  list(
    bundled = TRUE,
    targets_dir = targets_dir,
    h5_path = file.path(targets_dir, h5_name),
    matrix_path = file.path(targets_dir, "matrix.mtx.gz"),
    barcodes_path = file.path(targets_dir, "barcodes.tsv.gz"),
    features_path = file.path(targets_dir, "features.tsv.gz")
  )
}


.copy_file_to_gzip <- function(input_path, output_path, chunk_size = 1024L * 1024L) {
  input <- file(input_path, open = "rb")
  output <- gzfile(output_path, open = "wb")

  on.exit(close(input), add = TRUE)
  on.exit(close(output), add = TRUE)

  repeat {
    block <- readBin(input, what = "raw", n = as.integer(chunk_size))
    if (length(block) == 0L) break
    writeBin(block, output)
  }

  invisible(output_path)
}


.write_lines_gz <- function(lines, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(as.character(lines), con = con, sep = "\n", useBytes = TRUE)
  invisible(path)
}


.get_10x_export_parts <- function(
    data,
    barcode_column = "Name"
) {
  if (!is.list(data)) {
    stop("data must be a grandR-like list object.")
  }

  counts <- data$data$count
  coldata <- data$coldata
  gene_info <- data$gene.info

  if (is.null(counts)) {
    stop("Missing data$data$count.")
  }
  if (is.null(coldata) || !is.data.frame(coldata)) {
    stop("Missing or invalid data$coldata.")
  }
  if (is.null(gene_info) || !is.data.frame(gene_info)) {
    stop("Missing or invalid data$gene.info.")
  }

  counts <- Matrix::Matrix(counts, sparse = TRUE)

  matrix_gene_ids <- rownames(counts)
  matrix_cell_ids <- colnames(counts)

  if (is.null(matrix_gene_ids) || length(matrix_gene_ids) != nrow(counts)) {
    stop("data$data$count must have one row name per gene.")
  }
  if (is.null(matrix_cell_ids) || length(matrix_cell_ids) != ncol(counts)) {
    stop("data$data$count must have one column name per cell.")
  }

  # Match gene.info to the matrix row order.
  if (!"Gene" %in% names(gene_info)) {
    stop("data$gene.info must contain a 'Gene' column.")
  }

  gene_index <- match(matrix_gene_ids, as.character(gene_info$Gene))
  if (anyNA(gene_index)) {
    missing <- matrix_gene_ids[is.na(gene_index)]
    stop(
      "Some matrix genes are missing from data$gene.info: ",
      paste(utils::head(missing, 10L), collapse = ", ")
    )
  }

  gene_info <- gene_info[gene_index, , drop = FALSE]
  gene_ids <- as.character(gene_info$Gene)

  gene_names <- if ("Symbol" %in% names(gene_info)) {
    as.character(gene_info$Symbol)
  } else {
    gene_ids
  }

  missing_symbol <- is.na(gene_names) | !nzchar(gene_names)
  gene_names[missing_symbol] <- gene_ids[missing_symbol]

  # Match coldata to the matrix column order. Name is the correct default here
  # because it is identical to the cell IDs stored in the HDF5 profiles.
  if (!barcode_column %in% names(coldata)) {
    stop(
      "Column '", barcode_column, "' is missing from data$coldata. ",
      "Available columns: ", paste(names(coldata), collapse = ", ")
    )
  }

  if ("Name" %in% names(coldata)) {
    col_index <- match(matrix_cell_ids, as.character(coldata$Name))
    if (anyNA(col_index)) {
      missing <- matrix_cell_ids[is.na(col_index)]
      stop(
        "Some matrix cell names are missing from data$coldata$Name: ",
        paste(utils::head(missing, 10L), collapse = ", ")
      )
    }
    coldata <- coldata[col_index, , drop = FALSE]
  } else if (nrow(coldata) != ncol(counts)) {
    stop("data$coldata rows do not match matrix columns.")
  }

  barcodes <- as.character(coldata[[barcode_column]])

  if (length(barcodes) != ncol(counts) || anyNA(barcodes) || any(!nzchar(barcodes))) {
    stop("The selected barcode column must contain one non-empty value per matrix column.")
  }
  if (anyDuplicated(barcodes)) {
    stop("The selected barcode column contains duplicate values.")
  }

  list(
    counts = counts,
    barcodes = barcodes,
    features = data.frame(
      gene_id = gene_ids,
      gene_name = gene_names,
      feature_type = rep("Gene Expression", length(gene_ids)),
      stringsAsFactors = FALSE
    )
  )
}


#' Write a grandR simulation as a 10x-compatible reference
#'
#' The export uses the actual grandR object layout:
#' \code{data$data$count}, \code{data$coldata}, and \code{data$gene.info}.
#'
#' @param data A grandR object.
#' @param output_dir Destination directory.
#' @param slot Retained for backward compatibility; ignored.
#' @param overwrite Whether existing files may be replaced.
#' @param barcode_column Column in \code{data$coldata} written to
#'   \code{barcodes.tsv.gz}. The default \code{"Name"} matches the matrix and
#'   HDF5 cell identifiers. Use \code{"Barcode"} only when raw barcodes are
#'   explicitly required.
#'
#' @return Named paths to matrix, barcode, and feature files.
#' @details
#' The exported gene order and cell identifiers must match the likelihood file
#' used for fitting. The default `barcode_column = "Name"` preserves full cell
#' identifiers. Do not replace them by raw barcode substrings unless all other
#' files use those exact same identifiers.
#'
#' @seealso `SimulateReads()`, `merge_targets_dirs()`
#'
#' @examples
#' sim <- simulateBursting(
#'   kms = 1, Bs = 2, ds = 0.5, t.end = 1, no.cells = 8,
#'   gene.info = "gene_A", condition = "sim", ntr.sim = FALSE,
#'   parallel = FALSE, seed = 1, on_gene_error = "stop"
#' )
#' work <- tempfile("buRst_reference_")
#' files <- write_10x_reference(sim, output_dir = work)
#' file.exists(unlist(files))
#' unlink(work, recursive = TRUE)
#'
#' @export
write_10x_reference <- function(
    data,
    output_dir,
    slot = "count",
    overwrite = TRUE,
    barcode_column = "Name"
) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Need Matrix.")

  output_dir <- as.character(output_dir)[1L]
  if (!nzchar(output_dir)) stop("output_dir must not be empty.")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- list(
    matrix = file.path(output_dir, "matrix.mtx.gz"),
    barcodes = file.path(output_dir, "barcodes.tsv.gz"),
    features = file.path(output_dir, "features.tsv.gz")
  )

  existing <- vapply(paths, file.exists, logical(1))
  if (any(existing) && !isTRUE(overwrite)) {
    stop(
      "Output file(s) already exist: ",
      paste(unlist(paths[existing]), collapse = ", ")
    )
  }

  parts <- .get_10x_export_parts(
    data = data,
    barcode_column = barcode_column
  )

  # Let Matrix write the Matrix Market format. Then gzip the resulting file.
  # This avoids all Matrix-version-specific assumptions about summary() or slots.
  temporary_mtx <- tempfile(
    pattern = "buRst_matrix_",
    tmpdir = output_dir,
    fileext = ".mtx"
  )
  on.exit(unlink(temporary_mtx), add = TRUE)

  Matrix::writeMM(parts$counts, file = temporary_mtx)
  .copy_file_to_gzip(temporary_mtx, paths$matrix)

  .write_lines_gz(parts$barcodes, paths$barcodes)

  feature_lines <- paste(
    parts$features$gene_id,
    parts$features$gene_name,
    parts$features$feature_type,
    sep = "\t"
  )
  .write_lines_gz(feature_lines, paths$features)

  paths
}


#' Retrieve files created by SimulateReads
#'
#' @param data A grandR object returned by \code{SimulateReads()}.
#'
#' @return A named list of paths, or \code{NULL} when no manifest is attached.
#' @details
#' Count-only simulations do not create an output manifest. After `SimulateReads()`
#' or a read-evidence simulation, this list records the actual exported paths;
#' it does not copy or recreate the files.
#'
#' @seealso `SimulateReads()`, `simulateBursting()`
#'
#' @examples
#' sim <- simulateBursting(
#'   kms = 1, Bs = 2, ds = 0.5, t.end = 1, no.cells = 8,
#'   gene.info = "gene_A", condition = "sim", ntr.sim = FALSE,
#'   parallel = FALSE, seed = 1, on_gene_error = "stop"
#' )
#' simulation_files(sim)  # NULL: no read-evidence export has been requested.
#'
#' @export
simulation_files <- function(data) {
  attr(data, "buRst_simulation_files", exact = TRUE)
}


# -------------------------------------------------------------------
# Rcpp wrappers
# -------------------------------------------------------------------

#' Simulate newly synthesized RNA counts via Gillespie
#'
#' @param n Number of cells.
#' @param km Burst frequency.
#' @param B Mean burst size.
#' @param d Degradation rate.
#' @param t.end Labeling duration.
#'
#' @return Integer-like numeric vector of length n.
#' @examples
#' set.seed(1)
#' x <- gmRNA_new(n = 100, km = 1, B = 2, d = 0.5, t.end = 1)
#' head(x)
#' mean(x)
#'
#' @export
gmRNA_new <- function(n, km, B, d, t.end) {
  cpp_gmRNA_new(n, km, B, d, t.end)
}


#' Simulate old RNA survival counts via Gillespie
#'
#' @param n Number of cells.
#' @param d Degradation rate.
#' @param t.end Duration.
#' @param init Integer vector of initial old counts at t=0.
#'
#' @return Integer-like numeric vector of length n.
#' @examples
#' set.seed(1)
#' x <- gmRNA_old(n = 100, d = 0.5, t.end = 1, init = rep(4L, 100))
#' range(x)  # Surviving counts cannot exceed their initial counts.
#'
#' @export
gmRNA_old <- function(n, d, t.end, init) {
  cpp_gmRNA_old(n = n, d = d, tend = t.end, init = init)
}


#' Simulate joint total/new/old counts via Gillespie
#'
#' @param n Number of cells.
#' @param km Burst frequency.
#' @param B Mean burst size.
#' @param d Degradation rate.
#' @param t.end Duration.
#' @param init Integer vector of initial old counts at t=0.
#'
#' @return Numeric matrix with columns total, new, and old.
#' @examples
#' set.seed(1)
#' x <- gmRNA_total(n = 100, km = 1, B = 2, d = 0.5,
#'                  t.end = 1, init = rep(4L, 100))
#' head(x)
#' stopifnot(all(x[, 1] == x[, 2] + x[, 3]))
#'
#' @export
gmRNA_total <- function(n, km, B, d, t.end, init) {
  cpp_gmRNA_joint(
    n = n,
    km = km,
    B = B,
    d = d,
    t_end = t.end,
    init = init
  )
}


# -------------------------------------------------------------------
# grandR conversion
# -------------------------------------------------------------------

#' Convert simulated new/old matrices into a grandR object
#'
#' @param new Matrix of genes by cells containing new counts.
#' @param old Matrix of genes by cells containing old counts.
#' @param df_vals Data frame with columns kms, Bs, and ds.
#' @param timepoints One labeling duration or one value per cell.
#' @param condition Condition label.
#' @param gene.info Data frame with Gene/Symbol columns, or a character vector.
#' @param tp.label Timepoint label used in truth-column names.
#' @param prefix Prefix stored in the grandR object.
#'
#' @return A grandR object.
#' @details
#' Rows are genes and columns are cells. The corresponding row of `df_vals`
#' contains that gene's true parameters. The returned `count` slot is `new + old`;
#' the `ntr` slot is the new-to-total ratio, with zero for zero total counts.
#' True parameter columns use names such as `true_B.t0`.
#'
#' @seealso `simulateBursting()`, `write_10x_reference()`
#'
#' @examples
#' new <- matrix(c(0, 1, 2, 1, 0, 1), nrow = 2, byrow = TRUE)
#' old <- matrix(c(2, 1, 0, 3, 2, 1), nrow = 2, byrow = TRUE)
#' truth <- data.frame(kms = c(1, 2), Bs = c(2, 1), ds = c(0.5, 0.5))
#' sim <- simulation_to_grandR(new, old, truth,
#'                             gene.info = c("gene_A", "gene_B"))
#' grandR::GetMatrix(sim, mode.slot = "count")
#'
#' @export
simulation_to_grandR <- function(
    new,
    old,
    df_vals,
    timepoints = 2,
    condition = "sim",
    gene.info = paste("gene", seq_len(nrow(new))),
    tp.label = 0,
    prefix = "Simulated"
) {
  if (!requireNamespace("grandR", quietly = TRUE)) stop("Need grandR.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Need Matrix.")

  .validate_count_matrices(new, old)

  n_genes <- nrow(new)
  n_cells <- ncol(new)

  if (!is.data.frame(df_vals) || !all(c("kms", "Bs", "ds") %in% names(df_vals))) {
    stop("df_vals must be a data.frame containing kms, Bs, and ds.")
  }

  if (nrow(df_vals) != n_genes) {
    stop("df_vals must contain exactly one row per gene.")
  }

  gene.info <- .normalize_gene_info(gene.info, n_genes)

  if (length(timepoints) == 1L) {
    timepoints <- rep(as.numeric(timepoints), n_cells)
  } else if (length(timepoints) == n_cells) {
    timepoints <- as.numeric(timepoints)
  } else {
    stop("timepoints must have length 1 or ncol(new).")
  }

  barcodes <- .make_simulation_barcodes(n_cells)
  time_labels <- ifelse(timepoints == 0, "no4sU", paste0(timepoints, "h"))

  sample_names <- paste(
    condition,
    time_labels,
    "A",
    tp.label,
    barcodes,
    sep = "."
  )

  coldata <- grandR::MakeColdata(
    sample_names,
    design = c(
      grandR::Design$Condition,
      grandR::Design$dur.4sU,
      grandR::Design$Replicate,
      grandR::Design$hpi,
      grandR::Design$Barcode
    ),
    rownames = TRUE
  )
  coldata$no4sU <- timepoints == 0

  new <- Matrix::Matrix(new, sparse = TRUE)
  old <- Matrix::Matrix(old, sparse = TRUE)
  count <- new + old

  # New counts are non-zero only where total counts are positive, so the
  # NTR can be calculated on the sparse non-zero entries without creating a
  # dense genes-by-cells matrix. Use dgCMatrix slots directly: depending on
  # the Matrix version, summary(ntr) may be a matrix-like object without
  # `$i`/`$j` accessors.
  ntr <- methods::as(new, "dgCMatrix")
  count_c <- methods::as(count, "dgCMatrix")

  if (length(ntr@x) > 0L) {
    row_index <- ntr@i + 1L
    col_index <- rep.int(
      seq_len(ncol(ntr)),
      diff(ntr@p)
    )

    denominators <- as.numeric(
      count_c[cbind(row_index, col_index)]
    )

    if (length(denominators) != length(ntr@x) ||
        any(!is.finite(denominators)) ||
        any(denominators <= 0)) {
      stop(
        "Could not calculate sparse NTR values: invalid total-count ",
        "denominators at non-zero new-count positions."
      )
    }

    ntr@x <- ntr@x / denominators
  }

  data_slots <- list(
    count = count,
    ntr = ntr,
    new = new,
    old = old
  )

  data_slots <- lapply(data_slots, function(x) {
    rownames(x) <- gene.info$Gene
    colnames(x) <- sample_names
    x
  })

  gene.info[[paste0("true_B.t", tp.label)]] <- as.numeric(df_vals$Bs)
  gene.info[[paste0("true_km.t", tp.label)]] <- as.numeric(df_vals$kms)
  gene.info[[paste0("true_d.t", tp.label)]] <- as.numeric(df_vals$ds)

  out <- grandR::grandR(
    prefix = prefix,
    gene.info = gene.info,
    slots = data_slots,
    coldata = coldata,
    metadata = list(Description = "Simulated bursting data")
  )

  grandR::DefaultSlot(out) <- "count"
  out
}


# -------------------------------------------------------------------
# Distribution simulators
# -------------------------------------------------------------------

#' Simulate old RNA counts after a labeling interval
#'
#' @param t.end Duration.
#' @param km Burst frequency used for steady-state initialization.
#' @param d Degradation rate.
#' @param B Mean burst size used for steady-state initialization.
#' @param n Number of cells.
#' @param steady Whether to initialize from steady state.
#' @param init Initial old counts when \code{steady=FALSE}.
#' @param seed Optional seed.
#'
#' @return Numeric vector of length n.
#' @examples
#' x <- simulate_distr_old(t.end = 1, d = 0.5, n = 100,
#'                         steady = FALSE, init = rep(4L, 100), seed = 1)
#' mean(x)
#'
#' @export
simulate_distr_old <- function(
    t.end = 2,
    km = 3,
    d = log(2) / 4,
    B = 10,
    n = 10000,
    steady = TRUE,
    init = NULL,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  if (isTRUE(steady)) {
    init_df <- compute.distr.new(
      B = B,
      d = d,
      km = km,
      t.end = Inf,
      type = "steady"
    )
    init <- sample(init_df$k, n, prob = init_df$prob, replace = TRUE)
  } else {
    if (is.null(init)) stop("steady=FALSE requires init.")
    if (length(init) != n) stop("init must have length n.")
  }

  gmRNA_old(n = n, d = d, t.end = t.end, init = init)
}


#' Simulate new RNA counts after a labeling interval
#'
#' @param t.end Duration.
#' @param km Burst frequency.
#' @param d Degradation rate.
#' @param B Mean burst size.
#' @param n Number of cells.
#' @param seed Optional seed.
#'
#' @return Numeric vector of length n.
#' @examples
#' x <- simulate_distr_new(t.end = 1, km = 1, B = 2, d = 0.5,
#'                         n = 100, seed = 1)
#' summary(x)
#'
#' @export
simulate_distr_new <- function(
    t.end = 2,
    km = 3,
    d = log(2) / 4,
    B = 10,
    n = 10000,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  gmRNA_new(n = n, km = km, B = B, d = d, t.end = t.end)
}


#' Simulate total, new, and old counts
#'
#' @param km Burst frequency.
#' @param d Degradation rate.
#' @param B Mean burst size.
#' @param n Number of cells.
#' @param t.end Duration.
#' @param steady Whether to initialize from steady state.
#' @param init Initial old counts when \code{steady=FALSE}.
#' @param seed Optional seed.
#'
#' @return A list with total, new, and old vectors.
#' @examples
#' x <- simulate_distr_total(km = 1, B = 2, d = 0.5,
#'                           n = 100, t.end = 1, seed = 1)
#' head(data.frame(x))
#' stopifnot(all(x$total == x$new + x$old))
#'
#' @export
simulate_distr_total <- function(
    km = 3,
    d = log(2) / 4,
    B = 10,
    n = 10000,
    t.end,
    steady = TRUE,
    init = NULL,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  if (isTRUE(steady)) {
    init_df <- compute.distr.new(
      B = B,
      d = d,
      km = km,
      t.end = Inf,
      type = "steady"
    )
    init <- sample(init_df$k, n, prob = init_df$prob, replace = TRUE)
  } else {
    if (is.null(init)) stop("steady=FALSE requires init.")
    if (length(init) != n) stop("init must have length n.")
  }

  sim <- gmRNA_total(
    n = n,
    km = km,
    B = B,
    d = d,
    t.end = t.end,
    init = init
  )

  list(
    total = sim[, 1L],
    new = sim[, 2L],
    old = sim[, 3L]
  )
}


# -------------------------------------------------------------------
# Main count simulator
# -------------------------------------------------------------------

#' Simulate bursting counts and optionally simulate read evidence
#'
#' @param kms,ds,Bs Numeric vectors of per-gene kinetic parameters.
#' @param t.end Labeling duration.
#' @param no.cells Number of cells.
#' @param ntr.sim Whether to call \code{SimulateReads()}.
#' @param steady Whether to initialize at steady state.
#' @param init Initial old counts for a non-steady-state simulation.
#' @param condition Condition label.
#' @param gene.info Gene identifiers or Gene/Symbol data frame.
#' @param tp.label Truth-column suffix.
#' @param seed Optional seed.
#' @param burstMCMC Whether to write burstMCMC evidence.
#' @param fit.beta Whether to fit a Beta approximation for each cell and store
#'   alpha and beta slots.
#' @param parallel Whether to parallelize across genes.
#' @param ncores Number of workers.
#' @param on_gene_error Either \code{"zero"} to preserve the old fallback or
#'   \code{"stop"} to stop on the first failed gene.
#' @param prefix Prefix stored in the grandR object.
#' @param ... Passed to \code{SimulateReads()} when \code{ntr.sim=TRUE}.
#'
#' @return A grandR object.
#' @details
#' Vectors `kms`, `Bs` and `ds` must have the same length and order as the genes.
#' `kms` and `ds` are rates in the reciprocal time unit used for `t.end`.
#' When `steady = FALSE`, supply a genes-by-cells matrix or one count vector per
#' gene as `init`; these are the molecules present at the start of the interval.
#'
#' With `ntr.sim = FALSE`, the result contains simulated `count`, `new`, `old`
#' and `ntr` slots without creating labeling-evidence files. With `ntr.sim = TRUE`,
#' read-evidence settings and export paths can be passed through `...` to
#' `SimulateReads()`. Use `simulation_files()` to retrieve generated paths.
#' `on_gene_error = "stop"` is useful when validating an example, to avoid
#' replacing a failed simulation by zeros.
#'
#' @seealso `SimulateReads()`, `simulation_files()`, `simulate_distr_total()`
#'
#' @examples
#' # Count simulation only: no HDF5 files or read evidence are written.
#' sim <- simulateBursting(
#'   kms = c(1, 2), Bs = c(2, 1), ds = c(0.5, 0.5),
#'   gene.info = c("gene_A", "gene_B"), t.end = 1, no.cells = 40,
#'   condition = "control", ntr.sim = FALSE, parallel = FALSE,
#'   seed = 1, on_gene_error = "stop"
#' )
#' head(grandR::Coldata(sim))
#' grandR::GetMatrix(sim, mode.slot = "count")[, 1:4]
#'
#' @export
simulateBursting <- function(
    kms,
    ds,
    Bs,
    t.end,
    no.cells,
    ntr.sim = TRUE,
    steady = TRUE,
    init = NULL,
    condition = "Simulation",
    gene.info = paste("gene", seq_along(kms)),
    tp.label = 0,
    seed = NULL,
    burstMCMC = FALSE,
    fit.beta = FALSE,
    parallel = FALSE,
    ncores = NULL,
    ...,
    on_gene_error = c("zero", "stop"),
    prefix = "Simulated"
) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Need Matrix.")
  if (!requireNamespace("grandR", quietly = TRUE)) stop("Need grandR.")
  if (!requireNamespace("parallel", quietly = TRUE)) stop("Need parallel.")

  on_gene_error <- match.arg(on_gene_error)

  kms <- as.numeric(kms)
  ds <- as.numeric(ds)
  Bs <- as.numeric(Bs)

  G <- length(kms)

  if (G < 1L) stop("At least one gene must be simulated.")
  if (length(ds) != G || length(Bs) != G) {
    stop("kms, ds, and Bs must have identical lengths.")
  }
  if (any(!is.finite(kms)) || any(!is.finite(ds)) || any(!is.finite(Bs))) {
    stop("kms, ds, and Bs must contain only finite values.")
  }
  if (any(kms <= 0) || any(ds <= 0) || any(Bs <= 0)) {
    stop("kms, ds, and Bs must all be positive.")
  }

  no.cells <- as.integer(no.cells[1L])
  if (!is.finite(no.cells) || no.cells < 1L) {
    stop("no.cells must be a positive integer.")
  }

  if (!is.null(seed)) {
    old_kind <- RNGkind()
    on.exit(RNGkind(old_kind[1L], old_kind[2L], old_kind[3L]), add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
  }

  if (!isTRUE(steady)) {
    if (is.null(init)) stop("steady=FALSE requires init.")

    if (is.matrix(init) || inherits(init, "Matrix")) {
      if (nrow(init) != G || ncol(init) != no.cells) {
        stop("init matrix must be genes x cells.")
      }
    } else if (is.list(init)) {
      if (length(init) != G) {
        stop("init list must have one element per gene.")
      }
      if (!all(vapply(init, length, integer(1L)) == no.cells)) {
        stop("Every init[[g]] must have length no.cells.")
      }
    } else {
      stop("init must be a genes-by-cells matrix or a list with one vector per gene.")
    }
  }

  worker <- function(i) {
    init_i <- NULL
    if (!isTRUE(steady)) {
      init_i <- if (is.matrix(init) || inherits(init, "Matrix")) {
        as.numeric(init[i, ])
      } else {
        as.numeric(init[[i]])
      }
    }

    tryCatch(
      simulate_distr_total(
        B = Bs[i],
        d = ds[i],
        km = kms[i],
        t.end = t.end,
        n = no.cells,
        steady = steady,
        init = init_i
      ),
      error = function(e) {
        if (identical(on_gene_error, "stop")) {
          stop(
            sprintf("simulate_distr_total failed for gene %d: %s", i, conditionMessage(e)),
            call. = FALSE
          )
        }

        warning(
          sprintf(
            "simulate_distr_total failed for gene %d; replacing it with zeros: %s",
            i,
            conditionMessage(e)
          ),
          call. = FALSE
        )

        list(
          new = rep(0, no.cells),
          old = rep(0, no.cells),
          total = rep(0, no.cells)
        )
      }
    )
  }

  idx <- seq_len(G)
  ncores_use <- .resolve_ncores(ncores)

  if (!isTRUE(parallel) || G <= 1L || ncores_use <= 1L) {
    res <- lapply(idx, worker)
  } else if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(ncores_use)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    if (!is.null(seed)) {
      parallel::clusterSetRNGStream(cl, iseed = seed)
    }

    parallel::clusterExport(
      cl,
      varlist = c(
        "steady", "init", "no.cells", "t.end", "kms", "ds", "Bs",
        "simulate_distr_total", "worker", "on_gene_error"
      ),
      envir = environment()
    )

    res <- parallel::parLapply(cl, idx, worker)
  } else {
    res <- parallel::mclapply(
      idx,
      worker,
      mc.cores = ncores_use,
      mc.set.seed = !is.null(seed)
    )
  }

  old <- do.call(rbind, lapply(res, `[[`, "old"))
  new <- do.call(rbind, lapply(res, `[[`, "new"))

  df_vals <- data.frame(
    kms = kms,
    Bs = Bs,
    ds = ds,
    stringsAsFactors = FALSE
  )

  data <- simulation_to_grandR(
    new = Matrix::Matrix(new, sparse = TRUE),
    old = Matrix::Matrix(old, sparse = TRUE),
    df_vals = df_vals,
    timepoints = t.end,
    condition = condition,
    gene.info = gene.info,
    tp.label = tp.label,
    prefix = prefix
  )

  if (isTRUE(ntr.sim)) {
    data <- SimulateReads(
      data = data,
      burstMCMC = burstMCMC,
      fit.beta = fit.beta,
      seed = seed,
      parallel = parallel,
      ncores = ncores_use,
      ...
    )
  }

  grandR::Coldata(data, "steady") <- rep(
    isTRUE(steady),
    nrow(grandR::Coldata(data))
  )

  data
}


# -------------------------------------------------------------------
# Read simulation and bundled output
# -------------------------------------------------------------------

#' Simulate GRAND-SLAM-like evidence and write discrete profiles
#'
#' When \code{output_prefix} or \code{targets_dir} is supplied, a self-contained
#' target directory is created. For example, \code{output_prefix="./out/sim"}
#' creates \code{./out/sim.targets/} containing the HDF5 file and the standard
#' 10x matrix, barcode, and feature files.
#'
#' @param data A grandR object with count and ntr slots.
#' @param u.content Mean uridine fraction per read.
#' @param u.content.sd Standard deviation of uridine fraction.
#' @param read.length Simulated read length.
#' @param p.old Background conversion probability.
#' @param p.new Conversion probability for new RNA.
#' @param p.new.fit Conversion probability used during fitting.
#' @param seed Optional seed.
#' @param R_clip Quantization range.
#' @param burstMCMC Whether to write burstMCMC evidence.
#' @param fit.beta Whether to call \code{fit.ntr(..., beta.approx=TRUE)} and
#'   store fitted alpha and beta slots.
#' @param h5_path Legacy/direct HDF5 path. When bundled output is requested,
#'   only its basename is used inside the target directory.
#' @param overwrite_h5 Whether an existing HDF5 file may be replaced.
#' @param compression_level HDF5 compression level.
#' @param chunk_target Approximate HDF5 chunk target.
#' @param burstMC_file Optional burstMCMC output path.
#' @param parallel Whether to parallelize gene computation.
#' @param ncores Number of workers.
#' @param verbose_every Progress interval in sequential mode.
#' @param output_prefix Prefix used to create \code{<prefix>.targets}.
#' @param targets_dir Explicit target-directory path; overrides output_prefix.
#' @param h5_filename HDF5 filename used for bundled output.
#' @param write_10x Whether to write matrix.mtx.gz, barcodes.tsv.gz, and
#'   features.tsv.gz.
#' @param overwrite_10x Whether existing 10x files may be replaced.
#'
#' @return The input grandR object with ntr_true, fitted ntr, and an attached
#'   output manifest retrievable with \code{simulation_files()}.
#' @details
#' Use bundled output for later HDF5 fitting: `output_prefix` or `targets_dir`
#' keeps the discrete profiles and their matching count/reference files together.
#' Fitting inputs must retain the same gene and cell identifiers across files.
#' `fit.beta = TRUE` additionally creates the `alpha` and `beta` slots used by
#' `fit_gene_beta()`. All output paths in the example are temporary.
#'
#' @seealso `simulateBursting()`, `simulation_files()`, `write_10x_reference()`
#'
#' @examples
#' sim <- simulateBursting(
#'   kms = 1, Bs = 2, ds = 0.5, t.end = 1, no.cells = 8,
#'   gene.info = "gene_A", condition = "sim", ntr.sim = FALSE,
#'   parallel = FALSE, seed = 1, on_gene_error = "stop"
#' )
#' work <- tempfile("buRst_reads_")
#' dir.create(work)
#' with_reads <- SimulateReads(
#'   sim, output_prefix = file.path(work, "sim"),
#'   parallel = FALSE, seed = 1, verbose_every = 0
#' )
#' simulation_files(with_reads)
#' unlink(work, recursive = TRUE)
#'
#' @export
SimulateReads <- function(
    data,
    u.content = 0.25,
    u.content.sd = 0.05,
    read.length = 100,
    p.old = 1e-4,
    p.new = 0.04,
    p.new.fit = p.new,
    seed = NULL,
    R_clip = 60,
    burstMCMC = FALSE,
    fit.beta = FALSE,
    h5_path = "./sim_profiles.h5",
    overwrite_h5 = TRUE,
    compression_level = 4,
    chunk_target = 1e6,
    burstMC_file = paste0(data$prefix, ".txt"),
    parallel = FALSE,
    ncores = max(1L, parallel::detectCores() - 1L),
    verbose_every = 50,
    output_prefix = NULL,
    targets_dir = NULL,
    h5_filename = "4sU.Binom.discrete.h5",
    write_10x = !is.null(output_prefix) || !is.null(targets_dir),
    overwrite_10x = TRUE
) {
  if (!requireNamespace("grandR", quietly = TRUE)) stop("Need grandR.")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Need Matrix.")
  if (!requireNamespace("rhdf5", quietly = TRUE)) stop("Need rhdf5.")
  if (!requireNamespace("parallel", quietly = TRUE)) stop("Need parallel.")

  if (!is.finite(u.content) || u.content <= 0 || u.content >= 1) {
    stop("u.content must lie strictly between 0 and 1.")
  }
  if (!is.finite(u.content.sd) || u.content.sd <= 0) {
    stop("u.content.sd must be positive.")
  }
  if (u.content.sd^2 >= u.content * (1 - u.content)) {
    stop("u.content.sd is too large for the requested beta distribution.")
  }
  if (!is.finite(read.length) || read.length < 1) {
    stop("read.length must be positive.")
  }

  read.length <- as.integer(read.length)
  ncores <- .resolve_ncores(ncores)

  h5_path_was_supplied <- !missing(h5_path)
  h5_path_input <- if (
    !h5_path_was_supplied &&
    (!is.null(output_prefix) || !is.null(targets_dir) || isTRUE(write_10x))
  ) {
    NULL
  } else {
    h5_path
  }

  paths <- .resolve_targets_paths(
    data = data,
    h5_path = h5_path_input,
    output_prefix = output_prefix,
    targets_dir = targets_dir,
    h5_filename = h5_filename,
    write_10x = write_10x
  )

  if (isTRUE(paths$bundled)) {
    dir.create(paths$targets_dir, recursive = TRUE, showWarnings = FALSE)
  } else {
    dir.create(dirname(paths$h5_path), recursive = TRUE, showWarnings = FALSE)
  }

  h5_path <- paths$h5_path

  ntrs <- grandR::GetMatrix(data, mode.slot = "ntr")
  counts <- grandR::GetMatrix(data, mode.slot = "count")

  G <- nrow(counts)
  C <- ncol(counts)

  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("The count matrix must have gene and cell names.")
  }

  if (isTRUE(overwrite_h5) && file.exists(h5_path)) {
    unlink(h5_path)
  } else if (file.exists(h5_path)) {
    stop("HDF5 output already exists and overwrite_h5=FALSE: ", h5_path)
  }

  if (!rhdf5::h5createFile(h5_path)) {
    stop("Could not create HDF5 file: ", h5_path)
  }

  rhdf5::h5createGroup(h5_path, "genes")
  rhdf5::h5createGroup(h5_path, "genes_meta")
  rhdf5::h5createGroup(h5_path, "meta")
  rhdf5::h5write(colnames(counts), h5_path, "/meta/cell_names")
  rhdf5::h5write(rownames(counts), h5_path, "/meta/gene_names")

  shape1 <- u.content * (
    u.content * (1 - u.content) / u.content.sd^2 - 1
  )
  shape2 <- (1 - u.content) / u.content * shape1

  if (!is.null(seed)) {
    old_kind <- RNGkind()
    on.exit(RNGkind(old_kind[1L], old_kind[2L], old_kind[3L]), add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
  }

  gene_compute <- function(g_idx) {
    gene <- rownames(counts)[g_idx]
    counts_gene <- as.integer(counts[g_idx, ])
    ntr_gene <- as.numeric(ntrs[g_idx, ])

    q_list_gene <- vector("list", C)
    ntr_hat_gene <- rep(0, C)
    alpha_gene <- rep(1, C)
    beta_gene <- rep(1, C)
    burstMC_gene <- if (isTRUE(burstMCMC)) character(0) else NULL

    for (c_idx in seq_len(C)) {
      reads <- counts_gene[c_idx]
      ntr0 <- ntr_gene[c_idx]
      barcode <- colnames(counts)[c_idx]

      if (is.na(ntr0)) {
        q_list_gene[[c_idx]] <- integer(0)

        if (isTRUE(burstMCMC)) {
          burstMC_gene <- c(
            burstMC_gene,
            paste(c(gene, barcode, reads), collapse = "\t")
          )
        }
        next
      }

      if (reads <= 0L) {
        ntr_hat_gene[c_idx] <- 0
        q_list_gene[[c_idx]] <- integer(0)

        if (isTRUE(burstMCMC)) {
          burstMC_gene <- c(
            burstMC_gene,
            paste(c(gene, barcode, reads), collapse = "\t")
          )
        }
        next
      }

      u_per_read <- stats::rbinom(
        reads,
        size = read.length,
        prob = stats::rbeta(reads, shape1, shape2)
      )

      positive_u <- u_per_read[u_per_read > 0L]

      if (length(positive_u) == 0L) {
        # No informative U positions. Keep the cell uninformative rather than
        # constructing an invalid empty mix matrix.
        ntr_hat_gene[c_idx] <- 0
        q_list_gene[[c_idx]] <- rep.int(0L, reads + 1L)

        if (isTRUE(burstMCMC)) {
          burstMC_gene <- c(
            burstMC_gene,
            paste(c(gene, barcode, reads), collapse = "\t")
          )
        }
        next
      }

      us <- table(positive_u)
      u_names <- as.integer(names(us))
      u_histo <- integer(max(u_names))
      u_histo[u_names] <- as.integer(us)

      para_true <- .burst_model_par(
        ntr = ntr0,
        p.err = p.old,
        p.conv = p.new
      )
      
      mixmat <- .burst_create_mix_matrix(
        n.vector = u_histo,
        par = para_true
      )

      para_fit <- .burst_model_par(
        ntr = ntr0,
        p.err = p.old,
        p.conv = p.new.fit
      )
      if (isTRUE(fit.beta)) {
        
        ntr_fit <- .burst_fit_ntr(
          mixmat = mixmat,
          par = para_fit,
          beta.approx = TRUE
        )
        
        ntr_beta <- .extract_ntr_beta_fit(
          ntr_fit,
          gene = gene,
          barcode = barcode
        )
        
      } else {
        
        ntr_fit <- .burst_fit_ntr(
          mixmat = mixmat,
          par = para_fit,
          beta.approx = FALSE
        )
        
        ntr_beta <- list(
          ntr = as.numeric(ntr_fit)[1L],
          alpha = NA_real_,
          beta = NA_real_
        )
      }

      N <- as.integer(reads)
      cell <- prep_cell_ll01(mixmat, para_fit)
      ell_exact <- mode_from_exact_discrete(
        N,
        cell$ll0,
        cell$ll1,
        cell$c
      )$ell_exact
      q <- quantize_ell_u16_max0(ell_exact, R_clip = R_clip)

      ntr_hat_gene[c_idx] <- ntr_beta$ntr

      if (isTRUE(fit.beta)) {
        alpha_gene[c_idx] <- ntr_beta$alpha
        beta_gene[c_idx] <- ntr_beta$beta
      }
      q_list_gene[[c_idx]] <- q

      if (isTRUE(burstMCMC)) {
        evidence <- character(0)
        mm <- as.matrix(unclass(mixmat))
        ks <- .burst_get_mixmat_k(mixmat)
        ns <- .burst_get_mixmat_n(mixmat)

        for (m in seq_along(ns)) {
          for (b in seq_along(ks)) {
            count <- mm[b, m]
            if (count > 0) {
              evidence <- c(
                evidence,
                rep(paste0(ks[b], "/", ns[m]), count)
              )
            }
          }
        }

        burstMC_gene <- c(
          burstMC_gene,
          paste(c(gene, barcode, reads, evidence), collapse = "\t")
        )
      }
    }

    list(
      g_idx = g_idx,
      gene = gene,
      counts_gene = counts_gene,
      q_list_gene = q_list_gene,
      ntr_hat_gene = ntr_hat_gene,
      alpha_gene = alpha_gene,
      beta_gene = beta_gene,
      burstMC_gene = burstMC_gene
    )
  }

  idx <- seq_len(G)
  use_parallel <- isTRUE(parallel) &&
    .Platform$OS.type != "windows" &&
    ncores > 1L &&
    G > 1L

  results <- if (use_parallel) {
    parallel::mclapply(
      idx,
      gene_compute,
      mc.cores = ncores,
      mc.set.seed = !is.null(seed)
    )
  } else {
    lapply(idx, function(i) {
      if (verbose_every > 0L && (i %% verbose_every == 0L || i == G)) {
        message(sprintf("SimulateReads: %d / %d genes", i, G))
      }
      gene_compute(i)
    })
  }

  # Validate and write one result at a time. Use [[ ]] rather than `$`
  # so this remains safe for list-like return objects and avoids S4 index
  # evaluation errors.
  for (result_index in seq_along(results)) {
    result <- results[[result_index]]

    gene_result <- .sim_result_field(result, "gene", result_index)
    q_result <- .sim_result_field(result, "q_list_gene", result_index)
    counts_result <- .sim_result_field(result, "counts_gene", result_index)

    .simulate_reads_stage(
      write_one_gene_profile_from_disc(
        h5_path = h5_path,
        gene = gene_result,
        q_list_gene = q_result,
        counts_gene = counts_result,
        compression_level = compression_level,
        chunk_target = chunk_target
      ),
      paste0("writing the HDF5 profile for gene '", gene_result, "'")
    )
  }
  ntr_hat_mat <- matrix(
    0,
    nrow = G,
    ncol = C,
    dimnames = list(rownames(counts), colnames(counts))
  )

  alpha_mat <- matrix(
    1,
    nrow = G,
    ncol = C,
    dimnames = list(rownames(counts), colnames(counts))
  )

  beta_mat <- matrix(
    1,
    nrow = G,
    ncol = C,
    dimnames = list(rownames(counts), colnames(counts))
  )

  .simulate_reads_stage({
    for (result_index in seq_along(results)) {
      result <- results[[result_index]]

      g_idx_result <- as.integer(
        .sim_result_field(
          result,
          "g_idx",
          result_index
        )
      )[1L]

      ntr_result <- as.numeric(
        .sim_result_field(
          result,
          "ntr_hat_gene",
          result_index
        )
      )

      if (!is.finite(g_idx_result) ||
          g_idx_result < 1L ||
          g_idx_result > G) {
        stop(
          "Invalid g_idx in result ",
          result_index,
          ": ",
          g_idx_result,
          call. = FALSE
        )
      }

      if (length(ntr_result) != C) {
        stop(
          "Invalid ntr_hat_gene length in result ",
          result_index,
          ": expected ",
          C,
          ", got ",
          length(ntr_result),
          call. = FALSE
        )
      }

      ntr_hat_mat[g_idx_result, ] <- ntr_result

      if (isTRUE(fit.beta)) {
        alpha_result <- as.numeric(
          .sim_result_field(
            result,
            "alpha_gene",
            result_index
          )
        )

        beta_result <- as.numeric(
          .sim_result_field(
            result,
            "beta_gene",
            result_index
          )
        )

        if (length(alpha_result) != C) {
          stop(
            "Invalid alpha_gene length in result ",
            result_index,
            ": expected ",
            C,
            ", got ",
            length(alpha_result),
            call. = FALSE
          )
        }

        if (length(beta_result) != C) {
          stop(
            "Invalid beta_gene length in result ",
            result_index,
            ": expected ",
            C,
            ", got ",
            length(beta_result),
            call. = FALSE
          )
        }

        alpha_mat[g_idx_result, ] <- alpha_result
        beta_mat[g_idx_result, ] <- beta_result
      }
    }
  }, "collecting fitted NTR and beta approximation values")
  ntr_sim <- .simulate_reads_stage(
    Matrix::Matrix(ntr_hat_mat, sparse = TRUE),
    "converting fitted NTR values to a sparse matrix"
  )

  alpha_sim <- NULL
  beta_sim <- NULL

  if (isTRUE(fit.beta)) {
    alpha_sim <- .simulate_reads_stage(
      Matrix::Matrix(alpha_mat, sparse = TRUE),
      "converting fitted alpha values to a sparse matrix"
    )

    beta_sim <- .simulate_reads_stage(
      Matrix::Matrix(beta_mat, sparse = TRUE),
      "converting fitted beta values to a sparse matrix"
    )
  }

  if (isTRUE(burstMCMC)) {
    if (missing(burstMC_file) && isTRUE(paths$bundled)) {
      burstMC_file <- file.path(paths$targets_dir, "burstMCMC.txt")
    }

    burst_lines <- unlist(
      lapply(
        seq_along(results),
        function(i) .sim_result_field(results[[i]], "burstMC_gene", i)
      ),
      use.names = FALSE
    )
    burst_lines <- burst_lines[
      !is.na(burst_lines) & nzchar(burst_lines)
    ]

    if (length(burst_lines) > 0L) {
      dir.create(dirname(burstMC_file), recursive = TRUE, showWarnings = FALSE)
      if (isTRUE(overwrite_h5) && file.exists(burstMC_file)) {
        unlink(burstMC_file)
      }
      writeLines(burst_lines, con = burstMC_file, sep = "\n")
    }
  } else {
    burstMC_file <- NULL
  }

  data <- .simulate_reads_stage(
    grandR::AddSlot(
      data,
      "ntr_true",
      as.matrix(grandR::GetMatrix(data, mode.slot = "ntr"))
    ),
    "storing the true NTR slot"
  )
  data <- .simulate_reads_stage(
    suppressWarnings(grandR::AddSlot(data, "ntr", ntr_sim)),
    "replacing the fitted NTR slot"
  )
  if (isTRUE(fit.beta)) {
    data <- .simulate_reads_stage(
      suppressWarnings(grandR::AddSlot(data, "alpha", alpha_sim)),
      "storing fitted beta-posterior alpha values"
    )
    data <- .simulate_reads_stage(
      suppressWarnings(grandR::AddSlot(data, "beta", beta_sim)),
      "storing fitted beta-posterior beta values"
    )
  }

  tenx_paths <- NULL
  if (isTRUE(write_10x)) {
    tenx_paths <- .simulate_reads_stage(
      write_10x_reference(
        data = data,
        output_dir = paths$targets_dir,
        slot = "count",
        overwrite = overwrite_10x
      ),
      "writing the 10x matrix, barcodes, and features"
    )
  }

  manifest <- list(
    targets_dir = paths$targets_dir,
    h5 = normalizePath(h5_path, winslash = "/", mustWork = FALSE),
    matrix = if (!is.null(tenx_paths)) {
      normalizePath(tenx_paths$matrix, winslash = "/", mustWork = FALSE)
    } else {
      NULL
    },
    barcodes = if (!is.null(tenx_paths)) {
      normalizePath(tenx_paths$barcodes, winslash = "/", mustWork = FALSE)
    } else {
      NULL
    },
    features = if (!is.null(tenx_paths)) {
      normalizePath(tenx_paths$features, winslash = "/", mustWork = FALSE)
    } else {
      NULL
    },
    burstMCMC = if (!is.null(burstMC_file)) {
      normalizePath(burstMC_file, winslash = "/", mustWork = FALSE)
    } else {
      NULL
    }
  )

  if (!is.null(manifest$targets_dir)) {
    manifest$targets_dir <- normalizePath(
      manifest$targets_dir,
      winslash = "/",
      mustWork = FALSE
    )
  }

  attr(data, "buRst_simulation_files") <- manifest
  data
}
