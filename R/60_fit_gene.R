# fit_gene.R
# Fit bursting parameters for one gene under one condition.

# Helpers -----------------------------------------------------

.require_fit_namespaces <- function(require_hdf5 = FALSE) {
  packages <- c("grandR", "MASS", "numDeriv")
  if (isTRUE(require_hdf5)) packages <- c(packages, "rhdf5")

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

.resolve_gene_id <- function(data, gene) {
  gi <- grandR::GeneInfo(data)
  is_ensg <- grepl("^ENSG\\d+", gene)

  ensg <- if (is_ensg) {
    gene
  } else {
    gi$Gene[match(gene, gi$Symbol)]
  }

  if (length(ensg) != 1L || is.na(ensg) || !nzchar(ensg)) {
    stop("Could not map gene symbol to ENSG: ", gene, call. = FALSE)
  }

  symbol <- gi$Symbol[match(ensg, gi$Gene)]
  if (length(symbol) != 1L || is.na(symbol) || !nzchar(symbol)) {
    symbol <- if (is_ensg) NA_character_ else gene
  }

  list(
    gene_info = gi,
    is_ensg = is_ensg,
    ensg = ensg,
    symbol = symbol
  )
}

.get_gene_vector <- function(data, mode.slot, ensg, columns, integer = FALSE) {
  mat <- grandR::GetMatrix(
    data,
    mode.slot = mode.slot,
    genes = ensg,
    columns = columns
  )

  x <- as.numeric(mat)
  if (isTRUE(integer)) x <- as.integer(x)

  cell_names <- colnames(mat)
  if (is.null(cell_names) || length(cell_names) != length(x)) {
    cell_names <- as.character(columns)
  }
  names(x) <- cell_names

  x
}

.prepare_gene_data <- function(
    data,
    gene,
    fit.cond,
    steady,
    ref.cond,
    t.diff,
    min_cells = 5L
) {
  gene_id <- .resolve_gene_id(data, gene)
  cd <- grandR::Coldata(data)

  if (!"Condition" %in% names(cd)) {
    stop("Coldata(data) must contain a 'Condition' column.", call. = FALSE)
  }

  ok_mask <- rep(TRUE, nrow(cd))
  if ("no4sU" %in% names(cd)) {
    ok_mask <- !as.logical(cd$no4sU)
    ok_mask[is.na(ok_mask)] <- FALSE
  }

  fit_mask <- cd$Condition == fit.cond & ok_mask
  fit_mask[is.na(fit_mask)] <- FALSE
  fit_columns <- grandR::Columns(data, fit_mask, reorder = TRUE)

  if (length(fit_columns) < min_cells) {
    stop(
      "Too few fit cells for gene ", gene,
      " in condition ", fit.cond,
      " (", length(fit_columns), " found; need at least ", min_cells, ").",
      call. = FALSE
    )
  }

  counts_fit <- .get_gene_vector(
    data = data,
    mode.slot = "count",
    ensg = gene_id$ensg,
    columns = fit_columns,
    integer = FALSE
  )

  total_ref <- NULL
  ref_columns <- NULL

  if (!isTRUE(steady)) {
    if (is.null(ref.cond)) {
      stop("steady = FALSE requires ref.cond.", call. = FALSE)
    }
    if (is.null(t.diff)) {
      stop("steady = FALSE requires t.diff.", call. = FALSE)
    }

    ref_mask <- cd$Condition == ref.cond & ok_mask
    ref_mask[is.na(ref_mask)] <- FALSE
    ref_columns <- grandR::Columns(data, ref_mask, reorder = TRUE)

    if (length(ref_columns) < min_cells) {
      stop(
        "Too few reference cells for gene ", gene,
        " in condition ", ref.cond,
        " (", length(ref_columns), " found; need at least ", min_cells, ").",
        call. = FALSE
      )
    }

    total_ref <- .get_gene_vector(
      data = data,
      mode.slot = "count",
      ensg = gene_id$ensg,
      columns = ref_columns,
      integer = TRUE
    )
  }

  c(
    gene_id,
    list(
      fit_columns = fit_columns,
      fit_names = names(counts_fit),
      counts_fit = counts_fit,
      ref_columns = ref_columns,
      total_ref = total_ref
    )
  )
}

.get_initial_new_counts <- function(data, ensg, fit_columns, counts_fit) {
  new_fit <- tryCatch(
    .get_gene_vector(
      data = data,
      mode.slot = "new.count",
      ensg = ensg,
      columns = fit_columns,
      integer = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(new_fit)) {
    ntr_fit <- .get_gene_vector(
      data = data,
      mode.slot = "ntr",
      ensg = ensg,
      columns = fit_columns,
      integer = FALSE
    )
    new_fit <- as.integer(round(ntr_fit * counts_fit))
    names(new_fit) <- names(counts_fit)
  }

  new_fit
}

.information_flag <- function(x, type = c("total", "new")) {
  type <- match.arg(type)
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(if (type == "total") "too_few_total_counts" else "too_few_new_counts")
  }

  total <- sum(x, na.rm = TRUE)
  frac_nonzero <- mean(x > 0, na.rm = TRUE)
  n_distinct <- length(unique(x))

  if (type == "total") {
    if (total < 50) return("too_few_total_counts")
    if (n_distinct <= 3L) return("degenerate_total_distribution")
    if (frac_nonzero < 0.05) return("too_sparse_total_counts")
  } else {
    if (total < 10) return("too_few_new_counts")
    if (n_distinct <= 3L) return("degenerate_new_distribution")
    if (frac_nonzero < 0.03) return("too_sparse_new_counts")
  }

  "ok"
}

.validate_information <- function(counts_fit, new_fit) {
  total_flag <- .information_flag(counts_fit, "total")
  if (!identical(total_flag, "ok")) {
    stop("Insufficient total-count information: ", total_flag, ".", call. = FALSE)
  }

  burst_flag <- .information_flag(new_fit, "new")
  if (!identical(burst_flag, "ok")) {
    stop("Insufficient labeling information: ", burst_flag, ".", call. = FALSE)
  }

  invisible(list(total = total_flag, new = burst_flag))
}

.initial_fit_parameters <- function(
    counts_fit,
    new_fit,
    total_ref,
    t.label,
    steady,
    t.diff
) {
  # The labeling fraction is count-weighted, but the reference comparison
  # uses per-cell means so unequal population sizes do not shift the start.
  ntr_sum <- sum(new_fit, na.rm = TRUE) / pmax(1, sum(counts_fit, na.rm = TRUE))

  d_init_raw <- fit_d(
    ntr = ntr_sum,
    total = mean(counts_fit, na.rm = TRUE),
    t = t.label,
    t0 = if (isTRUE(steady)) NULL else t.diff,
    f0 = if (isTRUE(steady)) NULL else mean(total_ref, na.rm = TRUE)
  )

  d_init <- if (is.finite(d_init_raw) && d_init_raw > 0) {
    d_init_raw
  } else {
    log(2) / t.label
  }

  d_init <- round(d_init, 5L)
  d_init <- pmax(
    log(2) / 50,
    pmin(d_init, log(2) / (5 / 60))
  )

  v_init <- stats::var(new_fit, na.rm = TRUE)
  if (!is.finite(v_init)) v_init <- 0

  moment_init <- initial_params_noss(
    c(
      d = d_init,
      m = mean(new_fit, na.rm = TRUE),
      v = v_init,
      m_0 = 0,
      v_0 = 0
    ),
    t.end = t.label
  )

  B0_try <- unname(moment_init[1])
  km0_try <- unname(moment_init[2])

  bad_moment_init <- (
    !is.finite(B0_try) || !is.finite(km0_try) ||
      B0_try <= 0 || km0_try <= 0 ||
      B0_try < 1e-3 || B0_try > 100 ||
      km0_try < 1e-3 || km0_try > 20 ||
      B0_try * km0_try < 1e-2
  )

  if (bad_moment_init) {
    lambda0 <- max(mean(new_fit, na.rm = TRUE) * d_init, 1e-2)
    B0 <- 0.1
    km0 <- lambda0 / B0
  } else {
    B0 <- round(B0_try, 3L)
    km0 <- round(km0_try, 3L)
  }

  B0 <- min(max(B0, 1e-3), 100)
  km0 <- min(max(km0, 1e-3), 10)
  lambda0 <- B0 * km0

  mv0 <- m_v_total(B = B0, km = km0, d = d_init, t.end = Inf)
  mu0_tot <- as.numeric(mv0[1])
  sd0_tot <- as.numeric(mv0[3])

  if (!is.finite(mu0_tot)) mu0_tot <- max(counts_fit, na.rm = TRUE)
  if (!is.finite(sd0_tot) || sd0_tot < 0) sd0_tot <- 0

  observed_max <- max(c(counts_fit, total_ref), na.rm = TRUE)
  if (!is.finite(observed_max)) observed_max <- 0

  # Keep the conservative support used by the original fit_gene().
  max_observed_support <- if (isTRUE(steady)) {
    observed_max
  } else {
    max(observed_max + 1L, ceiling(mu0_tot + 10 * sd0_tot))
  }

  N0_fit <- max(
    max_observed_support + 1L,
    as.integer(ceiling(mu0_tot + 10 * sd0_tot))
  )
  N0_fit <- min(as.integer(N0_fit), 20000L)

  list(
    d = d_init,
    B = B0,
    km = km0,
    lambda = lambda0,
    N0 = N0_fit,
    par = c(
      log_lambda = log(lambda0),
      logB = log(B0),
      logd = log(d_init)
    )
  )
}

.make_reparameterized_objective <- function(loglik, N0_fit, pass_N0 = TRUE) {
  force(loglik)
  force(N0_fit)
  force(pass_N0)

  function(par) {
    if (length(par) != 3L || any(!is.finite(par))) return(-Inf)

    lambda <- exp(par[1])
    B <- exp(par[2])
    d <- exp(par[3])
    km <- lambda / B

    if (
      any(!is.finite(c(lambda, B, d, km))) ||
      any(c(lambda, B, d, km) <= 1e-12) ||
      any(c(lambda, B, d, km) > 1e3)
    ) {
      return(-Inf)
    }

    value <- if (isTRUE(pass_N0)) {
      loglik(B, km, d, N0_fit)
    } else {
      loglik(B, km, d)
    }

    if (length(value) != 1L || !is.finite(value)) -Inf else value
  }
}

.make_positive_definite <- function(x, eigen_floor = 1e-10) {
  x <- 0.5 * (x + t(x))
  eig <- eigen(x, symmetric = TRUE)
  eig$values[eig$values < eigen_floor] <- eigen_floor
  out <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
  dimnames(out) <- dimnames(x)
  out
}

.fit_likelihood <- function(
    loglik,
    initial,
    ndraw,
    seed,
    pass_N0 = TRUE,
    ridge_cor_threshold = -0.9
) {
  if (!is.numeric(ndraw) || length(ndraw) != 1L || ndraw < 1) {
    stop("ndraw must be a positive integer.", call. = FALSE)
  }
  ndraw <- as.integer(ndraw)

  objective <- .make_reparameterized_objective(
    loglik = loglik,
    N0_fit = initial$N0,
    pass_N0 = pass_N0
  )

  initial_value <- objective(initial$par)
  if (!is.finite(initial_value)) {
    stop("The likelihood is not finite at the initial parameter values.", call. = FALSE)
  }

  fit <- stats::optim(
    par = initial$par,
    fn = objective,
    method = "Nelder-Mead",
    control = list(
      fnscale = -1,
      reltol = 1e-6,
      maxit = 1000
    )
  )

  if (length(fit$par) != 3L || any(!is.finite(fit$par))) {
    stop("Optimization returned invalid parameter estimates.", call. = FALSE)
  }

  hessian <- numDeriv::hessian(objective, fit$par)
  dimnames(hessian) <- list(
    c("log_lambda.B", "logB.B", "logd"),
    c("log_lambda.B", "logB.B", "logd")
  )
  fit$hessian <- hessian

  H <- 0.5 * (hessian + t(hessian))
  Sigma <- tryCatch(solve(-H), error = function(e) NULL)
  if (is.null(Sigma) || any(!is.finite(Sigma))) {
    stop("Could not invert the observed information matrix (bad Hessian).", call. = FALSE)
  }
  dimnames(Sigma) <- dimnames(H)
  Sigma_pd <- .make_positive_definite(Sigma)

  set.seed(seed)
  draws_par <- MASS::mvrnorm(
    n = ndraw,
    mu = fit$par,
    Sigma = Sigma_pd
  )
  if (is.null(dim(draws_par))) {
    draws_par <- matrix(draws_par, nrow = 1L)
  }

  lambda_s <- exp(draws_par[, 1])
  B_s <- exp(draws_par[, 2])
  d_s <- exp(draws_par[, 3])
  km_s <- lambda_s / B_s

  ok <- (
    is.finite(lambda_s) & is.finite(B_s) &
      is.finite(d_s) & is.finite(km_s) &
      lambda_s > 0 & B_s > 0 & d_s > 0 & km_s > 0
  )

  if (!any(ok)) {
    stop("No valid Laplace draws were produced.", call. = FALSE)
  }

  samples <- list(
    B = B_s[ok],
    km = km_s[ok],
    d = d_s[ok],
    lambda = lambda_s[ok]
  )

  q3 <- function(x) {
    stats::quantile(
      x,
      probs = c(0.025, 0.5, 0.975),
      na.rm = TRUE,
      names = FALSE
    )
  }

  CI <- lapply(samples[c("km", "B", "d", "lambda")], q3)

  threshold <- abs(ridge_cor_threshold)
  if (!is.finite(threshold) || threshold <= 0 || threshold >= 1) threshold <- 0.90

  diagnostic <- identifiability_diagnostic(
    fit,
    corr_warn = threshold
  )

  list(
    fit = fit,
    samples = samples,
    CI = CI,
    diagnostic = diagnostic,
    estimates = c(
      lambda = unname(exp(fit$par[1])),
      B = unname(exp(fit$par[2])),
      d = unname(exp(fit$par[3])),
      km = unname(exp(fit$par[1] - fit$par[2]))
    )
  )
}

.log_interval_width <- function(ci) {
  if (length(ci) == 3L && all(is.finite(ci)) && all(ci > 0)) {
    log(ci[3]) - log(ci[1])
  } else {
    NA_real_
  }
}


.collapse_pair <- function(x) {
  if (length(x) != 2L || anyNA(x)) NA_character_ else paste(x, collapse = " -- ")
}

.diagnostic_columns <- function(diagnostic) {
  out <- list(
    id_flag = NA_character_,
    id_ridge_score = NA_real_,
    id_condition_number = NA_real_,
    id_max_corr_fit = NA_real_,
    id_worst_pair_fit = NA_character_,
    id_worst_corr_fit = NA_real_,
    id_max_corr_bio = NA_real_,
    id_worst_pair_bio = NA_character_,
    id_worst_corr_bio = NA_real_,
    id_diag_score = NA_real_
  )

  if (!is.null(diagnostic) && isTRUE(diagnostic$ok)) {
    out$id_flag <- diagnostic$flag
    out$id_ridge_score <- diagnostic$ridge_score
    out$id_condition_number <- diagnostic$condition_number
    out$id_max_corr_fit <- diagnostic$max_corr_fit
    out$id_worst_pair_fit <- .collapse_pair(diagnostic$worst_pair_fit)
    out$id_worst_corr_fit <- diagnostic$worst_corr_fit
    out$id_max_corr_bio <- diagnostic$max_corr_bio
    out$id_worst_pair_bio <- .collapse_pair(diagnostic$worst_pair_bio)
    out$id_worst_corr_bio <- diagnostic$worst_corr_bio
    out$id_diag_score <- diagnostic$diagnostic_score
  }

  out
}

.build_fit_output <- function(
    prepared,
    fit.cond,
    steady,
    ref.cond,
    t.label,
    t.diff,
    fitted
) {
  est <- fitted$estimates
  CI <- fitted$CI
  id <- .diagnostic_columns(fitted$diagnostic)

  data.frame(
    gene = prepared$symbol,
    ENSG = prepared$ensg,
    fit_cond = fit.cond,
    steady = steady,
    ref_cond = if (is.null(ref.cond)) NA_character_ else ref.cond,
    t_label = t.label,
    t_diff = if (is.null(t.diff)) NA_real_ else t.diff,

    km_hat = unname(est["km"]),
    B_hat = unname(est["B"]),
    d_hat = unname(est["d"]),
    lambda_hat = unname(est["lambda"]),

    km_lo = CI$km[1], km_med = CI$km[2], km_hi = CI$km[3],
    B_lo = CI$B[1], B_med = CI$B[2], B_hi = CI$B[3],
    d_lo = CI$d[1], d_med = CI$d[2], d_hi = CI$d[3],
    lambda_lo = CI$lambda[1],
    lambda_med = CI$lambda[2],
    lambda_hi = CI$lambda[3],

    id_flag = id$id_flag,
    id_ridge_score = id$id_ridge_score,
    id_condition_number = id$id_condition_number,
    id_max_corr_fit = id$id_max_corr_fit,
    id_worst_pair_fit = id$id_worst_pair_fit,
    id_worst_corr_fit = id$id_worst_corr_fit,
    id_max_corr_bio = id$id_max_corr_bio,
    id_worst_pair_bio = id$id_worst_pair_bio,
    id_worst_corr_bio = id$id_worst_corr_bio,
    id_diag_score = id$id_diag_score,

    log_width_B = .log_interval_width(CI$B),
    log_width_km = .log_interval_width(CI$km),
    log_width_d = .log_interval_width(CI$d),
    value = fitted$fit$value,
    conv = fitted$fit$convergence,
    stringsAsFactors = FALSE
  )
}

.add_truth_and_pit <- function(
    out,
    gene_info,
    ensg,
    samples,
    truth_suffix,
    add_pit
) {
  out$truth_B <- NA_real_
  out$truth_km <- NA_real_
  out$truth_d <- NA_real_
  out$truth_lambda <- NA_real_
  out$pit_B <- NA_real_
  out$pit_km <- NA_real_
  out$pit_d <- NA_real_
  out$pit_lambda <- NA_real_

  if (!isTRUE(add_pit) || is.null(truth_suffix)) return(out)

  truth_columns <- c(
    B = paste0("true_B.", truth_suffix),
    km = paste0("true_km.", truth_suffix),
    d = paste0("true_d.", truth_suffix)
  )

  if (!all(truth_columns %in% names(gene_info))) return(out)

  row_i <- match(ensg, gene_info$Gene)
  if (is.na(row_i)) return(out)

  truth <- c(
    B = as.numeric(gene_info[row_i, truth_columns["B"]]),
    km = as.numeric(gene_info[row_i, truth_columns["km"]]),
    d = as.numeric(gene_info[row_i, truth_columns["d"]])
  )
  truth <- c(truth, lambda = unname(truth["B"] * truth["km"]))

  out$truth_B <- unname(truth["B"])
  out$truth_km <- unname(truth["km"])
  out$truth_d <- unname(truth["d"])
  out$truth_lambda <- unname(truth["lambda"])

  if (all(is.finite(truth))) {
    pit <- tryCatch(
      pit_from_samples(samples, truth),
      error = function(e) NULL
    )

    if (!is.null(pit)) {
      out$pit_B <- unname(pit[["B"]])
      out$pit_km <- unname(pit[["km"]])
      out$pit_d <- unname(pit[["d"]])
      out$pit_lambda <- unname(pit[["lambda"]])
    }
  }

  out
}

.finalize_fit_output <- function(
    prepared,
    fitted,
    fit.cond,
    steady,
    ref.cond,
    t.label,
    t.diff,
    truth_suffix,
    add_pit
) {
  out <- .build_fit_output(
    prepared = prepared,
    fit.cond = fit.cond,
    steady = steady,
    ref.cond = ref.cond,
    t.label = t.label,
    t.diff = t.diff,
    fitted = fitted
  )

  out <- .add_truth_and_pit(
    out = out,
    gene_info = prepared$gene_info,
    ensg = prepared$ensg,
    samples = fitted$samples,
    truth_suffix = truth_suffix,
    add_pit = add_pit
  )

  attr(out, "fit") <- fitted$fit
  attr(out, "samples") <- fitted$samples
  attr(out, "CI") <- fitted$CI
  out
}

.run_gene_fit <- function(
    prepared,
    loglik,
    new_fit,
    fit.cond,
    t.label,
    steady,
    ref.cond,
    t.diff,
    ndraw,
    seed,
    truth_suffix,
    add_pit,
    ridge_cor_threshold,
    pass_N0,
    check_information = TRUE
) {
  if (isTRUE(check_information)) {
    .validate_information(prepared$counts_fit, new_fit)
  }

  initial <- .initial_fit_parameters(
    counts_fit = prepared$counts_fit,
    new_fit = new_fit,
    total_ref = prepared$total_ref,
    t.label = t.label,
    steady = steady,
    t.diff = t.diff
  )

  fitted <- .fit_likelihood(
    loglik = loglik,
    initial = initial,
    ndraw = ndraw,
    seed = seed,
    pass_N0 = pass_N0,
    ridge_cor_threshold = ridge_cor_threshold
  )

  .finalize_fit_output(
    prepared = prepared,
    fitted = fitted,
    fit.cond = fit.cond,
    steady = steady,
    ref.cond = ref.cond,
    t.label = t.label,
    t.diff = t.diff,
    truth_suffix = truth_suffix,
    add_pit = add_pit
  )
}


# Fitting functions -----------------------------------------------------

#' Fit bursting kinetics from posterior gene profiles
#'
#' Fits one gene in one condition using posterior gene profiles stored in HDF5.
#' The optimization is performed in `(log(lambda), log(B), log(d))` space,
#' where `lambda = B * km`, followed by a Laplace approximation.
#'
#' @param data A grandR object.
#' @param gene Gene symbol or ENSG identifier.
#' @param fit.cond Condition to fit; must match `Coldata(data)$Condition`.
#' @param t.label Labeling duration.
#' @param steady Whether to use the steady-state model.
#' @param ref.cond Reference condition used when `steady = FALSE`.
#' @param t.diff Modeled reference-to-observation interval when `steady = FALSE`.
#'   Must equal `t.label`; the reference represents RNA at labeling onset.
#' @param h5_path Path to the HDF5 gene-profile file.
#' @param reference Optional preloaded 10x reference. When `NULL`, it is loaded
#'   from `dirname(h5_path)`. Batch runners should pass a shared reference.
#' @param ndraw Number of Laplace draws.
#' @param seed Random-number seed used for Laplace draws.
#' @param truth_suffix Optional suffix for truth columns in `GeneInfo(data)`.
#' @param add_pit Whether to append truth and PIT columns when available.
#' @param ridge_cor_threshold Absolute correlation threshold used for the
#'   identifiability warning; the historical negative default is accepted.
#'
#' @return A one-row data frame. The fitted optimizer object, posterior draws,
#'   and intervals are attached as attributes `fit`, `samples`, and `CI`.
#' @export
fit_gene <- function(
    data,
    gene,
    fit.cond,
    t.label = 2,
    steady = TRUE,
    ref.cond = NULL,
    t.diff = NULL,
    h5_path,
    ndraw = 50000,
    seed = 1,
    truth_suffix = "t0",
    add_pit = TRUE,
    ridge_cor_threshold = -0.9,
    reference = NULL
) {
  .require_fit_namespaces(require_hdf5 = TRUE)

  prepared <- .prepare_gene_data(
    data = data,
    gene = gene,
    fit.cond = fit.cond,
    steady = steady,
    ref.cond = ref.cond,
    t.diff = t.diff
  )

  new_fit <- .get_initial_new_counts(
    data = data,
    ensg = prepared$ensg,
    fit_columns = prepared$fit_columns,
    counts_fit = prepared$counts_fit
  )

  if (is.null(reference)) {
    reference <- load_10x_reference(dirname(h5_path))
  }

  gp_fit <- read_gene_profiles_subset(
    h5_path = h5_path,
    gene = prepared$ensg,
    subset_names = prepared$fit_names,
    reference = reference
  )

  loglik <- make_loglik_from_gp(
    gp = gp_fit,
    counts = prepared$counts_fit,
    t.label = t.label,
    steady = steady,
    t.diff = t.diff,
    init = prepared$total_ref
  )

  .run_gene_fit(
    prepared = prepared,
    loglik = loglik,
    new_fit = new_fit,
    fit.cond = fit.cond,
    t.label = t.label,
    steady = steady,
    ref.cond = ref.cond,
    t.diff = t.diff,
    ndraw = ndraw,
    seed = seed,
    truth_suffix = truth_suffix,
    add_pit = add_pit,
    ridge_cor_threshold = ridge_cor_threshold,
    pass_N0 = TRUE,
    check_information = TRUE
  )
}


#' Fit bursting kinetics from MAP new-count estimates
#'
#' Uses `new.count` when available and otherwise reconstructs integer new counts
#' from `ntr * count`. `h5_path` is retained for backward compatibility but is
#' not used by this likelihood.
#'
#' @inheritParams fit_gene
#' @return See [fit_gene()].
#' @export
fit_gene_MAP <- function(
    data,
    gene,
    fit.cond,
    t.label = 2,
    steady = TRUE,
    ref.cond = NULL,
    t.diff = NULL,
    h5_path = NULL,
    ndraw = 50000,
    seed = 1,
    truth_suffix = "t0",
    add_pit = TRUE,
    ridge_cor_threshold = -0.9
) {
  .require_fit_namespaces(require_hdf5 = FALSE)
  invisible(h5_path)

  prepared <- .prepare_gene_data(
    data = data,
    gene = gene,
    fit.cond = fit.cond,
    steady = steady,
    ref.cond = ref.cond,
    t.diff = t.diff
  )

  new_fit <- .get_initial_new_counts(
    data = data,
    ensg = prepared$ensg,
    fit_columns = prepared$fit_columns,
    counts_fit = prepared$counts_fit
  )

  old_fit <- as.integer(pmax(prepared$counts_fit - new_fit, 0))
  names(old_fit) <- names(prepared$counts_fit)

  loglik <- make_loglik_from_truth(
    new_counts = new_fit,
    old_counts = old_fit,
    total_counts = as.integer(prepared$counts_fit),
    t.label = t.label,
    steady = steady,
    t.diff = t.diff,
    init = prepared$total_ref
  )

  .run_gene_fit(
    prepared = prepared,
    loglik = loglik,
    new_fit = new_fit,
    fit.cond = fit.cond,
    t.label = t.label,
    steady = steady,
    ref.cond = ref.cond,
    t.diff = t.diff,
    ndraw = ndraw,
    seed = seed,
    truth_suffix = truth_suffix,
    add_pit = add_pit,
    ridge_cor_threshold = ridge_cor_threshold,
    pass_N0 = TRUE,
    check_information = TRUE
  )
}


#' Fit bursting kinetics from exact new and old counts
#'
#' Intended mainly for simulations or diagnostics where exact new and old RNA
#' counts are available in the `new` and `old` grandR slots.
#'
#' @inheritParams fit_gene
#' @param h5_path Retained for backward compatibility; not used.
#' @return See [fit_gene()].
#' @export
fit_gene_truth <- function(
    data,
    gene,
    fit.cond,
    t.label = 2,
    steady = TRUE,
    ref.cond = NULL,
    t.diff = NULL,
    ndraw = 50000,
    seed = 1,
    h5_path = NULL,
    truth_suffix = "t0",
    add_pit = TRUE,
    ridge_cor_threshold = -0.9
) {
  .require_fit_namespaces(require_hdf5 = FALSE)
  invisible(h5_path)

  prepared <- .prepare_gene_data(
    data = data,
    gene = gene,
    fit.cond = fit.cond,
    steady = steady,
    ref.cond = ref.cond,
    t.diff = t.diff
  )

  new_truth <- .get_gene_vector(
    data = data,
    mode.slot = "new",
    ensg = prepared$ensg,
    columns = prepared$fit_columns,
    integer = TRUE
  )
  old_truth <- .get_gene_vector(
    data = data,
    mode.slot = "old",
    ensg = prepared$ensg,
    columns = prepared$fit_columns,
    integer = TRUE
  )

  # Preserve the original behavior: initialization may use new.count/ntr even
  # though the likelihood itself uses exact truth counts.
  new_fit <- .get_initial_new_counts(
    data = data,
    ensg = prepared$ensg,
    fit_columns = prepared$fit_columns,
    counts_fit = prepared$counts_fit
  )

  loglik <- make_loglik_from_truth(
    new_counts = new_truth,
    old_counts = old_truth,
    total_counts = as.integer(prepared$counts_fit),
    t.label = t.label,
    steady = steady,
    t.diff = t.diff,
    init = prepared$total_ref
  )

  .run_gene_fit(
    prepared = prepared,
    loglik = loglik,
    new_fit = new_fit,
    fit.cond = fit.cond,
    t.label = t.label,
    steady = steady,
    ref.cond = ref.cond,
    t.diff = t.diff,
    ndraw = ndraw,
    seed = seed,
    truth_suffix = truth_suffix,
    add_pit = add_pit,
    ridge_cor_threshold = ridge_cor_threshold,
    pass_N0 = TRUE,
    # The old truth function did not reject low-information simulated genes.
    check_information = FALSE
  )
}


# Identifiability diagnostic ---------------------------------------------------

identifiability_diagnostic <- function(
    fit,
    corr_warn = 0.90,
    corr_bad = 0.97,
    ridge_warn = 0.01,
    ridge_bad = 0.001,
    d_corr_warn = 0.95
) {
  out <- list(
    ok = FALSE,
    flag = "bad_hessian",
    info_mat = NULL,
    eigenvalues = NA_real_,
    ridge_score = NA_real_,
    condition_number = NA_real_,
    cov_fit = NULL,
    cor_fit = NULL,
    max_corr_fit = NA_real_,
    worst_pair_fit = NA_character_,
    worst_corr_fit = NA_real_,
    cov_bio = NULL,
    cor_bio = NULL,
    max_corr_bio = NA_real_,
    worst_pair_bio = NA_character_,
    worst_corr_bio = NA_real_,
    diagnostic_score = NA_real_
  )

  if (is.null(fit$hessian)) {
    out$flag <- "no_hessian"
    out$ok <- TRUE
    return(out)
  }

  info <- -fit$hessian
  info <- 0.5 * (info + t(info))
  out$info_mat <- info

  if (any(!is.finite(info))) {
    out$ok <- TRUE
    return(out)
  }

  eigenvalues <- tryCatch(
    eigen(info, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, nrow(info))
  )
  out$eigenvalues <- eigenvalues

  if (any(!is.finite(eigenvalues)) || any(eigenvalues <= 0)) {
    out$ok <- TRUE
    return(out)
  }

  out$ridge_score <- min(eigenvalues) / max(eigenvalues)
  out$condition_number <- max(eigenvalues) / min(eigenvalues)

  covariance <- tryCatch(solve(info), error = function(e) NULL)
  if (is.null(covariance) || any(!is.finite(covariance))) {
    out$ok <- TRUE
    return(out)
  }
  covariance <- 0.5 * (covariance + t(covariance))
  out$cov_fit <- covariance

  correlation <- tryCatch(stats::cov2cor(covariance), error = function(e) NULL)
  if (is.null(correlation) || any(!is.finite(correlation))) {
    out$ok <- TRUE
    return(out)
  }
  out$cor_fit <- correlation

  fit_summary <- .largest_off_diagonal_correlation(correlation)
  out$max_corr_fit <- fit_summary$absolute
  out$worst_pair_fit <- fit_summary$pair
  out$worst_corr_fit <- fit_summary$signed

  needed <- c("log_lambda.B", "logB.B", "logd")
  rn <- rownames(covariance)

  if (!is.null(rn) && all(needed %in% rn)) {
    vl <- covariance["log_lambda.B", "log_lambda.B"]
    vB <- covariance["logB.B", "logB.B"]
    vd <- covariance["logd", "logd"]
    c_lB <- covariance["log_lambda.B", "logB.B"]
    c_ld <- covariance["log_lambda.B", "logd"]
    c_Bd <- covariance["logB.B", "logd"]

    covariance_bio <- matrix(
      NA_real_,
      nrow = 3L,
      ncol = 3L,
      dimnames = list(
        c("logB", "logkm", "logd"),
        c("logB", "logkm", "logd")
      )
    )

    covariance_bio["logB", "logB"] <- vB
    covariance_bio["logkm", "logkm"] <- vl + vB - 2 * c_lB
    covariance_bio["logd", "logd"] <- vd
    covariance_bio["logB", "logkm"] <- c_lB - vB
    covariance_bio["logkm", "logB"] <- c_lB - vB
    covariance_bio["logB", "logd"] <- c_Bd
    covariance_bio["logd", "logB"] <- c_Bd
    covariance_bio["logkm", "logd"] <- c_ld - c_Bd
    covariance_bio["logd", "logkm"] <- c_ld - c_Bd

    out$cov_bio <- covariance_bio

    correlation_bio <- tryCatch(
      stats::cov2cor(covariance_bio),
      error = function(e) NULL
    )

    if (!is.null(correlation_bio) && all(is.finite(correlation_bio))) {
      out$cor_bio <- correlation_bio
      bio_summary <- .largest_off_diagonal_correlation(correlation_bio)
      out$max_corr_bio <- bio_summary$absolute
      out$worst_pair_bio <- bio_summary$pair
      out$worst_corr_bio <- bio_summary$signed
    }
  }

  corr_pool <- c(out$max_corr_fit, out$max_corr_bio)
  corr_pool <- corr_pool[is.finite(corr_pool)]
  max_corr_use <- if (length(corr_pool) > 0L) max(corr_pool) else NA_real_

  score_parts <- c(
    if (is.finite(out$ridge_score)) -log10(out$ridge_score) else NA_real_,
    if (is.finite(max_corr_use)) 1 / pmax(1 - max_corr_use, 1e-8) else NA_real_
  )
  score_parts <- score_parts[is.finite(score_parts)]
  out$diagnostic_score <- if (length(score_parts) > 0L) max(score_parts) else NA_real_

  fit_is_lambda_d <- (
    length(out$worst_pair_fit) == 2L &&
      setequal(out$worst_pair_fit, c("log_lambda.B", "logd"))
  )
  strong_fit_d <- (
    fit_is_lambda_d && is.finite(out$worst_corr_fit) &&
      abs(out$worst_corr_fit) >= d_corr_warn
  )

  bio_is_B_km <- (
    length(out$worst_pair_bio) == 2L &&
      setequal(out$worst_pair_bio, c("logB", "logkm"))
  )
  strong_bio_Bkm <- (
    bio_is_B_km && is.finite(out$worst_corr_bio) &&
      abs(out$worst_corr_bio) >= corr_bad
  )
  moderate_bio_Bkm <- (
    bio_is_B_km && is.finite(out$worst_corr_bio) &&
      abs(out$worst_corr_bio) >= corr_warn
  )

  severe_ridge <- is.finite(out$ridge_score) && out$ridge_score < ridge_bad
  moderate_ridge <- is.finite(out$ridge_score) && out$ridge_score < ridge_warn

  out$flag <- if (severe_ridge) {
    "strong_tradeoff"
  } else if (strong_fit_d) {
    "d_related_tradeoff"
  } else if (strong_bio_Bkm) {
    "strong_tradeoff_Bkm"
  } else if (moderate_bio_Bkm) {
    "B_km_tradeoff"
  } else if (moderate_ridge) {
    "moderate_tradeoff"
  } else {
    "ok"
  }

  out$ok <- TRUE
  out
}

.largest_off_diagonal_correlation <- function(x) {
  if (!is.matrix(x) || nrow(x) < 2L || ncol(x) < 2L) {
    return(list(absolute = NA_real_, signed = NA_real_, pair = NA_character_))
  }

  upper <- upper.tri(x)
  values <- x[upper]
  finite <- is.finite(values)

  if (!any(finite)) {
    return(list(absolute = NA_real_, signed = NA_real_, pair = NA_character_))
  }

  positions <- which(upper, arr.ind = TRUE)
  positions <- positions[finite, , drop = FALSE]
  values <- values[finite]
  index <- which.max(abs(values))
  ij <- positions[index, ]

  row_names <- rownames(x)
  col_names <- colnames(x)
  if (is.null(row_names)) row_names <- as.character(seq_len(nrow(x)))
  if (is.null(col_names)) col_names <- as.character(seq_len(ncol(x)))

  list(
    absolute = abs(values[index]),
    signed = values[index],
    pair = c(row_names[ij[1]], col_names[ij[2]])
  )
}


