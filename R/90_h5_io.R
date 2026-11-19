# HDF5 utilities for storing discrete likelihood profiles

#' Write one gene's discrete likelihood profile to an HDF5 file
#'
#' This function stores the quantized discrete log-likelihood profile
#' for a single gene across all cells. The likelihood vectors are provided
#' as a list of integer vectors (one per cell), representing the quantized
#' log-likelihood over \eqn{k = 0, \dots, N}.
#'
#' The vectors are flattened and written to an HDF5 dataset
#' \code{/genes/<gene>}. Metadata describing the per-cell lengths and
#' molecule counts are written to
#' \code{/genes_meta/<gene>/counts_reads} and
#' \code{/genes_meta/<gene>/lens_reads_plus1}.
#'
#' This layout allows efficient streaming writes during simulation
#' and fast reconstruction of per-cell likelihood vectors later.
#'
#' @param h5_path Path to the HDF5 file.
#' @param gene Gene identifier (used as dataset name).
#' @param q_list_gene List of integer vectors containing quantized
#'   likelihood values for each cell.
#' @param counts_gene Integer vector of molecule counts per cell.
#' @param compression_level HDF5 compression level (0–9).
#' @param chunk_target Approximate chunk size for HDF5 storage.
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
write_one_gene_profile_from_disc <- function(h5_path, gene, q_list_gene, counts_gene,
                                             compression_level = 4,
                                             chunk_target = 1e6) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) stop("Need rhdf5.")

  dset <- paste0("/genes/", gene)

  counts_gene <- as.integer(counts_gene)
  lens <- ifelse(counts_gene == 0L, 0L, counts_gene + 1L)

  # validate lengths
  ok <- vapply(seq_along(q_list_gene), function(i) {
    q <- q_list_gene[[i]]
    if (lens[i] == 0L) length(q) == 0L else length(q) == lens[i]
  }, logical(1))
  if (!all(ok)) {
    bad <- which(!ok)[1]
    stop(sprintf(
      "Length mismatch for gene '%s', cell index %d: counts=%d => expected %d, got %d",
      gene, bad, counts_gene[bad], lens[bad], length(q_list_gene[[bad]])
    ))
  }

  flat <- unlist(q_list_gene, use.names = FALSE)
  if (length(flat) == 0L) flat <- as.integer(0)
  flat <- as.integer(flat)

  ch <- min(length(flat), as.integer(chunk_target))
  if (ch < 1L) ch <- 1L

  # delete if present (FAST: no h5ls)
  try(rhdf5::h5delete(h5_path, dset), silent = TRUE)

  rhdf5::h5createDataset(
    file = h5_path,
    dataset = dset,
    dims = length(flat),
    storage.mode = "integer",
    chunk = c(ch),
    level = compression_level
  )
  rhdf5::h5write(flat, h5_path, dset)

  # per-gene meta datasets
  meta_grp <- paste0("/genes_meta/", gene)
  try(rhdf5::h5createGroup(h5_path, meta_grp), silent = TRUE)

  cpath <- paste0(meta_grp, "/counts_reads")
  lpath <- paste0(meta_grp, "/lens_reads_plus1")

  try(rhdf5::h5delete(h5_path, cpath), silent = TRUE)
  try(rhdf5::h5delete(h5_path, lpath), silent = TRUE)

  rhdf5::h5write(as.integer(counts_gene), h5_path, cpath)
  rhdf5::h5write(as.integer(lens),       h5_path, lpath)

  invisible(TRUE)
}

#' Compute PIT values for posterior samples against named truths
#'
#' @param samples_list Named list of numeric vectors (e.g. list(B=..., km=..., d=...)).
#' @param truth_named Named numeric vector with names matching samples_list entries.
#' @return Named numeric vector of PIT values in `[0,1]`.
#' @keywords internal
pit_from_samples <- function(samples_list, truth_named) {
  if (is.null(samples_list) || is.null(truth_named)) return(setNames(rep(NA_real_, length(truth_named)), names(truth_named)))
  sapply(names(truth_named), function(nm) {
    s <- samples_list[[nm]]
    if (is.null(s) || !is.numeric(s) || !is.finite(truth_named[[nm]])) return(NA_real_)
    mean(s <= truth_named[[nm]], na.rm = TRUE)
  })
}
safe_solve_cov <- function(hess) {
  Hneg <- -hess
  Sigma <- tryCatch(solve(Hneg), error = function(e) NULL)

  if (!is.null(Sigma)) {
    ev <- tryCatch(eigen(Sigma, symmetric = TRUE)$values, error = function(e) NULL)
    if (!is.null(ev) && all(is.finite(ev)) && all(ev > 0)) return(Sigma)
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Need Matrix package for nearPD fallback when covariance is not PD.")
  }

  Sigma0 <- tryCatch(solve(Hneg), error = function(e) diag(1, nrow(Hneg)))
  as.matrix(Matrix::nearPD(Sigma0, corr = FALSE)$mat)
}

draw_ci_from_laplace <- function(fit, ndraw = 100000, seed = 1, return_samples = TRUE) {
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Need MASS for mvrnorm().")

  mu <- fit$par
  Sigma <- safe_solve_cov(fit$hessian)

  set.seed(seed)
  S <- MASS::mvrnorm(ndraw, mu = mu, Sigma = Sigma)

  samples <- list(
    lambda = exp(S[, 1]),
    km     = exp(S[, 2]),
    d     = exp(S[, 3])
    # d      = rep(d_fixed, ndraw)
  )
  samples$B <- samples$lambda / samples$km

  CI <- list(
    lambda = quantile(samples$lambda, c(0.025, 0.5, 0.975), names = FALSE),
    km     = quantile(samples$km,     c(0.025, 0.5, 0.975), names = FALSE),
    d      = quantile(samples$d,     c(0.025, 0.5, 0.975), names = FALSE),
    B      = quantile(samples$B,      c(0.025, 0.5, 0.975), names = FALSE)
  )

  if (return_samples) list(CI = CI, samples = samples) else list(CI = CI)
}
