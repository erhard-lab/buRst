# Fit bursting kinetics using per-cell Beta(alpha, beta) uncertainty

# alpha and beta describe a Beta approximation to the per-cell NTR posterior.
# For a cell with total count N, this file converts that approximation into
# an N+1-vector over possible new counts k = 0,...,N.
#
# IMPORTANT:
# make_loglik_from_gp() expects the same profile convention as the HDF5 route:
# a relative negative-log-likelihood penalty whose minimum is zero.
# It does NOT expect a normalized probability vector.

.get_gene_uncertainty_vector <- function(
    data,
    mode.slot,
    ensg,
    columns
) {
  column_names <- as.character(columns)

  # Prefer the standard grandR accessor.
  out <- tryCatch(
    grandR::GetMatrix(
      data,
      mode.slot = mode.slot,
      genes = ensg,
      columns = columns
    ),
    error = function(e) NULL
  )

  if (!is.null(out)) {
    if (nrow(out) != 1L) {
      stop(
        "Expected exactly one row for gene ", ensg,
        " in slot '", mode.slot, "', got ", nrow(out), ".",
        call. = FALSE
      )
    }

    values <- as.numeric(out[1L, , drop = TRUE])
    names(values) <- colnames(out)

    missing <- setdiff(column_names, names(values))
    if (length(missing) > 0L) {
      stop(
        "Cells missing from uncertainty slot '", mode.slot, "': ",
        paste(utils::head(missing, 10L), collapse = ", "),
        call. = FALSE
      )
    }

    return(values[column_names])
  }

  # Fallback for slots stored directly in data$data.
  mat <- data$data[[mode.slot]]

  if (is.null(mat)) {
    stop(
      "Could not find uncertainty slot '", mode.slot,
      "' through grandR::GetMatrix() or data$data[['", mode.slot, "']].",
      call. = FALSE
    )
  }

  gene_index <- match(ensg, rownames(mat))
  if (is.na(gene_index)) {
    stop(
      "Gene ", ensg, " not found in uncertainty slot '", mode.slot, "'.",
      call. = FALSE
    )
  }

  column_index <- match(column_names, colnames(mat))
  if (anyNA(column_index)) {
    missing <- column_names[is.na(column_index)]
    stop(
      "Cells missing from uncertainty slot '", mode.slot, "': ",
      paste(utils::head(missing, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  values <- as.numeric(
    mat[gene_index, column_index, drop = TRUE]
  )
  names(values) <- column_names
  values
}


.beta_ntr_to_count_profile <- function(
    total_count,
    alpha,
    beta,
    R_clip = 60,
    probability_floor = 1e-300
) {
  total_count <- as.integer(round(total_count))

  if (
    length(total_count) != 1L ||
    is.na(total_count) ||
    total_count < 0L
  ) {
    stop(
      "total_count must be one non-negative integer.",
      call. = FALSE
    )
  }

  # A zero-count cell has no possible new-count uncertainty to contribute.
  if (total_count == 0L) {
    return(numeric(0))
  }

  if (
    length(alpha) != 1L ||
    length(beta) != 1L ||
    !is.finite(alpha) ||
    !is.finite(beta) ||
    alpha <= 0 ||
    beta <= 0
  ) {
    stop(
      "Positive-count cells require finite, strictly positive alpha and beta.",
      call. = FALSE
    )
  }

  k <- 0:total_count

  # Integrate the Beta approximation over the interval that rounds NTR*N to k.
  lower <- pmax(0, (k - 0.5) / total_count)
  upper <- pmin(1, (k + 0.5) / total_count)

  lower[1L] <- 0
  upper[length(upper)] <- 1

  mass <- stats::pbeta(
    upper,
    shape1 = alpha,
    shape2 = beta
  ) - stats::pbeta(
    lower,
    shape1 = alpha,
    shape2 = beta
  )

  mass[!is.finite(mass)] <- 0
  mass <- pmax(mass, probability_floor)

  if (!any(is.finite(mass) & mass > 0)) {
    stop(
      "Could not construct a count profile from Beta(",
      signif(alpha, 5), ", ", signif(beta, 5), ").",
      call. = FALSE
    )
  }

  # Convert posterior mass to a relative negative-log profile.
  #
  # best-supported k: 0
  # less-supported k: positive penalty
  log_mass <- log(mass)
  profile <- max(log_mass) - log_mass

  profile[!is.finite(profile)] <- R_clip
  profile <- pmin(pmax(profile, 0), R_clip)

  # Guard against floating-point drift.
  profile <- profile - min(profile)
  profile
}


make_gp_from_beta <- function(
    counts,
    alpha,
    beta,
    R_clip = 60
) {
  if (is.null(names(counts))) {
    stop("counts must be named by cell.", call. = FALSE)
  }

  cell_names <- names(counts)
  alpha <- alpha[cell_names]
  beta <- beta[cell_names]

  if (length(alpha) != length(counts) ||
      length(beta) != length(counts)) {
    stop(
      "alpha and beta must align with all fitted cells.",
      call. = FALSE
    )
  }

  profiles <- Map(
    function(total_count, a, b, cell_name) {
      total_count <- as.integer(round(total_count))

      if (total_count == 0L) {
        return(numeric(0))
      }

      if (
        !is.finite(a) || a <= 0 ||
        !is.finite(b) || b <= 0
      ) {
        stop(
          "Positive-count cell '", cell_name,
          "' has invalid alpha/beta: alpha=", a,
          ", beta=", b, ".",
          call. = FALSE
        )
      }

      .beta_ntr_to_count_profile(
        total_count = total_count,
        alpha = a,
        beta = b,
        R_clip = R_clip
      )
    },
    total_count = counts,
    a = alpha,
    b = beta,
    cell_name = cell_names
  )

  names(profiles) <- cell_names
  profiles
}


#' Fit bursting kinetics using a Beta approximation to per-cell NTR uncertainty
#'
#' This is analogous to fit_gene(), but does not read an HDF5 profile.
#' Each positive-count cell's Beta(alpha,beta) approximation is converted into
#' a relative negative-log profile over possible new counts.
#'
#' @param alpha_slot Slot containing Beta shape1 values.
#' @param beta_slot Slot containing Beta shape2 values.
#' @param R_clip Maximum relative negative-log penalty, matching the HDF5 route.
#' @inheritParams fit_gene
#'
#' @return See fit_gene().
#' @export
fit_gene_beta <- function(
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
    ridge_cor_threshold = -0.9,
    alpha_slot = "alpha",
    beta_slot = "beta",
    R_clip = 60
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

  alpha_fit <- .get_gene_uncertainty_vector(
    data = data,
    mode.slot = alpha_slot,
    ensg = prepared$ensg,
    columns = prepared$fit_columns
  )

  beta_fit <- .get_gene_uncertainty_vector(
    data = data,
    mode.slot = beta_slot,
    ensg = prepared$ensg,
    columns = prepared$fit_columns
  )

  informative <- prepared$counts_fit > 0

  invalid_alpha <- informative & (
    !is.finite(alpha_fit) |
      alpha_fit <= 0
  )
  invalid_beta <- informative & (
    !is.finite(beta_fit) |
      beta_fit <= 0
  )

  if (any(invalid_alpha) || any(invalid_beta)) {
    bad <- names(prepared$counts_fit)[invalid_alpha | invalid_beta]

    stop(
      "Invalid Beta approximation for ",
      length(bad), " positive-count cells of gene ", gene,
      " in condition ", fit.cond, ". Examples: ",
      paste(utils::head(bad, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  gp_fit <- make_gp_from_beta(
    counts = prepared$counts_fit,
    alpha = alpha_fit,
    beta = beta_fit,
    R_clip = R_clip
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
