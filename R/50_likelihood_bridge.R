# Likelihoods for read-level profiles and observed component counts.
# PMFs are normalized on their full FFT support. Extracting entries for
# observed counts must never condition the model on the sample maximum.

.validate_rna_counts <- function(x, name = "counts") {
  if (!is.numeric(x) || length(x) == 0L || any(!is.finite(x)) ||
      any(x < 0) || any(abs(x - round(x)) > 1e-8) ||
      any(x >= .Machine$integer.max)) {
    stop(name, " must contain finite, non-negative integer counts.", call. = FALSE)
  }
  as.integer(round(x))
}

.validate_labeling_interval <- function(t.label, steady, t.diff = NULL) {
  if (length(t.label) != 1L || !is.numeric(t.label) ||
      !is.finite(t.label) || t.label <= 0) {
    stop("t.label must be one finite positive duration.", call. = FALSE)
  }
  if (length(steady) != 1L || !is.logical(steady) || is.na(steady)) {
    stop("steady must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(steady)) return(invisible(t.label))
  if (is.null(t.diff)) t.diff <- t.label
  if (length(t.diff) != 1L || !is.numeric(t.diff) ||
      !is.finite(t.diff) || t.diff <= 0) {
    stop("t.diff must be one finite positive duration.", call. = FALSE)
  }
  tolerance <- 1e-8 * max(1, abs(t.label), abs(t.diff))
  if (abs(t.diff - t.label) > tolerance) {
    stop(
      "Non-steady-state inference currently requires t.diff = t.label. ",
      "The reference represents the population at labeling onset; ",
      "synthesis before labeling is not modeled.",
      call. = FALSE
    )
  }
  invisible(t.diff)
}

.build_kinetic_pmfs <- function(B, km, d, N0, Nmax, t.label,
                                steady = TRUE, t.diff = NULL, init = NULL) {
  .validate_labeling_interval(t.label, steady, t.diff)
  if (length(N0) != 1L || !is.finite(N0) || N0 < 1) {
    stop("N0 must be one finite positive support size.", call. = FALSE)
  }
  # The FFT must cover observed counts even if the initial support was capped.
  N0 <- max(ceiling(N0), Nmax + 1L)

  if (isTRUE(steady)) {
    init <- burst_pmf_new(
      t_end = Inf, km = km, d = d, B = B, N = N0, eps_tail = 1e-8
    )
  } else if (is.null(init)) {
    stop("steady = FALSE requires an initial RNA population.", call. = FALSE)
  }

  old_df <- burst_pmf_old(
    t_end = t.label, d = d, steady = FALSE,
    init = init, N = N0, eps_tail = 1e-8
  )
  new_df <- burst_pmf_new(
    t_end = t.label, km = km, d = d, B = B, N = N0, eps_tail = 1e-8
  )

  old_df$prob <- .safe_norm(old_df$prob)
  new_df$prob <- .safe_norm(new_df$prob)
  total_prob <- .safe_norm(.convolve_pmf_fft(old_df$prob, new_df$prob))
  total_df <- data.frame(k = seq_along(total_prob) - 1L, prob = total_prob)

  # These are probability subsets, not distributions to renormalize.
  list(
    p_new = .pmf_to_vec(new_df, Nmax),
    p_old = .pmf_to_vec(old_df, Nmax),
    p_tot = .pmf_to_vec(total_df, Nmax)
  )
}

cond_from_uncond_log <- function(N, p_new, p_old, p_tot_full) {
  # Normalize the joint allocation in log space. p_tot_full is retained
  # for compatibility; the denominator is computed from the components.
  index <- seq_len(N + 1L)
  log_weights <- log(p_new[index]) + log(rev(p_old[index]))
  normalizer <- .logsumexp_vec(log_weights)
  if (!is.finite(normalizer)) return(rep(0, N + 1L))
  exp(log_weights - normalizer)
}

build_Emat_grouped_from_gp <- function(gp, counts) {
  counts <- .validate_rna_counts(counts)
  if (!is.list(gp) || length(gp) != length(counts)) {
    stop("gp must contain one profile per cell.", call. = FALSE)
  }
  idx_pos <- which(counts > 0L)
  uniqN <- sort(unique(counts[idx_pos]))
  Emat <- lapply(uniqN, function(Ni) {
    indices <- which(counts == Ni)
    rows <- lapply(indices, function(j) {
      profile <- as.numeric(gp[[j]])
      if (length(profile) != Ni + 1L || any(!is.finite(profile))) {
        stop(
          "Cell ", j, " requires a finite profile of length ", Ni + 1L, ".",
          call. = FALSE
        )
      }
      evidence <- -profile
      evidence - max(evidence)
    })
    do.call(rbind, rows)
  })
  names(Emat) <- as.character(uniqN)
  list(Emat = Emat, uniqN = uniqN, idx_pos = idx_pos, N_vec = counts)
}

make_loglik_from_gp <- function(gp, counts, t.label = 2,
                                steady = TRUE, t.diff = NULL, init = NULL) {
  if (!requireNamespace("matrixStats", quietly = TRUE)) stop("Need matrixStats.")
  .validate_labeling_interval(t.label, steady, t.diff)
  counts <- .validate_rna_counts(counts)
  grouped <- build_Emat_grouped_from_gp(gp, counts)
  maxN <- max(counts)
  n_zero <- sum(counts == 0L)

  # Restore the original buffer-and-column-update strategy. R may still
  # copy on modification; this avoids sweep()'s recycling/permutation path.
  Mbuf <- lapply(grouped$Emat, function(E) matrix(0, nrow(E), ncol(E)))

  f <- function(B, km, d, N0) {
    kin <- .build_kinetic_pmfs(
      B, km, d, N0, maxN, t.label, steady, t.diff, init
    )
    # Evaluate the joint likelihood directly. This is algebraically the
    # total-count likelihood times the conditional measurement term, without
    # division, artificial conditional floors, or sample-maximum conditioning.
    log_new <- log(kin$p_new)
    log_old <- log(kin$p_old)
    ll <- 0
    if (n_zero > 0L) ll <- n_zero * (log_new[1L] + log_old[1L])

    for (g in seq_along(grouped$uniqN)) {
      Ni <- grouped$uniqN[g]
      index <- seq_len(Ni + 1L)
      log_joint <- log_new[index] + rev(log_old[index])
      E <- grouped$Emat[[g]]
      M <- Mbuf[[g]]
      for (k in seq_along(log_joint)) {
        M[, k] <- E[, k] + log_joint[k]
      }
      ll <- ll + sum(matrixStats::rowLogSumExps(M))
    }
    ll
  }
  attr(f, "Emat") <- grouped$Emat
  attr(f, "uniqN") <- grouped$uniqN
  attr(f, "N_vec") <- counts
  f
}

make_loglik_from_truth <- function(new_counts, old_counts, total_counts = NULL,
                                   t.label = 2, steady = TRUE, t.diff = NULL,
                                   init = NULL, eps = 1e-300, check_total = TRUE) {
  .validate_labeling_interval(t.label, steady, t.diff)
  new_counts <- .validate_rna_counts(new_counts, "new_counts")
  old_counts <- .validate_rna_counts(old_counts, "old_counts")
  if (length(new_counts) != length(old_counts)) {
    stop("new_counts and old_counts must have the same length.", call. = FALSE)
  }
  if (length(eps) != 1L || !is.finite(eps) || eps < 0 || eps >= 1) {
    stop("eps must be a finite probability floor in [0, 1).", call. = FALSE)
  }
  if (is.null(total_counts)) total_counts <- new_counts + old_counts
  total_counts <- .validate_rna_counts(total_counts, "total_counts")
  if (length(total_counts) != length(new_counts)) {
    stop("total_counts must have the same length as component counts.", call. = FALSE)
  }
  if (isTRUE(check_total) && any(total_counts != new_counts + old_counts)) {
    stop("total_counts != new + old.", call. = FALSE)
  }
  Nmax <- max(c(new_counts, old_counts, total_counts))

  f <- function(B, km, d, N0, use_total = FALSE) {
    kin <- .build_kinetic_pmfs(
      B, km, d, N0, Nmax, t.label, steady, t.diff, init
    )
    ll <- sum(log(pmax(kin$p_new[new_counts + 1L], eps))) +
      sum(log(pmax(kin$p_old[old_counts + 1L], eps)))
    if (!isTRUE(use_total)) return(ll)
    # Diagnostic only: total counts are already determined by the components.
    ll_total <- sum(log(pmax(kin$p_tot[total_counts + 1L], eps)))
    c(ll = ll, ll_with_total = ll + ll_total)
  }
  attr(f, "new_counts") <- new_counts
  attr(f, "old_counts") <- old_counts
  attr(f, "total_counts") <- total_counts
  attr(f, "Nmax") <- Nmax
  f
}
