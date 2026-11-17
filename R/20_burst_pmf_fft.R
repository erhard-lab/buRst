#' Bursting model PMFs via FFT
#'
#' Implements stable FFT-based computation of discrete RNA count distributions
#' under the bursting model, with adaptive support (tail control).
#'
#' @name burst_pmf_fft
NULL

# Invert one log-PGF on its full Fourier grid. Signed FFT coefficients are
# clipped as in the original implementation. Record the correction separately;
# a signed numerical Fourier approximation is not itself a valid PMF.
.fft_pmf_from_log <- function(log_pgf) {
  if (length(log_pgf) == 0L || any(is.na(log_pgf)) ||
      any(!is.finite(log_pgf) & Re(log_pgf) != -Inf)) {
    stop("Invalid log-PGF values during FFT evaluation.", call. = FALSE)
  }
  shift <- max(Re(log_pgf))
  if (!is.finite(shift)) {
    stop("The log-PGF has no finite scaling point.", call. = FALSE)
  }
  coefficients <- Re(fft(exp(log_pgf - shift), inverse = TRUE)) /
    length(log_pgf)
  if (any(!is.finite(coefficients))) {
    stop("Non-finite coefficients during FFT inversion.", call. = FALSE)
  }
  negative_mass <- sum(pmax(-coefficients, 0))
  positive_mass <- sum(pmax(coefficients, 0))
  probabilities <- .safe_norm(pmax(coefficients, 0))
  attr(probabilities, "fft_negative_mass") <- negative_mass
  attr(probabilities, "fft_negative_fraction") <- negative_mass / positive_mass
  attr(probabilities, "fft_min_coefficient") <- min(coefficients)
  probabilities
}

.choose_N_min <- function(N, m, se_d) {
  N_min <- as.integer(ceiling(m + 4 * se_d) + 1L)
  if (is.null(N)) return(N_min)
  max(as.integer(N), N_min)
}

.ensure_tail_N <- function(make_df,
                           N0 = 64L,
                           grow = 1.5,
                           eps = 1e-10,
                           tail_bins = 10L,
                           maxN = 2^18,
                           max_tries = 4L) {
  N <- as.integer(N0)
  tries <- 0L

  repeat {
    df <- make_df(N)

    p <- df$prob
    p[!is.finite(p)] <- 0
    p <- pmax(p, 0)
    sp <- sum(p)

    if (!is.finite(sp) || sp <= 0) {
      tries <- tries + 1L
      if (tries >= max_tries || N >= maxN) {
        attr(df, "N_used") <- N
        attr(df, "tail_proxy") <- NA_real_
        warning("PMF normalization failed; returning last df (check make_df).")
        return(list(df = df, N = N, tail_proxy = NA_real_))
      }
      N <- as.integer(ceiling(N * grow))
      next
    }

    p <- p / sp
    tb <- min(as.integer(tail_bins), length(p))
    tail_proxy <- sum(tail(p, tb))

    if (is.finite(tail_proxy) && tail_proxy < eps) {
      df$prob <- p
      attr(df, "N_used") <- N
      attr(df, "tail_proxy") <- tail_proxy
      attr(df, "tail_converged") <- tail_proxy < eps
      return(list(df = df, N = N, tail_proxy = tail_proxy))
    }

    tries <- tries + 1L
    if (tries >= max_tries || N >= maxN) {
      df$prob <- p
      attr(df, "N_used") <- N
      attr(df, "tail_proxy") <- tail_proxy
      attr(df, "tail_converged") <- tail_proxy < eps
      return(list(df = df, N = N, tail_proxy = tail_proxy))
    }

    N <- as.integer(ceiling(N * grow))
  }
}

.convolve_pmf_fft <- function(p1, p2) {
  n1 <- length(p1)
  n2 <- length(p2)
  fft_len <- 2^(ceiling(log2(n1 + n2 - 1)))
  f1 <- fft(c(p1, rep(0, fft_len - n1)))
  f2 <- fft(c(p2, rep(0, fft_len - n2)))
  out <- Re(fft(f1 * f2, inverse = TRUE)) / fft_len
  out[1:(n1 + n2 - 1)]
}

#' PMF of newly synthesized RNA under bursting
#'
#' @param t_end Labeling time.
#' @param km Burst initiation rate.
#' @param d Degradation rate.
#' @param B Mean burst size.
#' @param N Optional FFT support size.
#' @param eps_tail Tail-mass proxy threshold.
#' @param tail_bins Number of bins used for tail proxy.
#' @param maxN Maximum value of N to consider.
#' @param ... Additional arguments passed to ...
#' @return data.frame with columns k and prob.
#' @export
burst_pmf_new <- function(t_end = 2, km, d, B,
                          N = NULL,
                          eps_tail = 1e-12,
                          tail_bins = 10L,
                          maxN = 2^20,
                          ...) {
  m_v  <- m_v_total(B = B, km = km, d = d, t.end = t_end)
  se_d <- as.numeric(m_v[3])
  m    <- as.numeric(m_v[1])

  N_start <- .choose_N_min(N, m, se_d)

  make_df <- function(Nloc) {
    k <- 0:(Nloc - 1L)

    logA <- sapply(
      exp(-2 * pi * 1i * k / Nloc) - 1,
      phi,
      t.end = t_end, km = km, d = d, B = B
    )

    probabilities <- .fft_pmf_from_log(logA)
    df <- data.frame(k = k, prob = as.numeric(probabilities))
    for (field in c("fft_negative_mass", "fft_negative_fraction", "fft_min_coefficient")) {
      attr(df, field) <- attr(probabilities, field)
    }
    df
  }

  res <- .ensure_tail_N(
    make_df,
    N0 = N_start,
    eps = eps_tail,
    tail_bins = tail_bins,
    maxN = maxN
  )

  df <- res$df
  df$type <- "fft_new"
  df
}

#' PMF of old RNA after decay for time t_end
#'
#' Old RNA corresponds to the survival/decay of molecules present at t=0,
#' with no synthesis during [0, t_end].
#'
#' @param t_end Time horizon.
#' @param d Degradation rate.
#' @param steady If TRUE, uses steady-state init from bursting model.
#' @param init If steady=FALSE: initial distribution at t=0 (data.frame k/prob or raw samples).
#' @param km,B Only used if steady=TRUE to generate the steady-state init distribution.
#' @param N ...
#' @param eps_tail ...
#' @param tail_bins ...
#' @param maxN ...
#' @param max_tries Maximum number of attempts used when extending or fitting the PMF support.
#' @param ... Additional arguments passed to ...
#' @return data.frame with columns k and prob.
#' @export
burst_pmf_old <- function(t_end, d,
                          steady = TRUE,
                          init = NULL,
                          km = NULL,
                          B = NULL,
                          N = NULL,
                          eps_tail = 1e-10,
                          tail_bins = 10L,
                          maxN = 2^18,
                          max_tries = 4L,
                          ...) {

  # init at t=0
  if (isTRUE(steady)) {
    if (is.null(km) || is.null(B)) stop("steady=TRUE requires km and B to form steady-state init.")
    init <- burst_pmf_new(t_end = Inf, km = km, d = d, B = B, N = N, eps_tail = eps_tail, ...)
  } else {
    if (is.null(init)) stop("steady=FALSE requires init.")
    if (!is.data.frame(init)) {
      init <- as.data.frame(table(init))
      names(init) <- c("k", "prob")
      init$k <- as.numeric(as.character(init$k))
      init$prob <- init$prob / sum(init$prob)
    } else {
      stopifnot(all(c("k", "prob") %in% names(init)))
      init$k <- as.numeric(init$k)
      init$prob <- as.numeric(init$prob)
      init$prob <- init$prob / sum(init$prob)
    }
  }

  m0  <- sum(init$k * init$prob)
  ms0 <- sum((init$k^2) * init$prob)

  m_v  <- m_v_total(B = 0, km = 0, d = d, t.end = t_end, m_0 = m0, m_s_0 = ms0)
  se_d <- as.numeric(m_v[3])
  m    <- as.numeric(m_v[1])

  N_start <- .choose_N_min(N, m, se_d)

  make_df <- function(Nloc) {
    k <- 0:(Nloc - 1L)

    logA <- sapply(
      exp(-2 * pi * 1i * k / Nloc),
      G,
      t.end = t_end, d = d, val = init$k, prob = init$prob
    )

    probabilities <- .fft_pmf_from_log(logA)
    df <- data.frame(k = k, prob = as.numeric(probabilities))
    for (field in c("fft_negative_mass", "fft_negative_fraction", "fft_min_coefficient")) {
      attr(df, field) <- attr(probabilities, field)
    }
    df
  }

  res <- .ensure_tail_N(
    make_df,
    N0 = N_start,
    eps = eps_tail,
    tail_bins = tail_bins,
    maxN = maxN,
    max_tries = max_tries
  )

  df <- res$df
  df$type <- "fft_old"
  df
}

#' PMF of total RNA at t_end (old + new)
#'
#' @param t_end Time horizon.
#' @param km Burst initiation rate.
#' @param d Degradation rate.
#' @param B Mean burst size.
#' @param steady If TRUE, old pool init is steady-state.
#' @param init If steady=FALSE, initial distribution for old pool at t=0.
#' @param N Initial or fixed truncation/support size.
#' @param eps_tail Tail-probability tolerance used to determine whether the support is sufficiently large.
#' @param tail_bins Number of bins at the upper end of the support used to assess remaining tail mass.
#' @param maxN Maximum allowed support size.
#' @param ... Additional arguments passed to lower-level PMF routines.
#' @return data.frame(k, prob, type).
#' @export
burst_pmf_total <- function(t_end, km, d, B,
                            steady = TRUE,
                            init = NULL,
                            N = NULL,
                            eps_tail = 1e-12,
                            tail_bins = 10L,
                            maxN = 2^20,
                            ...) {
  # old pmf (decay)
  if (steady) {
    old_df <- burst_pmf_old(
      t_end = t_end, d = d, steady = TRUE, km = km, B = B,
      N = N, eps_tail = eps_tail, tail_bins = tail_bins, maxN = maxN, ...
    )
  } else {
    old_df <- burst_pmf_old(
      t_end = t_end, d = d, steady = FALSE, init = init,
      N = N, eps_tail = eps_tail, tail_bins = tail_bins, maxN = maxN, ...
    )
  }

  # new pmf
  new_df <- burst_pmf_new(
    t_end = t_end, km = km, d = d, B = B,
    N = N, eps_tail = eps_tail, tail_bins = tail_bins, maxN = maxN, ...
  )

  probs <- .safe_norm(.convolve_pmf_fft(old_df$prob, new_df$prob))

  df <- data.frame(
    k = seq_along(probs) - 1L,
    prob = probs,
    type = "fft_total"
  )

  attr(df, "tail_old") <- attr(old_df, "tail_proxy")
  attr(df, "tail_new") <- attr(new_df, "tail_proxy")
  df
}

# Backward-compatible wrappers (optional)
#' @keywords internal
#' @noRd
compute.distr.new <- function(t.end = 2, km, d, B, ...) burst_pmf_new(t_end = t.end, km = km, d = d, B = B, ...)
#' @keywords internal
#' @noRd
compute.distr.old <- function(t.end, km, d, B, steady = TRUE, init = NULL, ...) {
  if (steady) burst_pmf_old(t_end = t.end, d = d, steady = TRUE, km = km, B = B, ...)
  else burst_pmf_old(t_end = t.end, d = d, steady = FALSE, init = init, ...)
}
#' @keywords internal
#' @noRd
compute.distr.total <- function(t.end, km, d, B, steady = TRUE, init = NULL, ...) burst_pmf_total(t_end = t.end, km = km, d = d, B = B, steady = steady, init = init, ...)
