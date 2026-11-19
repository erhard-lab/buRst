#' Bursting model moment and generating-function helpers
#'
#' These functions implement analytic components used to compute RNA count
#' distributions under a bursty transcription model with first-order degradation.
#'
#' @keywords internal
#' @noRd
NULL

# φ(u): factorial-cumulant generating function used for FFT PMFs (new RNA)
phi <- function(u, t.end, km, d, B) {
  if (is.infinite(t.end) && t.end > 0) {
    return(-(km / d) * log(1 - B * u))
  }
  gp <- 1/(B + 1)
  M  <- function(u) gp / (1 - (1 + u) * (1 - gp))
  int.fun <- function(s) M(u * exp(-d * s)) - 1
  km * (
    integrate(function(s) Re(int.fun(s)), 0, t.end)$value +
      1i * integrate(function(s) Im(int.fun(s)), 0, t.end)$value
  )
}

# G(x): log PGF for old RNA given an initial distribution at t=0
G <- function(x, t.end, d, val, prob) {
  G0  <- (x - 1) * exp(-d * t.end) + 1
  use <- is.finite(prob) & prob > 0
  val <- val[use]
  prob <- prob[use]
  log_terms <- val * log(G0)
  log_terms[val == 0] <- 0
  .logsumexp_vec(log_terms + log(prob))
}

#' Estimate initial bursting parameters from non-steady-state data
#'
#' @param params Numeric vector c(d, m, v, m_0, v_0).
#' @param t.end Final time point.
#' @param B Optional known burst size.
#' @param km Optional known burst frequency.
#' @return Named numeric vector c(B=..., km=...).
#' @details
#' The moments describe counts per cell in the reference and final populations.
#' This moment-based result is an initialization, not a likelihood fit or an
#' uncertainty estimate. Supply time and rate values in consistent units.
#'
#' @examples
#' # d, final mean/variance, and reference mean/variance.
#' params <- c(d = 0.5, m = 5, v = 15, m_0 = 4, v_0 = 12)
#' initial_params_noss(params, t.end = 1, B = 2)
#'
#' @export
initial_params_noss <- function(params, t.end, B = NULL, km = NULL) {
  d   <- params[1]
  m   <- params[2]
  v   <- params[3]
  m_0 <- params[4]
  v_0 <- params[5]

  opt <- optimize(
    f = params_est,
    interval = c(0, 1e10),
    d = d, t.end = t.end, m = m, v = v, m_0 = m_0, v_0 = v_0
  )

  if (is.null(B)) B <- as.numeric(opt$minimum)

  if (is.na(opt$objective)) {
    return(c(B = NA_real_, km = NA_real_))
  }

  if (!is.null(km) && is.null(B)) {
    B <- (d * (m * exp(d * t.end) - m_0)) / (km * (exp(d * t.end) - 1))
  }
  if (is.null(km)) {
    km <- (d * (m * exp(d * t.end) - m_0)) / (B * (exp(d * t.end) - 1))
  }

  if (is.na(B)  || B  <= 0) B  <- 1e-6
  if (is.na(km) || km <= 0) km <- 1e-6

  c(B = B, km = km)
}

#' Estimate an initial degradation rate
#'
#' @param ntr Aggregate new-to-total RNA ratio in `[0, 1]`.
#' @param total Mean total RNA count in the fitted population.
#' @param t Positive labeling duration.
#' @param t0 Modeled reference-to-observation duration; must equal t.
#' @param f0 Mean reference RNA count for a non-steady-state initialization.
#'   total and f0 must use the same per-cell scale, not unequal population sums.
#' @return Initial degradation rate; NA when a positive reference is unavailable.
fit_d <- function(ntr, total, t, t0 = NULL, f0 = NULL) {
  .validate_labeling_interval(t, steady = is.null(f0), t.diff = t0)
  if (length(ntr) != 1L || !is.finite(ntr) || ntr < 0 || ntr > 1) {
    stop("ntr must be one finite ratio in [0, 1].", call. = FALSE)
  }
  if (length(total) != 1L || !is.finite(total) || total < 0) {
    stop("total must be one finite non-negative per-cell count.", call. = FALSE)
  }
  if (is.null(f0)) return(-log1p(-ntr) / t)
  if (length(f0) != 1L || !is.finite(f0) || f0 <= 0) return(NA_real_)
  survival <- min(total * (1 - ntr) / f0, 1)
  -log(survival) / t
}

#' Moments for total RNA count under bursting (and optional non-steady init moments)
#'
#' @param B Mean burst size.
#' @param km Burst initiation rate.
#' @param d Degradation rate.
#' @param t.end Time horizon.
#' @param m_0 Initial mean at t=0 (optional).
#' @param m_s_0 Initial second moment at t=0 (optional).
#' @return Named numeric vector of moments and derived quantities.
#' @keywords internal
#' @noRd
m_v_total <- function(B, km, d, t.end, m_0 = 0, m_s_0 = 0) {
  ps  <- 1/(B + 1)
  B_s <- (1 - ps)/(ps^2) + B^2

  m <- (B * km)/d + (m_0 + (-B * km)/d) * exp(-d * t.end)

  c1 <- ((B^2 * km^2))/d^2 +
    km * ((B * (1/2 - 2 * m_0))/d - B_s/(2 * d)) +
    m_s_0 - m_0

  m_s <- -(2 * (B^2) * (km^2) * exp(-d * t.end))/d^2 +
    (B^2 * km^2)/d^2 +
    (2 * B * km * m_0 * exp(-d * t.end))/d -
    (B * km * exp(-d * t.end))/d +
    (B * km)/(2 * d) +
    c1 * exp(-2 * d * t.end) +
    (km * B_s)/(2 * d) +
    m_0 * exp(-d * t.end)

  var_d <- m_s - m^2
  se_d  <- sqrt(var_d)

  # legacy extra fields kept for compatibility
  m_t <- (m_s/m)^3
  m_q <- (m_t^4 * m^4)/(m_s^6)

  return(c("E" = m, "var" = var_d, "se" = se_d, "E2" = m_s, "E3" = m_t, "E4" = m_q))
}

m_v_new <- function(B, km, d, t.end, t_0 = 0, m_0 = 0, m_s_0 = 0) {
  ps   <- 1/(B + 1)
  Es   <- B
  vars <- (1 - ps)/(ps^2)
  B_s  <- vars + Es^2

  c  <- (m_0 + (-B * km)/d)/(exp(-d * t_0))
  c1 <- (m_s_0 + (-B * km * m_0 - 0.5 * d * m_0 - 0.5 * km * B_s)/(d))/(exp(-2 * d * t_0))

  m   <- (B * km)/d + c * exp(-d * t.end)
  m_s <- ((B_s * km)/(2 * d)) + ((B * km * m)/d) + (m/2) + c1 * exp(-2 * d * t.end)
  var_d <- m_s - m^2
  se_d  <- sqrt(abs(var_d))

  c(E = m, `E^2` = m_s, var = var_d, se = se_d, c = c, c1 = c1)
}

params_est <- function(B, d, t.end, m_0, v_0, m, v) {
  m_s_0 <- v_0 + m_0^2
  m_s   <- v + m^2

  ps  <- 1/(B + 1)
  B_s <- (1 - ps)/(ps^2) + B^2

  km <- (d * (m * exp(d * t.end) - m_0))/(B * (exp(d * t.end) - 1))

  c1 <- ((B^2 * km^2))/d^2 +
    km * ((B * (1/2 - 2 * m_0))/d - B_s/(2 * d)) +
    m_s_0 - m_0

  ((-(2 * (B^2) * (km^2) * exp(-d * t.end))/d^2 +
      (B^2 * km^2)/d^2 +
      (2 * B * km * m_0 * exp(-d * t.end))/d -
      (B * km * exp(-d * t.end))/d +
      (B * km)/(2 * d) +
      c1 * exp(-2 * d * t.end) +
      (km * B_s)/(2 * d) +
      m_0 * exp(-d * t.end)) - m_s)^2
}
