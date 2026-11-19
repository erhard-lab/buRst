# Shared numerical utilities.

#' buRst: Transcriptional Bursting Inference from Metabolic RNA Labeling
#'
#' Infer burst frequency, mean burst size and RNA degradation rates from
#' single-cell metabolic RNA-labeling measurements. Supports steady-state and
#' transient models, uncertainty summaries, simulation and condition comparisons.
#'
#' @details
#' Start with the `bursting_same_mean` vignette for the forward model,
#' `inference` for experimental fitting, and `simulation` for controlled examples.
#' The full inference route integrates discrete labeling likelihoods and uses a
#' local Laplace approximation for kinetic-parameter uncertainty.
#'
#' @seealso `run_kinetics_batch()`, `simulateBursting()`,
#'   `compare_posteriors_gate_star_kmB()`
#'
#' @importFrom stats dbinom fft integrate optimize quantile setNames uniroot
#' @importFrom utils head read.delim tail
#' @keywords internal
"_PACKAGE"

.logsumexp_vec <- function(x) {
  if (length(x) == 0L) return(-Inf)
  largest <- max(Re(x))
  if (is.na(largest)) return(NA_real_)
  if (!is.finite(largest)) return(largest)
  largest + log(sum(exp(x - largest)))
}

.logspace_add_vec <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  out[!is.finite(m)] <- m[!is.finite(m)]
  out
}

.safe_norm <- function(p, eps = 0) {
  p <- Re(p)
  if (length(p) == 0L || any(!is.finite(p))) {
    stop("A PMF must contain finite probabilities.", call. = FALSE)
  }
  if (length(eps) != 1L || !is.finite(eps) || eps < 0 || eps >= 1) {
    stop("eps must be a finite probability floor in [0, 1).", call. = FALSE)
  }
  tolerance <- 1e-10 * max(1, max(abs(p)))
  if (any(p < -tolerance)) {
    stop("A PMF contains materially negative probabilities.", call. = FALSE)
  }
  # Remove FFT round-off without shifting every coefficient or inventing mass.
  p <- pmax(p, eps)
  total <- sum(p)
  if (!is.finite(total) || total <= 0) {
    stop("A PMF must have finite positive total mass.", call. = FALSE)
  }
  p / total
}

.pmf_to_vec <- function(df, N) {
  # Extract entries only: omitted mass remains omitted, not redistributed.
  out <- numeric(N + 1L)
  keep <- df$k >= 0 & df$k <= N
  out[df$k[keep] + 1L] <- df$prob[keep]
  out
}

.log_width <- function(x, probs = c(0.05, 0.95), eps = 1e-12) {
  q <- quantile(x + eps, probs = probs, na.rm = TRUE)
  as.numeric(log2(q[2] / q[1]))
}
