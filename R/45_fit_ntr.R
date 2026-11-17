# Internal replacement for grandR:::fit.ntr()
#
# Uses the binomial-mixture likelihood functions already
# implemented in buRst.
#
# Returns:
#   beta.approx = FALSE:
#       numeric NTR estimate
#
#   beta.approx = TRUE:
#       named numeric vector:
#       ntr, alpha, beta, integral
#
# The integral is on the log scale.

.burst_fit_ntr <- function(
    mixmat,
    par,
    beta.approx = FALSE,
    conversion.reads = FALSE,
    plot = FALSE,
    prior = c(1, 1)
) {
  
  # --------------------------------------------------
  # 1. Validate parameters
  # --------------------------------------------------
  
  if (
    length(prior) != 2L ||
    any(!is.finite(prior)) ||
    any(prior <= 0)
  ) {
    stop("prior must contain two positive finite values.")
  }
  
  if (
    !is.finite(par$p.err) ||
    !is.finite(par$p.conv) ||
    par$p.err < 0 ||
    par$p.err > 1 ||
    par$p.conv < 0 ||
    par$p.conv > 1
  ) {
    stop("Invalid conversion probabilities.")
  }
  
  # --------------------------------------------------
  # 2. Prepare likelihood vectors
  # --------------------------------------------------
  
  d <- prep_cell_ll01(mixmat, par)
  
  ll0 <- d$ll0
  ll1 <- d$ll1
  counts <- d$c
  
  # Binomial mixture log-likelihood
  loglik <- function(theta) {
    
    ell_theta_from_ll01(
      theta = theta,
      ll0 = ll0,
      ll1 = ll1,
      c = counts
    )
  }
  
  # Log posterior, including the beta prior
  logpost <- function(theta) {
    
    loglik(theta) +
      stats::dbeta(
        theta,
        shape1 = prior[1],
        shape2 = prior[2],
        log = TRUE
      )
  }
  
  # --------------------------------------------------
  # 3. Maximum posterior estimate
  # --------------------------------------------------
  
  opt <- stats::optimize(
    f = logpost,
    interval = c(0, 1),
    maximum = TRUE
  )
  
  ntr_hat <- opt$maximum
  logpost_max <- opt$objective
  
  # Explicit boundary checks, as in grandR
  ll_left <- logpost(0)
  ll_right <- logpost(1)
  
  if (ll_left > logpost_max) {
    
    ntr_hat <- 0
    logpost_max <- ll_left
    
  } else if (ll_right > logpost_max) {
    
    ntr_hat <- 1
    logpost_max <- ll_right
  }
  
  # --------------------------------------------------
  # 4. Conversion reads
  # --------------------------------------------------
  
  conversion_count <- NULL
  
  if (isTRUE(conversion.reads)) {
    
    zero_row <- match("0", rownames(mixmat))
    
    if (is.na(zero_row)) {
      conversion_count <- sum(mixmat)
    } else {
      conversion_count <-
        sum(mixmat) - sum(mixmat[zero_row, , drop = FALSE])
    }
  }
  
  # --------------------------------------------------
  # 5. Return NTR only if beta fit is disabled
  # --------------------------------------------------
  
  if (!isTRUE(beta.approx)) {
    
    if (isTRUE(conversion.reads)) {
      return(c(
        ntr = ntr_hat,
        conversion.reads = conversion_count
      ))
    }
    
    return(ntr_hat)
  }
  
  # --------------------------------------------------
  # 6. Determine posterior support
  # --------------------------------------------------
  
  # Match grandR's relative log-posterior threshold.
  threshold <- log(0.001)
  
  relative_logpost <- function(theta) {
    logpost(theta) - logpost_max
  }
  
  left <- if (
    relative_logpost(0) > threshold
  ) {
    
    0
    
  } else {
    
    stats::uniroot(
      function(theta) {
        relative_logpost(theta) - threshold
      },
      interval = c(0, ntr_hat)
    )$root
  }
  
  right <- if (
    relative_logpost(1) > threshold
  ) {
    
    1
    
  } else {
    
    stats::uniroot(
      function(theta) {
        relative_logpost(theta) - threshold
      },
      interval = c(ntr_hat, 1)
    )$root
  }
  
  # --------------------------------------------------
  # 7. Evaluate posterior on a 100-point grid
  # --------------------------------------------------
  
  x <- seq(
    from = left,
    to = right,
    length.out = 100L
  )
  
  log_density <- vapply(
    x,
    logpost,
    numeric(1)
  )
  
  # Stabilize exponentiation.
  log_density <- log_density - max(log_density)
  
  density <- exp(log_density)
  
  # --------------------------------------------------
  # 8. Construct numerical posterior CDF
  # --------------------------------------------------
  
  # Trapezoidal integration
  dx <- diff(x)
  
  interval_mass <- (
    density[-length(density)] +
      density[-1L]
  ) * dx / 2
  
  cumulative <- c(
    0,
    cumsum(interval_mass)
  )
  
  total_mass <- tail(cumulative, 1)
  
  if (!is.finite(total_mass) || total_mass <= 0) {
    stop("Could not normalize NTR posterior.")
  }
  
  posterior_cdf <- cumulative / total_mass
  
  # --------------------------------------------------
  # 9. Fit beta distribution to posterior CDF
  # --------------------------------------------------
  
  beta_objective <- function(shapes) {
    
    predicted <- stats::pbeta(
      x,
      shape1 = shapes[1],
      shape2 = shapes[2]
    )
    
    sum((predicted - posterior_cdf)^2)
  }
  
  # Positive shape parameters, matching grandR's
  # constrained optimization approach.
  beta_fit <- stats::constrOptim(
    theta = c(alpha = 3, beta = 3),
    f = beta_objective,
    grad = NULL,
    ui = diag(2),
    ci = c(0, 0)
  )
  
  shapes <- beta_fit$par
  
  names(shapes) <- c("alpha", "beta")
  
  # --------------------------------------------------
  # 10. Log-integrated likelihood
  #     Original grandR integration method
  # --------------------------------------------------
  
  # Evaluate the binomial-mixture log-likelihood
  # at the same grid points used for the beta fit.
  
  log_lik_grid <- vapply(
    x,
    loglik,
    numeric(1)
  )
  
  # Match the original grandR implementation:
  # lsse(logpost) + log(x[2] - x[1])
  #
  # Use the existing buRst log-sum-exp helper.
  
  log_integral <- .logsumexp_vec(log_lik_grid) +
    log(x[2] - x[1])
  
  # --------------------------------------------------
  # 11. Optional diagnostic plot
  # --------------------------------------------------
  
  if (isTRUE(plot)) {
    
    graphics::plot(
      x,
      posterior_cdf,
      type = "l",
      xlab = "NTR",
      ylab = "Cumulative probability",
      main = "NTR posterior and beta approximation"
    )
    
    graphics::lines(
      x,
      stats::pbeta(
        x,
        shape1 = shapes[1],
        shape2 = shapes[2]
      ),
      col = "red"
    )
    
    graphics::legend(
      "topleft",
      legend = c(
        "Numerical posterior",
        "Beta approximation"
      ),
      col = c("black", "red"),
      lty = 1
    )
  }
  
  # --------------------------------------------------
  # 12. Return results
  # --------------------------------------------------
  
  result <- c(
    ntr = ntr_hat,
    alpha = unname(shapes[1]),
    beta = unname(shapes[2])
  )
  
  if (isTRUE(conversion.reads)) {
    result <- c(
      result,
      conversion.reads = conversion_count
    )
  }
  
  c(
    result,
    integral = log_integral
  )
}