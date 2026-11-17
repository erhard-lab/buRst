#' Measurement likelihood helpers for quantized profiles
#'
#' @keywords internal
#' @noRd
NULL

extract_ll_vectors <- function(ak, an, p.err, p.conv) {
  v1_list <- list()
  v2_list <- list()

  idx <- 1
  for (n in an) {
    hak <- ak[ak <= n]
    v1_list[[idx]] <- dbinom(hak, n, prob = p.err, log = TRUE)
    v2_list[[idx]] <- dbinom(hak, n, prob = p.conv, log = TRUE)
    idx <- idx + 1
  }

  list(ll0 = unlist(v1_list), ll1 = unlist(v2_list))
}

extract_from_mixmat <- function(mixmat) {
  an <- .burst_get_mixmat_n(mixmat)
  ak <- .burst_get_mixmat_k(mixmat)
  counts_vec <- c()
  for (n in an) {
    hak <- ak[ak <= n]
    counts_vec <- c(counts_vec, as.numeric(unclass(mixmat[hak, n])))
  }
  list(ak = ak, an = an, counts = counts_vec)
}

prep_cell_ll01 <- function(mixmat, par) {
  d <- extract_from_mixmat(mixmat)
  ak <- d$ak
  an <- d$an
  c  <- d$counts

  lls <- extract_ll_vectors(ak, an, p.err = par$p.err, p.conv = par$p.conv)

  list(
    ll0 = as.numeric(lls$ll0),
    ll1 = as.numeric(lls$ll1),
    c   = as.numeric(c)
  )
}

ell_theta_from_ll01 <- function(theta, ll0, ll1, c) {
  if (theta <= 0) return(sum(c * ll0))
  if (theta >= 1) return(sum(c * ll1))
  a <- ll0 + log1p(-theta)
  b <- ll1 + log(theta)
  sum(c * .logspace_add_vec(a, b))
}

ell_grid_from_ll01 <- function(theta_vec, ll0, ll1, c) {
  vapply(theta_vec, ell_theta_from_ll01, numeric(1), ll0 = ll0, ll1 = ll1, c = c)
}

mode_from_exact_discrete <- function(N, ll0, ll1, c) {
  theta_grid <- (0:N) / N
  ell_exact  <- ell_grid_from_ll01(theta_grid, ll0, ll1, c)
  ell_exact[!is.finite(ell_exact)] <- -Inf
  k0 <- which.max(ell_exact) - 1L
  list(k0 = as.integer(k0), ell_exact = ell_exact)
}
