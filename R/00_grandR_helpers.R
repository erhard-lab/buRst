# Internal helpers for the binomial mixture model.
# These functions are not exported.

.burst_model_par <- function(
    ntr = 0.1,
    p.err = 4e-4,
    p.conv = 0.05,
    p.mconv = 0.05,
    shape = 5
) {
  structure(
    list(
      ntr = ntr,
      p.err = p.err,
      p.conv = p.conv,
      p.mconv = p.mconv,
      shape = shape,
      logLik = NA_real_
    ),
    class = "grmodel"
  )
}


.burst_rbinommix <- function(n, size, par) {
  
  n1 <- if (par$ntr <= 0) {
    round(-par$ntr)
  } else {
    stats::rbinom(
      1,
      size = n,
      prob = par$ntr
    )
  }
  
  c(
    stats::rbinom(
      n1,
      size = size,
      prob = par$p.conv
    ),
    stats::rbinom(
      n - n1,
      size = size,
      prob = par$p.err
    )
  )
}


.burst_create_mix_matrix <- function(n.vector, par) {
  
  n <- length(n.vector)
  
  if (n == 0L) {
    stop("n.vector must not be empty.")
  }
  
  m <- matrix(
    0,
    nrow = n + 1L,
    ncol = n,
    dimnames = list(
      as.character(0:n),
      as.character(seq_len(n))
    )
  )
  
  for (i in seq_len(n)) {
    
    if (n.vector[i] > 0) {
      
      draws <- .burst_rbinommix(
        n = n.vector[i],
        size = i,
        par = par
      )
      
      tab <- table(draws)
      
      m[names(tab), i] <- as.numeric(tab)
    }
  }
  
  m <- m[
    rowSums(m) > 0,
    colSums(m) > 0,
    drop = FALSE
  ]
  
  class(m) <- "MixMatrix"
  
  m
}


.burst_get_mixmat_k <- function(mixmat) {
  as.integer(rownames(mixmat))
}


.burst_get_mixmat_n <- function(mixmat) {
  as.integer(colnames(mixmat))
}