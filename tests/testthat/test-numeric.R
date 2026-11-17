test_that("full-support normalization preserves true zero probabilities", {
  expect_equal(.safe_norm(c(0, 2, 2)), c(0, 0.5, 0.5))
  expect_equal(.safe_norm(c(-1e-15, 0.4, 0.6)), c(0, 0.4, 0.6))
  expect_error(.safe_norm(c(0, 0)), "positive total mass")
  expect_error(.safe_norm(c(NA, 1)), "finite probabilities")
  expect_error(.safe_norm(c(-0.1, 1.1)), "negative probabilities")
})

test_that("extracting observed PMF entries never redistributes tail mass", {
  pmf <- data.frame(k = 0:2, prob = c(0.2, 0.2, 0.6))
  expect_equal(.pmf_to_vec(pmf, 1), c(0.2, 0.2))
  expect_equal(sum(.pmf_to_vec(pmf, 1)), 0.4)
})

test_that("log-sum-exp handles empty and impossible support", {
  expect_equal(.logsumexp_vec(numeric(0)), -Inf)
  expect_equal(.logsumexp_vec(c(-Inf, -Inf)), -Inf)
  expect_equal(.logsumexp_vec(c(-1000, -1001)), -1000 + log1p(exp(-1)))
})

test_that("a zero initial count has a defined PGF at zero", {
  expect_equal(G(0 + 0i, t.end = 0, d = 1, val = c(0, 1), prob = c(0.5, 0.5)),
               log(0.5) + 0i)
})

test_that("FFT convolution agrees with a direct finite convolution", {
  actual <- .convolve_pmf_fft(c(0.2, 0.8), c(0.3, 0.4, 0.3))
  expect_equal(actual, c(0.06, 0.32, 0.38, 0.24), tolerance = 1e-12)
})

test_that("FFT inversion preserves a point mass and rejects invalid PGFs", {
  pmf <- .fft_pmf_from_log(rep(0 + 0i, 16))
  
  expect_equal(
    as.numeric(pmf),
    c(1, rep(0, 15)),
    tolerance = 1e-12
  )
  
  expect_equal(attr(pmf, "fft_min_coefficient"), 0)
  expect_equal(attr(pmf, "fft_negative_fraction"), 0)
  expect_equal(attr(pmf, "fft_negative_mass"), 0)
  
  expect_error(
    .fft_pmf_from_log(c(NA_complex_, 0 + 0i)),
    "Invalid log-PGF"
  )
  expect_error(
    .fft_pmf_from_log(c(-Inf, -Inf)),
    "finite scaling point"
  )
})

test_that("adaptive support reports whether the tail criterion was met", {
  point_mass <- function(N) data.frame(k = seq_len(N) - 1L,
                                      prob = c(1, rep(0, N - 1L)))
  passed <- .ensure_tail_N(point_mass, N0 = 64, eps = 1e-8)
  expect_true(attr(passed$df, "tail_converged"))
  uniform <- function(N) data.frame(k = seq_len(N) - 1L, prob = rep(1 / N, N))
  stopped <- .ensure_tail_N(uniform, N0 = 64, eps = 1e-8, max_tries = 1)
  expect_false(attr(stopped$df, "tail_converged"))
})
