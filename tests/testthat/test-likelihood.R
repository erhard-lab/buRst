test_that("flat profiles retain the unconditional stationary count likelihood", {
  counts <- c(0, 0, 1, 2)
  profiles <- lapply(counts, function(n) if (n == 0) numeric(0) else rep(0, n + 1))
  B <- 2
  km <- 0.7
  d <- 0.3
  likelihood <- make_loglik_from_gp(profiles, counts, t.label = 2)
  expected <- sum(stats::dnbinom(counts, size = km / d, mu = B * km / d, log = TRUE))
  expect_equal(likelihood(B, km, d, N0 = 512), expected, tolerance = 1e-6)
})

test_that("a profile likelihood equals direct marginalization of binomial survival", {
  # No synthesis: old RNA starts at four and survives with probability 1/2.
  # All probability above the largest observed count (two) must remain omitted.
  counts <- c(0, 1, 2)
  profiles <- list(numeric(0), c(0.4, 1), c(0.7, 0.2, 1.2))
  likelihood <- make_loglik_from_gp(
    profiles, counts, t.label = 2, steady = FALSE, t.diff = 2, init = rep(4, 5)
  )
  # Profiles are centered at their best-supported allocation.
  expected <- sum(stats::dbinom(counts, size = 4, prob = 0.5, log = TRUE)) -
    (0.7 - 0.2)
  expect_equal(likelihood(B = 1, km = 0, d = log(2) / 2, N0 = 64),
               expected, tolerance = 1e-7)
})

test_that("component likelihood retains probability above the observed maximum", {
  likelihood <- make_loglik_from_truth(
    new_counts = c(0, 0), old_counts = c(0, 1),
    t.label = 2, steady = FALSE, t.diff = 2, init = rep(4, 5)
  )
  expected <- sum(stats::dbinom(c(0, 1), size = 4, prob = 0.5, log = TRUE))
  expect_equal(likelihood(B = 1, km = 0, d = log(2) / 2, N0 = 64),
               expected, tolerance = 1e-7)
})

test_that("uninformative zero-count cells are valid likelihood observations", {
  likelihood <- make_loglik_from_gp(list(numeric(0), numeric(0)), c(0, 0))
  expected <- 2 * stats::dnbinom(0, size = 0.7 / 0.3, mu = 2 * 0.7 / 0.3, log = TRUE)
  expect_equal(likelihood(B = 2, km = 0.7, d = 0.3, N0 = 512),
               expected, tolerance = 1e-6)
})

test_that("impossible allocations are not replaced by uniform probabilities", {
  expect_equal(cond_from_uncond_log(1, c(1, 0), c(1, 0), c(1, 0)), c(0, 0))
})

test_that("invalid counts and unsupported intervals fail before fitting", {
  expect_error(make_loglik_from_truth(c(0, 0.5), c(1, 1)), "integer counts")
  expect_error(make_loglik_from_truth(c(0, -1), c(1, 1)), "integer counts")
  expect_error(make_loglik_from_gp(list(c(0, 0)), 1, steady = FALSE,
                                   t.label = 2, t.diff = 4, init = 1),
               "requires t.diff = t.label", fixed = TRUE)
  expect_error(make_loglik_from_truth(0, 1, steady = FALSE,
                                      t.label = 2, t.diff = 4, init = 1),
               "requires t.diff = t.label", fixed = TRUE)
})
