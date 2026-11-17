test_that("unequal reference cell numbers do not change initial parameters", {
  counts <- c(4, 8, 12, 16)
  new <- c(1, 2, 3, 4)
  reference <- c(8, 12, 16, 20)
  initial <- .initial_fit_parameters(counts, new, reference, 2, FALSE, 2)
  repeated <- .initial_fit_parameters(counts, new, rep(reference, 7), 2, FALSE, 2)
  expect_equal(initial, repeated)
  expect_equal(initial$d, round(-log(mean(counts - new) / mean(reference)) / 2, 5))
})

test_that("steady and non-steady degradation starts use the expected ratios", {
  expect_equal(fit_d(0.25, 10, 2), -log(0.75) / 2)
  expect_equal(fit_d(0.25, 10, 2, t0 = 2, f0 = 14), -log(7.5 / 14) / 2)
  expect_equal(fit_d(0, 10, 2), 0)
  expect_equal(fit_d(1, 10, 2), Inf)
  expect_true(is.na(fit_d(0.25, 10, 2, t0 = 2, f0 = 0)))
  expect_error(fit_d(0.25, 10, 2, t0 = 4, f0 = 14), "requires t.diff = t.label", fixed = TRUE)
})
