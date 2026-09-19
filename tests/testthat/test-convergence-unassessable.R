# CONV-001: a convergence window in which every monitored scalar is constant
# must read as "cannot assess", not as "converged".
#
# R-hat and ESS are NA exactly when a chain is constant. Reducing an all-NA
# vector with `max(na.rm = TRUE)` yields -Inf and `min(na.rm = TRUE)` yields
# Inf, and both satisfy the stopping rules (`maxRhat <= threshold`,
# `minEss >= threshold`). So the worst mixing outcome there is -- a sampler
# that has not moved at all -- used to be reported as the best one, precisely
# where the user is relying on automatic stopping rather than reading traces.

test_that(".MinOrNA / .MaxOrNA report NA rather than a passing sentinel", {
  expect_identical(.MinOrNA(c(NA_real_, NA_real_)), NA_real_)
  expect_identical(.MaxOrNA(c(NA_real_, NA_real_)), NA_real_)
  expect_identical(.MinOrNA(numeric(0)), NA_real_)
  expect_identical(.MaxOrNA(numeric(0)), NA_real_)

  # Partially observed vectors still reduce normally, ignoring the NAs.
  expect_identical(.MinOrNA(c(NA, 5, 9)), 5)
  expect_identical(.MaxOrNA(c(NA, 5, 9)), 9)
})

test_that("an all-constant window does not report convergence", {
  paramNames <- c("log_posterior", "log_likelihood",
                  "tree_length", "rate_log_sd")
  n <- 60L

  # A frozen chain: every column constant across every retained draw.
  stuck <- matrix(rep(c(-100, -120, 2.5, 0.8), each = n), nrow = n,
                  dimnames = list(NULL, paramNames))
  runs <- list(list(samples = stuck, saved_idx = n),
               list(samples = stuck, saved_idx = n))

  res <- .CheckConvergence(runs, paramNames,
                           mcmc = list(minEss = 200, maxRhat = 1.01))

  expect_false(res$converged)
  expect_true(is.na(res$minEss))
  expect_true(is.na(res$maxRhat))
})

test_that("a mixing window is still assessed normally", {
  paramNames <- c("log_posterior", "log_likelihood",
                  "tree_length", "rate_log_sd")
  n <- 200L
  set.seed(1L)

  draws <- function() {
    matrix(c(rnorm(n, -100), rnorm(n, -120), rnorm(n, 2.5, 0.1),
             rnorm(n, 0.8, 0.05)), nrow = n,
           dimnames = list(NULL, paramNames))
  }
  runs <- list(list(samples = draws(), saved_idx = n),
               list(samples = draws(), saved_idx = n))

  res <- .CheckConvergence(runs, paramNames,
                           mcmc = list(minEss = 10, maxRhat = 1.5))

  # The point is that a real window still yields finite diagnostics -- the
  # guard must not turn every window into NA.
  expect_true(is.finite(res$minEss))
  expect_true(is.finite(res$maxRhat))
  expect_true(res$converged)
})
