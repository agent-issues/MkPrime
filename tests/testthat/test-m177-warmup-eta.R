# Tests for M-177: Warmup ETA anchor fix

# Regression: .CheckStabilisation() returns nStableConsecutive correctly
# and the anchor-based ETA does not slide forward on every batch.
# Requires MkPrime.d2 (development companion package) to run.

skip_if_not_installed("MkPrime.d2")

test_that(".CheckStabilisation returns 0 when not enough data", {
  history <- rnorm(15)  # < 2 * windowSize (20)
  result <- MkPrime:::.CheckStabilisation(history, 0L)
  expect_false(result$stable)
  expect_identical(result$nStableConsecutive, 0L)
})

test_that(".CheckStabilisation increments counter for stable windows", {
  set.seed(1)
  # 20 draws from stationary normal — both windows should have similar means
  history <- rnorm(20, mean = 5, sd = 0.1)
  result <- MkPrime:::.CheckStabilisation(history, 0L)
  expect_false(result$stable)
  expect_gte(result$nStableConsecutive, 1L)
})

test_that(".CheckStabilisation resets to 0 for clearly shifted windows", {
  # prev window near 0, recent window near 100 → large z, reset to 0
  history <- c(rnorm(10, mean = 0, sd = 1), rnorm(10, mean = 100, sd = 1))
  result <- MkPrime:::.CheckStabilisation(history, 2L)
  expect_false(result$stable)
  expect_identical(result$nStableConsecutive, 0L)
})

test_that(".CheckStabilisation declares stable after 3 consecutive passes", {
  set.seed(42)
  history <- rnorm(20, mean = 5, sd = 0.05)
  result <- MkPrime:::.CheckStabilisation(history, 2L)
  expect_true(result$stable)
  expect_identical(result$nStableConsecutive, 3L)
})
