# Tests for M-135: adaptive thinning from observed autocorrelation

# --- Unit tests for .AdaptThinning() ---

test_that(".AdaptThinning returns currentThin when too few rows", {
  mat <- matrix(rnorm(30 * 3), 30, 3,
                dimnames = list(NULL, c("log_posterior", "tree_length",
                                        "rate_log_sd")))
  expect_equal(MkPrime:::.AdaptThinning(mat, 20L, 15L), 20L)
})

test_that(".AdaptThinning returns nMoves floor for iid samples", {
  set.seed(4782)
  n <- 200L
  mat <- matrix(rnorm(n * 3), n, 3,
                dimnames = list(NULL, c("log_posterior", "tree_length",
                                        "rate_log_sd")))
  # iid: ACT ≈ 1 sample → ACT_iter ≈ currentThin → newThin ≈ currentThin * log(2)

  # With nMoves = 15 and currentThin = 15, log(2)*15 ≈ 10.4 < 15 → returns 15
  result <- MkPrime:::.AdaptThinning(mat, 15L, 15L)
  expect_equal(result, 15L)
})

test_that(".AdaptThinning increases thin for autocorrelated samples", {
  set.seed(3194)
  n <- 300L
  nMoves <- 10L
  currentThin <- 10L
  # AR(1) with rho = 0.95: ACT ≈ (1+0.95)/(1-0.95) = 39 samples
  # ACT_iter = 39 * 10 = 390; newThin ≈ 390 * log(2) ≈ 270
  x <- numeric(n)
  x[1] <- rnorm(1)
  for (i in 2:n) x[i] <- 0.95 * x[i - 1] + rnorm(1) * sqrt(1 - 0.95^2)

  mat <- matrix(NA_real_, n, 3,
                dimnames = list(NULL, c("log_posterior", "tree_length",
                                        "rate_log_sd")))
  mat[, 1] <- x
  mat[, 2] <- x + rnorm(n, sd = 0.1)
  mat[, 3] <- rnorm(n)  # iid — should not drive the result

  result <- MkPrime:::.AdaptThinning(mat, currentThin, nMoves)
  expect_gt(result, currentThin * 5L)  # substantially larger than default
  expect_lte(result, 50L * nMoves)     # within cap
})

test_that(".AdaptThinning respects the 50x cap", {
  set.seed(8816)
  n <- 200L
  nMoves <- 10L
  currentThin <- 10L
  # Extreme autocorrelation: near-constant with tiny drift
  x <- cumsum(rnorm(n, sd = 0.001))
  mat <- matrix(x, n, 1, dimnames = list(NULL, c("tree_length")))
  result <- MkPrime:::.AdaptThinning(mat, currentThin, nMoves)
  expect_lte(result, 50L * nMoves)
})

test_that(".AdaptThinning returns currentThin for constant input", {
  n <- 100L
  mat <- matrix(5.0, n, 2,
                dimnames = list(NULL, c("log_posterior", "tree_length")))
  result <- MkPrime:::.AdaptThinning(mat, 20L, 15L)
  expect_equal(result, 20L)
})

test_that(".AdaptThinning excludes kPrime and branch columns", {
  set.seed(6021)
  n <- 200L
  # Autocorrelated kPrime should NOT drive thinning; iid scalars should
  x <- numeric(n)
  x[1] <- rnorm(1)
  for (i in 2:n) x[i] <- 0.99 * x[i - 1] + rnorm(1) * 0.1
  mat <- matrix(rnorm(n * 4), n, 4,
                dimnames = list(NULL, c("log_posterior", "tree_length",
                                        "kPrime_1", "br_1")))
  mat[, 3] <- x  # kPrime: highly autocorrelated
  mat[, 4] <- x  # branch: highly autocorrelated
  result <- MkPrime:::.AdaptThinning(mat, 15L, 15L)
  # Should be near the floor (driven by iid scalars), not by kPrime/br
  expect_lte(result, 25L)
})
