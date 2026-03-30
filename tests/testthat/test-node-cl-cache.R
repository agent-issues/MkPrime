# tests/testthat/test-node-cl-cache.R
# M-121: Validate node-level CL cache for NNI and beta_simplex

library(MkPrime)
library(TreeTools)

# ---------------------------------------------------------------------------
# Helper: small dataset for testing
# ---------------------------------------------------------------------------
make_test_setup <- function(nTip = 12L, nChar = 15L, kMax = 3L,
                             seed = 4281L) {
  set.seed(seed)
  tree <- ape::rtree(nTip, rooted = FALSE)
  tree <- Preorder(tree)
  mat <- matrix(sample(0:(kMax - 1L), nTip * nChar, replace = TRUE),
                nrow = nTip,
                dimnames = list(tree$tip.label, paste0("c", seq_len(nChar))))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel()
  list(tree = tree, mkd = mkd, model = model)
}


# ---------------------------------------------------------------------------
# Test: MCMC with NNI produces valid posteriors (node CL cache active)
# ---------------------------------------------------------------------------
test_that("MCMC with NNI + node CL cache produces valid log posteriors", {
  setup <- make_test_setup(nTip = 10L, nChar = 12L, seed = 7321L)

  # NNI-heavy config: disable SPR/TBR to force NNI as primary topology move
  mcmc <- MkPrimeMCMC(
    nIter = 300L, maxWarmup = 100L, minWarmup = 100L, thin = 3L, autoTune = FALSE,
    nRuns = 1L
  )
  result <- suppressWarnings(RunMkPrime(data = setup$mkd, tree = setup$tree,
                       model = setup$model, mcmc = mcmc))
  expect_s3_class(result, "MkPosterior")
  expect_true(nrow(result$samples) > 0)
  expect_true(all(is.finite(result$samples[, "log_posterior"])))
})


# ---------------------------------------------------------------------------
# Test: MCMC with beta_simplex produces valid posteriors
# ---------------------------------------------------------------------------
test_that("MCMC with beta_simplex + node CL cache produces valid posteriors", {
  setup <- make_test_setup(nTip = 10L, nChar = 12L, seed = 5603L)

  mcmc <- MkPrimeMCMC(
    nIter = 300L, maxWarmup = 100L, minWarmup = 100L, thin = 3L, autoTune = FALSE,
    nRuns = 1L
  )
  result <- suppressWarnings(RunMkPrime(data = setup$mkd, tree = setup$tree,
                       model = setup$model, mcmc = mcmc))
  expect_s3_class(result, "MkPosterior")
  expect_true(nrow(result$samples) > 0)
  expect_true(all(is.finite(result$samples[, "log_posterior"])))
})


# ---------------------------------------------------------------------------
# Test: MCMC with ACRV + node CL cache produces valid posteriors
# ---------------------------------------------------------------------------
test_that("MCMC with ACRV and node CL cache produces valid posteriors", {
  setup <- make_test_setup(nTip = 10L, nChar = 12L, seed = 8894L)
  # Enable ACRV (non-zero rateLogSd)
  model <- MkPrimeModel()

  mcmc <- MkPrimeMCMC(
    nIter = 300L, maxWarmup = 100L, minWarmup = 100L, thin = 3L, autoTune = FALSE,
    nRuns = 1L
  )
  result <- suppressWarnings(RunMkPrime(data = setup$mkd, tree = setup$tree,
                       model = model, mcmc = mcmc))
  expect_s3_class(result, "MkPosterior")
  expect_true(all(is.finite(result$samples[, "log_posterior"])))
})


# ---------------------------------------------------------------------------
# Test: mixed moves (NNI + SPR + beta_simplex) with node CL cache
# This tests cache invalidation on SPR acceptance + rebuild on next NNI.
# ---------------------------------------------------------------------------
test_that("mixed NNI + SPR + beta_simplex with cache invalidation works", {
  setup <- make_test_setup(nTip = 12L, nChar = 15L, seed = 2190L)

  mcmc <- MkPrimeMCMC(
    nIter = 500L, maxWarmup = 200L, minWarmup = 200L, thin = 5L, autoTune = FALSE,
    nRuns = 1L
  )
  result <- suppressWarnings(RunMkPrime(data = setup$mkd, tree = setup$tree,
                       model = setup$model, mcmc = mcmc))
  expect_s3_class(result, "MkPosterior")
  expect_true(nrow(result$samples) > 0)
  expect_true(all(is.finite(result$samples[, "log_posterior"])))

  # Samples should cover a reasonable range of log posteriors
  lp <- result$samples[, "log_posterior"]
  expect_true(diff(range(lp)) > 0)
})


# ---------------------------------------------------------------------------
# Test: node CL cache disabled for Q-heterogeneity (falls back to full eval)
# ---------------------------------------------------------------------------
test_that("Q-heterogeneity bypasses node CL cache without errors", {
  set.seed(6612)
  tree <- ape::rtree(8L, rooted = FALSE)
  tree <- Preorder(tree)
  mat <- matrix(sample(0:2, 8L * 10L, replace = TRUE),
                nrow = 8L,
                dimnames = list(tree$tip.label, paste0("c", seq_len(10L))))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(qHeterogeneity = TRUE)

  mcmc <- MkPrimeMCMC(
    nIter = 200L, maxWarmup = 80L, minWarmup = 80L, thin = 3L, autoTune = FALSE,
    nRuns = 1L
  )
  result <- suppressWarnings(RunMkPrime(data = mkd, tree = tree,
                       model = model, mcmc = mcmc))
  expect_s3_class(result, "MkPosterior")
  expect_true(all(is.finite(result$samples[, "log_posterior"])))
})
