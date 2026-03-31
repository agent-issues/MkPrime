# tests/testthat/test-node-cl-cache.R
# M-121: Validate node-level CL cache for NNI and beta_simplex

library("MkPrime")
library("TreeTools")

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


# ---------------------------------------------------------------------------
# Regression test: partial-CL ascertainment correction matches full eval.
#
# Previously, cache_total_loglik (used by partial-CL paths for NNI,
# beta_simplex, Dirichlet) skipped the ascertainment correction for
# transformational partitions (type 1) and also passed an empty rates
# vector when ACRV was off.  This inflated the cached log-likelihood
# relative to full-eval, causing all full-eval moves (tree_length,
# rate_loss, rate_log_sd, kPrime) to be systematically rejected
# ("wall pattern").
# ---------------------------------------------------------------------------
test_that("partial-CL path does not freeze tree_length (ascertainment fix)", {
  # All-transformational dataset: no neomorphic characters
  set.seed(3847)
  tree <- ape::rtree(10L, rooted = FALSE)
  tree <- Preorder(tree)
  nChar <- 20L
  mat <- matrix(sample(0:2, 10L * nChar, replace = TRUE),
                nrow = 10L,
                dimnames = list(tree$tip.label, paste0("c", seq_len(nChar))))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)  # no neomorphic → all transformational

  model <- MkPrimeModel()
  mcmc <- MkPrimeMCMC(
    nIter = 600L, maxWarmup = 150L, minWarmup = 150L, thin = 3L,
    autoTune = FALSE, nRuns = 1L
  )
  result <- suppressWarnings(RunMkPrime(data = mkd, tree = tree,
                                         model = model, mcmc = mcmc))

  tl <- result$samples[, "tree_length"]
  n_unique <- length(unique(tl))
  # With the bug, tree_length would be nearly completely frozen
  # (ascertainment gap ~8 LL units → acceptance rate ~0.03%).
  # Expect at least 5 distinct values from normal mixing;
  # the buggy code would typically produce 1-2.
  expect_true(n_unique >= 5,
              info = paste("tree_length has only", n_unique,
                           "unique values; expected >= 5"))
})


# ---------------------------------------------------------------------------
# Regression test M-145: slice sampler must invalidate node CL cache.
#
# Before the fix, slice_scalar_impl accepted without setting
# nodeCL.valid = false.  A subsequent NNI using partial-CL evaluation
# would read stale CLs (computed with the old parameter value), producing
# incorrect log-likelihoods.  The diagnostic drift counter catches this:
# every 100 iterations, do_move_impl compares state->logLik against a
# fresh full evaluation and increments diagDriftCount on mismatch.
# ---------------------------------------------------------------------------
test_that("slice sampler invalidates CL cache (M-145 regression)", {
  skip_if_not_installed("TreeSearch")
  dat <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  mkd <- suppressWarnings(MkPrimeData(dat))
  model <- MkPrimeModel()
  tree <- ape::rtree(length(dat), tip.label = names(dat))
  tree <- TreeTools::Preorder(tree)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)

  mcmcData <- MkPrime:::.InitMcmcData(mkd, model)
  state <- MkPrime:::.InitState(tree, mkd, model)
  chainState <- MkPrime:::.InitMcmcChain(state)
  fill_partition_cache(mcmcData, chainState)
  allocate_cl_workspace(mcmcData, chainState)

  hasNeo <- any(mkd$type == "neomorphic")
  nEdge <- nrow(tree$edge)
  mcmcCfg <- MkPrimeMCMC(
    nIter = 100L, minWarmup = 50L,
    gibbsSpr = FALSE, gibbsSubtreeSwap = FALSE
  )
  moves <- MkPrime:::.BuildMoves(nEdge, sum(mkd$type == "transformational"),
                                  hasNeo, mcmcCfg)
  moveTypeCodes <- vapply(
    moves, function(m) MkPrime:::.kMoveTypes[[m$name]], integer(1L)
  )
  moveWeights <- vapply(moves, `[[`, numeric(1), "weight")
  nMoves <- length(moves)
  transIdx <- which(mkd$type == "transformational")
  transIdx0 <- if (length(transIdx)) transIdx - 1L else integer(0)

  scaleTunings <- matrix(0.5, 1, nMoves)
  sliceParamCodes <- vapply(moves, function(m) m$sliceParamIdx %||% 0L,
                            integer(1L))
  sliceWidths <- matrix(1.0, 1, nMoves)
  jointRhos <- matrix(0.0, 1, nMoves)
  moveIntPars <- integer(nMoves)

  # Run enough iterations for the 100-iter drift check to fire multiple times
  result <- run_mcmc_batch_cpp(
    mcmcData, list(chainState), 1.0,
    moveTypeCodes, transIdx0, sliceParamCodes, moveWeights,
    scaleTunings, 10, 1L, moveIntPars, sliceWidths, jointRhos,
    500L, 1L, 500L, 10L,
    hasNeo, nEdge
  )

  # Before M-145 fix, drift > 0 because slice sampler left stale CL cache.
  expect_equal(result$diag_counters[["drift"]], 0L,
               info = "slice sampler should invalidate nodeCL after accepting")
})
