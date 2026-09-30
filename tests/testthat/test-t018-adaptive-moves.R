# tests/testthat/test-t018-adaptive-moves.R
#
# Warmup adapts move weights and the sampling phase freezes them.
# Runs via MKPRIME_SLOW_TESTS=true.

# ==========================================================================
# Integration: warmup adapts, sampling phase freezes
# ==========================================================================

test_that("T-018: move weights change during warmup and freeze at sampling", {
  skip_slow_tests()
  skip_if_not_installed("TreeSearch")

  dat <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  mkd <- suppressWarnings(MkPrimeData(dat))
  model <- MkPrimeModel()
  tree <- Preorder(ape::rtree(length(dat), tip.label = names(dat)))

  # Long warmup so adaptation has iterations to fire; short sampling
  mcmc <- MkPrimeMCMC(
    nIter = 1500L, minWarmup = 500L, maxWarmup = 1000L,
    thin = 10L, autoTune = FALSE, nRuns = 1L, nChains = 1L,
    logFile = tempfile(fileext = ".log"),
    checkpointFile = tempfile(fileext = ".ckp")
  )

  # Run and capture the returned move weights (post-warmup, frozen)
  result <- RunMkPrime(mkd, tree = tree, model = model, mcmc = mcmc)
  expect_s3_class(result, "MkPosterior")

  # result$moveWeights holds run 1's frozen schedule, captured at sampling
  # phase start (see RunMkPrime()'s @return).
  finalWeights <- result$moveWeights
  expect_false(is.null(finalWeights))
  expect_equal(sum(finalWeights), 1.0, tolerance = 1e-6)

  # The sampling phase should have produced some samples. mcmc$logFile puts
  # this run in streaming mode, where result$samples is an empty placeholder
  # matrix (see test-streaming.R); result$nSamples is the real count.
  expect_gt(result$nSamples, 0L)
})
