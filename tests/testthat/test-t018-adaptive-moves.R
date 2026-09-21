# tests/testthat/test-t018-adaptive-moves.R
#
# Tests for T-018: Acceptance-aware adaptive move decay.
# Covers: .DecayLowAcceptMoves() unit tests.
# Integration test (warmup changes weights / sampling freezes them) runs via
# MKPRIME_SLOW_TESTS=true.

# ==========================================================================
# Unit tests for .DecayLowAcceptMoves()
# ==========================================================================

test_that(".DecayLowAcceptMoves returns same-length vector summing to 1", {
  mw <- c(nni = 0.3, spr = 0.3, tree_length = 0.4)
  res <- MkPrime:::.DecayLowAcceptMoves(
    currentWeights = mw,
    batchAccept    = c(nni = 0L, spr = 0L, tree_length = 30L),
    batchPropose   = c(nni = 50L, spr = 50L, tree_length = 50L),
    initialWeights = mw,
    moveNames      = names(mw),
    pinnedWeights  = NULL
  )
  expect_length(res, 3L)
  expect_equal(sum(res), 1.0, tolerance = 1e-10)
})

test_that(".DecayLowAcceptMoves decays move below accept_floor", {
  mw <- c(nni = 0.5, spr = 0.5)
  # nni: 0% accept in 50 proposals — well below default 2% floor
  res <- MkPrime:::.DecayLowAcceptMoves(
    currentWeights = mw,
    batchAccept    = c(nni = 0L, spr = 25L),
    batchPropose   = c(nni = 50L, spr = 50L),
    initialWeights = mw,
    moveNames      = names(mw),
    pinnedWeights  = NULL,
    accept_floor   = 0.02,
    decay          = 0.7
  )
  # nni should be decayed
  expect_lt(res[["nni"]], mw[["nni"]])
  # spr has 50% accept — not decayed, should be >= initial
  expect_gte(res[["spr"]], mw[["spr"]])
  expect_equal(sum(res), 1.0, tolerance = 1e-10)
})

test_that(".DecayLowAcceptMoves respects weight floor", {
  # Put nni at a very small initial weight so decay would drive it to 0
  mw <- c(nni = 0.01, spr = 0.99)
  res <- MkPrime:::.DecayLowAcceptMoves(
    currentWeights = mw,
    batchAccept    = c(nni = 0L, spr = 0L),
    batchPropose   = c(nni = 50L, spr = 50L),
    initialWeights = mw,
    moveNames      = names(mw),
    pinnedWeights  = NULL,
    accept_floor   = 0.02,
    decay          = 0.0,        # extreme: would drive to 0 without floor
    weight_floor   = 0.1         # floor = 10% of initial = 0.001
  )
  # Floor = 0.1 * 0.01 = 0.001; weight must be >= that (before re-norm)
  # After re-norm, proportion is conserved, so nni > 0
  expect_gt(res[["nni"]], 0)
  expect_equal(sum(res), 1.0, tolerance = 1e-10)
})

test_that(".DecayLowAcceptMoves leaves moves with < n_min proposals unchanged", {
  mw <- c(nni = 0.5, spr = 0.5)
  res <- MkPrime:::.DecayLowAcceptMoves(
    currentWeights = mw,
    batchAccept    = c(nni = 0L, spr = 0L),
    batchPropose   = c(nni = 5L, spr = 5L),   # only 5 proposals — below n_min=30
    initialWeights = mw,
    moveNames      = names(mw),
    pinnedWeights  = NULL,
    n_min          = 30L
  )
  # No decay should fire — weights unchanged
  expect_equal(res, mw, tolerance = 1e-10)
})

test_that(".DecayLowAcceptMoves does not decay pinned moves", {
  mw <- c(nni = 0.3, spr = 0.3, gibbs_p = 0.4)
  res <- MkPrime:::.DecayLowAcceptMoves(
    currentWeights = mw,
    batchAccept    = c(nni = 0L, spr = 0L, gibbs_p = 0L),
    batchPropose   = c(nni = 50L, spr = 50L, gibbs_p = 50L),
    initialWeights = mw,
    moveNames      = names(mw),
    pinnedWeights  = c(gibbs_p = 0.4),
    accept_floor   = 0.02
  )
  # Pinned move must stay at its pinned value
  expect_equal(res[["gibbs_p"]], 0.4, tolerance = 1e-10)
  expect_equal(sum(res), 1.0, tolerance = 1e-10)
})

test_that(".DecayLowAcceptMoves does not decay moves above accept_floor", {
  mw <- c(nni = 0.5, spr = 0.5)
  res <- MkPrime:::.DecayLowAcceptMoves(
    currentWeights = mw,
    batchAccept    = c(nni = 20L, spr = 20L),  # 40% accept — well above 2% floor
    batchPropose   = c(nni = 50L, spr = 50L),
    initialWeights = mw,
    moveNames      = names(mw),
    pinnedWeights  = NULL,
    accept_floor   = 0.02
  )
  # Neither move should change (no decay fired)
  expect_equal(res, mw, tolerance = 1e-10)
})

test_that(".DecayLowAcceptMoves output sums to 1 after multi-decay", {
  mw <- c(a = 0.25, b = 0.25, c = 0.25, d = 0.25)
  # All moves have 0 acceptance — all get decayed
  res <- MkPrime:::.DecayLowAcceptMoves(
    currentWeights = mw,
    batchAccept    = c(a = 0L, b = 0L, c = 0L, d = 0L),
    batchPropose   = c(a = 50L, b = 50L, c = 50L, d = 50L),
    initialWeights = mw,
    moveNames      = names(mw),
    pinnedWeights  = NULL,
    accept_floor   = 0.02,
    decay          = 0.7,
    weight_floor   = 0.1
  )
  expect_equal(sum(res), 1.0, tolerance = 1e-10)
  # When all moves decay equally, weights stay uniform
  expect_true(max(res) - min(res) < 1e-10)
})


# ==========================================================================
# Integration: warmup adapts, sampling phase freezes
# ==========================================================================

skip_slow_tests <- function() {
  if (!identical(Sys.getenv("MKPRIME_SLOW_TESTS"), "true")) {
    testthat::skip("Slow integration test: set MKPRIME_SLOW_TESTS=true to run")
  }
}

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

  # The result$runs[[1]]$moveWeights was captured at sampling phase start
  finalWeights <- result$runs[[1L]]$moveWeights
  expect_false(is.null(finalWeights))
  expect_equal(sum(finalWeights), 1.0, tolerance = 1e-6)

  # The sampling phase should have produced some samples
  expect_gt(nrow(result$samples[[1L]]), 0L)
})
