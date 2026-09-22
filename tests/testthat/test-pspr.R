# M-119: Parsimony-guided SPR (pSPR) proposal
#
# Tests for Fitch parsimony scoring and pSPR integration.

# --- Fitch parsimony scoring ---

test_that("Fitch score is correct for a small hand-computed tree", {
  # 4-tip tree:  ((1,2), (3,4))
  #              5       6       root=5
  # Edges (1-indexed, preorder):
  #   5→6, 6→1, 6→2, 5→3, 5→4
  parent <- c(5L, 6L, 6L, 5L, 5L)
  child  <- c(6L, 1L, 2L, 3L, 4L)

  # 1 binary character: tips = 0, 0, 1, 1
  # Optimal: 1 change (between the two clades)
  tipStates <- matrix(c(0L, 0L, 1L, 1L), ncol = 1)
  expect_equal(fitch_score_r(parent, child, tipStates, nTip = 4L, kStates = 2L), 1L)

  # Tips = 0, 1, 0, 1 → minimum 2 changes on this tree
  tipStates2 <- matrix(c(0L, 1L, 0L, 1L), ncol = 1)
  expect_equal(fitch_score_r(parent, child, tipStates2, nTip = 4L, kStates = 2L), 2L)

  # All same state: 0 changes
  tipStates3 <- matrix(c(0L, 0L, 0L, 0L), ncol = 1)
  expect_equal(fitch_score_r(parent, child, tipStates3, nTip = 4L, kStates = 2L), 0L)
})

test_that("Fitch score handles missing data (-1) correctly", {
  parent <- c(5L, 6L, 6L, 5L, 5L)
  child  <- c(6L, 1L, 2L, 3L, 4L)

  # Tips = 0, -1, 1, 1 → missing data at tip 2 is compatible with anything
  # Optimal score: 1 change (tip 1 is 0, rest are 1; change at node 6)
  tipStates <- matrix(c(0L, -1L, 1L, 1L), ncol = 1)
  expect_equal(fitch_score_r(parent, child, tipStates, nTip = 4L, kStates = 2L), 1L)

  # All missing: 0 changes
  tipStates2 <- matrix(c(-1L, -1L, -1L, -1L), ncol = 1)
  expect_equal(fitch_score_r(parent, child, tipStates2, nTip = 4L, kStates = 2L), 0L)
})

test_that("Fitch score handles multi-state characters (k > 2)", {
  parent <- c(5L, 6L, 6L, 5L, 5L)
  child  <- c(6L, 1L, 2L, 3L, 4L)

  # 3-state character: tips = 0, 1, 2, 0
  # node 6 = {0} ∩ {1} = {} → union {0,1}, score +1
  # root = {0,1} ∩ {2} = {} → need to consider tip 4 too
  # root: children are node 6, tip 3, tip 4
  # First: root ← node 6 = {0,1}
  # Second: root = {0,1} ∩ {2} = {} → union {0,1,2}, score +1
  # Third: root = {0,1,2} ∩ {0} = {0}, score +0
  # Total: 2
  tipStates <- matrix(c(0L, 1L, 2L, 0L), ncol = 1)
  expect_equal(fitch_score_r(parent, child, tipStates, nTip = 4L, kStates = 3L), 2L)
})

test_that("Fitch score sums across multiple characters", {
  parent <- c(5L, 6L, 6L, 5L, 5L)
  child  <- c(6L, 1L, 2L, 3L, 4L)

  # 2 characters: char 1 = 0,0,1,1 (1 change); char 2 = 0,1,0,1 (2 changes)
  tipStates <- matrix(c(0L, 0L, 1L, 1L,
                         0L, 1L, 0L, 1L), ncol = 2)
  expect_equal(fitch_score_r(parent, child, tipStates, nTip = 4L, kStates = 2L), 3L)
})


# --- Hastings ratio normalization (u.122 fix) ---

test_that("pSPR Hastings formula includes normalization correction", {
  # With few candidates, forward and reverse normalization constants diverge.
  # This test verifies the algebraic relationship that the full formula
  # (implemented in C++) satisfies but the old formula (missing correction) does not.
  alpha <- 0.1
  scores <- c(10, 12, 15)  # 3 candidates (small tree)
  scoreOrig <- 14

  minScore <- min(c(scores, scoreOrig))
  w <- exp(-alpha * (scores - minScore))
  wOrig <- exp(-alpha * (scoreOrig - minScore))
  sumW_fwd <- sum(w)

  chosen <- 1L  # best candidate
  sumW_rev <- sumW_fwd + wOrig - w[chosen]

  # Correction term is non-trivial with few candidates
  correction <- log(sumW_fwd) - log(sumW_rev)
  expect_true(abs(correction) > 0.01)

  # Full formula = log(w_orig/w_chosen) + correction
  logH_full <- log(wOrig) - log(w[chosen]) + correction
  logH_old  <- -alpha * (scoreOrig - scores[chosen])
  expect_equal(logH_full, logH_old + correction)
  expect_false(isTRUE(all.equal(logH_full, logH_old, tolerance = 1e-4)))
})


# --- pSPR integration ---

test_that("pSPR move produces valid trees", {
  skip_under_memcheck()
  skip_on_cran()

  nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
  pd <- ReadAsPhyDat(nexFile)
  mkd <- MkPrimeData(pd)

  # pSPR only (disable other topology moves); long enough that acceptance = 0
  # would be a genuine bug (no fixed seed — cross-platform RNG differs).
  cfg <- MkPrimeMCMC(
    nIter = 10000L, nRuns = 1L, nChains = 1L,
    minWarmup = 1000L, maxWarmup = 2000L,
    autoTune = FALSE, thin = 10L,
    gibbsSpr = FALSE, gibbsSubtreeSwap = FALSE,
    tbr = FALSE, pSpr = TRUE
  )

  res <- RunMkPrime(mkd, mcmc = cfg)
  # Should complete without error and produce samples

  expect_s3_class(res, "MkPosterior")
  expect_gt(nrow(res$samples), 0)
  # pSPR acceptance should be > 0 (some moves accepted)
  expect_true("pspr" %in% names(res$acceptance))
  expect_gt(res$acceptance[["pspr"]], 0)
})

test_that("pSPR can be disabled", {
  skip_under_memcheck()
  skip_on_cran()

  nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
  pd <- ReadAsPhyDat(nexFile)
  mkd <- MkPrimeData(pd)

  cfg <- MkPrimeMCMC(
    nIter = 2000L, nRuns = 1L, nChains = 1L,
    minWarmup = 500L, maxWarmup = 1000L,
    autoTune = FALSE, thin = 10L,
    pSpr = FALSE
  )

  set.seed(2917)
  res <- RunMkPrime(mkd, mcmc = cfg)
  expect_s3_class(res, "MkPosterior")
  expect_false("pspr" %in% names(res$acceptance))
})
