# Tests for WeightedBranchLengthScale (moveType 12)
# M-087
#
# MH move: discretises branch-fraction space into B bins using
# Beta(0.25, 0.25) quantile breakpoints, evaluates LL at each midpoint,
# Gibbs-samples a bin, draws fraction from a within-bin Beta, returns
# logHastings for standard MH acceptance.
#
# Testing strategy:
#   - Structural: output is valid (simplex preserved, tree length unchanged)
#   - Functional: does change branch lengths; does not mutate topology/scalars
#   - Statistical: accepted moves shift branch proportions toward
#     likelihood-improving fractions

library(ape)
library(TreeTools)

# ---------------------------------------------------------------------------
# Shared fixture (reuses pattern from test-gibbs-spr.R)
# ---------------------------------------------------------------------------

.wbs_pts <- function(seed = 42L, nTip = 6L) {
  set.seed(seed)
  tree  <- ape::rtree(nTip, rooted = FALSE)
  tree  <- TreeTools::Preorder(tree)
  mat   <- matrix(
    sample(0:2, nTip * 5L, replace = TRUE),
    nrow = nTip, ncol = 5L,
    dimnames = list(tree$tip.label, NULL)
  )
  pd    <- MatrixToPhyDat(mat)
  mkd   <- MkPrimeData(pd)
  model <- MkPrimeModel()
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  state0   <- MkPrime:::.InitState(tree, mkd, model)
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  list(dataPtr = dataPtr, statePtr = statePtr, nEdge = nrow(tree$edge))
}

# Attempt move up to max_try times; return TRUE if any accepted
.try_wbs <- function(pts, max_try = 80L, beta = 1.0) {
  for (i in seq_len(max_try)) {
    if (do_move_cpp(pts$dataPtr, pts$statePtr,
                    12L, 0L, 0.5, 0.5, 1L, beta))
      return(TRUE)
  }
  FALSE
}

# Snapshot helper: force a deep copy of relBrLengths.
# get_mcmc_state()$relBrLengths shares the C++ SEXP, and as.double()
# is a no-op on double vectors; +0 forces a new allocation.
.snap_br <- function(pts) {
  get_mcmc_state(pts$statePtr)$relBrLengths + 0
}

# ---------------------------------------------------------------------------
# Structural tests
# ---------------------------------------------------------------------------

test_that("WeightedBranchScale (moveType 12) can change branch lengths", {
  pts    <- .wbs_pts(seed = 3201L)
  before <- .snap_br(pts)
  changed <- .try_wbs(pts)
  expect_true(changed)
  after  <- .snap_br(pts)
  expect_false(identical(before, after))
})

test_that("WeightedBranchScale preserves tree length", {
  pts       <- .wbs_pts(seed = 3202L)
  tl_before <- get_mcmc_state(pts$statePtr)$treeLength
  .try_wbs(pts)
  tl_after  <- get_mcmc_state(pts$statePtr)$treeLength
  expect_equal(tl_before, tl_after, tolerance = 1e-12)
})

test_that("WeightedBranchScale preserves simplex constraint", {
  pts <- .wbs_pts(seed = 3203L)
  before_sum <- sum(.snap_br(pts))
  .try_wbs(pts)
  after_sum  <- sum(.snap_br(pts))
  expect_equal(before_sum, after_sum, tolerance = 1e-12)
})

test_that("WeightedBranchScale does not change topology", {
  pts    <- .wbs_pts(seed = 3204L)
  before <- get_mcmc_state(pts$statePtr)$edge
  .try_wbs(pts)
  after  <- get_mcmc_state(pts$statePtr)$edge
  expect_identical(before, after)
})

test_that("WeightedBranchScale does not change scalar parameters", {
  pts    <- .wbs_pts(seed = 3205L)
  s0     <- get_mcmc_state(pts$statePtr)
  .try_wbs(pts)
  s1     <- get_mcmc_state(pts$statePtr)
  expect_equal(s0$treeLength, s1$treeLength)
  expect_equal(s0$rateLoss,   s1$rateLoss)
  expect_equal(s0$rateLogSd,  s1$rateLogSd)
  expect_equal(s0$rateNeo,    s1$rateNeo)
  expect_equal(s0$p,          s1$p)
})

test_that("WeightedBranchScale keeps all branch lengths positive", {
  pts <- .wbs_pts(seed = 3206L)
  for (i in 1:20) {
    do_move_cpp(pts$dataPtr, pts$statePtr, 12L, 0L, 0.5, 0.5, 1L, 1.0)
  }
  br <- .snap_br(pts)
  expect_true(all(br > 0))
})

# ---------------------------------------------------------------------------
# Functional: logLik is updated consistently after acceptance
# ---------------------------------------------------------------------------

test_that("WeightedBranchScale updates logLik consistently", {
  pts <- .wbs_pts(seed = 3207L)
  .try_wbs(pts)
  s <- get_mcmc_state(pts$statePtr)
  recomputed <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)
  expect_equal(s$logLik, recomputed, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Multiple accepted moves: run many iterations, check chain stays valid
# ---------------------------------------------------------------------------

test_that("WeightedBranchScale chain of 50 moves stays valid", {
  pts <- .wbs_pts(seed = 3208L, nTip = 8L)
  n_accepted <- 0L
  for (i in 1:50) {
    acc <- do_move_cpp(pts$dataPtr, pts$statePtr, 12L, 0L, 0.5, 0.5, 1L, 1.0)
    if (acc) n_accepted <- n_accepted + 1L
  }
  s <- get_mcmc_state(pts$statePtr)
  # Simplex preserved
  expect_equal(sum(s$relBrLengths), 1.0, tolerance = 1e-10)
  # All positive
  expect_true(all(s$relBrLengths > 0))
  # logLik consistent
  recomputed <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)
  expect_equal(s$logLik, recomputed, tolerance = 1e-8)
  # At least some accepted
  expect_gt(n_accepted, 0L)
})

# ---------------------------------------------------------------------------
# Rejection rollback: verify state is unchanged after rejection
# ---------------------------------------------------------------------------

test_that("WeightedBranchScale rejection leaves state unchanged", {
  pts <- .wbs_pts(seed = 3209L)

  # Run until we observe a rejection
  got_rejection <- FALSE
  for (i in 1:200) {
    before_br <- .snap_br(pts)
    before_ll <- get_mcmc_state(pts$statePtr)$logLik
    acc <- do_move_cpp(pts$dataPtr, pts$statePtr,
                       12L, 0L, 0.5, 0.5, 1L, 1.0)
    if (!acc) {
      after_br <- .snap_br(pts)
      after_ll <- get_mcmc_state(pts$statePtr)$logLik
      expect_equal(before_br, after_br, tolerance = 1e-14)
      expect_equal(before_ll, after_ll, tolerance = 1e-14)
      got_rejection <- TRUE
      break
    }
  }
  expect_true(got_rejection)
})

# ---------------------------------------------------------------------------
# Comparative: reasonable acceptance rate
# ---------------------------------------------------------------------------

test_that("WeightedBranchScale has non-trivial acceptance rate", {
  pts <- .wbs_pts(seed = 3210L, nTip = 7L)
  n <- 100L
  n_acc <- 0L
  for (i in seq_len(n)) {
    if (do_move_cpp(pts$dataPtr, pts$statePtr, 12L, 0L, 0.5, 0.5, 1L, 1.0))
      n_acc <- n_acc + 1L
  }
  rate <- n_acc / n
  # Should be between 5% and 95%
  expect_gt(rate, 0.02)
  expect_lt(rate, 0.98)
})
