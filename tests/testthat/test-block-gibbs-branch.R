# Tests for block Gibbs branch-length sweep (moveType 15)
# M-054 reframed
#
# Random-permutation-scan MH-within-Gibbs over all edge pairs.
# Each pair uses the same bin-based approximate conditional as
# WeightedBranchLengthScale (M-087), accepted/rejected independently.
#
# Testing strategy:
#   - Structural: simplex preserved, tree length unchanged, topology unchanged
#   - Functional: does change branch lengths, logLik updated consistently
#   - Statistical: reasonable within-sweep acceptance (not 0%, not 100%)
#   - Scheduler: dim field plumbed through correctly

library(ape)
library(TreeTools)


# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

.bgb_pts <- function(seed = 7711L, nTip = 6L) {
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
  set_branch_bins(dataPtr, 10L)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  list(dataPtr = dataPtr, statePtr = statePtr, nEdge = nrow(tree$edge))
}

# Attempt block sweep via do_move_cpp (moveType 15)
.try_bgb <- function(pts, max_try = 20L, beta = 1.0) {
  for (i in seq_len(max_try)) {
    if (do_move_cpp(pts$dataPtr, pts$statePtr,
                    15L, 0L, 0.5, 0.5, 1L, beta))
      return(TRUE)
  }
  FALSE
}

# Deep copy of relBrLengths
.snap_br <- function(pts) {
  get_mcmc_state(pts$statePtr)$relBrLengths + 0
}


# ---------------------------------------------------------------------------
# Structural tests
# ---------------------------------------------------------------------------

test_that("BlockGibbsBranch (moveType 15) can change branch lengths", {
  pts    <- .bgb_pts(seed = 7712L)
  before <- .snap_br(pts)
  changed <- .try_bgb(pts)
  expect_true(changed)
  after  <- .snap_br(pts)
  expect_false(identical(before, after))
})

test_that("BlockGibbsBranch preserves tree length", {
  pts       <- .bgb_pts(seed = 7713L)
  tl_before <- get_mcmc_state(pts$statePtr)$treeLength
  .try_bgb(pts)
  tl_after  <- get_mcmc_state(pts$statePtr)$treeLength
  expect_equal(tl_before, tl_after, tolerance = 1e-12)
})

test_that("BlockGibbsBranch preserves simplex sum", {
  pts <- .bgb_pts(seed = 7714L)
  before_sum <- sum(.snap_br(pts))
  .try_bgb(pts)
  after_sum  <- sum(.snap_br(pts))
  expect_equal(before_sum, after_sum, tolerance = 1e-10)
})

test_that("BlockGibbsBranch does not change topology", {
  pts    <- .bgb_pts(seed = 7715L)
  before <- get_mcmc_state(pts$statePtr)$edge
  .try_bgb(pts)
  after  <- get_mcmc_state(pts$statePtr)$edge
  expect_identical(before, after)
})

test_that("BlockGibbsBranch does not change scalar parameters", {
  pts <- .bgb_pts(seed = 7716L)
  s0  <- get_mcmc_state(pts$statePtr)
  .try_bgb(pts)
  s1  <- get_mcmc_state(pts$statePtr)
  expect_equal(s0$treeLength, s1$treeLength)
  expect_equal(s0$rateLoss,   s1$rateLoss)
  expect_equal(s0$rateLogSd,  s1$rateLogSd)
  expect_equal(s0$rateNeo,    s1$rateNeo)
  expect_equal(s0$p,          s1$p)
})

test_that("BlockGibbsBranch keeps all branch lengths positive", {
  pts <- .bgb_pts(seed = 7717L)
  for (i in 1:10) {
    do_move_cpp(pts$dataPtr, pts$statePtr, 15L, 0L, 0.5, 0.5, 1L, 1.0)
  }
  br <- .snap_br(pts)
  expect_true(all(br > 0))
})


# ---------------------------------------------------------------------------
# Functional: logLik updated consistently
# ---------------------------------------------------------------------------

test_that("BlockGibbsBranch updates logLik consistently", {
  pts <- .bgb_pts(seed = 7718L)
  .try_bgb(pts)
  s <- get_mcmc_state(pts$statePtr)
  recomputed <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)
  expect_equal(s$logLik, recomputed, tolerance = 1e-8)
})

test_that("BlockGibbsBranch chain of 20 sweeps stays valid", {
  pts <- .bgb_pts(seed = 7719L, nTip = 8L)
  n_accepted <- 0L
  for (i in 1:20) {
    acc <- do_move_cpp(pts$dataPtr, pts$statePtr,
                       15L, 0L, 0.5, 0.5, 1L, 1.0)
    if (acc) n_accepted <- n_accepted + 1L
  }
  s <- get_mcmc_state(pts$statePtr)
  expect_equal(sum(s$relBrLengths), 1.0, tolerance = 1e-10)
  expect_true(all(s$relBrLengths > 0))
  recomputed <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)
  expect_equal(s$logLik, recomputed, tolerance = 1e-8)
  expect_gt(n_accepted, 0L)
})


# ---------------------------------------------------------------------------
# Tempering: works with heated chains (beta < 1)
# ---------------------------------------------------------------------------

test_that("BlockGibbsBranch works with heated chain (beta = 0.3)", {
  pts <- .bgb_pts(seed = 7720L)
  changed <- .try_bgb(pts, beta = 0.3)
  expect_true(changed)
  s <- get_mcmc_state(pts$statePtr)
  expect_true(all(s$relBrLengths > 0))
  expect_equal(sum(s$relBrLengths), 1.0, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# Scheduler dim field
# ---------------------------------------------------------------------------

test_that(".BuildMoves sets dim = nEdge for block_gibbs_branch", {
  tree  <- ape::rtree(8L, rooted = FALSE)
  tree  <- TreeTools::Preorder(tree)
  nEdge <- nrow(tree$edge)
  mcmc  <- MkPrimeMCMC(nIter = 100L, blockGibbsBranch = TRUE)
  moves <- MkPrime:::.BuildMoves(nEdge, nTrans = 3L, hasNeo = TRUE, mcmc)
  bgb   <- Filter(function(m) m$name == "block_gibbs_branch", moves)
  expect_length(bgb, 1L)
  expect_equal(bgb[[1]]$dim, nEdge)
})

test_that(".BuildMoves sets dim = 1 for all standard moves", {
  tree  <- ape::rtree(6L, rooted = FALSE)
  tree  <- TreeTools::Preorder(tree)
  nEdge <- nrow(tree$edge)
  mcmc  <- MkPrimeMCMC(nIter = 100L)
  moves <- MkPrime:::.BuildMoves(nEdge, nTrans = 3L, hasNeo = TRUE, mcmc)
  dims  <- vapply(moves, function(m) m$dim %||% 1L, integer(1L))
  expect_true(all(dims == 1L))
})


# ---------------------------------------------------------------------------
# .AdaptMoveWeights dim-adjusted scoring
# ---------------------------------------------------------------------------

test_that(".AdaptMoveWeights uses dim in score computation", {
  # Two moves: A (dim=1, cheap, moderate acceptance) vs B (dim=10, expensive,
  # similar acceptance). Without dim adjustment, A dominates. With dim
  # adjustment, B should score competitively.
  nMoves <- 2L
  moveNames <- c("moveA", "moveB")
  currentWeights <- c(0.5, 0.5)
  names(currentWeights) <- moveNames

  acceptCount <- c(moveA = 30L, moveB = 25L)
  proposeCount <- c(moveA = 100L, moveB = 100L)
  # Move B costs 10x as much wall time
  moveTimeNs <- c(moveA = 1e8, moveB = 1e9)

  # Without dim: A should dominate (same accept rate, 10x cheaper)
  w_nodim <- MkPrime:::.AdaptMoveWeights(
    currentWeights, acceptCount, proposeCount, moveTimeNs, moveNames,
    moveDim = c(1L, 1L), pinnedWeights = NULL, warmupProgress = 0.8
  )
  expect_gt(w_nodim[["moveA"]], w_nodim[["moveB"]])

  # With dim=10 for B: B should gain weight relative to no-dim case
  w_dim <- MkPrime:::.AdaptMoveWeights(
    currentWeights, acceptCount, proposeCount, moveTimeNs, moveNames,
    moveDim = c(1L, 10L), pinnedWeights = NULL, warmupProgress = 0.8
  )
  # B's relative share should increase when dim is accounted for
  ratio_nodim <- w_nodim[["moveB"]] / w_nodim[["moveA"]]
  ratio_dim   <- w_dim[["moveB"]] / w_dim[["moveA"]]
  expect_gt(ratio_dim, ratio_nodim)
})
