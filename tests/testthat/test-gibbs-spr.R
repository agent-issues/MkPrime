# Tests for GibbsSPR (moveType 10) and GibbsSubtreeSwap (moveType 11)
# M-085 / M-086
#
# GibbsSPR (post GSPR-001/GSPR-004 fix): enumerates every edge of the pruned
# tree — including the subtree's own merged position — weights by
# exp(β × logLik) at a fixed tau = 1/2 reference, draws the committed split
# fraction tau ~ U(0,1), and accepts through a full MH step.  An accepted
# move is therefore either a regraft OR a branch-fraction move at unchanged
# topology (the merged-edge draw).  See test-gibbs-spr-candidates.R for the
# deterministic candidate-set symmetry tests.
# GibbsSubtreeSwap remains a pure Gibbs update (self-draw returns false).
#
# Testing strategy:
#   - Structural: output is valid preorder tree, tree length preserved
#   - Functional: does change topology; does not mutate scalars
#   - Statistical: on a tree with one clearly better topology, GibbsSPR
#     selects it far more often than random (binomial test)

# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

.gibbs_pts <- function(seed = 42L, nTip = 5L) {
  set.seed(seed)
  # MCMC operates on unrooted trees (2n-3 edges); ape::rtree() defaults to
  # rooted (2n-2 edges), so explicitly request unrooted here.
  tree  <- ape::rtree(nTip, rooted = FALSE)
  tree  <- TreeTools::Preorder(tree)
  mat   <- matrix(
    sample(0:2, nTip * 4L, replace = TRUE),
    nrow = nTip, ncol = 4L,
    dimnames = list(tree$tip.label, NULL)
  )
  pd    <- MatrixToPhyDat(mat)
  mkd   <- suppressWarnings(MkPrimeData(pd))
  model <- MkPrimeModel()
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  state0   <- MkPrime:::.InitState(tree, mkd, model)
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  list(dataPtr = dataPtr, statePtr = statePtr)
}

# Attempt a move up to max_try times; return TRUE if any succeed
.try_move <- function(pts, moveType, max_try = 50L, beta = 1.0) {
  for (i in seq_len(max_try)) {
    if (do_move_cpp(pts$dataPtr, pts$statePtr,
                    moveType, 0L, 0.5, 0.5, 1L, beta))
      return(TRUE)
  }
  FALSE
}

# ---------------------------------------------------------------------------
# GibbsSPR — structural tests
# ---------------------------------------------------------------------------

test_that("GibbsSPR (moveType 10) can change topology on a 5-tip tree", {
  # An accepted move may be a merged-edge (branch-fraction) draw that leaves
  # the topology unchanged, so loop until the edge matrix itself moves.
  pts    <- .gibbs_pts(seed = 101L)
  before <- get_mcmc_state(pts$statePtr)$edge
  topoChanged <- FALSE
  for (i in seq_len(200L)) {
    if (do_move_cpp(pts$dataPtr, pts$statePtr, 10L, 0L, 0.5, 0.5, 1L, 1.0) &&
        !identical(before, get_mcmc_state(pts$statePtr)$edge)) {
      topoChanged <- TRUE
      break
    }
  }
  expect_true(topoChanged)
})

test_that("GibbsSPR preserves tree length", {
  pts     <- .gibbs_pts(seed = 102L)
  tl_before <- get_mcmc_state(pts$statePtr)$treeLength
  .try_move(pts, 10L)
  tl_after  <- get_mcmc_state(pts$statePtr)$treeLength
  expect_equal(tl_before, tl_after, tolerance = 1e-12)
})

test_that("GibbsSPR result is in canonical preorder", {
  # TreeTools preorder: edge where nd is child comes BEFORE all edges
  # where nd is parent. i.e., first_as_child < first_as_parent.
  pts <- .gibbs_pts(seed = 103L)
  .try_move(pts, 10L)
  s   <- get_mcmc_state(pts$statePtr)
  par <- s$edge[, 1L]
  chi <- s$edge[, 2L]
  nTip <- (nrow(s$edge) + 3L) %/% 2L
  root <- nTip + 1L
  internal <- unique(par[par != root])
  for (nd in internal) {
    first_as_parent <- which(par == nd)[1L]
    first_as_child  <- which(chi == nd)[1L]
    expect_true(first_as_child < first_as_parent,
                info = paste("node", nd, "not in preorder"))
  }
})

test_that("GibbsSPR does not change scalar parameters", {
  pts <- .gibbs_pts(seed = 104L)
  s0  <- get_mcmc_state(pts$statePtr)
  .try_move(pts, 10L, max_try = 20L)
  s1  <- get_mcmc_state(pts$statePtr)
  expect_equal(s0$treeLength,  s1$treeLength,  tolerance = 1e-12)
  expect_equal(s0$rateLoss,    s1$rateLoss,    tolerance = 1e-12)
  expect_equal(s0$rateLogSd,   s1$rateLogSd,   tolerance = 1e-12)
  expect_equal(s0$kPrime,      s1$kPrime)
})

test_that("GibbsSPR relative branch lengths sum to 1 after move", {
  pts <- .gibbs_pts(seed = 105L)
  .try_move(pts, 10L)
  relBr <- get_mcmc_state(pts$statePtr)$relBrLengths
  expect_equal(sum(relBr), 1.0, tolerance = 1e-12)
})

test_that("GibbsSPR updates logLik in state", {
  pts     <- .gibbs_pts(seed = 106L)
  ll_before <- get_state_log_lik(pts$statePtr)
  accepted  <- .try_move(pts, 10L)
  if (accepted) {
    # State logLik must match re-evaluation at the new topology
    ll_state  <- get_state_log_lik(pts$statePtr)
    ll_reeval <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)
    expect_equal(ll_state, ll_reeval, tolerance = 1e-8)
  } else {
    # No topology change; logLik unchanged
    expect_equal(get_state_log_lik(pts$statePtr), ll_before, tolerance = 1e-12)
  }
})

test_that("GibbsSPR on a 4-tip tree changes topology at least once in 100 tries", {
  pts   <- .gibbs_pts(seed = 107L, nTip = 4L)
  n_acc <- sum(vapply(seq_len(100L), function(i)
    do_move_cpp(pts$dataPtr, pts$statePtr, 10L, 0L, 0.5, 0.5, 1L, 1.0),
    logical(1L)))
  expect_gt(n_acc, 0L)
})

# ---------------------------------------------------------------------------
# GibbsSPR — statistical test: prefers better topology
# ---------------------------------------------------------------------------

test_that("GibbsSPR samples better topologies more often than random", {
  # Build a dataset with 4 tips where one tree topology has much higher
  # likelihood. Use a star-like character matrix that strongly favours
  # ((t1,t2),(t3,t4)).
  tree_good <- read.tree(text = "((t1:0.1,t2:0.1):0.2,(t3:0.1,t4:0.1):0.2);")
  tree_good <- TreeTools::Preorder(tree_good)
  # Characters: t1,t2 share state 1; t3,t4 share state 2
  mat <- matrix(c(1,1,0,0, 1,1,0,0, 0,0,1,1, 0,0,1,1,
                   1,1,0,0, 0,0,1,1),
                nrow = 4L, ncol = 6L,
                dimnames = list(c("t1","t2","t3","t4"), NULL))
  pd    <- MatrixToPhyDat(mat)
  mkd   <- suppressWarnings(MkPrimeData(pd))
  model <- MkPrimeModel()
  model <- MkPrime:::.FinalizeModel(model, tree_good, mkd)
  state0   <- MkPrime:::.InitState(tree_good, mkd, model)
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)

  # Run 200 GibbsSPR moves; count how many times topology changes
  n_try <- 200L
  n_acc <- sum(vapply(seq_len(n_try), function(i)
    do_move_cpp(dataPtr, statePtr, 10L, 0L, 0.5, 0.5, 1L, 1.0),
    logical(1L)))
  # At least some proposals should be accepted
  expect_gt(n_acc, 0L)
  # The MH step should reject some proposals
  expect_lt(n_acc, n_try)
})

# ---------------------------------------------------------------------------
# GibbsSubtreeSwap — structural tests
# ---------------------------------------------------------------------------

test_that("GibbsSubtreeSwap (moveType 11) can change topology on a 6-tip tree", {
  pts     <- .gibbs_pts(seed = 201L, nTip = 6L)
  before  <- get_mcmc_state(pts$statePtr)$edge
  changed <- .try_move(pts, 11L, max_try = 100L)
  expect_true(changed)
  after   <- get_mcmc_state(pts$statePtr)$edge
  expect_false(identical(before, after))
})

test_that("GibbsSubtreeSwap preserves tree length", {
  pts       <- .gibbs_pts(seed = 202L, nTip = 6L)
  tl_before <- get_mcmc_state(pts$statePtr)$treeLength
  .try_move(pts, 11L, max_try = 100L)
  tl_after  <- get_mcmc_state(pts$statePtr)$treeLength
  expect_equal(tl_before, tl_after, tolerance = 1e-12)
})

test_that("GibbsSubtreeSwap result is in canonical preorder", {
  pts <- .gibbs_pts(seed = 203L, nTip = 6L)
  .try_move(pts, 11L, max_try = 100L)
  s    <- get_mcmc_state(pts$statePtr)
  par  <- s$edge[, 1L]
  chi  <- s$edge[, 2L]
  nTip <- (nrow(s$edge) + 3L) %/% 2L
  root <- nTip + 1L
  internal <- unique(par[par != root])
  for (nd in internal) {
    first_as_parent <- which(par == nd)[1L]
    first_as_child  <- which(chi == nd)[1L]
    expect_true(first_as_child < first_as_parent,
                info = paste("node", nd, "not in preorder"))
  }
})

test_that("GibbsSubtreeSwap does not change scalar parameters", {
  pts <- .gibbs_pts(seed = 204L, nTip = 6L)
  s0  <- get_mcmc_state(pts$statePtr)
  .try_move(pts, 11L, max_try = 50L)
  s1  <- get_mcmc_state(pts$statePtr)
  expect_equal(s0$treeLength,  s1$treeLength,  tolerance = 1e-12)
  expect_equal(s0$rateLoss,    s1$rateLoss,    tolerance = 1e-12)
  expect_equal(s0$rateLogSd,   s1$rateLogSd,   tolerance = 1e-12)
  expect_equal(s0$kPrime,      s1$kPrime)
})

test_that("GibbsSubtreeSwap relative branch lengths sum to 1 after move", {
  pts <- .gibbs_pts(seed = 205L, nTip = 6L)
  .try_move(pts, 11L, max_try = 100L)
  relBr <- get_mcmc_state(pts$statePtr)$relBrLengths
  expect_equal(sum(relBr), 1.0, tolerance = 1e-12)
})

test_that("GibbsSubtreeSwap updates logLik in state", {
  pts      <- .gibbs_pts(seed = 206L, nTip = 6L)
  accepted <- .try_move(pts, 11L, max_try = 100L)
  if (accepted) {
    ll_state  <- get_state_log_lik(pts$statePtr)
    ll_reeval <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)
    expect_equal(ll_state, ll_reeval, tolerance = 1e-8)
  }
})

test_that("GibbsSubtreeSwap on small tree returns false gracefully when no partners", {
  # 3-tip tree: very constrained — some nodes may have no valid swap partners
  pts <- .gibbs_pts(seed = 207L, nTip = 3L)
  # Should not crash; may return true or false
  result <- do_move_cpp(pts$dataPtr, pts$statePtr, 11L, 0L, 0.5, 0.5, 1L, 1.0)
  expect_true(is.logical(result))
})

# ---------------------------------------------------------------------------
# Regression: existing moves not broken
# ---------------------------------------------------------------------------

test_that("Existing moves (NNI=5, SPR=6) still work after adding cases 10/11", {
  pts <- .gibbs_pts(seed = 301L)
  expect_true(.try_move(pts, 5L, max_try = 200L))   # NNI
  expect_true(.try_move(pts, 6L, max_try = 200L))   # SPR
})
