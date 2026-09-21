# Tests for WeightedSubtreeSwap (moveType 14)  — M-089
#
# WeightedSubtreeSwap extends GibbsSubtreeSwap by marginalising over branch
# fractions at each candidate swap partner.  Samples partner ∝ marginal
# weight, then samples bin and fraction, accepts/rejects via MH.
#
# Testing strategy:
#   - Structural: valid preorder output, tree length preserved, simplex
#     preserved, topology can change
#   - Functional: scalar parameters unchanged, logLik consistent with
#     recomputation, all branch lengths positive
#   - Chain: multi-move chain stays valid
#   - Acceptance: non-trivial acceptance rate
#   - Rejection: state unchanged on rejection

# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

.wss_pts <- function(seed = 42L, nTip = 6L) {
  set.seed(seed)
  tree  <- ape::rtree(nTip, rooted = FALSE)
  tree  <- Preorder(tree)
  nEdge <- nrow(tree$edge)
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
  list(dataPtr = dataPtr, statePtr = statePtr, nEdge = nEdge)
}

.try_wss <- function(pts, max_try = 100L, beta = 1.0) {
  for (i in seq_len(max_try)) {
    if (do_move_cpp(pts$dataPtr, pts$statePtr,
                    14L, 0L, 0.5, 0.5, 1L, beta))
      return(TRUE)
  }
  FALSE
}

.snap_br <- function(pts)
  get_mcmc_state(pts$statePtr)$relBrLengths + 0

# ---------------------------------------------------------------------------
# Structural tests
# ---------------------------------------------------------------------------

test_that("WeightedSubtreeSwap (moveType 14) can change topology", {
  pts    <- .wss_pts(seed = 301L)
  before <- get_mcmc_state(pts$statePtr)$edge
  changed <- .try_wss(pts)
  expect_true(changed)
  after  <- get_mcmc_state(pts$statePtr)$edge
  expect_false(identical(before, after))
})

test_that("WeightedSubtreeSwap preserves tree length", {
  pts      <- .wss_pts(seed = 302L)
  tl_before <- get_mcmc_state(pts$statePtr)$treeLength
  .try_wss(pts)
  tl_after  <- get_mcmc_state(pts$statePtr)$treeLength
  expect_equal(tl_before, tl_after, tolerance = 1e-12)
})

test_that("WeightedSubtreeSwap result is in canonical preorder", {
  pts <- .wss_pts(seed = 303L)
  .try_wss(pts)
  s   <- get_mcmc_state(pts$statePtr)
  edge <- s$edge
  # Check non-root internal nodes: first_as_child < first_as_parent
  root_nodes <- setdiff(unique(edge[, 1]), unique(edge[, 2]))
  internal   <- setdiff(unique(edge[, 1]), root_nodes)
  for (nd in internal) {
    first_child  <- min(which(edge[, 2] == nd))
    parent_rows  <- which(edge[, 1] == nd)
    if (length(parent_rows) == 0) next
    first_parent <- min(parent_rows)
    expect_true(first_child < first_parent,
                label = paste("node", nd, "preorder violated"))
  }
})

test_that("WeightedSubtreeSwap preserves simplex constraint (sum = 1)", {
  pts <- .wss_pts(seed = 304L)
  .try_wss(pts)
  br <- get_mcmc_state(pts$statePtr)$relBrLengths
  expect_equal(sum(br), 1.0, tolerance = 1e-12)
})

test_that("WeightedSubtreeSwap does not change scalar parameters", {
  pts    <- .wss_pts(seed = 305L)
  before <- get_mcmc_state(pts$statePtr)
  .try_wss(pts)
  after  <- get_mcmc_state(pts$statePtr)
  expect_equal(after$treeLength, before$treeLength, tolerance = 1e-14)
  expect_equal(after$rateLoss,   before$rateLoss,   tolerance = 1e-14)
  expect_equal(after$rateLogSd,  before$rateLogSd,  tolerance = 1e-14)
  expect_equal(after$rateNeo,    before$rateNeo,     tolerance = 1e-14)
  expect_equal(after$p,          before$p,           tolerance = 1e-14)
  expect_identical(as.integer(after$kPrime), as.integer(before$kPrime))
})

# ---------------------------------------------------------------------------
# Functional tests
# ---------------------------------------------------------------------------

test_that("WeightedSubtreeSwap updates logLik consistently", {
  pts <- .wss_pts(seed = 306L)
  .try_wss(pts)
  stored <- get_mcmc_state(pts$statePtr)$logLik
  recomp <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)
  expect_equal(stored, recomp, tolerance = 1e-8)
})

test_that("WeightedSubtreeSwap keeps all branch lengths positive", {
  pts <- .wss_pts(seed = 307L)
  .try_wss(pts)
  br <- get_mcmc_state(pts$statePtr)$relBrLengths
  expect_true(all(br > 0))
})

# ---------------------------------------------------------------------------
# Chain validity
# ---------------------------------------------------------------------------

test_that("chain of 30 WeightedSubtreeSwap moves stays valid", {
  pts   <- .wss_pts(seed = 308L, nTip = 7L)
  n_acc <- 0L
  for (i in seq_len(30L)) {
    acc <- do_move_cpp(pts$dataPtr, pts$statePtr,
                       14L, 0L, 0.5, 0.5, 1L, 1.0)
    if (acc) n_acc <- n_acc + 1L
  }
  s <- get_mcmc_state(pts$statePtr)
  expect_equal(sum(s$relBrLengths), 1.0, tolerance = 1e-12)
  expect_true(all(s$relBrLengths > 0))
  expect_equal(s$logLik,
               eval_full_loglik_cpp(pts$dataPtr, pts$statePtr),
               tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Acceptance rate
# ---------------------------------------------------------------------------

test_that("WeightedSubtreeSwap has non-trivial acceptance rate", {
  pts   <- .wss_pts(seed = 309L, nTip = 7L)
  n_acc <- 0L
  nTry  <- 80L
  for (i in seq_len(nTry)) {
    if (do_move_cpp(pts$dataPtr, pts$statePtr,
                    14L, 0L, 0.5, 0.5, 1L, 1.0))
      n_acc <- n_acc + 1L
  }
  expect_gt(n_acc, 0L)
  expect_lt(n_acc, nTry)
})

# ---------------------------------------------------------------------------
# Rejection leaves state unchanged
# ---------------------------------------------------------------------------

test_that("WeightedSubtreeSwap rejection leaves state unchanged", {
  pts <- .wss_pts(seed = 310L)
  found_rejection <- FALSE
  for (i in seq_len(200L)) {
    br_pre   <- .snap_br(pts)
    edge_pre <- get_mcmc_state(pts$statePtr)$edge
    ll_pre   <- get_mcmc_state(pts$statePtr)$logLik
    acc <- do_move_cpp(pts$dataPtr, pts$statePtr,
                       14L, 0L, 0.5, 0.5, 1L, 1.0)
    if (!acc) {
      br_post   <- .snap_br(pts)
      edge_post <- get_mcmc_state(pts$statePtr)$edge
      ll_post   <- get_mcmc_state(pts$statePtr)$logLik
      expect_equal(br_post, br_pre, tolerance = 1e-15)
      expect_identical(edge_post, edge_pre)
      expect_equal(ll_post, ll_pre, tolerance = 1e-15)
      found_rejection <- TRUE
      break
    }
  }
  expect_true(found_rejection, label = "should find at least one rejection")
})
