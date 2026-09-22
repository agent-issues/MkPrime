# Tests for WeightedSPR (moveType 13)  — M-088
#
# WeightedSPR extends GibbsSPR by marginalising over branch fractions at each
# candidate reattachment position.  It samples topology ∝ marginal weight,
# then samples a bin and fraction within the chosen topology, and
# accepts/rejects via MH.
#
# Testing strategy:
#   - Structural: valid preorder output, tree length preserved, simplex
#     preserved, topology can change
#   - Functional: scalar parameters unchanged, logLik consistent with
#     recomputation, all branch lengths positive
#   - Chain: multi-move chain stays valid
#   - Acceptance: non-trivial acceptance rate

# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

.wspr_pts <- function(seed = 42L, nTip = 6L) {
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

# Attempt moveType up to max_try times; return TRUE on first success
.try_wspr <- function(pts, max_try = 100L, beta = 1.0) {
  for (i in seq_len(max_try)) {
    if (do_move_cpp(pts$dataPtr, pts$statePtr,
                    13L, 0L, 0.5, 0.5, 1L, beta))
      return(TRUE)
  }
  FALSE
}

# Deep-copy relBrLengths (avoid SEXP aliasing; see M-087 note)
.snap_br <- function(pts)
  get_mcmc_state(pts$statePtr)$relBrLengths + 0

# ---------------------------------------------------------------------------
# Structural tests
# ---------------------------------------------------------------------------

test_that("WeightedSPR (moveType 13) can change topology", {
  pts    <- .wspr_pts(seed = 201L)
  before <- get_mcmc_state(pts$statePtr)$edge
  changed <- .try_wspr(pts)
  expect_true(changed)
  after  <- get_mcmc_state(pts$statePtr)$edge
  expect_false(identical(before, after))
})

test_that("WeightedSPR preserves tree length", {
  pts      <- .wspr_pts(seed = 202L)
  tl_before <- get_mcmc_state(pts$statePtr)$treeLength
  .try_wspr(pts)
  tl_after  <- get_mcmc_state(pts$statePtr)$treeLength
  expect_equal(tl_before, tl_after, tolerance = 1e-12)
})

test_that("WeightedSPR result is in canonical preorder", {
  # TreeTools preorder: edge where nd is child comes BEFORE all edges
  # where nd is parent. i.e., first_as_child < first_as_parent.
  # Root (nTip + 1) is excluded: it never appears as a child.
  pts <- .wspr_pts(seed = 203L)
  .try_wspr(pts)
  s   <- get_mcmc_state(pts$statePtr)
  edge <- s$edge
  internal <- setdiff(unique(edge[, 1]), unique(edge[, 2]))
  internal <- setdiff(unique(edge[, 1]), internal)  # parents that ARE children
  for (nd in internal) {
    first_child  <- min(which(edge[, 2] == nd))
    parent_rows  <- which(edge[, 1] == nd)
    if (length(parent_rows) == 0) next
    first_parent <- min(parent_rows)
    expect_true(first_child < first_parent,
                label = paste("node", nd, "preorder violated"))
  }
})

test_that("WeightedSPR preserves simplex constraint (sum = 1)", {
  pts <- .wspr_pts(seed = 204L)
  .try_wspr(pts)
  br <- get_mcmc_state(pts$statePtr)$relBrLengths
  expect_equal(sum(br), 1.0, tolerance = 1e-12)
})

test_that("WeightedSPR does not change scalar parameters", {
  pts    <- .wspr_pts(seed = 205L)
  before <- get_mcmc_state(pts$statePtr)
  .try_wspr(pts)
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

test_that("WeightedSPR updates logLik consistently with recomputation", {
  pts <- .wspr_pts(seed = 206L)
  .try_wspr(pts)
  stored <- get_mcmc_state(pts$statePtr)$logLik
  recomp <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)
  expect_equal(stored, recomp, tolerance = 1e-8)
})

test_that("WeightedSPR keeps all branch lengths positive", {
  pts <- .wspr_pts(seed = 207L)
  .try_wspr(pts)
  br <- get_mcmc_state(pts$statePtr)$relBrLengths
  expect_true(all(br > 0))
})

# ---------------------------------------------------------------------------
# Chain validity
# ---------------------------------------------------------------------------

test_that("chain of 30 WeightedSPR moves stays valid", {
  pts   <- .wspr_pts(seed = 208L, nTip = 7L)
  n_acc <- 0L
  for (i in seq_len(30L)) {
    acc <- do_move_cpp(pts$dataPtr, pts$statePtr,
                       13L, 0L, 0.5, 0.5, 1L, 1.0)
    if (acc) n_acc <- n_acc + 1L
  }
  s <- get_mcmc_state(pts$statePtr)
  expect_equal(sum(s$relBrLengths), 1.0, tolerance = 1e-12)
  expect_true(all(s$relBrLengths > 0))
  # logLik consistent with recomputation
  expect_equal(s$logLik,
               eval_full_loglik_cpp(pts$dataPtr, pts$statePtr),
               tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Acceptance rate
# ---------------------------------------------------------------------------

test_that("WeightedSPR has non-trivial acceptance rate", {
  pts   <- .wspr_pts(seed = 209L, nTip = 7L)
  n_acc <- 0L
  nTry  <- 80L
  for (i in seq_len(nTry)) {
    if (do_move_cpp(pts$dataPtr, pts$statePtr,
                    13L, 0L, 0.5, 0.5, 1L, 1.0))
      n_acc <- n_acc + 1L
  }
  # At least 1 acceptance (2%) and not all (< 100%)
  expect_gt(n_acc, 0L)
  expect_lt(n_acc, nTry)
})

# ---------------------------------------------------------------------------
# Rejection leaves state unchanged
# ---------------------------------------------------------------------------

test_that("WeightedSPR rejection leaves state unchanged", {
  pts <- .wspr_pts(seed = 210L)
  # Run until we get at least one rejection
  snap_edge <- get_mcmc_state(pts$statePtr)$edge
  snap_br   <- .snap_br(pts)
  snap_ll   <- get_mcmc_state(pts$statePtr)$logLik
  found_rejection <- FALSE
  for (i in seq_len(200L)) {
    br_pre   <- .snap_br(pts)
    edge_pre <- get_mcmc_state(pts$statePtr)$edge
    ll_pre   <- get_mcmc_state(pts$statePtr)$logLik
    acc <- do_move_cpp(pts$dataPtr, pts$statePtr,
                       13L, 0L, 0.5, 0.5, 1L, 1.0)
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
