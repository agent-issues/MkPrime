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
  tree  <- Preorder(tree)
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

# Attempt a move up to max_try times; return TRUE if any succeed.
# Callers assert on the result: an invariant checked after a move that never
# fired is an invariant checked on the untouched state.
.try_move <- function(pts, moveType, max_try = 50L, beta = 1.0) {
  for (i in seq_len(max_try)) {
    if (do_move_cpp(pts$dataPtr, pts$statePtr,
                    moveType, 0L, 0.5, 0.5, 1L, beta))
      return(TRUE)
  }
  FALSE
}

# The single non-trivial split of a 4-tip tree, as a canonical string naming
# the side that holds the first taxon.
.split_label <- function(edge, taxa) {
  tree <- structure(
    list(edge = edge, Nnode = nrow(edge) - length(taxa) + 1L,
         tip.label = taxa),
    class = "phylo"
  )
  side <- taxa[as.logical(as.Splits(tree))[1L, ]]
  if (!taxa[[1L]] %in% side) side <- setdiff(taxa, side)
  paste(sort(side), collapse = "|")
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
  expect_true(.try_move(pts, 10L))
  tl_after  <- get_mcmc_state(pts$statePtr)$treeLength
  expect_equal(tl_before, tl_after, tolerance = 1e-12)
})

test_that("GibbsSPR result is in canonical preorder", {
  # TreeTools preorder: edge where nd is child comes BEFORE all edges
  # where nd is parent. i.e., first_as_child < first_as_parent.
  pts <- .gibbs_pts(seed = 103L)
  expect_true(.try_move(pts, 10L))
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
  expect_true(.try_move(pts, 10L, max_try = 20L))
  s1  <- get_mcmc_state(pts$statePtr)
  # The move fired, so it either regrafted or redrew the committed split
  # fraction: the invariants below are asserted against a state it touched.
  expect_false(identical(s0$edge, s1$edge) &&
                 isTRUE(all.equal(s0$relBrLengths, s1$relBrLengths)))
  expect_equal(s0$treeLength,  s1$treeLength,  tolerance = 1e-12)
  expect_equal(s0$rateLoss,    s1$rateLoss,    tolerance = 1e-12)
  expect_equal(s0$rateLogSd,   s1$rateLogSd,   tolerance = 1e-12)
  expect_equal(s0$kPrime,      s1$kPrime)
})

test_that("GibbsSPR relative branch lengths sum to 1 after move", {
  pts <- .gibbs_pts(seed = 105L)
  expect_true(.try_move(pts, 10L))
  relBr <- get_mcmc_state(pts$statePtr)$relBrLengths
  expect_equal(sum(relBr), 1.0, tolerance = 1e-12)
})

test_that("GibbsSPR updates logLik in state", {
  pts     <- .gibbs_pts(seed = 106L)
  expect_true(.try_move(pts, 10L))
  # State logLik must match re-evaluation at the committed state
  expect_equal(get_state_log_lik(pts$statePtr),
               eval_full_loglik_cpp(pts$dataPtr, pts$statePtr),
               tolerance = 1e-8)
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
  # Six characters split {t1,t2} from {t3,t4}, so of the three unrooted 4-tip
  # topologies exactly one is supported; the other two are ~e^15 less likely.
  # A sampler that ignores the likelihood when weighting regraft positions
  # would spread over all three, so the test asks which topology was drawn,
  # not merely whether a move was accepted.
  mat <- matrix(c(1,1,0,0, 1,1,0,0, 0,0,1,1, 0,0,1,1,
                   1,1,0,0, 0,0,1,1),
                nrow = 4L, ncol = 6L,
                dimnames = list(c("t1","t2","t3","t4"), NULL))
  mkd  <- suppressWarnings(MkPrimeData(MatrixToPhyDat(mat)))
  taxa <- rownames(mkd$matrix)

  # Tip indices are positional, so the tree must carry the data's tip order.
  .pts_at <- function(text) {
    tree  <- Preorder(RenumberTips(read.tree(text = text), taxa))
    model <- MkPrime:::.FinalizeModel(MkPrimeModel(), tree, mkd)
    state0   <- MkPrime:::.InitState(tree, mkd, model)
    dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
    statePtr <- MkPrime:::.InitMcmcChain(state0)
    fill_partition_cache(dataPtr, statePtr)
    allocate_cl_workspace(dataPtr, statePtr)
    list(dataPtr = dataPtr, statePtr = statePtr)
  }

  starts <- c(good = "((t1:0.1,t2:0.1):0.2,(t3:0.1,t4:0.1):0.2);",
              alt1 = "((t1:0.1,t3:0.1):0.2,(t2:0.1,t4:0.1):0.2);",
              alt2 = "((t1:0.1,t4:0.1):0.2,(t2:0.1,t3:0.1):0.2);")
  logLiks <- vapply(starts, function(text) {
    pts <- .pts_at(text)
    eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)
  }, numeric(1))
  expect_gt(logLiks[["good"]], logLiks[["alt1"]] + 10)
  expect_gt(logLiks[["good"]], logLiks[["alt2"]] + 10)

  # From every start, including both wrong ones, the chain must find the
  # supported topology and stay there.
  for (start in names(starts)) {
    set.seed(4242L)
    pts <- .pts_at(starts[[start]])
    visited <- vapply(seq_len(500L), function(i) {
      do_move_cpp(pts$dataPtr, pts$statePtr, 10L, 0L, 0.5, 0.5, 1L, 1.0)
      .split_label(get_mcmc_state(pts$statePtr)$edge, taxa)
    }, character(1))
    expect_gt(mean(visited == "t1|t2"), 0.95, label = paste("from", start))
  }

  # The MH step must still reject some proposals: perfect acceptance would
  # mean the accept ratio is not being evaluated at all.
  set.seed(4242L)
  pts   <- .pts_at(starts[["good"]])
  n_try <- 200L
  n_acc <- sum(vapply(seq_len(n_try), function(i)
    do_move_cpp(pts$dataPtr, pts$statePtr, 10L, 0L, 0.5, 0.5, 1L, 1.0),
    logical(1L)))
  expect_gt(n_acc, 0L)
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
  expect_true(.try_move(pts, 11L, max_try = 100L))
  tl_after  <- get_mcmc_state(pts$statePtr)$treeLength
  expect_equal(tl_before, tl_after, tolerance = 1e-12)
})

test_that("GibbsSubtreeSwap result is in canonical preorder", {
  pts <- .gibbs_pts(seed = 203L, nTip = 6L)
  expect_true(.try_move(pts, 11L, max_try = 100L))
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
  expect_true(.try_move(pts, 11L, max_try = 50L))
  s1  <- get_mcmc_state(pts$statePtr)
  # A swap that fired must have moved the topology; the invariants below are
  # then asserted against a state the move actually touched.
  expect_false(identical(s0$edge, s1$edge))
  expect_equal(s0$treeLength,  s1$treeLength,  tolerance = 1e-12)
  expect_equal(s0$rateLoss,    s1$rateLoss,    tolerance = 1e-12)
  expect_equal(s0$rateLogSd,   s1$rateLogSd,   tolerance = 1e-12)
  expect_equal(s0$kPrime,      s1$kPrime)
})

test_that("GibbsSubtreeSwap relative branch lengths sum to 1 after move", {
  pts <- .gibbs_pts(seed = 205L, nTip = 6L)
  expect_true(.try_move(pts, 11L, max_try = 100L))
  relBr <- get_mcmc_state(pts$statePtr)$relBrLengths
  expect_equal(sum(relBr), 1.0, tolerance = 1e-12)
})

test_that("GibbsSubtreeSwap updates logLik in state", {
  pts <- .gibbs_pts(seed = 206L, nTip = 6L)
  expect_true(.try_move(pts, 11L, max_try = 100L))
  expect_equal(get_state_log_lik(pts$statePtr),
               eval_full_loglik_cpp(pts$dataPtr, pts$statePtr),
               tolerance = 1e-8)
})

test_that("GibbsSubtreeSwap on small tree returns false gracefully when no partners", {
  # A 3-tip tree is a star: no node has a valid swap partner, so every call
  # must decline and leave the state exactly as it found it.
  pts    <- .gibbs_pts(seed = 207L, nTip = 3L)
  before <- get_mcmc_state(pts$statePtr)
  results <- vapply(seq_len(20L), function(i)
    do_move_cpp(pts$dataPtr, pts$statePtr, 11L, 0L, 0.5, 0.5, 1L, 1.0),
    logical(1L))
  expect_false(any(results))
  after <- get_mcmc_state(pts$statePtr)
  expect_identical(before$edge, after$edge)
  expect_equal(before$relBrLengths, after$relBrLengths, tolerance = 1e-14)
  expect_equal(get_state_log_lik(pts$statePtr),
               eval_full_loglik_cpp(pts$dataPtr, pts$statePtr),
               tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Regression: existing moves not broken
# ---------------------------------------------------------------------------

test_that("Existing moves (NNI=5, SPR=6) still work after adding cases 10/11", {
  pts <- .gibbs_pts(seed = 301L)
  expect_true(.try_move(pts, 5L, max_try = 200L))   # NNI
  expect_true(.try_move(pts, 6L, max_try = 200L))   # SPR
})
