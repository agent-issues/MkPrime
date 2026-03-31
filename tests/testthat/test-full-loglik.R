# Tests for compute_full_loglik / compute_full_loglik_at (M-083)
#
# eval_full_loglik_cpp(dataPtr, statePtr)        — current-topology path
# eval_full_loglik_at_cpp(dataPtr, statePtr,     — arbitrary-topology path
#                         parent, child, edgeLen)

library("ape")
library("TreeTools")

# ---------------------------------------------------------------------------
# Shared fixture helpers
# ---------------------------------------------------------------------------

.make_trans_tree <- function() {
  read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
}

.make_trans_pd <- function() {
  mat <- matrix(
    c(0, 1, 0, 1,
      1, 1, 0, 0,
      0, 1, 2, 0),
    nrow = 4, ncol = 3,
    dimnames = list(paste0("t", 1:4), NULL)
  )
  MatrixToPhyDat(mat)
}

.make_ptrs <- function(with_workspace = TRUE) {
  tree  <- .make_trans_tree()
  pd    <- .make_trans_pd()
  mkd   <- MkPrimeData(pd)
  model <- MkPrimeModel()
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  # RunMkPrime reorders to canonical preorder before .InitState; replicate that here.
  tree <- TreeTools::Preorder(tree)

  state0   <- MkPrime:::.InitState(tree, mkd, model)
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  if (with_workspace) allocate_cl_workspace(dataPtr, statePtr)

  s <- get_mcmc_state(statePtr)
  list(dataPtr  = dataPtr,
       statePtr = statePtr,
       parent   = s$edge[, 1L],
       child    = s$edge[, 2L],
       relBr    = s$relBrLengths,
       treeLen  = s$treeLength)
}


# ---------------------------------------------------------------------------
# Test: current-topology path matches cached logLik
# ---------------------------------------------------------------------------

test_that("eval_full_loglik_cpp matches cached logLik (no workspace)", {
  pts    <- .make_ptrs(with_workspace = FALSE)
  cached <- get_state_log_lik(pts$statePtr)
  evaled <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)

  expect_equal(evaled, cached, tolerance = 1e-10)
  expect_true(is.finite(evaled))
})

test_that("eval_full_loglik_cpp matches cached logLik (with workspace)", {
  pts    <- .make_ptrs(with_workspace = TRUE)
  cached <- get_state_log_lik(pts$statePtr)
  evaled <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)

  expect_equal(evaled, cached, tolerance = 1e-10)
  expect_true(is.finite(evaled))
})


# ---------------------------------------------------------------------------
# Test: eval_full_loglik_at_cpp with same topology == eval_full_loglik_cpp
# ---------------------------------------------------------------------------

test_that("eval_full_loglik_at_cpp with current topology == eval_full_loglik_cpp", {
  pts     <- .make_ptrs(with_workspace = TRUE)
  edgeLen <- pts$treeLen * pts$relBr

  at_val   <- eval_full_loglik_at_cpp(pts$dataPtr, pts$statePtr,
                                      pts$parent, pts$child, edgeLen)
  full_val <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)

  expect_equal(at_val, full_val, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# Test: scaled edge lengths produce a different likelihood
# ---------------------------------------------------------------------------

test_that("eval_full_loglik_at_cpp detects different edge lengths", {
  pts        <- .make_ptrs(with_workspace = TRUE)
  edgeLen    <- pts$treeLen * pts$relBr
  edgeLenMod <- edgeLen * 2.0

  orig_val <- eval_full_loglik_at_cpp(pts$dataPtr, pts$statePtr,
                                      pts$parent, pts$child, edgeLen)
  mod_val  <- eval_full_loglik_at_cpp(pts$dataPtr, pts$statePtr,
                                      pts$parent, pts$child, edgeLenMod)

  expect_true(is.finite(orig_val))
  expect_true(is.finite(mod_val))
  expect_false(isTRUE(all.equal(orig_val, mod_val, tolerance = 1e-6)))
})


# ---------------------------------------------------------------------------
# Test: NNI-proposed topology evaluates without error
# ---------------------------------------------------------------------------

test_that("eval_full_loglik_at_cpp with NNI topology gives finite result", {
  set.seed(4231)
  pts  <- .make_ptrs(with_workspace = TRUE)
  edge <- cbind(pts$parent, pts$child)

  nni_res <- nni_proposal(edge, nTip = 4L,
                          treeLength   = pts$treeLen,
                          relBrLengths = pts$relBr)

  skip_if(is.infinite(nni_res$logHastings) && nni_res$logHastings < 0,
          "No valid NNI on this tree")

  nni_par <- nni_res$edge[, 1L]
  nni_ch  <- nni_res$edge[, 2L]
  nni_el  <- pts$treeLen * nni_res$rel_br_lengths

  nni_ll <- eval_full_loglik_at_cpp(pts$dataPtr, pts$statePtr,
                                    nni_par, nni_ch, nni_el)
  expect_true(is.finite(nni_ll))
})


# ---------------------------------------------------------------------------
# Test: calling _at_ does NOT mutate state's cached logLik
# ---------------------------------------------------------------------------

test_that("eval_full_loglik_at_cpp does not mutate state logLik", {
  pts     <- .make_ptrs(with_workspace = TRUE)
  before  <- get_state_log_lik(pts$statePtr)
  edgeLen <- pts$treeLen * pts$relBr * 3.0

  eval_full_loglik_at_cpp(pts$dataPtr, pts$statePtr,
                          pts$parent, pts$child, edgeLen)

  after <- get_state_log_lik(pts$statePtr)
  expect_equal(before, after, tolerance = 1e-14)
})


# ---------------------------------------------------------------------------
# Test: workspace vs. no-workspace give identical results
# ---------------------------------------------------------------------------

test_that("workspace and no-workspace paths give identical results", {
  pts_ws  <- .make_ptrs(with_workspace = TRUE)
  pts_nws <- .make_ptrs(with_workspace = FALSE)
  edgeLen <- pts_ws$treeLen * pts_ws$relBr * 1.5

  ws_val  <- eval_full_loglik_at_cpp(pts_ws$dataPtr,  pts_ws$statePtr,
                                     pts_ws$parent,  pts_ws$child, edgeLen)
  nws_val <- eval_full_loglik_at_cpp(pts_nws$dataPtr, pts_nws$statePtr,
                                     pts_nws$parent, pts_nws$child, edgeLen)

  expect_equal(ws_val, nws_val, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# Test: repeated calls are deterministic (no side effects)
# ---------------------------------------------------------------------------

test_that("eval_full_loglik_at_cpp is deterministic across calls", {
  pts     <- .make_ptrs(with_workspace = TRUE)
  edgeLen <- pts$treeLen * pts$relBr

  val1 <- eval_full_loglik_at_cpp(pts$dataPtr, pts$statePtr,
                                   pts$parent, pts$child, edgeLen)
  val2 <- eval_full_loglik_at_cpp(pts$dataPtr, pts$statePtr,
                                   pts$parent, pts$child, edgeLen)

  expect_equal(val1, val2, tolerance = 1e-14)
})
