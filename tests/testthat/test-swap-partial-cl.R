# M-111: Validate partial CL evaluation for Gibbs subtree swap
#
# For multiple tree sizes and model types, compare evaluate_swap_candidate()
# (partial CL, O(depth×C)) against compute_full_loglik_at() (full pruning,
# O(N×C)) for every valid swap partner of every node.

# Shared fixture
.swap_fixture <- function(seed, nTip, nChar, maxState,
                          model = MkPrimeModel()) {
  set.seed(seed)
  tree <- Preorder(rtree(nTip, rooted = FALSE))
  mat  <- matrix(sample(0:maxState, nTip * nChar, replace = TRUE),
                 nrow = nTip, ncol = nChar,
                 dimnames = list(tree$tip.label, NULL))
  pd   <- suppressWarnings(MatrixToPhyDat(mat))
  mkd  <- suppressWarnings(MkPrimeData(pd))
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  state0   <- MkPrime:::.InitState(tree, mkd, model)
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  list(tree = tree, dataPtr = dataPtr, statePtr = statePtr)
}

# Helper: validate all swap candidates for every node
.validate_all_swaps <- function(fix, tol = 1e-8) {
  maxDiff <- 0
  nCand   <- 0L
  for (nodeA in fix$tree$edge[, 2]) {
    res <- validate_swap_partial_cl(fix$dataPtr, fix$statePtr, nodeA)
    if (nrow(res) == 0) next
    nCand   <- nCand + nrow(res)
    maxDiff <- max(maxDiff, max(abs(res$ll_partial - res$ll_full)))
  }
  expect_lt(maxDiff, tol,
            label = sprintf("maxDiff (n=%d candidates)", nCand))
  invisible(nCand)
}

test_that("partial CL matches full eval: 8-tip transformational", {
  fix <- .swap_fixture(seed = 3847, nTip = 8, nChar = 10, maxState = 3)
  .validate_all_swaps(fix)
})

test_that("partial CL matches full eval: 8-tip binary (neomorphic)", {
  fix <- .swap_fixture(seed = 3847, nTip = 8, nChar = 10, maxState = 1)
  .validate_all_swaps(fix)
})

test_that("partial CL matches full eval: 20-tip transformational", {
  fix <- .swap_fixture(seed = 5603, nTip = 20, nChar = 15, maxState = 4)
  .validate_all_swaps(fix)
})

test_that("partial CL matches full eval: 20-tip with ACRV", {
  fix <- .swap_fixture(seed = 5603, nTip = 20, nChar = 15, maxState = 4,
                       model = MkPrimeModel(nCat = 4L))
  .validate_all_swaps(fix)
})

test_that("partial CL matches full eval: 54-tip tree", {
  fix <- .swap_fixture(seed = 1842, nTip = 54, nChar = 20, maxState = 5)
  .validate_all_swaps(fix)
})
