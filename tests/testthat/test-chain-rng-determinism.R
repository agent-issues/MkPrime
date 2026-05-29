# T-017-IIa: ChainRng determinism and coverage gate tests
#
# Test 1: Same-seed reproducibility
#   Two calls to run_mcmc_batch_cpp with identical R seeds must produce
#   identical scalar trajectories.
#
# Test 2: Coverage gate — exactly one R-RNG draw consumed
#   run_mcmc_batch_cpp draws exactly ONE R::unif_rand() value (the base seed
#   for the per-chain Weyl-mix stream).  All other random draws go through
#   ChainRng.  Verified by checking that the RNG position advances by
#   exactly 1 during a batch call.

library("ape")
library("TreeTools")

# ---------------------------------------------------------------------------
# Shared fixture: minimal MkPrime MCMC setup (mirrors test-cache-aware-scheduling.R)
# ---------------------------------------------------------------------------

.crng_setup <- function(seed = 42L, nTip = 8L, nChar = 4L) {
  set.seed(seed)
  tree <- ape::rtree(nTip, br = runif)
  tree$edge.length <- abs(tree$edge.length) + 0.01
  tree <- ape::unroot(tree)
  tree <- TreeTools::Preorder(tree)

  tipStates <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                      nrow = nTip, ncol = nChar,
                      dimnames = list(tree$tip.label, NULL))
  pd  <- TreeTools::MatrixToPhyDat(tipStates)
  mkd <- suppressWarnings(MkPrimeData(pd, knownStates = rep(2L, nChar)))

  model  <- MkPrimeModel(coding = "variable", nCat = 1L)
  model  <- MkPrime:::.FinalizeModel(model, tree, mkd)
  state0 <- MkPrime:::.InitState(tree, mkd, model)
  state0$log_prior <- MkPrime:::LogPrior(state0, model, mkd)

  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)

  nEdge    <- nrow(tree$edge)
  # Simple move set: scale tree_length (0), beta_simplex (4), NNI (5)
  moveTypes   <- c(0L, 4L, 5L)
  moveWeights <- c(1.0, 2.0, 2.0)
  nMoves      <- length(moveTypes)
  scaleTunings <- matrix(0.5, nrow = 1L, ncol = nMoves)
  sliceWidths  <- matrix(1.0, nrow = 1L, ncol = nMoves)
  jointRhos    <- matrix(0.0, nrow = 1L, ncol = nMoves)

  list(
    dataPtr   = dataPtr,
    statePtr  = statePtr,
    nEdge     = nEdge,
    moveTypes = moveTypes, moveWeights = moveWeights, nMoves = nMoves,
    scaleTunings = scaleTunings, sliceWidths = sliceWidths,
    jointRhos = jointRhos
  )
}

.run_batch <- function(s, statePtr, nBatch = 30L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  run_mcmc_batch_cpp(
    s$dataPtr, list(statePtr), 1.0,
    s$moveTypes, integer(0), integer(0), s$moveWeights,
    s$scaleTunings, 0.5, 1L, integer(s$nMoves),
    s$sliceWidths, s$jointRhos,
    nBatch = nBatch, startIter = 0L, warmup = nBatch, thin = 1L,
    hasNeo = FALSE, nEdge = s$nEdge
  )
}

# ---------------------------------------------------------------------------
# Test 1: Same-seed reproducibility
# ---------------------------------------------------------------------------

test_that("T-017-IIa: same seed produces identical scalar trajectories", {
  # Two independent state objects built from the same initial conditions
  s1 <- .crng_setup(seed = 1701L)
  s2 <- .crng_setup(seed = 1701L)

  # Reset states (already fresh from .crng_setup)
  r1 <- .run_batch(s1, s1$statePtr, nBatch = 40L, seed = 999L)
  r2 <- .run_batch(s2, s2$statePtr, nBatch = 40L, seed = 999L)

  # Scalar accept/propose counts must be identical
  expect_equal(r1$accept_counts,  r2$accept_counts,
               info = "Same seed must yield identical accept counts")
  expect_equal(r1$propose_counts, r2$propose_counts,
               info = "Same seed must yield identical propose counts")
})


# ---------------------------------------------------------------------------
# Test 2: Coverage gate — exactly one R-RNG draw consumed per batch call
# ---------------------------------------------------------------------------

test_that("T-017-IIa: run_mcmc_batch_cpp consumes exactly one R-RNG draw", {
  s <- .crng_setup(seed = 1702L)

  # Establish the expected value of the SECOND draw from R's RNG with this seed.
  set.seed(12345L)
  runif(1)           # draw 1 — consumed by run_mcmc_batch_cpp as base_seed
  ref2 <- runif(1)   # draw 2 — must be the next available draw after the batch

  # Reset seed to same state and run the batch call
  set.seed(12345L)
  .run_batch(s, s$statePtr, nBatch = 50L)  # seed arg already set above

  # The next R-RNG draw must equal ref2 (draw 2 from the same seed sequence)
  next_draw <- runif(1)
  expect_equal(next_draw, ref2,
               info = paste0(
                 "run_mcmc_batch_cpp must consume exactly one R-RNG draw. ",
                 "Expected next draw = ", ref2,
                 ", got ", next_draw
               ))
})
