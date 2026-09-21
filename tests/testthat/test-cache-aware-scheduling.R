# M-159: Cache-aware move scheduling tests

test_that("cacheBonus increases partial-CL move selection when cache is valid", {
  set.seed(4917)
  tree <- ape::rtree(10, br = runif)
  tree$edge.length <- abs(tree$edge.length) + 0.01
  tree <- ape::unroot(tree)
  nTip <- ape::Ntip(tree)

  nChar <- 4
  kStates <- 2L
  tipStates <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                      nrow = nTip, ncol = nChar,
                      dimnames = list(tree$tip.label, NULL))
  pd <- MatrixToPhyDat(tipStates)
  mkd <- MkPrimeData(pd, knownStates = rep(kStates, nChar))

  model <- MkPrimeModel(coding = "variable", nCat = 1L)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  tree0 <- Preorder(tree)

  state0 <- MkPrime:::.InitState(tree0, mkd, model)
  state0$log_prior <- MkPrime:::LogPrior(state0, model, mkd)

  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)

  # Run two short batches: one with cacheBonus=1, one with cacheBonus=10
  run_batch <- function(bonus) {
    statePtr <- MkPrime:::.InitMcmcChain(state0)
    fill_partition_cache(dataPtr, statePtr)
    allocate_cl_workspace(dataPtr, statePtr)

    # Simple move set: NNI (5) and tree_length scale (0)
    moveTypes <- c(5L, 0L)
    moveWeights <- c(0.5, 0.5)
    nMoves <- length(moveTypes)

    scaleTunings <- matrix(0.5, nrow = 1, ncol = nMoves)
    sliceWidths  <- matrix(1.0, nrow = 1, ncol = nMoves)
    jointRhos    <- matrix(0.0, nrow = 1, ncol = nMoves)

    run_mcmc_batch_cpp(
      dataPtr, list(statePtr), 1.0,
      moveTypes, integer(0), integer(0), moveWeights,
      scaleTunings, 10.0, 1L, integer(nMoves),
      sliceWidths, jointRhos,
      nBatch = 2000L, startIter = 0L, warmup = 2000L, thin = 1L,
      hasNeo = FALSE, nEdge = nrow(tree0$edge),
      cacheBonus = bonus
    )
  }

  res_base  <- run_batch(1.0)
  res_boost <- run_batch(10.0)

  # NNI is move index 1 (0-indexed in C++, 1-indexed in R matrix)
  nni_base  <- res_base$propose_counts[1, 1]
  nni_boost <- res_boost$propose_counts[1, 1]

  # With cacheBonus=10, NNI should be proposed substantially more often
  # because once the cache is populated by an NNI, subsequent selections
  # will use boosted weights favoring NNI.
  expect_gt(nni_boost, nni_base)

  # Cache hits should be positive when bonus is active
  expect_gt(res_boost$cache_hits, 0)
  # No cache hits when bonus is 1.0 (haveCacheBoost is false)
  expect_equal(res_base$cache_hits, 0L)
})

test_that("cacheBonus = 1 reproduces default behaviour", {
  set.seed(6213)
  tree <- ape::rtree(8, br = runif)
  tree$edge.length <- abs(tree$edge.length) + 0.01
  tree <- ape::unroot(tree)
  nTip <- ape::Ntip(tree)

  nChar <- 3
  kStates <- 2L
  tipStates <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                      nrow = nTip, ncol = nChar,
                      dimnames = list(tree$tip.label, NULL))
  pd <- MatrixToPhyDat(tipStates)
  mkd <- MkPrimeData(pd, knownStates = rep(kStates, nChar))

  model <- MkPrimeModel(coding = "variable", nCat = 1L)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  tree0 <- Preorder(tree)

  state0 <- MkPrime:::.InitState(tree0, mkd, model)
  state0$log_prior <- MkPrime:::LogPrior(state0, model, mkd)

  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)

  moveTypes <- c(5L, 0L)
  moveWeights <- c(0.5, 0.5)
  nMoves <- length(moveTypes)

  scaleTunings <- matrix(0.5, nrow = 1, ncol = nMoves)
  sliceWidths  <- matrix(1.0, nrow = 1, ncol = nMoves)
  jointRhos    <- matrix(0.0, nrow = 1, ncol = nMoves)

  res <- run_mcmc_batch_cpp(
    dataPtr, list(statePtr), 1.0,
    moveTypes, integer(0), integer(0), moveWeights,
    scaleTunings, 10.0, 1L, integer(nMoves),
    sliceWidths, jointRhos,
    nBatch = 500L, startIter = 0L, warmup = 500L, thin = 1L,
    hasNeo = FALSE, nEdge = nrow(tree0$edge),
    cacheBonus = 1.0
  )

  # With cacheBonus=1, cache_hits should be 0 (boost disabled)
  expect_equal(res$cache_hits, 0L)
  # Total proposals should equal nBatch * nChains
  total_propose <- sum(res$propose_counts)
  expect_equal(total_propose, 500L)
})

test_that("MkPrimeMCMC stores cacheBonus", {
  mcmc <- suppressWarnings(MkPrimeMCMC(nIter = 100, cacheBonus = 3))
  expect_equal(mcmc$cacheBonus, 3)

  mcmc_default <- suppressWarnings(MkPrimeMCMC(nIter = 100))
  expect_equal(mcmc_default$cacheBonus, 5)
})

test_that("MkPrimeMCMC rejects invalid cacheBonus", {
  expect_error(suppressWarnings(MkPrimeMCMC(nIter = 100, cacheBonus = 0.5)),
               "cacheBonus")
  expect_error(suppressWarnings(MkPrimeMCMC(nIter = 100, cacheBonus = -1)),
               "cacheBonus")
})
