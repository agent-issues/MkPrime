# Under marginal_k, the committed `state->logLik` must be the marginal-over-k
# likelihood, and `state->partLogLik` must never be read as if it were.
#
# `fill_partition_cache` computes `partLogLik` with kPrime pinned to kObs, so its
# entries are FIXED-kPrime partition sums. Nothing under marginal_k refreshes
# them. Any `hasPLC` fast path that reads one therefore commits a fixed-kPrime
# total as the marginal logLik -- 5.2 nats above the marginal truth on the 8-tip
# matrix below, which is an offset on every other move's MH ratio and on the
# reported log_likelihood trace. The sign is data-dependent: inflation costs
# efficiency through spurious rejections, deflation biases the posterior. The
# assertions are two-sided for that reason.
#
# The check is the general one -- committed logLik equals a cold recompute --
# and lives in R so that it costs nothing on the hot path.

.SliceCohTree <- function() {
  Preorder(ape::read.tree(text = paste0(
    "(((t1:0.05,t2:0.07):0.04,(t3:0.06,t4:0.05):0.03):0.05,",
    "((t5:0.04,t6:0.06):0.05,(t7:0.05,t8:0.04):0.06):0.04);"
  )))
}

# Two binary characters (neomorphic -- rate_loss, and hence slice paramIdx 1,
# exists only when hasNeo) plus four multistate ones (transformational, the
# characters marginal_k integrates k' over). Every marginal_k fixture in the
# repository is trans-only, which is why paramIdx 1 has no other coverage.
.SliceCohData <- function() {
  tips <- paste0("t", 1:8)
  mat <- matrix(
    c(0, 0, 1, 1, 0, 1, 0, 1,
      0, 1, 1, 0, 1, 0, 0, 1,
      0, 1, 2, 0, 1, 2, 0, 1,
      2, 2, 0, 1, 1, 0, 2, 0,
      1, 0, 0, 2, 2, 1, 0, 1,
      0, 2, 1, 1, 0, 0, 2, 2),
    nrow = 8L, ncol = 6L, dimnames = list(tips, NULL)
  )
  MkPrimeData(MatrixToPhyDat(mat), neomorphic = 1:2)
}

.SliceCohFixture <- function(likelihoodMode) {
  tree  <- .SliceCohTree()
  mkd   <- .SliceCohData()
  model <- MkPrimeModel(kPrimePrior = "geometric",
                        likelihoodMode = likelihoodMode,
                        coding = "none", relabel = FALSE)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)

  state0 <- MkPrime:::.InitState(tree, mkd, model)
  state0$tree_length    <- sum(tree$edge.length)
  state0$rel_br_lengths <- tree$edge.length / state0$tree_length

  fixture <- list(tree = tree, mkd = mkd, model = model, state0 = state0)
  fixture$dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  fixture$statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(fixture$dataPtr, fixture$statePtr)
  allocate_cl_workspace(fixture$dataPtr, fixture$statePtr)
  fixture
}

# The cold value comes from a SECOND state built from scratch at the same
# parameters, because `eval_full_loglik_cpp` reads the charLL cache that the move
# under test has just written -- which is the thing that must not be trusted.
.SliceCohCold <- function(fixture, s) {
  state0 <- fixture$state0
  state0$tree_length    <- s$treeLength
  state0$rel_br_lengths <- s$relBrLengths
  state0$rate_loss      <- s$rateLoss
  state0$rate_log_sd    <- s$rateLogSd
  state0$p              <- s$p
  if (!is.null(state0$rate_neo))   state0$rate_neo   <- s$rateNeo
  if (!is.null(state0$beta_scale)) state0$beta_scale <- s$betaScale

  dataPtr  <- MkPrime:::.InitMcmcData(fixture$mkd, fixture$model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  get_mcmc_state(statePtr)$logLik
}

.SliceCohRun <- function(fixture, sliceParamCodes, nBatch = 20L) {
  nMoves <- length(sliceParamCodes)
  run_mcmc_batch_cpp(
    dataPtr = fixture$dataPtr, stateXPtrs = list(fixture$statePtr), betas = 1.0,
    moveTypeCodes = rep(19L, nMoves),
    transIdxCpp = which(fixture$mkd$type == "transformational") - 1L,
    sliceParamCodes = as.integer(sliceParamCodes),
    moveWeights = rep(1.0, nMoves),
    chainScaleTunings = matrix(0.5, 1L, nMoves), chainBsmpTunings = 10,
    chainIntWalkWins = 1L, moveIntParams = integer(nMoves),
    sliceWidths = matrix(1.0, 1L, nMoves), jointRhos = matrix(0.0, 1L, nMoves),
    nBatch = nBatch, startIter = 1L, warmup = nBatch, thin = nBatch,
    hasNeo = TRUE, nEdge = nrow(fixture$tree$edge)
  )
  get_mcmc_state(fixture$statePtr)
}


test_that("marginal_k slice_rate_loss commits a marginal logLik", {
  set.seed(101)
  fixture <- .SliceCohFixture("marginal_k")
  s <- .SliceCohRun(fixture, 1L)

  # A fixed-kPrime commit reads -84.1 here against a cold marginal -89.3.
  expect_equal(s$logLik, .SliceCohCold(fixture, s), tolerance = 1e-8)
  # The slice must also have moved, or the identity above is trivially true.
  expect_false(isTRUE(all.equal(s$rateLoss, fixture$state0$rate_loss)))
})


test_that("marginal_k slice ignores a partLogLik refilled by direct dispatch", {
  set.seed(202)
  fixture <- .SliceCohFixture("marginal_k")

  # `.BuildMoves` drops the k'-moves under marginal_k, but `do_move_cpp`
  # dispatches them regardless, and the k' sweep repopulates partLogLik with
  # fixed-kPrime sums. An empty cache at init is therefore not on its own enough:
  # the slice sites need their own exclusion. The sweep leaves logLik fixed-k
  # (here -83.9 against a marginal -89.2); an accepted marginal slice must heal
  # that rather than build on it.
  expect_true(do_move_cpp(fixture$dataPtr, fixture$statePtr, 25L, 0L,
                          0.5, 0.5, 2L, 1.0))
  s <- .SliceCohRun(fixture, 1L)

  expect_equal(s$logLik, .SliceCohCold(fixture, s), tolerance = 1e-8)
})


test_that("sampled_k keeps the partitioned slice fast path coherent", {
  # Control: sampled_k must be untouched by any of this. It does not pin the
  # clear's marginal_k condition -- an empty cache under sampled_k simply routes
  # the accept to the full-recompute branch, which sums to the same value.
  set.seed(303)
  fixture <- .SliceCohFixture("sampled_k")
  s <- .SliceCohRun(fixture, 1L)

  expect_equal(s$logLik, .SliceCohCold(fixture, s), tolerance = 1e-8)
})


test_that("marginal_k slice leaves the charLL cache cold when it exhausts", {
  # The shrink loop's trial evaluations refill the per-(char, k') cache. When
  # no trial is accepted the scalar reverts, so the cache no longer describes
  # the state; mh_logit_p and gibbs_p_marginal would read it without a rebuild.
  # A stored log-prior far above anything the target can reach puts every
  # trial below the slice, so the loop always exhausts.
  set.seed(404)
  fixture <- .SliceCohFixture("marginal_k")
  state0 <- fixture$state0
  state0$log_prior <- 1e6
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(fixture$dataPtr, statePtr)
  allocate_cl_workspace(fixture$dataPtr, statePtr)
  expect_true(get_marginal_cache_state(statePtr)$charLLReady)

  expect_false(do_move_cpp(fixture$dataPtr, statePtr, 19L, 0L,
                           1.0, 10, 1L, 1.0))
  expect_equal(get_mcmc_state(statePtr)$treeLength, state0$tree_length)
  expect_false(get_marginal_cache_state(statePtr)$charLLReady)
})
