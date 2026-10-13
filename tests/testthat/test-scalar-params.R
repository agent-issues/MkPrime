# Scalar parameter registry (#441 S1): kScalarParams in src/mcmc.cpp.

.ScalarFixture <- function() {
  set.seed(4411L)
  ntax <- 7L
  mat <- matrix(sample.int(3, ntax * 10, replace = TRUE) - 1L,
                nrow = ntax, ncol = 10,
                dimnames = list(paste0("t", seq_len(ntax)), NULL))
  mat[, 1:4] <- pmin(mat[, 1:4], 1L)
  mkd <- MkPrimeData(TreeTools::MatrixToPhyDat(mat), neomorphic = 1:4)
  dataPtr <- prepare_mcmc_data(
    partitions_r = mkd$partitions, kObs_r = mkd$kObs,
    charTypes_r = mkd$type, hasNeo = TRUE, nCat = 4L,
    codingStr = "variable", relabelFlag = TRUE,
    treeLengthShape = 1.5, treeLengthRate = 1.0,
    rateLossMeanlog = 0.0, rateLossSdlog = 1.0,
    rateLogSdShape = 1.0, rateLogSdRate = 1.0,
    rateNeoMeanlog = 0.0, rateNeoSdlog = 2.0,
    kprimeHyperA = 1.0, kprimeHyperB = 1.0,
    kPriorLogseries = TRUE, kprimeLogseriesC = 0.7,
    kPriorBetaGeometric = FALSE, qHeterogeneity = FALSE,
    nBetaCat = 4L, betaScaleShape = 1.0, betaScaleRate = 1.0,
    kPriorEmpiricalGeometric = FALSE, empLogBody = numeric(0)
  )
  tree <- TreeTools::Preorder(
    ape::rtree(ntax, br = function(n) runif(n, 0.05, 0.3))
  )
  list(mkd = mkd, dataPtr = dataPtr, tree = tree)
}

.ScalarState <- function(fx) {
  edgeLen <- fx$tree$edge.length
  Init <- function(logLik, logPrior) {
    init_mcmc_state(
      fx$tree$edge[, 1], fx$tree$edge[, 2], edgeLen / sum(edgeLen),
      sum(edgeLen), rateLoss = 0.8, rateLogSd = 0.5, rateNeo = 1.2,
      p = 0.5, kPrime = as.integer(fx$mkd$kObs), logLik = logLik,
      logPrior = logPrior
    )
  }
  draft <- Init(0, 0)
  Init(eval_full_loglik_cpp(fx$dataPtr, draft),
       eval_log_prior_cpp(fx$dataPtr, draft))
}

.ScalarTrace <- function(fx, moveType, charIdx, seed, n = 40L) {
  statePtr <- .ScalarState(fx)
  set.seed(seed)
  accepted <- vapply(seq_len(n), function(i) {
    do_move_cpp(fx$dataPtr, statePtr, moveType = moveType, charIdx = charIdx,
                scaleTuning = 0.6, betaSimplexTuning = 1.0,
                intWalkWindow = 1L, beta = 1.0)
  }, logical(1))
  list(accepted = accepted, state = get_mcmc_state(statePtr),
       logLik = get_state_log_lik(statePtr),
       fresh = eval_full_loglik_cpp(fx$dataPtr, statePtr))
}

test_that("R and C++ registries list the same scalars in the same order", {
  expect_identical(MkPrime:::.kScalarParams, scalar_param_names())
  expect_identical(MkPrime:::.ScalarParamIndex("rate_neo"), 3L)
  expect_error(MkPrime:::.ScalarParamIndex("rate_nao"), "not a registered")
})

test_that("generic scale_param move reproduces each dedicated scale move", {
  fx <- .ScalarFixture()
  dedicated <- c(tree_length = 0L, rate_loss = 1L, rate_log_sd = 2L,
                 rate_neo = 3L, beta_scale = 16L)
  for (param in names(dedicated)) {
    idx <- MkPrime:::.ScalarParamIndex(param)
    seed <- 100L + idx
    old <- .ScalarTrace(fx, dedicated[[param]], 0L, seed)
    new <- .ScalarTrace(fx, 40L, idx, seed)
    expect_true(any(old$accepted), label = param)
    expect_identical(new, old, label = param)
    expect_equal(new$logLik, new$fresh, tolerance = 1e-10, label = param)
  }
})

test_that("scale_param refuses an unregistered index", {
  fx <- .ScalarFixture()
  statePtr <- .ScalarState(fx)
  before <- get_mcmc_state(statePtr)
  for (idx in c(-1L, length(MkPrime:::.kScalarParams))) {
    expect_false(do_move_cpp(fx$dataPtr, statePtr, 40L, idx, 0.6, 1.0, 1L,
                             1.0))
  }
  expect_identical(get_mcmc_state(statePtr), before)
})

test_that(".ScaleParamMove builds a move dispatched by code 40", {
  move <- MkPrime:::.ScaleParamMove("rate_loss", weight = 2)
  expect_identical(MkPrime:::.MoveTypeCode(move), 40L)
  expect_identical(move$sliceParamIdx, 1L)
  expect_identical(move$name, "scale_param_rate_loss")
  expect_identical(MkPrime:::.MoveTypeCode(list(name = "rate_neo")), 3L)
})
