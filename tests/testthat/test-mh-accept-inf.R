# #145 F7-08: a chain started from a -Inf state gets logAlpha = +Inf for
# every proposal. Rejecting non-finite logAlpha outright trapped it there.

InfStartFixture <- function() {
  set.seed(145L)
  ntax <- 6L
  mat <- matrix(sample.int(2, ntax * 8, replace = TRUE) - 1L,
                nrow = ntax, ncol = 8,
                dimnames = list(paste0("t", 1:ntax), NULL))
  mkd <- MkPrimeData(TreeTools::MatrixToPhyDat(mat), neomorphic = 1:3)
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

InfState <- function(fx, kPrime, logLik) {
  edgeLen <- fx$tree$edge.length
  Init <- function(logPrior) {
    init_mcmc_state(
      fx$tree$edge[, 1], fx$tree$edge[, 2], edgeLen / sum(edgeLen),
      sum(edgeLen), rateLoss = 1, rateLogSd = 0.5, rateNeo = 1, p = 0.5,
      kPrime = as.integer(kPrime), logLik = logLik, logPrior = logPrior
    )
  }
  Init(eval_log_prior_cpp(fx$dataPtr, Init(0)))
}

test_that("MH move escapes a -Inf starting state", {
  fx <- InfStartFixture()
  statePtr <- InfState(fx, fx$mkd$kObs, -Inf)
  expect_equal(get_state_log_lik(statePtr), -Inf)

  accepted <- do_move_cpp(fx$dataPtr, statePtr, moveType = 0L, charIdx = 0L,
                          scaleTuning = 0.4, betaSimplexTuning = 1.0,
                          intWalkWindow = 1L, beta = 1.0)
  expect_true(accepted)
  expect_true(is.finite(get_state_log_lik(statePtr)))
  expect_equal(get_state_log_lik(statePtr),
               eval_full_loglik_cpp(fx$dataPtr, statePtr))
})

test_that("block Gibbs branch sweep escapes a -Inf starting state", {
  fx <- InfStartFixture()
  statePtr <- InfState(fx, fx$mkd$kObs, -Inf)

  set.seed(1L)
  do_move_cpp(fx$dataPtr, statePtr, moveType = 15L, charIdx = 0L,
              scaleTuning = 0.4, betaSimplexTuning = 1.0,
              intWalkWindow = 1L, beta = 1.0)
  expect_true(is.finite(get_state_log_lik(statePtr)))
  expect_equal(get_state_log_lik(statePtr),
               eval_full_loglik_cpp(fx$dataPtr, statePtr))
})

# #339 RT7-06: at beta = 0, beta * (newLL - curLL) is 0 * Inf = NaN when the
# current logLik is -Inf, so the +Inf escape above never fired.
test_that("MH move escapes a -Inf starting state at beta = 0", {
  fx <- InfStartFixture()
  statePtr <- InfState(fx, fx$mkd$kObs, -Inf)

  set.seed(339L)
  accepted <- vapply(seq_len(5L), function(i)
    do_move_cpp(fx$dataPtr, statePtr, moveType = 0L, charIdx = 0L,
                scaleTuning = 0.4, betaSimplexTuning = 1.0,
                intWalkWindow = 1L, beta = 0), logical(1))
  expect_true(any(accepted))
  expect_true(is.finite(get_state_log_lik(statePtr)))
  expect_equal(get_state_log_lik(statePtr),
               eval_full_loglik_cpp(fx$dataPtr, statePtr))
})

# #485: these moves weight candidates by exp(beta * logLik) before any MH
# step, so at beta = 0 the -Inf current state gave 0 * -Inf = NaN weights.
EscapeAtBetaZero <- function(moveType, tries = 10L) {
  fx <- InfStartFixture()
  statePtr <- InfState(fx, fx$mkd$kObs, -Inf)
  for (i in seq_len(tries)) {
    do_move_cpp(fx$dataPtr, statePtr, moveType = moveType, charIdx = 0L,
                scaleTuning = 1.0, betaSimplexTuning = 1.0,
                intWalkWindow = 1L, beta = 0)
    if (is.finite(get_state_log_lik(statePtr))) break
  }
  list(fx = fx, statePtr = statePtr)
}

test_that("subtree swaps and slice escape a -Inf state at beta = 0", {
  set.seed(485L)
  for (moveType in c(11L, 14L, 19L)) {
    run <- EscapeAtBetaZero(moveType)
    expect_true(is.finite(get_state_log_lik(run$statePtr)),
                info = paste("moveType", moveType))
    expect_equal(get_state_log_lik(run$statePtr),
                 eval_full_loglik_cpp(run$fx$dataPtr, run$statePtr),
                 info = paste("moveType", moveType))
  }
})

test_that("subtree swaps never enter a -Inf state at beta = 0", {
  fx <- InfStartFixture()
  logLik <- eval_full_loglik_cpp(fx$dataPtr, InfState(fx, fx$mkd$kObs, 0))
  statePtr <- InfState(fx, fx$mkd$kObs, logLik)

  set.seed(4850L)
  for (i in seq_len(40)) {
    do_move_cpp(fx$dataPtr, statePtr, moveType = c(11L, 14L)[i %% 2 + 1],
                charIdx = 0L, scaleTuning = 1.0, betaSimplexTuning = 1.0,
                intWalkWindow = 1L, beta = 0)
    expect_true(is.finite(get_state_log_lik(statePtr)))
  }
  expect_equal(get_state_log_lik(statePtr),
               eval_full_loglik_cpp(fx$dataPtr, statePtr))
})
