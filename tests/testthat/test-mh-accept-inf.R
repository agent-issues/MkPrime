# #145 F7-08: a chain started from a -Inf state gets logAlpha = +Inf for
# every proposal. Rejecting non-finite logAlpha outright trapped it there.

test_that("MH move escapes a -Inf starting state", {
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
  edgeLen <- tree$edge.length
  statePtr <- init_mcmc_state(
    tree$edge[, 1], tree$edge[, 2], edgeLen / sum(edgeLen), sum(edgeLen),
    rateLoss = 1, rateLogSd = 0.5, rateNeo = 1, p = 0.5,
    kPrime = as.integer(mkd$kObs), logLik = -Inf, logPrior = 0
  )
  statePtr <- init_mcmc_state(
    tree$edge[, 1], tree$edge[, 2], edgeLen / sum(edgeLen), sum(edgeLen),
    rateLoss = 1, rateLogSd = 0.5, rateNeo = 1, p = 0.5,
    kPrime = as.integer(mkd$kObs), logLik = -Inf,
    logPrior = eval_log_prior_cpp(dataPtr, statePtr)
  )
  expect_equal(get_state_log_lik(statePtr), -Inf)

  accepted <- do_move_cpp(dataPtr, statePtr, moveType = 0L, charIdx = 0L,
                          scaleTuning = 0.4, betaSimplexTuning = 1.0,
                          intWalkWindow = 1L, beta = 1.0)
  expect_true(accepted)
  expect_true(is.finite(get_state_log_lik(statePtr)))
  expect_equal(get_state_log_lik(statePtr),
               eval_full_loglik_cpp(dataPtr, statePtr))
})
