# Fixtures for test-extensions.R and test-extension-likelihood.R (#441).

ExtFixture <- function() {
  set.seed(4412L)
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

ExtState <- function(fx, extParams) {
  edgeLen <- fx$tree$edge.length
  Init <- function(logLik, logPrior) {
    init_mcmc_state(
      fx$tree$edge[, 1], fx$tree$edge[, 2], edgeLen / sum(edgeLen),
      sum(edgeLen), rateLoss = 0.8, rateLogSd = 0.5, rateNeo = 1.2,
      p = 0.5, kPrime = as.integer(fx$mkd$kObs), logLik = logLik,
      logPrior = logPrior, extParams = extParams
    )
  }
  draft <- Init(0, 0)
  Init(eval_full_loglik_cpp(fx$dataPtr, draft),
       eval_log_prior_cpp(fx$dataPtr, draft))
}

ExtRunData <- function() {
  set.seed(1)
  tree <- ape::rtree(8)
  tree$tip.label <- paste0("t", 1:8)
  mat <- matrix(sample(0:1, 8 * 20, TRUE), 8, 20,
                dimnames = list(tree$tip.label, NULL))
  list(tree = tree, pd = TreeTools::MatrixToPhyDat(mat))
}
