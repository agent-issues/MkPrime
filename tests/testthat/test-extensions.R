# Engine extensions, parameters and priors only (#441 S2a):
# inst/include/MkPrime/extension.h, src/extensions.cpp, R/extensions.R.
# toy_extension_cpp(): toy_sigma ~ LogNormal(0, 0.5), toy_rate ~ Gamma(2, 4).

ToyPrior <- function(x) {
  stats::dlnorm(x[[1]], 0, 0.5, log = TRUE) +
    stats::dgamma(x[[2]], shape = 2, rate = 4, log = TRUE)
}

Toy <- function() MkPrimeExtension(function() toy_extension_cpp())

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

test_that("an extension pointer describes itself", {
  expect_identical(extension_abi_version(), 1L)
  info <- extension_info_cpp(toy_extension_cpp())
  expect_identical(info$name, "toy")
  expect_identical(info$abi_version, 1L)
  expect_identical(info$params, c("toy_sigma", "toy_rate"))
  expect_identical(info$init, c(1, 0.5))

  ext <- Toy()
  expect_s3_class(ext, "MkPrimeExtension")
  expect_identical(ext$init, c(toy_sigma = 1, toy_rate = 0.5))
})

test_that("malformed extensions are refused", {
  expect_error(MkPrimeExtension(toy_extension_cpp()), "must be a function")
  expect_error(MkPrimeExtension(function() 1), "Not an MkPrimeExtension")
  expect_error(extension_info_cpp(NULL), "Not an MkPrimeExtension")
  expect_error(MkPrimeModel(extensions = list(Toy(), Toy())), "defined twice")
  expect_error(MkPrimeModel(extensions = list("toy")), "must be a list")
  fake <- structure(list(params = c("tree_length", "br_3")),
                    class = "MkPrimeExtension")
  expect_error(MkPrimeModel(extensions = list(fake)), "clash")
})

test_that("models without extensions are unchanged", {
  expect_identical(MkPrimeModel(extensions = NULL), MkPrimeModel())
  expect_identical(MkPrimeModel(extensions = list()), MkPrimeModel())
  expect_null(MkPrimeModel()$extensions)
  model <- MkPrimeModel(extensions = Toy())
  expect_length(model$extensions, 1L)
  expect_identical(MkPrime:::.ExtensionParamNames(model),
                   c("toy_sigma", "toy_rate"))
  expect_identical(MkPrime:::.WithExtensionMoves(list(1), MkPrimeModel()),
                   list(1))
  moves <- MkPrime:::.WithExtensionMoves(list(), model)
  expect_identical(vapply(moves, `[[`, "", "name"),
                   c("scale_param_toy_sigma", "scale_param_toy_rate"))
  expect_identical(vapply(moves, `[[`, 0L, "sliceParamIdx"),
                   length(MkPrime:::.kScalarParams) + 0:1)
})

test_that("the extension prior is the extension's own", {
  ptr <- toy_extension_cpp()
  for (x in list(c(1, 0.5), c(0.3, 1.7), c(2.4, 0.05))) {
    expect_equal(extension_log_prior_cpp(ptr, x), ToyPrior(x),
                 tolerance = 1e-12)
  }
  expect_identical(extension_log_prior_cpp(ptr, c(-1, 0.5)), -Inf)
  expect_error(extension_log_prior_cpp(ptr, 1), "has 2 parameters")
  model <- MkPrimeModel(extensions = Toy())
  expect_equal(MkPrime:::.ExtensionLogPrior(c(0.3, 1.7), model),
               ToyPrior(c(0.3, 1.7)), tolerance = 1e-12)
})

test_that("scale moves on extension parameters touch only the prior", {
  fx <- ExtFixture()
  base <- ExtState(fx, numeric(0))
  baseLP <- eval_log_prior_cpp(fx$dataPtr, base)
  attach_extensions_cpp(fx$dataPtr, list(toy_extension_cpp()))
  # A state whose extension block does not fit the attached extensions.
  expect_identical(eval_log_prior_cpp(fx$dataPtr, base), -Inf)

  # .InitMcmcChain() (stepping stone) carries the extension block.
  edgeLen <- fx$tree$edge.length
  rState <- list(tree = fx$tree, rel_br_lengths = edgeLen / sum(edgeLen),
                 tree_length = sum(edgeLen), rate_loss = 0.8,
                 rate_log_sd = 0.5, rate_neo = 1.2, p = 0.5,
                 kPrime = as.integer(fx$mkd$kObs), log_lik = 0, log_prior = 0,
                 ext_params = c(1, 0.5))
  expect_identical(get_mcmc_state(MkPrime:::.InitMcmcChain(rState))$extParams,
                   c(1, 0.5))

  statePtr <- ExtState(fx, c(1, 0.5))
  start <- get_mcmc_state(statePtr)
  expect_equal(start$logPrior - baseLP, ToyPrior(c(1, 0.5)),
               tolerance = 1e-10)
  nScalar <- length(MkPrime:::.kScalarParams)
  expect_false(do_move_cpp(fx$dataPtr, statePtr, 40L, nScalar + 2L, 0.6,
                           1.0, 1L, 1.0))

  set.seed(4413L)
  prev <- start
  nAccept <- 0L
  for (i in seq_len(60L)) {
    j <- (i %% 2L) + 1L
    accepted <- do_move_cpp(fx$dataPtr, statePtr, 40L, nScalar + j - 1L,
                            0.6, 1.0, 1L, 1.0)
    now <- get_mcmc_state(statePtr)
    nAccept <- nAccept + accepted
    expect_identical(now$extParams[-j], prev$extParams[-j])
    if (!accepted) expect_identical(now$extParams, prev$extParams)
    expect_identical(now$logLik, start$logLik)
    expect_identical(now[setdiff(names(now), c("extParams", "logPrior",
                                               "logPost"))],
                     start[setdiff(names(start), c("extParams", "logPrior",
                                                   "logPost"))])
    expect_equal(now$logPrior, eval_log_prior_cpp(fx$dataPtr, statePtr),
                 tolerance = 1e-10)
    expect_equal(now$logPrior - baseLP, ToyPrior(now$extParams),
                 tolerance = 1e-10)
    prev <- now
  }
  expect_gt(nAccept, 5L)
  expect_lt(nAccept, 60L)
})

ExtRunData <- function() {
  set.seed(1)
  tree <- ape::rtree(8)
  tree$tip.label <- paste0("t", 1:8)
  mat <- matrix(sample(0:1, 8 * 20, TRUE), 8, 20,
                dimnames = list(tree$tip.label, NULL))
  list(tree = tree, pd = TreeTools::MatrixToPhyDat(mat))
}

test_that("RunMkPrime samples extension parameters from their prior", {
  skip_on_cran()
  d <- ExtRunData()
  set.seed(4414L)
  res <- allow_warning(RunMkPrime(
    d$pd, d$tree,
    model = MkPrimeModel(extensions = list(Toy())),
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 20000L, thin = 5L,
                       minWarmup = 1000L, maxWarmup = 1000L,
                       autoTune = FALSE, maxTime = 120)
  ), "maxWarmup|stabilis")
  s <- res$samples
  expect_true(all(c("toy_sigma", "toy_rate") %in% colnames(s)))
  expect_true(all(s[, "toy_sigma"] > 0 & s[, "toy_rate"] > 0))
  expect_gt(length(unique(s[, "toy_rate"])), 20L)
  # Prior: E[log toy_sigma] = 0 (sd 0.5); E[toy_rate] = 0.5 (sd 0.35).
  expect_lt(abs(mean(log(s[, "toy_sigma"]))), 0.35)
  expect_lt(abs(mean(s[, "toy_rate"]) - 0.5), 0.25)
})

test_that("extension parameters survive tempering and checkpoint resume", {
  skip_on_cran()
  d <- ExtRunData()
  dir <- withr::local_tempdir()
  ckp <- file.path(dir, "run.ckp")
  model <- MkPrimeModel(extensions = list(Toy()))
  set.seed(4415L)
  res <- allow_warning(RunMkPrime(
    d$pd, d$tree, model = model,
    mcmc = MkPrimeMCMC(nRuns = 1L, nChains = 2L, nIter = 2000L, thin = 5L,
                       minWarmup = 500L, maxWarmup = 500L, autoTune = FALSE,
                       checkpointFile = ckp, maxTime = 60)
  ), "maxWarmup|stabilis")
  expect_true(all(c("toy_sigma", "toy_rate") %in% colnames(res$samples)))
  ckpt <- readRDS(ckp)
  expect_true(all(c("toy_sigma", "toy_rate") %in% ckpt$paramNames))
  for (ch in ckpt$runs[[1]]$chains) {
    expect_length(ch$ext_params, 2L)
    expect_true(all(ch$ext_params > 0))
    expect_false(identical(ch$ext_params, c(1, 0.5)))
  }

  resumed <- allow_warning(ResumeMkPrime(
    ckp, d$pd, d$tree,
    mcmc = list(progressFn = function(info) NULL, nIter = 3000L)
  ), "maxWarmup|stabilis")
  s <- resumed$samples
  expect_true(all(c("toy_sigma", "toy_rate") %in% colnames(s)))
  expect_gt(nrow(s), nrow(res$samples))
  expect_true(all(s[, "toy_sigma"] > 0 & s[, "toy_rate"] > 0))
})
