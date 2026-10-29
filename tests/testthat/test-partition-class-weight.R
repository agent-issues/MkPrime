# The Dirichlet(classRateConcentration) prior on class_w (#359) and the
# class-rate normaliser shared by every partition of one evaluation (#362).

test_that("classRateConcentration must be a positive finite scalar (#359)", {
  bad <- list(0, -1, NA_real_, Inf, c(1, 2), "1")
  for (alpha in bad) {
    expect_error(MkPrimeModel(classRateConcentration = alpha),
                 "classRateConcentration.*single positive finite number",
                 label = deparse(alpha))
  }

  set.seed(3)
  mat <- matrix(sample(0:2, 48, TRUE), 6,
                dimnames = list(paste0("t", 1:6), NULL))
  mkd <- MkPrimeData(MatrixToPhyDat(mat))
  for (alpha in bad) {
    model <- MkPrimeModel()
    model$classRateConcentration <- alpha
    expect_error(.FinalizeModel(model, NULL, mkd),
                 "classRateConcentration.*single positive finite number",
                 label = deparse(alpha))
  }

  # A model list saved before the field existed takes the default.
  model <- MkPrimeModel()
  model$classRateConcentration <- NULL
  expect_identical(.FinalizeModel(model, NULL, mkd)$classRateConcentration, 1)
})


.ClassWeightPriors <- function(alpha, unlink) {
  set.seed(3)
  nTip <- 8L
  nChar <- 12L
  mat <- matrix(sample(0:2, nTip * nChar, TRUE), nTip,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  tree <- Preorder(NJTree(pd, edgeLengths = TRUE))
  tree$edge.length <- pmax(tree$edge.length, 0.05)
  model <- .FinalizeModel(MkPrimeModel(classRateConcentration = alpha),
                          NULL, mkd)
  spec <- .ValidatePartitionArgs(rep(1:3, length.out = mkd$nChar), unlink,
                                 mkd)
  mkdK <- mkd
  mkdK$partitions <- .BuildPartitions(mkd, spec$partition)
  set.seed(1)
  sP <- .InitStatePartitioned(tree, mkdK, model, spec)
  set.seed(1)
  sL <- .InitState(tree, mkd, model)
  run <- .InitRun(tree, mkdK, model, MkPrimeMCMC(nChains = 1L),
                  list(list(name = "x", weight = 1)), partitionSpec = spec)
  dp <- .InitMcmcData(mkdK, model)
  sp <- .DeserialiseChain(run$chains[[1]])
  fill_partition_cache(dp, sp)
  list(legacy = sL$log_prior, r = sP$log_prior,
       cpp = get_mcmc_state(sp)$logPrior,
       wrapper = eval_log_prior_partitioned_cpp(
         dp, sp, sP$class_rate_log_sd, sP$class_w, sP$class_rate, 1),
       classW = sP$class_w)
}

test_that("linked rate multiplier adds no class_w Dirichlet term (#359)", {
  for (alpha in c(0.4, 1, 5)) {
    lp <- .ClassWeightPriors(alpha, character(0))
    expect_equal(lp$r, lp$legacy, tolerance = 1e-10)
    expect_equal(lp$cpp, lp$legacy, tolerance = 1e-10)
    expect_equal(lp$wrapper, lp$legacy, tolerance = 1e-10)
  }
})

test_that("unlinked rate multiplier keeps the class_w Dirichlet term", {
  for (alpha in c(0.4, 1, 5)) {
    lp <- .ClassWeightPriors(alpha, "ratemultiplier")
    w <- lp$classW
    K <- length(w)
    logDir <- lgamma(K * alpha) - K * lgamma(alpha) +
      (alpha - 1) * sum(log(w))
    expect_equal(lp$r, lp$legacy + logDir, tolerance = 1e-10)
    expect_equal(lp$cpp, lp$r, tolerance = 1e-10)
    expect_equal(lp$wrapper, lp$r, tolerance = 1e-10)
  }
})


test_that("evaluators share one class-rate normaliser exactly (#362)", {
  set.seed(20261004L)
  nTip <- 7L
  mat <- matrix(sample.int(2, nTip * 15, replace = TRUE) - 1L, nTip,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  for (j in seq_len(ncol(mat))) {
    if (length(unique(mat[, j])) < 2L) mat[1, j] <- 1L - mat[1, j]
  }
  mat[, 10:15] <- sample.int(3, nTip * 6, replace = TRUE) - 1L
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = c(1L, 2L, 5L, 9L))
  expect_gt(sum(mkd$type == "neomorphic"), 0)
  expect_gt(sum(mkd$type == "transformational"), 0)
  tree <- Preorder(NJTree(pd, edgeLengths = TRUE))
  tree$edge.length <- pmax(tree$edge.length, 0.05)

  # Classes differ in their neo:trans mix, so the normaliser is not 1.
  partition <- c(1L, 1L, 2L, 3L, 1L, rep(2:3, length.out = mkd$nChar - 5L))
  chain <- .PartitionedChain(tree, mkd, classW = c(0.2, 0.5, 0.3),
                             classRateLogSd = c(0.3, 0.6, 0.9),
                             partition = partition, rateNeo = 1.7)
  st <- get_mcmc_state(chain$statePtr)
  expect_true(is.finite(st$logLik))
  expect_identical(.PartitionedLogLik(chain), st$logLik)
  expect_identical(eval_full_loglik_cpp(chain$dataPtr, chain$statePtr),
                   st$logLik)
})


test_that("dirichlet_class_w_alpha, not beta_simplex, sets the class-weight step (#76)", {
  set.seed(76)
  nTip <- 7L
  mat <- matrix(sample.int(3, nTip * 12, replace = TRUE) - 1L, nTip,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  tree <- Preorder(NJTree(pd, edgeLengths = TRUE))
  tree$edge.length <- pmax(tree$edge.length, 0.05)

  MaxStep <- function(classAlpha, branchAlpha) {
    chain <- .PartitionedChain(tree, mkd, classW = c(0.2, 0.5, 0.3),
                               classRateLogSd = c(0.3, 0.6, 0.9),
                               partition = rep(1:3, length.out = mkd$nChar))
    steps <- vapply(seq_len(200), function(i) {
      before <- get_mcmc_state(chain$statePtr)$classW + 0
      do_move_cpp(chain$dataPtr, chain$statePtr, 32L, 0L,
                  classAlpha, branchAlpha, 3L, 1)
      max(abs(get_mcmc_state(chain$statePtr)$classW - before))
    }, double(1))
    max(steps)
  }
  # A concentration of 1e6 moves each weight by about 1e-3.
  expect_lt(MaxStep(classAlpha = 1e6, branchAlpha = 1), 0.01)
  expect_gt(MaxStep(classAlpha = 1, branchAlpha = 1e6), 0.01)

  moves <- list(list(name = "dirichlet_simplex_class_w",
                     type = "dirichlet_simplex_class_w"),
                list(name = "branch_lengths", type = "beta_simplex"))
  tuning <- list(dirichlet_class_w_alpha = 7, beta_simplex = 10)
  expect_identical(.BuildScaleTuningMatrix(list(tuning), moves)[1, 1], 7)

  # Every proposal rejected: the step shrinks, so alpha rises.
  adapted <- .AdaptTuning(
    tuning,
    acceptCount = c(dirichlet_simplex_class_w = 0L, branch_lengths = 0L),
    proposeCount = c(dirichlet_simplex_class_w = 50L, branch_lengths = 0L),
    moves
  )
  expect_gt(adapted$dirichlet_class_w_alpha, 7)
  expect_identical(adapted$beta_simplex, 10)
})
