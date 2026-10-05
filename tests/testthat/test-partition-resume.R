# A resumed K-class chain must keep scoring each character under its own
# class (#239), honour classRateConcentration (#240) and apply the
# partition-aware guards the fresh path applies (#242).

.PartitionResumeFixture <- function() {
  set.seed(11)
  nTip <- 9L
  nChar <- 20L
  mat <- matrix(sample(0:2, nTip * nChar, TRUE), nTip,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1, j] <- (mat[1, j] + 1L) %% 3L
  }
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  tree <- Preorder(NJTree(pd, edgeLengths = TRUE))
  tree$edge.length <- pmax(tree$edge.length, 0.05)
  partition <- c(rep(1L, 10L), rep(2L, mkd$nChar - 10L))
  mkdK <- mkd
  mkdK$partitions <- .BuildPartitions(mkd, partition)
  list(pd = pd, mkd = mkd, mkdK = mkdK, tree = tree, partition = partition)
}

.PartitionResumeMcmc <- function(ckp, lg, nIter) {
  MkPrimeMCMC(nRuns = 1L, nChains = 1L, nIter = nIter, thin = 5L,
              maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
              checkEvery = 100L, nCore = 1L, maxTime = 60,
              checkpointFile = ckp, logFile = lg)
}

.RecomputeLogLik <- function(chain, mkd, model) {
  dp <- .InitMcmcData(mkd, model)
  eval_full_loglik_cpp(dp, .DeserialiseChain(chain))
}

test_that("resumed K-class chains keep their partition (#239)", {
  skip_on_cran()
  fx <- .PartitionResumeFixture()
  td <- tempfile("mkp_part_resume_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  unlinkArg <- c("shape", "ratemultiplier")

  ckp <- file.path(td, "a.ckp")
  lg <- file.path(td, "a.log")
  set.seed(2)
  suppressWarnings(suppressMessages(RunMkPrime(
    fx$pd, fx$tree, mcmc = .PartitionResumeMcmc(ckp, lg, 300L),
    partition = fx$partition, unlink = unlinkArg)))
  model <- readRDS(ckp)$model

  Copy <- function(stem) {
    newCkp <- file.path(td, paste0(stem, ".ckp"))
    newLog <- file.path(td, paste0(stem, ".log"))
    file.copy(lg, newLog)
    ck <- readRDS(ckp)
    ck$mcmc$checkpointFile <- newCkp
    ck$mcmc$logFile <- newLog
    saveRDS(ck, newCkp)
    c(ckp = newCkp, log = newLog)
  }
  Check <- function(file) {
    ck <- readRDS(file)
    expect_gt(ck$iter, 300L)
    chain <- ck$runs[[1]]$chains[[1]]
    expect_equal(chain$log_lik, .RecomputeLogLik(chain, fx$mkdK, model),
                 tolerance = 1e-8)
  }

  explicit <- Copy("b")
  suppressWarnings(suppressMessages(ResumeMkPrime(
    explicit[["ckp"]], fx$pd, fx$tree, mcmc = list(nIter = 400L))))
  Check(explicit[["ckp"]])

  auto <- Copy("c")
  suppressWarnings(suppressMessages(RunMkPrime(
    fx$pd, fx$tree,
    mcmc = .PartitionResumeMcmc(auto[["ckp"]], auto[["log"]], 400L),
    partition = fx$partition, unlink = unlinkArg)))
  Check(auto[["ckp"]])
})

test_that("resume applies the marginal_k partition guard (#242)", {
  skip_on_cran()
  fx <- .PartitionResumeFixture()
  td <- tempfile("mkp_part_resume_mk_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ckp <- file.path(td, "a.ckp")
  lg <- file.path(td, "a.log")
  set.seed(3)
  suppressWarnings(suppressMessages(RunMkPrime(
    fx$pd, fx$tree, mcmc = .PartitionResumeMcmc(ckp, lg, 200L),
    model = MkPrimeModel(kPrimePrior = "geometric"),
    partition = fx$partition, unlink = "shape")))
  expect_error(
    suppressWarnings(suppressMessages(ResumeMkPrime(
      ckp, fx$pd, fx$tree,
      model = MkPrimeModel(kPrimePrior = "geometric",
                           likelihoodMode = "marginal_k"),
      mcmc = list(nIter = 300L)))),
    "partition API")
})

test_that("a partitioned state rejects data built without its partition", {
  fx <- .PartitionResumeFixture()
  model <- .FinalizeModel(MkPrimeModel(), fx$tree, fx$mkd)
  spec <- .ValidatePartitionArgs(fx$partition, "ratemultiplier", fx$mkd)
  state <- .InitStatePartitioned(fx$tree, fx$mkdK, model, spec)
  chain <- list(
    edge = state$tree$edge, rel_br_lengths = state$rel_br_lengths,
    tree_length = state$tree_length, rate_loss = state$rate_loss,
    rate_log_sd = state$rate_log_sd, rate_neo = 1, p = state$p,
    kPrime = state$kPrime, log_lik = state$log_lik,
    log_prior = state$log_prior, class_rate_log_sd = state$class_rate_log_sd,
    class_w = state$class_w, class_rate = state$class_rate,
    nChar_c = state$nChar_c, eta_neo = 1)
  statePtr <- .DeserialiseChain(chain)
  expect_error(fill_partition_cache(.InitMcmcData(fx$mkd, model), statePtr),
               "2 classes but its data has 1")
  expect_silent(fill_partition_cache(.InitMcmcData(fx$mkdK, model), statePtr))
})

test_that("classRateConcentration reaches the C++ prior (#240)", {
  fx <- .PartitionResumeFixture()
  spec <- .ValidatePartitionArgs(fx$partition, "ratemultiplier", fx$mkd)
  for (alpha in c(0.2, 5)) {
    model <- .FinalizeModel(
      MkPrimeModel(classRateConcentration = alpha), fx$tree, fx$mkd)
    state <- .InitStatePartitioned(fx$tree, fx$mkdK, model, spec)
    state$class_w <- c(0.25, 0.75)
    state$class_rate <- .PartitionWToClassRate(state$class_w, state$nChar_c)
    dp <- .InitMcmcData(fx$mkdK, model)
    statePtr <- init_mcmc_state(
      state$tree$edge[, 1], state$tree$edge[, 2],
      state$rel_br_lengths, state$tree_length,
      state$rate_loss, state$rate_log_sd, 1.0, state$p,
      as.integer(state$kPrime), state$log_lik, 0,
      classRateLogSd = state$class_rate_log_sd, classW = state$class_w,
      classRate = state$class_rate, nCharPerClass = state$nChar_c)
    rPrior <- LogPrior(state, model, fx$mkdK)
    expect_equal(
      eval_log_prior_partitioned_cpp(dp, statePtr,
                                     classRateLogSd = state$class_rate_log_sd,
                                     classW = state$class_w,
                                     classRate = state$class_rate,
                                     etaNeo = 1),
      rPrior, tolerance = 1e-10)

    # The first MH ratio uses the C++ prior, whatever the state was given.
    fill_partition_cache(dp, statePtr)
    expect_equal(get_mcmc_state(statePtr)$logPrior, rPrior,
                 tolerance = 1e-10)
  }
})
