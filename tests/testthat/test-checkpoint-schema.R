# A checkpoint carries no schema version for its `model` and `mcmc`, and its
# chain state is serialised field by field: #55, #56, #108 and #225.

SchemaFixture <- function(model = MkPrimeModel(), nRuns = 1L, ...) {
  mat <- matrix(c(0, 1, 0, 1, 2, 0, 0, 1, 1, 2, 1, 0), 6, 2,
                dimnames = list(paste0("t", 1:6), NULL))
  pd <- MatrixToPhyDat(mat)
  tree <- ape::read.tree(
    text = "(((t1:0.1,t2:0.2):0.1,t3:0.1):0.1,((t4:0.1,t5:0.2):0.1,t6:0.3):0.1);")
  dir <- tempfile("ckp")
  dir.create(dir)
  logFile <- file.path(dir, "run.log")
  ckpFile <- file.path(dir, "run.ckp")
  set.seed(5500)
  allow_warning(
    RunMkPrime(pd, tree, model = model, nRuns = nRuns, nIter = 150L,
               thin = 10L,
               maxWarmup = 50L, minWarmup = 50L, autoTune = FALSE,
               logFile = logFile, checkpointFile = ckpFile, ...),
    "without stabilisation")
  list(pd = pd, tree = tree, dir = dir, logFile = logFile, ckpFile = ckpFile)
}

# Remove a field, and extend the run so that the resume has work to do.
StripField <- function(ckpFile, part, field) {
  ck <- readRDS(ckpFile)
  ck[[part]][[field]] <- NULL
  ck$mcmc$nIter <- 300L
  saveRDS(ck, ckpFile)
}


LastSample <- function(logFiles) {
  max(as.numeric(rownames(ReadMkLog(logFiles))))
}


test_that("Resume fills a missing mcmc$nCore instead of crashing (#225)", {
  # The unguarded test is reached only with more than one run.
  f <- SchemaFixture(nRuns = 2L, nCore = 1L)
  on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)
  StripField(f$ckpFile, "mcmc", "nCore")
  logFiles <- readRDS(f$ckpFile)$logFilePaths

  expect_warning(
    resumed <- ResumeMkPrime(f$ckpFile, f$pd),
    "nCore"
  )
  expect_s3_class(resumed, "MkPosterior")
  expect_gt(LastSample(logFiles), 150)
})


test_that("Resume restores the prior variant a checkpoint omits (#55)", {
  f <- SchemaFixture()
  on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)
  expect_identical(readRDS(f$ckpFile)$model$priorVariant, "unconditional")
  StripField(f$ckpFile, "model", "priorVariant")

  expect_warning(
    resumed <- ResumeMkPrime(f$ckpFile, f$pd),
    "priorVariant"
  )
  expect_identical(resumed$model$priorVariant, "unconditional")

  # The default follows the checkpoint's own k' prior.
  ck <- readRDS(f$ckpFile)
  ck$model$kPrimePrior <- "logseries"
  ck$model$priorVariant <- NULL
  expect_warning(migrated <- .MigrateCheckpoint(ck), "priorVariant")
  expect_identical(migrated$model$priorVariant, "conditional")
})


test_that("Chain serialisation round-trips every sampled field (#56)", {
  nEdge <- 9L
  edge <- ape::read.tree(
    text = "(((t1,t2),t3),((t4,t5),t6));")
  edge <- TreeTools::Preorder(edge)$edge
  for (partitioned in c(FALSE, TRUE)) {
    classArgs <- if (partitioned) {
      list(classRateLogSd = c(0.3, 0.7), classW = c(0.4, 0.6),
           classRate = c(0.8, 1.3), nCharPerClass = c(1L, 1L),
           etaNeo = 1.7, useHyperpriorOnSigma = TRUE, hyperTau = 0.9,
           classZ = c(-0.2, 0.4))
    } else {
      list()
    }
    ptr <- do.call(init_mcmc_state, c(list(
      edge[, 1], edge[, 2], rep(1 / nEdge, nEdge), 2.5, 1.1, 0.6, 1.4, 0.3,
      c(2L, 3L), -12.5, -3.25, betaScale = 0.8, kprimeAlpha = 3,
      kprimeBeta = 0.25), classArgs))
    Fields <- function(p) {
      s <- get_mcmc_state(p)
      s[!startsWith(names(s), "diag")]
    }
    expect_identical(Fields(.DeserialiseChain(.SerialiseChain(ptr))),
                     Fields(ptr))
  }
})


test_that("A fresh run starts from the model's k' hyperparameters (#56)", {
  f <- SchemaFixture(model = MkPrimeModel(kPrimePrior = "beta_geometric",
                                          kprimeAlpha = 3, kprimeBeta = 0.25))
  on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)
  mkd <- MkPrimeData(f$pd)
  model <- .FinalizeModel(
    MkPrimeModel(kPrimePrior = "beta_geometric", kprimeAlpha = 3,
                 kprimeBeta = 0.25), f$tree, mkd)
  run <- .InitRun(TreeTools::Preorder(f$tree), mkd, model,
                  MkPrimeMCMC(nChains = 1L), moves = list())
  expect_identical(run$chains[[1]]$kprime_alpha, 3)
  expect_identical(run$chains[[1]]$kprime_beta, 0.25)

  # ...and the checkpoint carries the sampled values, not 1.
  ch <- readRDS(f$ckpFile)$runs[[1]]$chains[[1]]
  expect_false(is.null(ch$kprime_alpha))
  expect_false(identical(c(ch$kprime_alpha, ch$kprime_beta), c(1, 1)))
})


test_that("Auto-resume honours the caller's run length (#108)", {
  f <- SchemaFixture()
  on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)

  # Same call with a longer nIter: previously the checkpoint's nIter won and
  # the resubmission did nothing.
  set.seed(5501)
  expect_warning(
    resumed <- RunMkPrime(f$pd, f$tree, nRuns = 1L, nIter = 300L,
                          thin = 7L, maxWarmup = 50L, minWarmup = 50L,
                          autoTune = FALSE, logFile = f$logFile,
                          checkpointFile = f$ckpFile),
    "thin"
  )
  expect_gt(LastSample(f$logFile), 150)
  expect_identical(readRDS(f$ckpFile)$mcmc$thin, 10L)
})


test_that(".ResumeMcmc() reports what it applies and what it keeps (#108)", {
  ck <- MkPrimeMCMC(nIter = 10000L, thin = 5L)
  local_mkp_verbosity(1)
  expect_message(
    merged <- .ResumeMcmc(ck, list(nIter = 200L, maxTime = 60)),
    "nIter"
  )
  expect_identical(merged$nIter, 200L)
  expect_identical(merged$maxTime, 60)

  # Unresolved defaults in a whole object are not requests.
  expect_no_warning(.ResumeMcmc(ck, MkPrimeMCMC(nIter = 10000L, thin = 5L)))
  expect_warning(.ResumeMcmc(ck, list(nRuns = 3L)), "nRuns")
  expect_identical(suppressWarnings(.ResumeMcmc(ck, list(nRuns = 3L)))$nRuns,
                   ck$nRuns)
})


test_that("MkPrimeRecover() finds per-run tree files (#108)", {
  f <- SchemaFixture()
  on.exit(unlink(f$dir, recursive = TRUE), add = TRUE)
  recovered <- MkPrimeRecover(f$logFile)
  expect_gt(length(recovered$trees), 0L)

  dir <- tempfile("ckp2")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  logFile <- file.path(dir, "multi.log")
  set.seed(5502)
  result <- RunMkPrime(f$pd, f$tree, nRuns = 2L, nIter = 150L, thin = 10L,
                       maxWarmup = 50L, minWarmup = 50L, autoTune = FALSE,
                       logFile = logFile, nCore = 1L)
  unlink(file.path(dir, "multi.ckp"))
  for (lf in list(logFile, result$logFile)) {
    recovered <- MkPrimeRecover(lf)
    expect_length(recovered$per_run, 2L)
    expect_gt(length(recovered$trees), 0L)
    expect_gt(length(recovered$per_run[[2]]$trees), 0L)
  }
})


test_that(".SaveCheckpoint() warns when it cannot replace the file (#108)", {
  target <- tempfile("ckp-dir")
  dir.create(target)
  writeLines("x", file.path(target, "occupied"))
  on.exit(unlink(c(target, paste0(target, ".tmp")), recursive = TRUE),
          add = TRUE)
  # Renaming a file over a non-empty directory fails on every platform.
  expect_warning(
    .SaveCheckpoint(list(), MkPrimeMCMC(), 1L, character(0), target),
    "Checkpoint not updated"
  )
})
