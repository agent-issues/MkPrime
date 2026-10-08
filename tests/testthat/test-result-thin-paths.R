# Results and recovery report each run's own thin and the files the run
# actually wrote (#319).

.ThinPathData <- function() {
  set.seed(319)
  tree <- ape::rtree(8)
  mat <- replicate(15, as.character(
    ape::rTraitDisc(tree, k = 3, rate = 2, states = c("0", "1", "2"))
  ))
  rownames(mat) <- tree$tip.label
  list(tree = tree, pd = MatrixToPhyDat(mat))
}

.ThinPathMcmc <- function(...) {
  do.call(MkPrimeMCMC, utils::modifyList(list(
    nIter = 400L, minWarmup = 100L, maxWarmup = 100L, autoTune = FALSE,
    nRuns = 1L, nChains = 1L, nCore = 1L, thin = 10L, maxTime = 60
  ), list(...)))
}

.RunThinPath <- function(data, ...) {
  allow_warning(RunMkPrime(data$pd, data$tree, mcmc = .ThinPathMcmc(...)),
                "stabilis")
}

test_that(".PostBurninData discards each run's trees at its own treeThin", {
  m <- matrix(seq_len(20), 10, 2, dimnames = list(NULL, c("a", "b")))
  posterior <- list(
    samples = rbind(m, m), mcmc = list(thin = 10L), treeThin = 10L,
    nRuns = 2L, burnin = 4L,
    per_run = list(
      list(samples = m, trees = as.list(1:10), thin = 10L, treeThin = 10L),
      list(samples = m, trees = as.list(1:5), thin = 10L, treeThin = 20L)
    )
  )
  pb <- MkPrime:::.PostBurninData(posterior)
  expect_equal(unlist(pb$per_run[[1]]$trees), 5:10)
  expect_equal(unlist(pb$per_run[[2]]$trees), 3:5)
})

test_that("results carry each run's adapted thin and treeThin", {
  skip_on_cran()
  data <- .ThinPathData()
  dir <- withr::local_tempdir()
  logFile <- file.path(dir, "thin.log")
  set.seed(3191)
  .RunThinPath(data, nRuns = 2L, logFile = logFile, treeFile = FALSE)

  ckpFile <- file.path(dir, "thin.ckp")
  cp <- readRDS(ckpFile)
  cp$runs[[2]]$thin <- 20L
  cp$runs[[2]]$treeThin <- 40L
  saveRDS(cp, ckpFile)
  result <- ResumeMkPrime(ckpFile, data$pd)

  expect_equal(result$thin, 10L)
  expect_equal(result$treeThin, 10L)
  expect_equal(result$per_run[[2]]$thin, 20L)
  expect_equal(result$per_run[[2]]$treeThin, 40L)
})

test_that("MkPrimeRecover(logFile) reports treeThin", {
  skip_on_cran()
  data <- .ThinPathData()
  dir <- withr::local_tempdir()
  logFile <- file.path(dir, "rec.log")
  set.seed(3192)
  .RunThinPath(data, logFile = logFile, treeThin = 30L,
               treeFile = file.path(dir, "rec_trees.nwk"))

  rec <- MkPrimeRecover(logFile)
  expect_equal(rec$thin, 10L)
  expect_equal(rec$treeThin, 30L)
  pb <- MkPrime:::.PostBurninData(rec, burnin = 6L)
  expect_length(pb$trees, length(rec$trees) - 2L)
})

test_that("overwrite = TRUE removes a discarded job's extra run files", {
  skip_on_cran()
  data <- .ThinPathData()
  dir <- withr::local_tempdir()
  logFile <- file.path(dir, "ow.log")
  treeFile <- file.path(dir, "ow.nwk")
  set.seed(3193)
  .RunThinPath(data, nRuns = 2L, nIter = 200L, logFile = logFile,
               treeFile = treeFile)
  stale <- file.path(dir, c("ow_1.log", "ow_2.log", "ow_1.nwk", "ow_2.nwk"))
  expect_true(all(file.exists(stale)))

  allow_warning(
    RunMkPrime(data$pd, data$tree, overwrite = TRUE, mcmc = .ThinPathMcmc(
      nIter = 200L, logFile = logFile, treeFile = treeFile
    )),
    "stabilis"
  )
  expect_false(any(file.exists(stale)))
  expect_equal(MkPrimeRecover(logFile)$logFile, logFile)
})

test_that("a resume that replaces a lost log records the replacement", {
  skip_on_cran()
  data <- .ThinPathData()
  dir <- withr::local_tempdir()
  logFile <- file.path(dir, "lost.log")
  ckpFile <- file.path(dir, "lost.ckp")
  set.seed(3194)
  .RunThinPath(data, logFile = logFile)
  expect_gt(readRDS(ckpFile)$runs[[1]]$tree_saved_idx, 0L)
  unlink(logFile)

  cancelFile <- file.path(dir, "cancel")
  file.create(cancelFile)
  allow_warning(
    ResumeMkPrime(ckpFile, data$pd,
                  mcmc = list(nIter = 600L, cancelFile = cancelFile)),
    "without drawing"
  )
  cp <- readRDS(ckpFile)
  expect_match(cp$logFilePaths, "mkp_resume_")
  expect_true(file.exists(cp$logFilePaths))
  expect_equal(cp$runs[[1]]$tree_saved_idx, cp$runs[[1]]$saved_idx)
})

test_that("MkPrimeRecover() returns a temp-log run's trees and runs", {
  skip_on_cran()
  data <- .ThinPathData()
  dir <- withr::local_tempdir()
  logFile <- file.path(dir, "tmp.log")
  ckpFile <- file.path(dir, "tmp.ckp")
  set.seed(3195)
  result <- .RunThinPath(data, nRuns = 2L, logFile = logFile)

  mkpEnv <- environment(RunMkPrime)$.mkp_env
  withr::defer(mkpEnv$recovery <- NULL)
  cp <- readRDS(ckpFile)
  mkpEnv$recovery <- list(
    logFiles = result$logFile, paramNames = cp$paramNames, model = cp$model,
    data = NULL, mcmc = cp$mcmc, isTempLog = TRUE,
    tempFiles = MkPrime:::.TempRunFiles(result$logFile, ckpFile, 2L),
    time = Sys.time()
  )
  rec <- MkPrimeRecover()
  expect_equal(rec$nRuns, 2L)
  expect_length(rec$per_run, 2L)
  expect_length(rec$trees, length(result$trees))
  expect_false(file.exists(ckpFile))
  expect_false(any(file.exists(result$logFile)))
})

test_that("resume finds per-run checkpoints whose name has regex characters", {
  skip_on_cran()
  data <- .ThinPathData()
  dir <- withr::local_tempdir()
  ckpFile <- file.path(dir, "run+1.ckp")
  set.seed(3196)
  .RunThinPath(data, nRuns = 2L, nIter = 200L,
               logFile = file.path(dir, "run+1.log"), checkpointFile = ckpFile)
  cp <- readRDS(ckpFile)
  perRun <- file.path(dir, c("run+1_1.ckp", "run+1_2.ckp"))
  for (i in 1:2) {
    MkPrime:::.SaveCheckpoint(
      list(cp$runs[[i]]), cp$mcmc, cp$iter, cp$paramNames, perRun[i],
      moveWeights = cp$moveWeights, phase = cp$phase, model = cp$model
    )
  }
  unlink(ckpFile)
  expect_s3_class(ResumeMkPrime(ckpFile, data$pd), "MkPosterior")

  # Run 1's state is not to be taken from run 2's checkpoint.
  unlink(c(ckpFile, perRun[1]))
  expect_error(
    allow_warning(ResumeMkPrime(ckpFile, data$pd), "missing state"),
    "1 of 2 runs"
  )
})
