# Tests for parallel run mode (M-095 / M-096)
#
# The orchestration logic is tested with future::plan("sequential"), which
# exercises the full .RunParallelRuns() code path in the current process
# (futures are evaluated synchronously, so package internals are available).
#
# The multisession test additionally validates cross-process execution; it
# requires the MkPrime package to be *installed* (not just load_all()'d),
# because multisession workers start fresh R processes and load the package
# from the library.  That test is skipped during development via load_all().

.is_mkprime_installed <- function() {
  # Returns TRUE only when MkPrime is properly installed (not just loaded via
  # devtools::load_all()).  pkgload::is_dev_package() distinguishes the two.
  if (!requireNamespace("pkgload", quietly = TRUE)) return(TRUE)
  !isTRUE(pkgload::is_dev_package("MkPrime"))
}

test_that("MkPrimeMCMC() accepts parallel and pollInterval", {
  mcmc <- MkPrimeMCMC(parallel = TRUE, pollInterval = 5L)
  expect_true(isTRUE(mcmc$parallel))
  expect_equal(mcmc$pollInterval, 5L)
})

test_that("MkPrimeMCMC() validates parallel and pollInterval", {
  expect_error(MkPrimeMCMC(parallel = "yes"),  "parallel")
  expect_error(MkPrimeMCMC(parallel = NA),     "parallel")
  expect_error(MkPrimeMCMC(pollInterval = 0L), "pollInterval")
  expect_error(MkPrimeMCMC(pollInterval = -1L),"pollInterval")
})

test_that("parallel = TRUE with nRuns = 1 falls back to sequential", {
  skip_if_not_installed("future")
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)

  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 400L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                        parallel = TRUE))
  expect_s3_class(result, "MkPosterior")
})

test_that("parallel mode auto-assigns logFile when logFile = NULL", {
  # Regression: .BuildResult() was using mcmc$logFile (NULL) not logFilePaths
  # to determine streaming mode, causing a crash on r$samples subscript.
  skip_if_not_installed("future")
  library(ape)
  library(future)

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan("sequential")

  # No logFile supplied — samples loaded into memory after temp-log cleanup
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 400L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                        parallel = TRUE, pollInterval = 1L))

  expect_s3_class(result, "MkPosterior")
  expect_true(nrow(result$samples) > 0L)
})

test_that("parallel orchestration (sequential plan) returns valid MkPosterior", {
  skip_if_not_installed("future")
  library(ape)
  library(future)

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  # Sequential plan: futures run synchronously in the current process,
  # so package internals are available.  Exercises the full orchestration
  # code path without requiring an installed package.
  future::plan("sequential")

  logFile <- tempfile(fileext = ".log")
  on.exit(unlink(c(logFile,
                   sub("\\.log$", "_1.log", logFile),
                   sub("\\.log$", "_2.log", logFile))), add = TRUE)

  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(
      nRuns       = 2L,
      nIter       = 400L,
      maxWarmup   = 200L,
      minWarmup   = 200L,
      autoTune    = FALSE,
      logFile     = logFile,
      parallel    = TRUE,
      pollInterval = 1L
    ))

  expect_s3_class(result, "MkPosterior")
  expect_true(!is.null(result$logFile))
  expect_true(result$nSamples > 0L)
  expect_false(is.null(result$stop_reason))

  samp <- ReadMkLog(result$logFile)
  expect_true(nrow(samp) > 0L)
  expect_true("log_posterior" %in% colnames(samp))
})

test_that("parallel mode saves checkpoint when checkpointFile is set", {
  # Regression: checkpoint save was inside sequential else-branch only;
  # parallel = TRUE + checkpointFile silently wrote nothing.
  skip_if_not_installed("future")
  library(ape)
  library(future)

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan("sequential")

  cp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(cp_file), add = TRUE)

  RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 400L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                        parallel = TRUE, pollInterval = 1L,
                        checkpointFile = cp_file))

  expect_true(file.exists(cp_file))
  cp <- readRDS(cp_file)
  expect_equal(length(cp$runs), 2L)
  for (r in cp$runs) {
    expect_true(is.finite(r$chains[[1]]$log_lik))
    expect_true(r$saved_idx > 0L)
  }
})

test_that("parallel orchestration (multisession, 2 workers) returns valid MkPosterior", {
  skip_if_not_installed("future")
  skip_if(!.is_mkprime_installed(),
          "MkPrime not installed — multisession workers need installed package")
  library(ape)
  library(future)

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan("multisession", workers = 2L)

  logFile <- tempfile(fileext = ".log")
  on.exit(unlink(c(logFile,
                   sub("\\.log$", "_1.log", logFile),
                   sub("\\.log$", "_2.log", logFile))), add = TRUE)

  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(
      nRuns       = 2L,
      nIter       = 400L,
      maxWarmup   = 200L,
      minWarmup   = 200L,
      autoTune    = FALSE,
      logFile     = logFile,
      parallel    = TRUE,
      pollInterval = 2L
    ))

  expect_s3_class(result, "MkPosterior")
  expect_true(result$nSamples > 0L)
  expect_false(is.null(result$stop_reason))

  samp <- ReadMkLog(result$logFile)
  expect_true(nrow(samp) > 0L)
  expect_true("log_posterior" %in% colnames(samp))
})
