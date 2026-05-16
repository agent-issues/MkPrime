# Tests for parallel run mode (M-095 / M-096)
#
# Parallel orchestration uses callr::r_bg() to spawn fresh R worker processes.
# Tests that exercise nCore > 1 therefore require:
#   (a) the callr package to be installed, and
#   (b) MkPrime itself to be *installed* (not just load_all()'d), because
#       callr workers start fresh R processes and load the package from the
#       library.
# Serial-path tests (nCore = 1) have no such constraint.

.is_mkprime_installed <- function() {
  # Returns TRUE only when MkPrime is properly installed (not just loaded via
  # devtools::load_all()).  pkgload::is_dev_package() distinguishes the two.
  if (!requireNamespace("pkgload", quietly = TRUE)) return(TRUE)
  !isTRUE(pkgload::is_dev_package("MkPrime"))
}

test_that("MkPrimeMCMC() accepts nCore and pollInterval", {
  mcmc <- suppressWarnings(MkPrimeMCMC(nCore = 2L, pollInterval = 5L))
  expect_equal(mcmc$nCore, 2L)
  expect_equal(mcmc$pollInterval, 5L)
})

test_that("MkPrimeMCMC() validates nCore and pollInterval", {
  expect_error(MkPrimeMCMC(nCore = 0L),                      "nCore")
  expect_error(MkPrimeMCMC(nCore = -1L),                     "nCore")
  expect_error(suppressWarnings(MkPrimeMCMC(nCore = "two")), "nCore")
  expect_error(MkPrimeMCMC(pollInterval = 0L),               "pollInterval")
  expect_error(MkPrimeMCMC(pollInterval = -1L),              "pollInterval")
})

test_that("nCore defaults to getOption('mc.cores', 1L)", {
  old <- getOption("mc.cores")
  on.exit(options(mc.cores = old), add = TRUE)

  options(mc.cores = 3L)
  mcmc3 <- suppressWarnings(MkPrimeMCMC())
  expect_equal(mcmc3$nCore, 3L)

  options(mc.cores = NULL)
  mcmc1 <- MkPrimeMCMC()
  expect_equal(mcmc1$nCore, 1L)
})

test_that("nCore = 1 runs serially", {
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)

  result <- suppressWarnings(RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 400L, maxWarmup = 200L,
                       minWarmup = 200L, autoTune = FALSE, nCore = 1L)))
  expect_s3_class(result, "MkPosterior")
})

test_that("nCore > 1 with nRuns = 1 silently runs serially", {
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)

  result <- suppressWarnings(RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 400L, maxWarmup = 200L,
                       minWarmup = 200L, autoTune = FALSE, nCore = 4L)))
  expect_s3_class(result, "MkPosterior")
})

test_that("parallel orchestration with nCore = 2 returns valid MkPosterior", {
  skip_if_not_installed("callr")
  skip_if_not(.is_mkprime_installed(),
              "MkPrime not installed — callr workers need installed package")
  library("ape")

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)

  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(
      nRuns       = 2L,
      nIter       = 400L,
      maxWarmup   = 200L,
      minWarmup   = 200L,
      autoTune    = FALSE,
      nCore       = 2L,
      pollInterval = 1L
    ))

  expect_s3_class(result, "MkPosterior")
  expect_true(result$nSamples > 0L)
})

test_that("mc.cores option triggers parallel mode", {
  skip_if_not_installed("callr")
  skip_if_not(.is_mkprime_installed(),
              "MkPrime not installed — callr workers need installed package")
  library("ape")

  old <- getOption("mc.cores")
  on.exit(options(mc.cores = old), add = TRUE)

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)

  options(mc.cores = 2L)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(
      nRuns       = 2L,
      nIter       = 400L,
      maxWarmup   = 200L,
      minWarmup   = 200L,
      autoTune    = FALSE,
      pollInterval = 1L
    ))

  expect_s3_class(result, "MkPosterior")
  expect_true(result$nSamples > 0L)
})

test_that("pool dispatch: nRuns > nCore launches in waves", {
  # nRuns = 4, nCore = 2 → at most 2 workers active at any time, rolling
  # launch as each finishes. All 4 runs must still complete and contribute
  # samples to the final MkPosterior.
  skip_if_not_installed("callr")
  skip_if_not(.is_mkprime_installed(),
              "MkPrime not installed — callr workers need installed package")
  library("ape")

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)

  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(
      nRuns        = 4L,
      nIter        = 300L,
      maxWarmup    = 150L,
      minWarmup    = 150L,
      autoTune     = FALSE,
      nCore        = 2L,
      pollInterval = 1L
    ))

  expect_s3_class(result, "MkPosterior")
  expect_true(result$nSamples > 0L)
  # All four runs must have run to completion and contributed samples.
  expect_equal(result$nRuns, 4L)
  expect_length(result$per_run, 4L)
  perRunRows <- vapply(result$per_run,
                       function(r) nrow(r$samples),
                       integer(1L))
  expect_true(all(perRunRows > 0L))
})

test_that("parallel mode auto-assigns logFile when logFile = NULL", {
  # Regression: .BuildResult() was using mcmc$logFile (NULL) not logFilePaths
  # to determine streaming mode, causing a crash on r$samples subscript.
  skip_if_not_installed("callr")
  skip_if_not(.is_mkprime_installed(),
              "MkPrime not installed — callr workers need installed package")
  library("ape")

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)

  # No logFile supplied — samples loaded into memory after temp-log cleanup
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 400L, maxWarmup = 200L,
                       minWarmup = 200L, autoTune = FALSE,
                       nCore = 2L, pollInterval = 1L))

  expect_s3_class(result, "MkPosterior")
  expect_true(nrow(result$samples) > 0L)
})

test_that("parallel mode saves checkpoint when checkpointFile is set", {
  # Regression: checkpoint save was inside sequential else-branch only;
  # nCore > 1 + checkpointFile silently wrote nothing.
  skip_if_not_installed("callr")
  skip_if_not(.is_mkprime_installed(),
              "MkPrime not installed — callr workers need installed package")
  library("ape")

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(cp_file), add = TRUE)

  suppressWarnings(RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 400L, maxWarmup = 200L,
                       minWarmup = 200L, autoTune = FALSE,
                       nCore = 2L, pollInterval = 1L,
                       checkpointFile = cp_file)))

  expect_true(file.exists(cp_file))
  cp <- readRDS(cp_file)
  expect_equal(length(cp$runs), 2L)
  for (r in cp$runs) {
    expect_true(is.finite(r$chains[[1]]$log_lik))
    expect_true(r$saved_idx > 0L)
  }
})
