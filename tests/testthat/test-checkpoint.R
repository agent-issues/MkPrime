# Tests for checkpointing (Phase 5: M-035)
skip_slow_tests()

test_that("Checkpoint file is written at check intervals", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(cp_file), add = TRUE)

  set.seed(8901)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                        checkEvery = 300L, checkpointFile = cp_file))

  # Checkpoint should exist

  expect_true(file.exists(cp_file))

  # Read and verify structure
  cp <- readRDS(cp_file)
  expect_equal(cp$version, 1L)
  expect_true(inherits(cp$timestamp, "POSIXct"))
  expect_true(is.list(cp$runs))
  expect_true(cp$iter > 0)
})


test_that("Checkpoint contains valid run state", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(cp_file), add = TRUE)

  set.seed(1450)
  RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 1000L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                        checkEvery = 300L, checkpointFile = cp_file))

  cp <- readRDS(cp_file)
  expect_equal(length(cp$runs), 2)

  # Each run should have chains with valid log_lik
  for (r in cp$runs) {
    expect_true(is.list(r$chains))
    expect_true(is.finite(r$chains[[1]]$log_lik))
    expect_true(r$saved_idx > 0)
  }
})


test_that("ResumeMkPrime continues from checkpoint", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(cp_file), add = TRUE)

  # Run with time limit to stop early and checkpoint
  set.seed(7766)
  result1 <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 5000L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                        checkEvery = 300L, checkpointFile = cp_file,
                        maxTime = 0.5))

  expect_true(file.exists(cp_file))
  n_samples_before <- nrow(result1$samples)

  # Resume
  result2 <- ResumeMkPrime(cp_file, pd, tree)

  # Should have more samples (or at least as many)
  expect_gte(nrow(result2$samples), n_samples_before)
  expect_s3_class(result2, "MkPosterior")
})


test_that("RunMkPrime auto-resumes from existing checkpoint", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(cp_file), add = TRUE)

  # Run with time limit to create a checkpoint
  set.seed(4821)
  mcmcConf <- MkPrimeMCMC(nRuns = 1L, nIter = 5000L, thin = 5L,
                           maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE, checkEvery = 300L,
                           checkpointFile = cp_file, maxTime = 0.5)
  result1 <- RunMkPrime(pd, tree, mcmc = mcmcConf)
  expect_true(file.exists(cp_file))
  n1 <- nrow(result1$samples)

  # Calling RunMkPrime again with same mcmc should auto-resume
  result2 <- RunMkPrime(pd, tree, mcmc = mcmcConf)
  expect_s3_class(result2, "MkPosterior")
  expect_gte(nrow(result2$samples), n1)
})


test_that("RunMkPrime overwrite = TRUE ignores existing checkpoint", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(cp_file), add = TRUE)

  # Create a checkpoint
  set.seed(3952)
  mcmcConf <- MkPrimeMCMC(nRuns = 1L, nIter = 1000L, thin = 5L,
                           maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE, checkEvery = 300L,
                           checkpointFile = cp_file)
  RunMkPrime(pd, tree, mcmc = mcmcConf)
  expect_true(file.exists(cp_file))

  # overwrite = TRUE starts fresh
  set.seed(6149)
  result <- RunMkPrime(pd, tree, mcmc = mcmcConf, overwrite = TRUE)
  expect_s3_class(result, "MkPosterior")
  expect_equal(nrow(result$samples), 160L)
})


test_that("Checkpoint without checkpointFile does nothing", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(6611)
  # No checkpointFile: should run normally
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                        checkEvery = 100L))

  expect_s3_class(result, "MkPosterior")
  expect_equal(nrow(result$samples), 60L)
})
