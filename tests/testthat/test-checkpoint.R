# Tests for checkpointing (Phase 5: M-035)

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
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000L, thin = 5L, warmup = 200L,
                        check_every = 300L, checkpoint_file = cp_file))

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
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 1000L, thin = 5L, warmup = 200L,
                        check_every = 300L, checkpoint_file = cp_file))

  cp <- readRDS(cp_file)
  expect_equal(length(cp$runs), 2)

  # Each run should have chains with valid log_lik
  for (r in cp$runs) {
    expect_true(is.list(r$chains))
    expect_true(is.finite(r$chains[[1]]$log_lik))
    expect_true(r$saved_idx > 0)
  }
})


test_that("resume_mkprime continues from checkpoint", {
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
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 5000L, thin = 5L, warmup = 200L,
                        check_every = 300L, checkpoint_file = cp_file,
                        max_time = 0.5))

  expect_true(file.exists(cp_file))
  n_samples_before <- nrow(result1$samples)

  # Resume
  result2 <- resume_mkprime(cp_file, pd, tree)

  # Should have more samples (or at least as many)
  expect_gte(nrow(result2$samples), n_samples_before)
  expect_s3_class(result2, "MkPosterior")
})


test_that("Checkpoint without checkpoint_file does nothing", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(6611)
  # No checkpoint_file: should run normally
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, warmup = 200L,
                        check_every = 100L))

  expect_s3_class(result, "MkPosterior")
  expect_equal(nrow(result$samples), 60L)
})
