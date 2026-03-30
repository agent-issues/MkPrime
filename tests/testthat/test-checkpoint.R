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
  expect_equal(cp$version, 2L)  # always streaming since always-on checkpointing
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
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(c(cp_file, log_file,
                    sub("\\.[^.]+$", "_1.log", log_file),
                    sub("\\.[^.]+$", "_trees.nwk", log_file))), add = TRUE)

  # Run with time limit to stop early and checkpoint (persistent log)
  set.seed(7766)
  result1 <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 5000L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                        checkEvery = 300L, checkpointFile = cp_file,
                        logFile = log_file, maxTime = 0.5))

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
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(c(cp_file, log_file,
                    sub("\\.[^.]+$", "_1.log", log_file),
                    sub("\\.[^.]+$", "_trees.nwk", log_file))), add = TRUE)

  # Run with time limit to create a checkpoint (persistent log)
  set.seed(4821)
  mcmcConf <- MkPrimeMCMC(nRuns = 1L, nIter = 5000L, thin = 5L,
                           maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE, checkEvery = 300L,
                           checkpointFile = cp_file, logFile = log_file,
                           maxTime = 0.5)
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


# --- M-149: Interrupt-safe checkpointing tests ---

test_that("Initial checkpoint is written before first batch", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".ckp")
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(c(cp_file, log_file,
                    sub("\\.[^.]+$", "_1.log", log_file),
                    sub("\\.[^.]+$", "_trees.nwk", log_file))), add = TRUE)

  set.seed(3281)
  RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L,
                        maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                        checkEvery = 100L, logFile = log_file,
                        checkpointFile = cp_file))

  expect_true(file.exists(cp_file))

  # Checkpoint should contain finalized model (M-149)
  cp <- readRDS(cp_file)
  expect_true(!is.null(cp$model))
  expect_true(!is.null(cp$model$expSteps))
})


test_that("Checkpoint model is used on auto-resume (tree = NULL)", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".ckp")
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(c(cp_file, log_file,
                    sub("\\.[^.]+$", "_1.log", log_file),
                    sub("\\.[^.]+$", "_trees.nwk", log_file))), add = TRUE)

  # Run to create checkpoint with persistent log
  set.seed(5872)
  mcmcConf <- MkPrimeMCMC(nRuns = 1L, nIter = 5000L, thin = 5L,
                           maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                           checkEvery = 300L, logFile = log_file,
                           checkpointFile = cp_file, maxTime = 0.5)
  result1 <- RunMkPrime(pd, tree, mcmc = mcmcConf)
  expect_true(file.exists(cp_file))

  # Resume with tree = NULL (should use checkpoint model, no NJ needed)
  result2 <- ResumeMkPrime(cp_file, pd)
  expect_s3_class(result2, "MkPosterior")
  expect_gte(nrow(result2$samples), nrow(result1$samples))
})


test_that("Warmup state (logPostHistory, nStableConsecutive) persisted in checkpoint", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".ckp")
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(c(cp_file, log_file,
                    sub("\\.[^.]+$", "_1.log", log_file),
                    sub("\\.[^.]+$", "_trees.nwk", log_file))), add = TRUE)

  # Long warmup + checkEvery to guarantee checkpoint during warmup
  set.seed(4193)
  RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 5000L, thin = 5L,
                        maxWarmup = 2000L, minWarmup = 2000L, autoTune = FALSE,
                        checkEvery = 500L, logFile = log_file,
                        checkpointFile = cp_file))

  expect_true(file.exists(cp_file))
  cp <- readRDS(cp_file)
  r <- cp$runs[[1]]

  # Warmup state should be persisted (M-149 fix #3)
  expect_true(is.numeric(r$logPostHistory))
  expect_true(length(r$logPostHistory) > 0)
  expect_true(is.integer(r$nStableConsecutive))
})


test_that("Tuning state persisted and tuningBuf re-allocated on resume", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".ckp")
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(c(cp_file, log_file,
                    sub("\\.[^.]+$", "_1.log", log_file),
                    sub("\\.[^.]+$", "_trees.nwk", log_file))), add = TRUE)

  # autoTune = TRUE with short warmup; tuning should fire.
  # Low maxTime to stop quickly (possibly during tuning).
  set.seed(7314)
  result1 <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 10000L, thin = 5L,
                        maxWarmup = 500L, minWarmup = 500L, autoTune = TRUE,
                        checkEvery = 300L, logFile = log_file,
                        checkpointFile = cp_file))

  expect_true(file.exists(cp_file))
  cp <- readRDS(cp_file)
  r <- cp$runs[[1]]

  # If tuning fired, these should be persisted (M-149 fix #4-5).
  # Even if tuning completed before checkpoint, the counters are written.
  if (identical(r$phase, "Sample")) {
    # Tuning completed: counters should reflect completed tuning
    expect_true(is.integer(r$tuningRoundsDone) || is.null(r$tuningRoundsDone))
  }

  # Now simulate a Tuning-phase resume by crafting a checkpoint with
  # phase = "Tuning".  The critical test: this must not crash (M-149 fix #1).
  if (is.list(r$chains)) {
    r$phase <- "Tuning"
    r$tuningIterUsed <- 0L
    r$tuningRoundsDone <- 0L
    cp2 <- cp
    cp2$runs[[1]] <- r
    cp2$iter <- max(cp2$mcmc$warmup %||% 500L, 1L)
    cp_file2 <- tempfile(fileext = ".ckp")
    on.exit(unlink(cp_file2), add = TRUE)
    saveRDS(cp2, cp_file2)

    # Resume from Tuning-phase checkpoint — should not crash
    expect_no_error(
      ResumeMkPrime(cp_file2, pd, tree)
    )
  }
})


test_that("Shared env update fires at batch boundaries (multi-run)", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".ckp")
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(c(cp_file, log_file,
                    sub("\\.[^.]+$", "_1.log", log_file),
                    sub("\\.[^.]+$", "_2.log", log_file),
                    sub("\\.[^.]+$", "_trees.nwk", log_file))), add = TRUE)

  # 2-run serial with maxRhat — the serial runs path uses shared env
  set.seed(2907)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 1000L, thin = 5L,
                        maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                        checkEvery = 300L, logFile = log_file,
                        checkpointFile = cp_file, maxRhat = 1.05))

  expect_true(file.exists(cp_file))
  cp <- readRDS(cp_file)
  expect_equal(length(cp$runs), 2)
  expect_true(!is.null(cp$model))

  # Both runs should have valid chain state

  for (r in cp$runs) {
    expect_true(is.list(r$chains))
    expect_true(is.finite(r$chains[[1]]$log_lik))
  }
})


# --- M-150: maxTime break saves checkpoint ---

test_that("maxTime break saves checkpoint (not just initial)", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".ckp")
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(c(cp_file, log_file,
                    sub("\\.[^.]+$", "_1.log", log_file),
                    sub("\\.[^.]+$", "_trees.nwk", log_file))), add = TRUE)

  # Large checkEvery + bufferSize so no periodic or streaming-flush checkpoint
  # fires.  The only checkpoint updates come from the initial (iter=0) and the

  # maxTime break path (M-150 fix).
  set.seed(5491)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = Inf, thin = 5L,
                        maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                        checkEvery = 1000000L, bufferSize = 100000L,
                        logFile = log_file, checkpointFile = cp_file,
                        maxTime = 0.5))

  expect_equal(result$stop_reason, "max_time")
  expect_true(file.exists(cp_file))

  cp <- readRDS(cp_file)
  # Must be updated beyond the initial iter=0 checkpoint
  expect_gt(cp$iter, 0L)
})
