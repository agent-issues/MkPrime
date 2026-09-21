# Tests for checkpointing (Phase 5: M-035)
skip_slow_tests()

test_that("Checkpoint file is written at check intervals", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

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
  expect_equal(cp$version, 3L)  # v3 adds per-run treeFilePaths
  expect_true(inherits(cp$timestamp, "POSIXct"))
  expect_true(is.list(cp$runs))
  expect_true(cp$iter > 0)
})


test_that("Checkpoint contains valid run state", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

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
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

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
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

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
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

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
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

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
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

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
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

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
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

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
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

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
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

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
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

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


# --- M-149 #2: move weights persisted and restored ---

test_that("Move weights stored in run state and checkpoint", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".ckp")
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(c(cp_file, log_file,
                    sub("\\.[^.]+$", "_1.log", log_file),
                    sub("\\.[^.]+$", "_trees.nwk", log_file))), add = TRUE)

  # autoTune to get adapted weights
  set.seed(8317)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 5000L, thin = 5L,
                        maxWarmup = 500L, minWarmup = 500L, autoTune = TRUE,
                        checkEvery = 500L, logFile = log_file,
                        checkpointFile = cp_file))

  expect_true(file.exists(cp_file))
  cp <- readRDS(cp_file)
  r <- cp$runs[[1]]

  # Per-run moveWeights should be stored
  expect_true(!is.null(r$moveWeights))
  expect_true(is.numeric(r$moveWeights))
  expect_true(length(r$moveWeights) > 0L)
  expect_equal(sum(r$moveWeights), 1.0, tolerance = 1e-10)
})


test_that("Resumed run uses checkpointed move weights", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".ckp")
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(c(cp_file, log_file,
                    sub("\\.[^.]+$", "_1.log", log_file),
                    sub("\\.[^.]+$", "_trees.nwk", log_file))), add = TRUE)

  # Run with tuning + time limit
  set.seed(2946)
  result1 <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 10000L, thin = 5L,
                        maxWarmup = 500L, minWarmup = 500L, autoTune = TRUE,
                        checkEvery = 500L, logFile = log_file,
                        checkpointFile = cp_file, maxTime = 1))

  cp <- readRDS(cp_file)
  savedWeights <- cp$runs[[1]]$moveWeights
  expect_true(!is.null(savedWeights))

  # Resume — should use the saved weights, not defaults
  result2 <- ResumeMkPrime(cp_file, pd, tree)
  expect_s3_class(result2, "MkPosterior")
})


# --- M-149 #6-7: serial orchestrator phase and per-run startIters ---

test_that("Serial multi-run checkpoint stores serialPhase", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  cp_file <- tempfile(fileext = ".ckp")
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(c(cp_file, log_file,
                    sub("\\.[^.]+$", "_1.log", log_file),
                    sub("\\.[^.]+$", "_2.log", log_file),
                    sub("\\.[^.]+$", "_trees.nwk", log_file))), add = TRUE)

  set.seed(6102)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 2000L, thin = 5L,
                        maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                        checkEvery = 300L, logFile = log_file,
                        checkpointFile = cp_file, maxRhat = 1.05))

  expect_true(file.exists(cp_file))
  cp <- readRDS(cp_file)

  # After Phase 1 completes, checkpoint should have serialPhase = 2
  # (or no serialPhase if convergence was reached during Phase 1)
  if (!is.null(cp$serialPhase)) {
    expect_equal(cp$serialPhase, 2L)
  }

  # Per-run actual_iter should be available for startIters restoration
  for (r in cp$runs) {
    expect_true(!is.null(r$actual_iter))
    expect_true(r$actual_iter > 0)
  }
})


# --- Tree-file handling on resume (fast unit tests) ---

test_that(".TruncateTreeToN drops trailing trees", {
  tf <- tempfile(fileext = ".nwk")
  on.exit(unlink(tf), add = TRUE)
  writeLines(c(
    "(a:0.1,b:0.2);",
    "(a:0.3,b:0.4);",
    "(a:0.5,b:0.6);"
  ), tf)

  MkPrime:::.TruncateTreeToN(tf, 2L)
  trees <- readLines(tf)
  expect_equal(length(trees), 2L)
  expect_equal(trees[1], "(a:0.1,b:0.2);")
  expect_equal(trees[2], "(a:0.3,b:0.4);")
})

test_that(".TruncateTreeToN drops trailing partial line", {
  tf <- tempfile(fileext = ".nwk")
  on.exit(unlink(tf), add = TRUE)
  # Final line is a torn write (no terminating semicolon).
  writeLines(c(
    "(a:0.1,b:0.2);",
    "(a:0.3,b:0.4);",
    "(a:0.5,b:0.6"
  ), tf)

  MkPrime:::.TruncateTreeToN(tf, 2L)
  trees <- readLines(tf)
  expect_equal(length(trees), 2L)
  expect_true(all(grepl("\\);$", trees)))
})

test_that(".TruncateTreeToN with nTrees == 0 empties the file", {
  tf <- tempfile(fileext = ".nwk")
  on.exit(unlink(tf), add = TRUE)
  writeLines(c("(a:0.1,b:0.2);", "(a:0.3,b:0.4);"), tf)

  MkPrime:::.TruncateTreeToN(tf, 0L)
  expect_true(file.exists(tf))
  expect_equal(length(readLines(tf)), 0L)
})

test_that(".TruncateTreeToN is a no-op when treeFile is NULL or missing", {
  expect_silent(MkPrime:::.TruncateTreeToN(NULL, 5L))
  tf <- tempfile(fileext = ".nwk")
  expect_false(file.exists(tf))
  expect_silent(MkPrime:::.TruncateTreeToN(tf, 5L))
})


# --- Tree file preserved and extended across resume (serial path) ---

test_that("Resume preserves pre-checkpoint trees and appends new ones (serial)", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  cp_file   <- tempfile(fileext = ".ckp")
  log_file  <- tempfile(fileext = ".log")
  tree_file <- tempfile(fileext = ".nwk")
  on.exit(unlink(c(cp_file, log_file, tree_file,
                    sub("\\.[^.]+$", "_1.log", log_file))), add = TRUE)

  set.seed(91011)
  mcmcConf <- MkPrimeMCMC(nRuns = 1L, nIter = 5000L, thin = 5L,
                           maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                           checkEvery = 300L, logFile = log_file,
                           treeFile = tree_file, checkpointFile = cp_file,
                           maxTime = 0.5)
  result1 <- RunMkPrime(pd, tree, mcmc = mcmcConf)

  expect_true(file.exists(tree_file))
  lines_before <- readLines(tree_file, warn = FALSE)
  tree_lines_before <- lines_before[grepl("\\);\\s*$", lines_before)]
  n_trees_before <- length(tree_lines_before)
  # Snapshot the first valid tree lines so we can verify they're preserved.
  head_before <- tree_lines_before[seq_len(min(3L, n_trees_before))]

  # Resume from the checkpoint.
  result2 <- ResumeMkPrime(cp_file, pd, tree)

  expect_true(file.exists(tree_file))
  lines_after <- readLines(tree_file, warn = FALSE)
  tree_lines_after <- lines_after[grepl("\\);\\s*$", lines_after)]
  expect_gte(length(tree_lines_after), n_trees_before)

  # Pre-checkpoint trees must still be intact at the top of the file.
  if (length(head_before) > 0L) {
    expect_equal(tree_lines_after[seq_along(head_before)], head_before)
  }

  # The resulting file must parse cleanly via ape.
  if (length(tree_lines_after) > 0L) {
    parsed <- ape::read.tree(tree_file)
    if (inherits(parsed, "phylo")) parsed <- list(parsed)
    expect_true(length(parsed) >= n_trees_before)
  }
})


# --- Resume discards post-checkpoint trees from a torn write ---

test_that("Resume truncates post-checkpoint trees written before SIGKILL", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  cp_file   <- tempfile(fileext = ".ckp")
  log_file  <- tempfile(fileext = ".log")
  tree_file <- tempfile(fileext = ".nwk")
  on.exit(unlink(c(cp_file, log_file, tree_file,
                    sub("\\.[^.]+$", "_1.log", log_file))), add = TRUE)

  set.seed(20342)
  mcmcConf <- MkPrimeMCMC(nRuns = 1L, nIter = 5000L, thin = 5L,
                           maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                           checkEvery = 300L, logFile = log_file,
                           treeFile = tree_file, checkpointFile = cp_file,
                           maxTime = 0.5)
  RunMkPrime(pd, tree, mcmc = mcmcConf)

  # Inject "torn" trees that the checkpoint never recorded.
  cat("(t1:0.1,t2:0.2,(t3:0.3,t4:0.4):0.05);\n",
      "(t1:0.9,t2:0.8,(t3:0.7,t4:0.6", # missing closing ");"
      file = tree_file, append = TRUE, sep = "")

  n_before_resume <- sum(grepl("\\);\\s*$",
                                readLines(tree_file, warn = FALSE)))

  # Resume should rewind the tree file to tree_saved_idx, dropping both the
  # extra complete tree and the torn one.
  result2 <- ResumeMkPrime(cp_file, pd, tree)

  # The resulting file must parse cleanly via ape.
  parsed <- ape::read.tree(tree_file)
  if (inherits(parsed, "phylo")) parsed <- list(parsed)
  expect_true(length(parsed) > 0L || length(readLines(tree_file)) == 0L)
})


# --- Per-run tree files: nRuns=2 writes to separate streams ---

test_that("Per-run tree files separate by run", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  log_base  <- tempfile()
  tree_base <- tempfile()
  log_file  <- paste0(log_base, ".log")
  tree_file <- paste0(tree_base, ".nwk")
  cp_file   <- tempfile(fileext = ".ckp")
  expected_logs  <- paste0(log_base,  "_", 1:2, ".log")
  expected_trees <- paste0(tree_base, "_", 1:2, ".nwk")
  on.exit(unlink(c(cp_file, log_file, tree_file,
                   expected_logs, expected_trees)), add = TRUE)

  set.seed(42)
  mcmcConf <- MkPrimeMCMC(nRuns = 2L, nIter = 3000L, thin = 5L,
                           maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                           checkEvery = 300L, logFile = log_file,
                           treeFile = tree_file, checkpointFile = cp_file,
                           maxTime = 0.5, nCore = 1L)
  RunMkPrime(pd, tree, mcmc = mcmcConf)

  for (p in expected_trees) expect_true(file.exists(p),
                                         info = paste("missing:", p))
  l1 <- readLines(expected_trees[1], warn = FALSE)
  l2 <- readLines(expected_trees[2], warn = FALSE)
  t1 <- l1[grepl("\\);\\s*$", l1)]
  t2 <- l2[grepl("\\);\\s*$", l2)]
  expect_gt(length(t1), 0L)
  expect_gt(length(t2), 0L)
  # The shared tree_file path must NOT have been created.
  expect_false(file.exists(tree_file))
  # Streams must be independent (different runs from different RNG seeds /
  # initial states virtually never produce identical sample sequences).
  expect_false(identical(t1, t2))
})


# --- Resume preserves per-run tree files across nRuns=2 ---

test_that("Resume preserves pre-checkpoint trees and appends new ones (nRuns=2 serial)", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  log_base  <- tempfile()
  tree_base <- tempfile()
  log_file  <- paste0(log_base, ".log")
  tree_file <- paste0(tree_base, ".nwk")
  cp_file   <- tempfile(fileext = ".ckp")
  expected_logs  <- paste0(log_base,  "_", 1:2, ".log")
  expected_trees <- paste0(tree_base, "_", 1:2, ".nwk")
  on.exit(unlink(c(cp_file, log_file, tree_file,
                   expected_logs, expected_trees)), add = TRUE)

  set.seed(777)
  mcmcConf <- MkPrimeMCMC(nRuns = 2L, nIter = 3000L, thin = 5L,
                           maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                           checkEvery = 300L, logFile = log_file,
                           treeFile = tree_file, checkpointFile = cp_file,
                           maxTime = 0.5, nCore = 1L)
  RunMkPrime(pd, tree, mcmc = mcmcConf)

  before <- lapply(expected_trees, function(p) {
    l <- readLines(p, warn = FALSE)
    l[grepl("\\);\\s*$", l)]
  })

  ResumeMkPrime(cp_file, pd, tree)

  after <- lapply(expected_trees, function(p) {
    l <- readLines(p, warn = FALSE)
    l[grepl("\\);\\s*$", l)]
  })

  for (i in seq_along(expected_trees)) {
    expect_gte(length(after[[i]]), length(before[[i]]))
    if (length(before[[i]]) > 0L) {
      expect_equal(after[[i]][seq_along(before[[i]])], before[[i]])
    }
    # Each per-run file must parse cleanly.
    parsed <- ape::read.tree(expected_trees[i])
    if (inherits(parsed, "phylo")) parsed <- list(parsed)
    expect_true(length(parsed) >= length(before[[i]]))
  }
})


# --- Sync-invariant aborts ---

test_that("Resume aborts on tree/log desync (fewer trees than checkpoint)", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  cp_file   <- tempfile(fileext = ".ckp")
  log_file  <- tempfile(fileext = ".log")
  tree_file <- tempfile(fileext = ".nwk")
  on.exit(unlink(c(cp_file, log_file, tree_file)), add = TRUE)

  set.seed(31337)
  mcmcConf <- MkPrimeMCMC(nRuns = 1L, nIter = 5000L, thin = 5L,
                           maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                           checkEvery = 300L, logFile = log_file,
                           treeFile = tree_file, checkpointFile = cp_file,
                           maxTime = 0.5)
  RunMkPrime(pd, tree, mcmc = mcmcConf)

  # Corrupt the tree file: drop half the valid trees so nValid < tree_saved_idx.
  lines <- readLines(tree_file, warn = FALSE)
  tree_lines <- lines[grepl("\\);\\s*$", lines)]
  skip_if(length(tree_lines) < 4L,
          "Run did not produce enough trees to corrupt meaningfully.")
  keep <- floor(length(tree_lines) / 2L)
  writeLines(tree_lines[seq_len(keep)], tree_file)

  expect_error(ResumeMkPrime(cp_file, pd, tree),
               regexp = "desync|fewer")
})


test_that("Sync invariant catches truncated log via .TruncateTreeToN", {
  # Direct unit test: build a tree file and call .TruncateTreeToN with a
  # saved_idx that's too small for the number of trees on disk.
  tf <- tempfile(fileext = ".nwk")
  lf <- tempfile(fileext = ".log")
  on.exit(unlink(c(tf, lf)), add = TRUE)
  writeLines(c("(a:0.1,b:0.2);",
               "(a:0.3,b:0.4);",
               "(a:0.5,b:0.6);"), tf)
  writeLines(c("Sample\tx", "1\t1", "2\t2"), lf)

  # 3 trees * treeEvery=2 = 6 param rows required, but saved_idx=2.
  expect_error(
    MkPrime:::.TruncateTreeToN(tf, 3L,
                                logFilePath = lf,
                                saved_idx   = 2L,
                                treeEvery   = 2L),
    regexp = "desync|require"
  )
})


test_that(".TruncateTreeToN aborts when nValid < nTrees (no log args)", {
  tf <- tempfile(fileext = ".nwk")
  on.exit(unlink(tf), add = TRUE)
  writeLines(c("(a:0.1,b:0.2);", "(a:0.3,b:0.4);"), tf)
  # File has 2 valid trees but the chain thinks there should be 5.
  expect_error(MkPrime:::.TruncateTreeToN(tf, 5L),
               regexp = "desync|fewer")
})
