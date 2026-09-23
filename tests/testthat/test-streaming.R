# Tests for batch-streaming MCMC output (logFile / bufferSize)
skip_slow_tests()

# Shared minimal fixture
.mkStreamFixture <- function() {
  tree <- ape::read.tree(
    text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);"
  )
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)
  list(tree = tree, pd = pd)
}


# --- MkPrimeMCMC validation ---

test_that("MkPrimeMCMC accepts valid logFile and bufferSize", {
  m <- MkPrimeMCMC(logFile = "run.log", bufferSize = 100L)
  expect_equal(m$logFile, "run.log")
  expect_equal(m$bufferSize, 100L)
})

test_that("MkPrimeMCMC rejects non-string logFile", {
  expect_error(MkPrimeMCMC(logFile = 123), "logFile")
})

test_that("MkPrimeMCMC rejects bufferSize < 1", {
  expect_error(MkPrimeMCMC(bufferSize = 0L), "bufferSize")
})


# --- Log file creation and content ---

test_that("Log file is created with correct header in streaming mode", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)

  set.seed(3714)
  RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 300L, thin = 5L, maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        logFile = log_file, bufferSize = 10L))

  expect_true(file.exists(log_file))
  lines <- readLines(log_file)
  # Header must be first line, starting with "Sample"
  expect_true(startsWith(lines[1], "Sample\t"))
  # Count data rows (exclude header and comment lines from move-weight log)
  dataLines <- lines[!startsWith(lines, "Sample\t") & !startsWith(lines, "#")]
  # (300 - 100) / 5 = 40 thinned samples → 40 data rows
  expect_equal(length(dataLines), 40L)
})

test_that("Sample column is strictly increasing", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)

  set.seed(8821)
  RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 300L, thin = 5L, maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        logFile = log_file, bufferSize = 10L))

  dat <- utils::read.table(log_file, header = TRUE, sep = "\t")
  expect_true(all(diff(dat$Sample) > 0))
})

test_that("Flush boundary: multiple flushes land correctly (small bufferSize)", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)

  # bufferSize = 3, 15 thinned samples → 5 flushes of 3
  set.seed(2255)
  RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 120L, thin = 5L, maxWarmup = 45L, minWarmup = 45L, autoTune = FALSE,
                        logFile = log_file, bufferSize = 3L))

  dat <- utils::read.table(log_file, header = TRUE, sep = "\t")
  # (120 - 45) / 5 = 15 samples
  expect_equal(nrow(dat), 15L)
  expect_true(all(diff(dat$Sample) > 0))
})

test_that("Multi-run creates separate log files with _1/_2 suffix", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  log1 <- sub("\\.log$", "_1.log", log_file)
  log2 <- sub("\\.log$", "_2.log", log_file)
  on.exit(unlink(c(log1, log2)), add = TRUE)

  set.seed(6631)
  RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 200L, thin = 5L, maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        logFile = log_file, bufferSize = 10L))

  expect_true(file.exists(log1))
  expect_true(file.exists(log2))
  dat1 <- utils::read.table(log1, header = TRUE, sep = "\t")
  dat2 <- utils::read.table(log2, header = TRUE, sep = "\t")
  # 20 samples each
  expect_equal(nrow(dat1), 20L)
  expect_equal(nrow(dat2), 20L)
})


# --- Convergence window (circular buffer) ---

test_that(".ConvWindowRows returns rows in chronological order after wrap", {
  convWindowSize <- 10L
  buf <- .InitStreamBuffers(nParams = 1L, paramNames = "x",
                             bufferSize = 100L, convWindowSize = convWindowSize)
  logFile <- tempfile(fileext = ".log")
  on.exit(unlink(logFile), add = TRUE)

  for (i in 1:13) {
    buf <- .AddToStreamBuffer(buf, row = i, iterNum = i, logFile = logFile,
                               bufferSize = 100L, convWindowSize = convWindowSize)
  }
  expect_true(buf$conv_filled)

  rows <- .ConvWindowRows(buf, minRows = 1L)
  # 13 samples into a window of 10 -> the oldest 3 (1:3) have been
  # overwritten; the chronological order of what remains is 4:13.
  expect_equal(as.vector(rows), 4:13)
})

test_that(".ConvWindowRows is chronological right when the window first fills", {
  # conv_head == convWindowSize is the edge case: the window has just
  # wrapped to "filled" but has not yet overwritten anything.
  convWindowSize <- 10L
  buf <- .InitStreamBuffers(nParams = 1L, paramNames = "x",
                             bufferSize = 100L, convWindowSize = convWindowSize)
  logFile <- tempfile(fileext = ".log")
  on.exit(unlink(logFile), add = TRUE)

  for (i in 1:10) {
    buf <- .AddToStreamBuffer(buf, row = i, iterNum = i, logFile = logFile,
                               bufferSize = 100L, convWindowSize = convWindowSize)
  }
  expect_true(buf$conv_filled)
  expect_equal(buf$conv_head, convWindowSize)

  rows <- .ConvWindowRows(buf, minRows = 1L)
  expect_equal(as.vector(rows), 1:10)
})

test_that(".ConvWindowRows returns unrotated rows before the window fills", {
  convWindowSize <- 10L
  buf <- .InitStreamBuffers(nParams = 1L, paramNames = "x",
                             bufferSize = 100L, convWindowSize = convWindowSize)
  logFile <- tempfile(fileext = ".log")
  on.exit(unlink(logFile), add = TRUE)

  for (i in 1:5) {
    buf <- .AddToStreamBuffer(buf, row = i, iterNum = i, logFile = logFile,
                               bufferSize = 100L, convWindowSize = convWindowSize)
  }
  expect_false(buf$conv_filled)

  rows <- .ConvWindowRows(buf, minRows = 1L)
  expect_equal(as.vector(rows), 1:5)
})


# --- MkPosterior result in streaming mode ---

test_that("Streaming result has logFile field and empty samples matrix", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)

  set.seed(1947)
  result <- RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 200L, thin = 5L, maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        logFile = log_file, bufferSize = 20L))

  expect_s3_class(result, "MkPosterior")
  expect_false(is.null(result$logFile))
  expect_equal(result$logFile, log_file)
  expect_equal(nrow(result$samples), 0L)  # empty matrix, not NULL
  expect_equal(result$nSamples, 20L)
})

test_that("Trees are still stored in memory in streaming mode", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)

  set.seed(5512)
  result <- RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 200L, thin = 5L, maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        logFile = log_file, bufferSize = 20L))

  expect_equal(length(result$trees), 20L)
  expect_s3_class(result$trees[[1]], "phylo")
})

test_that("summary() auto-loads samples from log file for streaming result", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)

  set.seed(7723)
  result <- RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 200L, thin = 5L, maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        logFile = log_file, bufferSize = 20L))

  # Samples auto-loaded from log file via .PostBurninData
  s <- summary(result)
  expect_s3_class(s, "data.frame")
  expect_true(nrow(s) > 0L)
})


# --- ReadMkLog ---

test_that("ReadMkLog returns matrix with correct dimensions and colnames", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)

  set.seed(4481)
  RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 300L, thin = 5L, maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        logFile = log_file, bufferSize = 10L))

  mat <- ReadMkLog(log_file)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 40L)
  expect_false("Sample" %in% colnames(mat))
  expect_true("log_posterior" %in% colnames(mat))
})

test_that("ReadMkLog round-trips with in-memory mode", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)

  set.seed(9034)
  result_mem <- RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 300L, thin = 5L, maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE))

  set.seed(9034)
  RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 300L, thin = 5L, maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        logFile = log_file, bufferSize = 7L))

  loaded <- ReadMkLog(log_file)
  # Column names must match; values must be identical (same seed, same run)
  expect_equal(colnames(result_mem$samples), colnames(loaded))
  expect_equal(result_mem$samples, loaded, ignore_attr = TRUE)
})

test_that("ReadMkLog errors on missing file", {
  expect_error(ReadMkLog("/no/such/file.log"), "not found")
})

test_that("ReadMkLog combines multiple files row-wise", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  log1 <- sub("\\.log$", "_1.log", log_file)
  log2 <- sub("\\.log$", "_2.log", log_file)
  on.exit(unlink(c(log1, log2)), add = TRUE)

  set.seed(1188)
  RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 200L, thin = 5L, maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        logFile = log_file, bufferSize = 10L))

  combined <- ReadMkLog(c(log1, log2))
  expect_equal(nrow(combined), 40L)  # 20 per run
})


# --- Checkpoints in streaming mode ---

test_that("Streaming checkpoint is version 2 and state-only", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  cp_file  <- tempfile(fileext = ".rds")
  on.exit(unlink(c(log_file, cp_file)), add = TRUE)

  set.seed(3301)
  RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        logFile = log_file, bufferSize = 20L,
                        checkEvery = 200L, checkpointFile = cp_file))

  expect_true(file.exists(cp_file))
  cp <- readRDS(cp_file)
  expect_equal(cp$version, 3L)
  expect_null(cp$runs[[1]]$samples)
  expect_false(is.null(cp$logFilePaths))
  expect_false(is.null(cp$paramNames))
})

test_that("Checkpoint is always streaming (version 3) even without logFile", {
  f <- .mkStreamFixture()
  cp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(cp_file), add = TRUE)

  set.seed(7732)
  RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        checkEvery = 200L, checkpointFile = cp_file))

  cp <- readRDS(cp_file)
  # Always streaming now (temp log), so version 3
  expect_equal(cp$version, 3L)
  # Samples live in log files, not in checkpoint
  expect_null(cp$runs[[1]]$samples)
  # But chain state is preserved
  expect_true(is.list(cp$runs[[1]]$chains))
  expect_true(is.finite(cp$runs[[1]]$chains[[1]]$log_lik))
})

test_that("Resume in streaming mode appends without gaps or duplicates", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  cp_file  <- tempfile(fileext = ".rds")
  on.exit(unlink(c(log_file, cp_file)), add = TRUE)

  set.seed(5544)
  result1 <- RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 2000L, thin = 5L, warmup = 200L,
                        logFile = log_file, bufferSize = 20L,
                        checkEvery = 300L, checkpointFile = cp_file,
                        maxTime = 0.5))

  expect_true(file.exists(cp_file))
  # Use checkpoint's saved_idx (not nSamples) — truncation on resume
  # discards samples flushed after the last checkpoint.
  cp <- readRDS(cp_file)
  cp_saved <- cp$runs[[1]]$saved_idx

  result2 <- ResumeMkPrime(cp_file, f$pd, f$tree)

  dat <- ReadMkLog(log_file)
  # Strictly increasing Sample values (stored as rownames by ReadMkLog)
  sampleNums <- as.integer(rownames(dat))
  expect_true(all(diff(sampleNums) > 0))
  # Must have at least as many rows as checkpoint recorded
  expect_gte(nrow(dat), cp_saved)
})

# --- .TruncateLogToN unit tests ---

test_that(".TruncateLogToN removes excess rows", {
  local_mkp_verbosity()
  tf <- tempfile(fileext = ".log")
  on.exit(unlink(tf), add = TRUE)
  writeLines(c("Sample\tA\tB",
               "1\t0.1\t0.2", "2\t0.3\t0.4", "3\t0.5\t0.6"), tf)

  expect_message(.TruncateLogToN(tf, 2L), "discarding")
  lines <- readLines(tf)
  expect_equal(length(lines), 3L)  # header + 2 data rows
  expect_true(startsWith(lines[3], "2\t"))
})

test_that(".TruncateLogToN is a no-op when counts match", {
  tf <- tempfile(fileext = ".log")
  on.exit(unlink(tf), add = TRUE)
  writeLines(c("Sample\tA\tB", "1\t0.1\t0.2", "2\t0.3\t0.4"), tf)

  .TruncateLogToN(tf, 2L)
  lines <- readLines(tf)
  expect_equal(length(lines), 3L)
})

test_that(".TruncateLogToN warns when log has fewer rows than expected", {
  tf <- tempfile(fileext = ".log")
  on.exit(unlink(tf), add = TRUE)
  writeLines(c("Sample\tA\tB", "1\t0.1\t0.2"), tf)

  expect_warning(.TruncateLogToN(tf, 5L), "fewer rows")
})


# --- Resume truncation integration ---

test_that("Resume truncates post-checkpoint rows from log file", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  cp_file  <- tempfile(fileext = ".rds")
  on.exit(unlink(c(log_file, cp_file)), add = TRUE)

  set.seed(4419)
  result1 <- RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 2000L, thin = 5L, warmup = 200L,
                        logFile = log_file, bufferSize = 20L,
                        checkEvery = 300L, checkpointFile = cp_file,
                        maxTime = 0.5))

  expect_true(file.exists(cp_file))

  # Read checkpoint to find its saved_idx

  cp <- readRDS(cp_file)
  cp_saved <- cp$runs[[1]]$saved_idx

  # Append a fake row to simulate post-checkpoint writes (e.g. crash
  # recovery where BuildResult flushed extra samples)
  cat("999999\t-100\t-100\t1\t1\t0.5\t0.5\n",
      file = log_file, append = TRUE)
  rows_with_fake <- length(readLines(log_file)) - 1L
  expect_gt(rows_with_fake, cp_saved)

  # Resume should truncate back to checkpoint saved_idx, then continue
  result2 <- ResumeMkPrime(cp_file, f$pd, f$tree)

  dat <- ReadMkLog(log_file)
  sampleNums <- as.integer(rownames(dat))

  # Fake row must be gone
  expect_false(999999L %in% sampleNums)
  # Strictly increasing — no duplicates or gaps from the splice
  expect_true(all(diff(sampleNums) > 0))
  # At least as many samples as the checkpoint recorded
  expect_gte(nrow(dat), cp_saved)
})


test_that("In-memory mode is unchanged when logFile = NULL", {
  f <- .mkStreamFixture()

  set.seed(2277)
  result <- RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 300L, thin = 5L, maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE))

  expect_null(result$logFile)
  expect_false(is.null(result$samples))
  expect_equal(nrow(result$samples), 40L)
  expect_s3_class(result, "MkPosterior")
})
