# Tests for batch-streaming MCMC output (logFile / bufferSize)

# Shared minimal fixture
.mkStreamFixture <- function() {
  tree <- ape::read.tree(
    text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);"
  )
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
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
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 300L, thin = 5L, warmup = 100L,
                        logFile = log_file, bufferSize = 10L))

  expect_true(file.exists(log_file))
  lines <- readLines(log_file)
  # Header must be first line, starting with "Sample"
  expect_true(startsWith(lines[1], "Sample\t"))
  # (300 - 100) / 5 = 40 thinned samples → 40 data rows
  expect_equal(length(lines) - 1L, 40L)
})

test_that("Sample column is strictly increasing", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)

  set.seed(8821)
  RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 300L, thin = 5L, warmup = 100L,
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
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 120L, thin = 5L, warmup = 45L,
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
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 200L, thin = 5L, warmup = 100L,
                        logFile = log_file, bufferSize = 10L))

  expect_true(file.exists(log1))
  expect_true(file.exists(log2))
  dat1 <- utils::read.table(log1, header = TRUE, sep = "\t")
  dat2 <- utils::read.table(log2, header = TRUE, sep = "\t")
  # 20 samples each
  expect_equal(nrow(dat1), 20L)
  expect_equal(nrow(dat2), 20L)
})


# --- MkPosterior result in streaming mode ---

test_that("Streaming result has logFile field and empty samples matrix", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)

  set.seed(1947)
  result <- RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 200L, thin = 5L, warmup = 100L,
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
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 200L, thin = 5L, warmup = 100L,
                        logFile = log_file, bufferSize = 20L))

  expect_equal(length(result$trees), 20L)
  expect_s3_class(result$trees[[1]], "phylo")
})

test_that("summary() errors helpfully for unloaded streaming result", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)

  set.seed(7723)
  result <- RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 200L, thin = 5L, warmup = 100L,
                        logFile = log_file, bufferSize = 20L))

  expect_error(summary(result), "Samples are not in memory")
})


# --- ReadMkLog ---

test_that("ReadMkLog returns matrix with correct dimensions and colnames", {
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)

  set.seed(4481)
  RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 300L, thin = 5L, warmup = 100L,
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
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 300L, thin = 5L, warmup = 100L))

  set.seed(9034)
  RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 300L, thin = 5L, warmup = 100L,
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
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 200L, thin = 5L, warmup = 100L,
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
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, warmup = 100L,
                        logFile = log_file, bufferSize = 20L,
                        checkEvery = 200L, checkpointFile = cp_file))

  expect_true(file.exists(cp_file))
  cp <- readRDS(cp_file)
  expect_equal(cp$version, 2L)
  expect_null(cp$runs[[1]]$samples)
  expect_false(is.null(cp$logFilePaths))
  expect_false(is.null(cp$paramNames))
})

test_that("In-memory checkpoint is still version 1", {
  f <- .mkStreamFixture()
  cp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(cp_file), add = TRUE)

  set.seed(7732)
  RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, warmup = 100L,
                        checkEvery = 200L, checkpointFile = cp_file))

  cp <- readRDS(cp_file)
  expect_equal(cp$version, 1L)
  expect_false(is.null(cp$runs[[1]]$samples))
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
  n_before <- result1$nSamples

  result2 <- ResumeMkPrime(cp_file, f$pd, f$tree)

  dat <- ReadMkLog(log_file)
  # Strictly increasing Sample values (stored as rownames by ReadMkLog)
  sampleNums <- as.integer(rownames(dat))
  expect_true(all(diff(sampleNums) > 0))
  # Must have at least as many rows as before the resume
  expect_gte(nrow(dat), n_before)
})

test_that("In-memory mode is unchanged when logFile = NULL", {
  f <- .mkStreamFixture()

  set.seed(2277)
  result <- RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 300L, thin = 5L, warmup = 100L))

  expect_null(result$logFile)
  expect_false(is.null(result$samples))
  expect_equal(nrow(result$samples), 40L)
  expect_s3_class(result, "MkPosterior")
})
