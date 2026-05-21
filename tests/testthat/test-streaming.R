# Tests for batch-streaming MCMC output (logFile / bufferSize)
skip_slow_tests()

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


# --- .FlushBuffer torn-write regression (STREAM-005) ----------------------

test_that(".FlushBuffer writes exactly one line per row with full columns", {
  tf <- tempfile(fileext = ".log")
  on.exit(unlink(tf), add = TRUE)
  # Tab-separated header so .FlushBuffer's appends slot in below it.
  paramNames <- paste0("p", seq_len(43L))   # mimic the 44-col scalar layout
  writeLines(paste(c("Sample", paramNames), collapse = "\t"), tf)

  nRows <- 25L
  nCols <- length(paramNames)
  buffer <- matrix(seq_len(nRows * nCols) + 0.5,
                   nrow = nRows, ncol = nCols,
                   dimnames = list(NULL, paramNames))
  iterNums <- seq_len(nRows) * 100L

  MkPrime:::.FlushBuffer(buffer, nRows, iterNums, tf)

  lines <- readLines(tf)
  # Header + nRows data lines, no torn trailing fragment
  expect_equal(length(lines), nRows + 1L)
  fieldCounts <- vapply(lines, function(l) length(strsplit(l, "\t",
                                                            fixed = TRUE)[[1]]),
                        integer(1L))
  # Every line, header included, must carry Sample + 43 params = 44 fields
  expect_true(all(fieldCounts == nCols + 1L),
              label = "no line has a truncated column count")
})

test_that(".FlushBuffer with nRows = 0 is a no-op", {
  tf <- tempfile(fileext = ".log")
  on.exit(unlink(tf), add = TRUE)
  writeLines("Sample\tA\tB", tf)

  buffer <- matrix(0, nrow = 5L, ncol = 2L, dimnames = list(NULL, c("A", "B")))
  MkPrime:::.FlushBuffer(buffer, 0L, integer(5L), tf)

  expect_equal(readLines(tf), "Sample\tA\tB")
})

test_that("Streaming run emits no torn log rows", {
  # Integration check: every data line in the log must have the same
  # column count as the header.  Catches torn writes from .FlushBuffer
  # that an earlier `cat(paste(lines, collapse='\n'), '\n', ...)` could
  # produce when a single write() split across an NFS boundary.
  f <- .mkStreamFixture()
  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)

  set.seed(9921)
  RunMkPrime(f$pd, f$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 400L, thin = 5L,
                        maxWarmup = 100L, minWarmup = 100L,
                        autoTune = FALSE,
                        logFile = log_file, bufferSize = 7L))

  lines <- readLines(log_file)
  lines <- lines[nzchar(lines)]
  isHeader <- startsWith(lines, "Sample\t")
  isComment <- startsWith(lines, "#")
  expect_true(any(isHeader))
  nFields <- length(strsplit(lines[which(isHeader)[1]], "\t",
                              fixed = TRUE)[[1]])
  dataLines <- lines[!isHeader & !isComment]
  expect_gt(length(dataLines), 0L)
  dataFields <- vapply(dataLines,
                       function(l) length(strsplit(l, "\t", fixed = TRUE)[[1]]),
                       integer(1L))
  expect_true(all(dataFields == nFields),
              label = "every data row matches header column count")
})


# --- z_samples / saved_idx alignment invariant (STREAM-005) ---------------

test_that("Streaming ecology run keeps z_samples aligned with log rows", {
  # Integration check: for an ecology-aware streaming run, the number of
  # z snapshots returned in $z_samples must equal the number of rows the
  # log file actually contains, otherwise RelabelEcology() can't run.
  set.seed(42L)
  tips <- paste0("t", 1:6)
  mat <- matrix(c(
    0, 1, 0, 1, 0, 1,
    0, 0, 1, 1, 0, 1,
    0, 1, 1, 0, 1, 0,
    0, 0, 1, 1, 2, 2
  ), nrow = 6L, ncol = 4L, byrow = FALSE,
  dimnames = list(tips, NULL))
  pd  <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, ecology = 4L)
  tree <- TreeTools::Preorder(ape::rtree(6L, tip.label = tips))

  log_file <- tempfile(fileext = ".log")
  on.exit(unlink(log_file), add = TRUE)

  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "none")
  mcmc  <- MkPrimeMCMC(nIter = 200L, nChains = 1L, nRuns = 1L,
                        thin = 5L, treeThin = 5L,
                        maxWarmup = 50L, minWarmup = 50L,
                        autoTune = FALSE,
                        logFile = log_file, bufferSize = 4L,
                        checkpointFile = NULL)
  res <- RunMkPrime(mkd, tree = tree, model = model, mcmc = mcmc)

  samp <- ReadMkLog(log_file)
  expect_gt(nrow(samp), 0L)
  expect_equal(length(res$z_samples), nrow(samp))
})

test_that("RelabelEcology aborts with diagnostics when z_samples misaligned", {
  # Default behaviour (no trimZSamples): abort with an informative message
  # that names the surplus count and suggests trimZSamples = "tail".
  samples <- matrix(c(2.5, 0.3, 0.4,
                       0.8, 0.2, 0.6),
                     nrow = 2L, byrow = TRUE,
                     dimnames = list(NULL, c("phi", "pi0", "theta_1")))
  zList <- list(
    matrix(c(0L, 1L, 2L), nrow = 3L, ncol = 1L),
    matrix(c(2L, 0L, 1L), nrow = 3L, ncol = 1L),
    matrix(c(1L, 1L, 0L), nrow = 3L, ncol = 1L)  # surplus
  )
  fake <- list(
    samples = samples,
    z_samples = zList,
    model = list(ecologyAware = TRUE, magnitudeMode = "global"),
    data  = list(refEcology = 0L, kEcology = 2L)
  )
  class(fake) <- "MkPosterior"

  # Must abort, not warn-and-trim
  expect_error(RelabelEcology(fake, magnitudeMode = "global"),
               "z_samples")
})

test_that("RelabelEcology trims tail surplus when trimZSamples = 'tail'", {
  # When the caller explicitly confirms the surplus entries are at the tail
  # (typical: streaming run interrupted mid-flush), RelabelEcology should
  # warn and keep only the first nSamples z entries.
  samples <- matrix(c(2.5, 0.3, 0.4,
                       0.8, 0.2, 0.6),
                     nrow = 2L, byrow = TRUE,
                     dimnames = list(NULL, c("phi", "pi0", "theta_1")))
  zList <- list(
    matrix(c(0L, 1L, 2L), nrow = 3L, ncol = 1L),  # paired with row 1
    matrix(c(2L, 0L, 1L), nrow = 3L, ncol = 1L),  # paired with row 2
    matrix(c(1L, 1L, 0L), nrow = 3L, ncol = 1L)   # surplus at tail
  )
  fake <- list(
    samples = samples,
    z_samples = zList,
    model = list(ecologyAware = TRUE, magnitudeMode = "global"),
    data  = list(refEcology = 0L, kEcology = 2L)
  )
  class(fake) <- "MkPosterior"

  expect_warning(
    out <- RelabelEcology(fake, magnitudeMode = "global",
                          trimZSamples = "tail"),
    "z_samples"
  )
  expect_equal(length(out$z_samples), 2L)
  expect_true(isTRUE(attr(out, "relabelled")))
})


# --- STREAM-006: stored tree edge.length must equal log tree_length ---
# Bug: brColStart formula in .RunMkPrime / .ResumeMkPrime omitted the ecology
# block (phi, pi0, theta_*) when ecologyAware=TRUE. The first nEco edges of
# every stored phylo had bogus values (kPrime integers scaled by tree_length);
# only the trailing edges were correct. Fix: derive brColStart from paramNames.
test_that("STREAM-006 stored tree edge.length sums to log tree_length (ecology-aware)", {
  set.seed(42)
  # Minimal ecology-aware fixture: 6 tips, 5 chars, 3 ecology states.
  mat <- matrix(c(
    0, 1, 0, 1, 0,
    1, 0, 1, 0, 1,
    0, 0, 1, 1, 0,
    1, 1, 0, 0, 1,
    0, 1, 1, 0, 0,
    1, 0, 0, 1, 1
  ), nrow = 6, byrow = TRUE,
  dimnames = list(paste0("t", 1:6), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  eco <- setNames(c(0L, 0L, 1L, 1L, 2L, 2L), paste0("t", 1:6))
  mkd <- MkPrimeData(pd, ecology = eco, neomorphic = integer(0))
  set.seed(7)
  tree <- TreeTools::Preorder(ape::rtree(6, tip.label = paste0("t", 1:6)))
  tree$edge.length <- rep(0.1, nrow(tree$edge))
  model <- MkPrimeModel(ecologyAware = TRUE, magnitudeMode = "global",
                        kPrimePrior = "geometric",
                        rho0Alpha = 2, rho0Beta = 2,
                        thetaAlpha = 2, thetaBeta = 2,
                        sigmaPhi = 0.5)
  logFile <- tempfile(fileext = ".log")
  ckpFile <- tempfile(fileext = ".ckp")
  mcmc <- MkPrimeMCMC(nIter = 400L, nChains = 1L, nRuns = 1L,
                      thin = 20L, treeThin = 20L,
                      minWarmup = 100L, maxWarmup = 200L,
                      logFile = logFile, checkpointFile = ckpFile)
  res <- RunMkPrime(mkd, tree = tree, model = model, mcmc = mcmc)
  samp <- ReadMkLog(logFile)
  expect_gte(length(res$trees), 2L)
  expect_equal(nrow(samp), length(res$trees))
  # Each stored phylo's sum(edge.length) must equal log tree_length to high
  # precision (Dirichlet branch proportions sum to 1).
  for (i in seq_along(res$trees)) {
    el_sum <- sum(res$trees[[i]]$edge.length)
    tl_log <- samp[i, "tree_length"]
    expect_equal(el_sum, tl_log, tolerance = 1e-6,
                 info = paste("sample i =", i))
  }
})
