# Tests for interrupt recovery (M-132)
#
# The tryCatch(interrupt = ...) mechanism can't be tested in-process
# (testthat runs synchronously), so we test the components:
# 1. MkPrimeRecover() with planted metadata (simulated interrupt)
# 2. .CleanupTempLogs() deletes files and resets state
# 3. .CleanupStaleTempLogs() handles stale temp files from prior runs
# 4. No-logFile runs produce in-memory samples (temp log transparent to user)
# 5. Explicit-logFile runs behave as before

library("ape")
library("TreeTools")

# Minimal test data
.mkTestData <- function() {
  pd <- MatrixToPhyDat(matrix(
    c(0L, 1L, 0L, 1L,
      0L, 0L, 1L, 1L),
    nrow = 4, dimnames = list(paste0("t", 1:4), paste0("c", 1:2))
  ))
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  list(pd = pd, tree = tree)
}


# --- MkPrimeRecover tests ---

test_that("MkPrimeRecover returns NULL when no interrupted run exists", {
  mkp_env <- environment(RunMkPrime)$.mkp_env
  mkp_env$recovery <- NULL
  expect_message(rec <- MkPrimeRecover(), "No interrupted run")
  expect_null(rec)
})

test_that("MkPrimeRecover retrieves partial results from planted metadata", {
  mkp_env <- environment(RunMkPrime)$.mkp_env

  # Write a fake temp log with header + 5 data rows
  tmpLog <- tempfile("mkp_test_", fileext = ".log")
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                   "rate_log_sd", "p")
  header <- paste(c("Sample", paramNames), collapse = "\t")
  set.seed(8421)
  rows <- vapply(seq_len(5), function(i) {
    paste(c(i * 100, rnorm(length(paramNames))), collapse = "\t")
  }, character(1))
  writeLines(c(header, rows), tmpLog)

  d <- .mkTestData()
  mcmc <- MkPrimeMCMC(nIter = 1000L, maxWarmup = 200L, minWarmup = 200L,
                       autoTune = FALSE)
  model <- MkPrimeModel()
  mkd <- MkPrimeData(d$pd)

  # Plant recovery metadata
  mkp_env$recovery <- list(
    logFiles   = tmpLog,
    paramNames = paramNames,
    model      = model,
    data       = mkd,
    mcmc       = mcmc,
    isTempLog  = TRUE,
    time       = Sys.time()
  )

  expect_message(rec <- MkPrimeRecover(), "Recovered 5 sample")
  expect_s3_class(rec, "MkPosterior")
  expect_equal(nrow(rec$samples), 5L)
  expect_true(isTRUE(rec$partial))
  expect_equal(rec$stop_reason, "interrupted")

  # Temp log should be cleaned up

  expect_false(file.exists(tmpLog))

  # Second call should return NULL
  expect_message(rec2 <- MkPrimeRecover(), "No interrupted run")
  expect_null(rec2)
})

test_that("MkPrimeRecover handles missing log files gracefully", {
  mkp_env <- environment(RunMkPrime)$.mkp_env
  mkp_env$recovery <- list(
    logFiles   = tempfile("gone_", fileext = ".log"),
    paramNames = "log_posterior",
    model      = MkPrimeModel(),
    data       = NULL,
    mcmc       = MkPrimeMCMC(nIter = 100L, maxWarmup = 50L, minWarmup = 50L,
                              autoTune = FALSE),
    isTempLog  = TRUE,
    time       = Sys.time()
  )
  expect_message(rec <- MkPrimeRecover(), "no longer exist")
  expect_null(rec)
  expect_null(mkp_env$recovery)
})

test_that("MkPrimeRecover handles empty log file", {
  mkp_env <- environment(RunMkPrime)$.mkp_env
  tmpLog <- tempfile("mkp_empty_", fileext = ".log")
  writeLines("Sample\tlog_posterior", tmpLog)  # header only

  mkp_env$recovery <- list(
    logFiles   = tmpLog,
    paramNames = "log_posterior",
    model      = MkPrimeModel(),
    data       = MkPrimeData(.mkTestData()$pd),
    mcmc       = MkPrimeMCMC(nIter = 100L, maxWarmup = 50L, minWarmup = 50L,
                              autoTune = FALSE),
    isTempLog  = TRUE,
    time       = Sys.time()
  )
  expect_message(rec <- MkPrimeRecover(), "no samples")
  expect_null(rec)
})


# --- Temp log cleanup tests ---

test_that(".CleanupTempLogs deletes files and resets state", {
  mkp_env <- environment(RunMkPrime)$.mkp_env
  f1 <- tempfile(fileext = ".log")
  f2 <- tempfile(fileext = ".log")
  writeLines("test", f1); writeLines("test", f2)
  mkp_env$active_temp_logs <- c(f1, f2)
  mkp_env$recovery <- list(logFiles = c(f1, f2))

  MkPrime:::.CleanupTempLogs(c(f1, f2))

  expect_false(file.exists(f1))
  expect_false(file.exists(f2))
  expect_null(mkp_env$active_temp_logs)
  expect_null(mkp_env$recovery)
})

test_that(".CleanupStaleTempLogs removes leftover files from prior run", {
  mkp_env <- environment(RunMkPrime)$.mkp_env
  stale <- tempfile("mkp_stale_", fileext = ".log")
  writeLines("stale data", stale)
  mkp_env$active_temp_logs <- stale
  mkp_env$recovery <- list(logFiles = stale)

  MkPrime:::.CleanupStaleTempLogs()

  expect_false(file.exists(stale))
  expect_null(mkp_env$active_temp_logs)
  expect_null(mkp_env$recovery)
})


# --- Integration: RunMkPrime with/without logFile ---

test_that("RunMkPrime without logFile produces in-memory samples", {
  d <- .mkTestData()
  mcmc <- MkPrimeMCMC(nIter = 500L, thin = 5L,
                       maxWarmup = 100L, minWarmup = 100L,
                       autoTune = FALSE, nRuns = 1L, nChains = 1L)
  set.seed(7723)
  res <- suppressWarnings(RunMkPrime(d$pd, d$tree, mcmc = mcmc))
  expect_s3_class(res, "MkPosterior")
  expect_true(nrow(res$samples) > 0)
  expect_null(res$logFile)
})

test_that("RunMkPrime with explicit logFile streams to disk", {
  d <- .mkTestData()
  logPath <- tempfile("mkp_explicit_", fileext = ".log")
  on.exit(unlink(logPath), add = TRUE)
  mcmc <- MkPrimeMCMC(nIter = 500L, thin = 5L,
                       maxWarmup = 100L, minWarmup = 100L,
                       autoTune = FALSE, nRuns = 1L, nChains = 1L,
                       logFile = logPath)
  set.seed(7723)
  res <- suppressWarnings(RunMkPrime(d$pd, d$tree, mcmc = mcmc))
  expect_s3_class(res, "MkPosterior")
  expect_true(file.exists(logPath))
  expect_equal(res$logFile, logPath)
  # Streaming mode: samples not in memory
  expect_equal(nrow(res$samples), 0L)
})

test_that("starting a new RunMkPrime cleans up stale temp logs", {
  mkp_env <- environment(RunMkPrime)$.mkp_env
  stale <- tempfile("mkp_stale_", fileext = ".log")
  writeLines("old data", stale)
  mkp_env$active_temp_logs <- stale
  mkp_env$recovery <- list(logFiles = stale)

  d <- .mkTestData()
  mcmc <- MkPrimeMCMC(nIter = 500L, thin = 5L,
                       maxWarmup = 100L, minWarmup = 100L,
                       autoTune = FALSE, nRuns = 1L, nChains = 1L)
  set.seed(3491)
  res <- suppressWarnings(RunMkPrime(d$pd, d$tree, mcmc = mcmc))

  # Stale file should be gone
  expect_false(file.exists(stale))
  expect_s3_class(res, "MkPosterior")
})


# --- MkPrimeRecover(logFile = ...) tests ---

# Helper: write a fake log file with header + n data rows
.writeFakeLog <- function(path, paramNames, n = 5, seed = 7392) {
  header <- paste(c("Sample", paramNames), collapse = "\t")
  set.seed(seed)
  rows <- vapply(seq_len(n), function(i) {
    paste(c(i * 100, rnorm(length(paramNames))), collapse = "\t")
  }, character(1))
  writeLines(c(header, rows), path)
}

test_that("MkPrimeRecover(logFile) reads single log file", {
  tmpLog <- tempfile("mkp_recov_", fileext = ".log")
  on.exit(unlink(tmpLog), add = TRUE)
  params <- c("log_posterior", "log_likelihood", "tree_length", "rate_log_sd")
  .writeFakeLog(tmpLog, params, n = 10)

  expect_message(rec <- MkPrimeRecover(logFile = tmpLog), "Recovered 10")
  expect_s3_class(rec, "MkPosterior")
  expect_equal(nrow(rec$samples), 10L)
  expect_true(isTRUE(rec$partial))
  expect_equal(rec$stop_reason, "recovered")
  # No checkpoint → model and mcmc are NULL
  expect_null(rec$model)
  expect_null(rec$mcmc)
})

test_that("MkPrimeRecover(logFile) discovers _1/_2 multi-run files", {
  tmpDir <- tempfile("mkp_multi_")
  dir.create(tmpDir)
  on.exit(unlink(tmpDir, recursive = TRUE), add = TRUE)
  basePath <- file.path(tmpDir, "test.log")
  params <- c("log_posterior", "log_likelihood", "tree_length")

  .writeFakeLog(file.path(tmpDir, "test_1.log"), params, n = 8, seed = 1122)
  .writeFakeLog(file.path(tmpDir, "test_2.log"), params, n = 6, seed = 3344)

  expect_message(rec <- MkPrimeRecover(logFile = basePath), "Recovered 14")
  expect_s3_class(rec, "MkPosterior")
  expect_equal(nrow(rec$samples), 14L)
  expect_equal(rec$nRuns, 2L)
  expect_length(rec$per_run, 2L)
  expect_equal(nrow(rec$per_run[[1]]$samples), 8L)
  expect_equal(nrow(rec$per_run[[2]]$samples), 6L)
})

test_that("MkPrimeRecover(logFile) extracts metadata from checkpoint", {
  tmpDir <- tempfile("mkp_ckp_")
  dir.create(tmpDir)
  on.exit(unlink(tmpDir, recursive = TRUE), add = TRUE)
  logPath <- file.path(tmpDir, "run.log")
  ckpPath <- file.path(tmpDir, "run.ckp")
  params <- c("log_posterior", "log_likelihood", "tree_length")

  .writeFakeLog(logPath, params, n = 5)

  # Write a minimal checkpoint
  mcmc <- MkPrimeMCMC(nIter = 5000L, thin = 10L,
                       maxWarmup = 100L, minWarmup = 100L,
                       autoTune = FALSE)
  model <- MkPrimeModel()
  saveRDS(list(mcmc = mcmc, model = model, iter = 500L,
               timestamp = Sys.time(), version = 2L),
          ckpPath)

  rec <- suppressMessages(MkPrimeRecover(logFile = logPath))
  expect_s3_class(rec, "MkPosterior")
  expect_false(is.null(rec$mcmc))
  expect_false(is.null(rec$model))
  expect_equal(rec$mcmc$nIter, 5000L)
})

test_that("MkPrimeRecover(logFile) returns NULL for missing file", {
  expect_message(
    rec <- MkPrimeRecover(logFile = tempfile("nonexistent")),
    "No log files found"
  )
  expect_null(rec)
})

test_that("MkPrimeRecover(logFile) returns NULL for empty log file", {
  tmpLog <- tempfile("mkp_empty_", fileext = ".log")
  on.exit(unlink(tmpLog), add = TRUE)
  writeLines("Sample\tlog_posterior", tmpLog)

  expect_message(
    rec <- MkPrimeRecover(logFile = tmpLog),
    "no samples"
  )
  expect_null(rec)
})

test_that("print.MkPosterior works with NULL data and mcmc", {
  params <- c("log_posterior", "log_likelihood", "tree_length")
  samples <- matrix(rnorm(15), nrow = 5, ncol = 3,
                    dimnames = list(NULL, params))
  result <- MkPosterior(
    samples = samples, trees = list(), acceptance = numeric(0),
    model = NULL, data = NULL, mcmc = NULL, warmup = 0L, tuning = NULL
  )
  result$partial <- TRUE
  result$nSamples <- 5L
  result$stop_reason <- "recovered"

  # Should not error (cli output goes to messages, not stdout)
  expect_no_error(print(result))
})
