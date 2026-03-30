# M-151: plot(), print(), .PostBurninData() for streaming multi-run results
#
# Streaming multi-run sets per_run[[i]]$samples = NULL because samples
# are written to disk. These tests verify that all MkPosterior methods
# handle this correctly.

# --- Helpers ---

# Create a fake streaming multi-run MkPosterior with NULL per-run samples
# and TSV log files on disk.
.FakeStreamingMultiRun <- function(nPerRun = 20L, nRuns = 2L) {
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                   "rate_log_sd", "p")

  set.seed(4621)
  logPaths <- character(nRuns)
  for (run in seq_len(nRuns)) {
    mat <- matrix(rnorm(nPerRun * length(paramNames)), nPerRun,
                  length(paramNames))
    colnames(mat) <- paramNames
    df <- data.frame(Sample = seq_len(nPerRun) * 10L, mat,
                     check.names = FALSE)
    logPaths[run] <- tempfile(fileext = ".log")
    utils::write.table(df, logPaths[run], sep = "\t", row.names = FALSE,
                       quote = FALSE)
  }

  combined <- do.call(rbind, lapply(logPaths, ReadMkLog))
  emptyMat <- matrix(numeric(0), nrow = 0L, ncol = length(paramNames),
                     dimnames = list(NULL, paramNames))

  mcmc <- suppressWarnings(
    MkPrimeMCMC(nIter = nPerRun * 10L, thin = 10L, maxWarmup = 50L)
  )
  mkd <- structure(list(nChar = 10L,
                        type = rep("transformational", 10)),
                   class = "MkPrimeData")

  result <- MkPosterior(
    samples    = emptyMat,
    trees      = list(),
    acceptance = c(nni = 0.3, spr = 0.1),
    model      = MkPrimeModel(),
    data       = mkd,
    mcmc       = mcmc,
    warmup     = 100L,
    tuning     = list()
  )
  result$nRuns    <- nRuns
  result$logFile  <- logPaths
  result$nSamples <- nPerRun * nRuns

  # Mimic .BuildResult streaming: per_run with NULL samples

  result$per_run <- lapply(seq_len(nRuns), function(i) {
    list(samples = NULL, trees = list(), acceptance = c(nni = 0.3, spr = 0.1),
         saved_idx = nPerRun)
  })

  list(result = result, logPaths = logPaths, combined = combined)
}


# --- Tests ---

test_that("M-151a: plot() works for streaming multi-run (NULL per-run samples)", {
  info <- .FakeStreamingMultiRun()
  on.exit(unlink(info$logPaths))

  # Should auto-load per-run samples from log files via .PostBurninData
  expect_no_error(plot(info$result))
})


test_that("M-151b: print() works for streaming multi-run (NULL per-run samples)", {
  info <- .FakeStreamingMultiRun()
  on.exit(unlink(info$logPaths))

  # Should not error or produce empty sample counts
  expect_no_error(print(info$result))

  # Verify the per-run count fallback logic directly:
  # nrow(NULL) is NULL, so saved_idx should be used
  perRunN <- nrow(info$result$per_run[[1]]$samples) %||%
             info$result$per_run[[1]]$saved_idx
  expect_equal(perRunN, 20L)
})


test_that("M-151c: .PostBurninData() with burnin > 0 for streaming multi-run", {
  info <- .FakeStreamingMultiRun(nPerRun = 30L)
  on.exit(unlink(info$logPaths))

  # Set a burnin
  info$result$burnin <- 5L

  # Should not crash — loads per-run samples from log files
  pb <- MkPrime:::.PostBurninData(info$result)

  # Combined samples should have (30 - 5) * 2 = 50 rows
  expect_equal(nrow(pb$samples), 50L)

  # Per-run samples should have 25 rows each
  expect_equal(nrow(pb$per_run[[1]]$samples), 25L)
  expect_equal(nrow(pb$per_run[[2]]$samples), 25L)
})


test_that("M-151: .PostBurninData() auto-loads per-run from log files (bi=0)", {
  info <- .FakeStreamingMultiRun()
  on.exit(unlink(info$logPaths))

  pb <- MkPrime:::.PostBurninData(info$result)

  # Combined samples should be loaded
  expect_gt(nrow(pb$samples), 0L)

  # Per-run samples should be loaded (no longer NULL)
  expect_false(is.null(pb$per_run[[1]]$samples))
  expect_false(is.null(pb$per_run[[2]]$samples))
  expect_equal(nrow(pb$per_run[[1]]$samples), 20L)
})


test_that("M-151: summary() works for streaming multi-run", {
  info <- .FakeStreamingMultiRun()
  on.exit(unlink(info$logPaths))

  s <- summary(info$result)
  expect_true(is.data.frame(s))
  expect_true("ESS" %in% names(s))
})


test_that("M-151: SetBurnin works for streaming multi-run (NULL per-run)", {
  info <- .FakeStreamingMultiRun()
  on.exit(unlink(info$logPaths))

  # Should use saved_idx fallback, not crash on nrow(NULL)
  result2 <- SetBurnin(info$result, 5)
  expect_equal(result2$burnin, 5L)
})


test_that("M-151: plot() falls back gracefully when per-run unresolvable", {
  info <- .FakeStreamingMultiRun()
  # Delete log files so per-run can't be loaded
  unlink(info$logPaths)
  info$result$logFile <- NULL

  # Should fall back to single-color trace (combined samples already loaded)
  info$result$samples <- info$combined
  expect_no_error(plot(info$result))
})
