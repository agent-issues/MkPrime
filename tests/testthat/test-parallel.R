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

test_that("pool dispatch: wave behaviour confirmed by launch times (PAR-004)", {
  # Regression: a bug that launched all nRuns workers simultaneously
  # (ignoring nCore) would dispatch runs 3 and 4 at the same time as
  # runs 1 and 2. With correct wave dispatch (nCore = 2, nRuns = 4),
  # runs 3 and 4 cannot launch until the polling loop detects that at
  # least one of runs 1/2 has finished. We assert this using
  # result$.par_launch_times — parent-side Sys.time() stamps recorded in
  # .RunParallelRuns at the moment each callr worker is dispatched.
  #
  # Why parent launch times rather than log-file mtimes?
  # NTFS mtime has 1-second resolution. On a fast machine all four short
  # MCMC runs can finish writing within one NTFS tick, making
  # min(mtimes[3:4]) == max(mtimes[1:2]) and the old expect_gte() pass
  # trivially even when all four workers were launched simultaneously.
  # Sys.time() in the parent process has ~10-16 ms resolution on Windows,
  # far finer than the ~1 s polling gap, so a false pass is impossible.
  #
  # Flakiness margin: pollInterval = 1L and each run takes ~1-2 s, so
  # wave 2 launches at least ~1 s after wave 1 — well above the 10-16 ms
  # clock resolution.
  #
  # TODO(PAR-004 followup): interrupt behaviour (Ctrl-C during nCore=2
  # run should kill workers and show the parallel-specific message added
  # under PAR-002) is not covered by automated tests — triggering SIGINT
  # reliably from testthat is platform-dependent. Verified manually.
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
      nIter        = 2000L,
      maxWarmup    = 1000L,
      minWarmup    = 1000L,
      autoTune     = FALSE,
      nCore        = 2L,
      pollInterval = 1L
    ))

  expect_s3_class(result, "MkPosterior")
  expect_equal(result$nRuns, 4L)

  # Retrieve parent-side launch timestamps (numeric seconds since epoch).
  t <- result$.par_launch_times
  expect_length(t, 4L)
  expect_false(any(is.na(t)),
               info = "All four workers must have been launched")
  # Wave 1 (runs 1, 2) launches first; wave 2 (runs 3, 4) can only launch
  # after the polling loop detects a wave-1 finish, so the earliest wave-2
  # launch must be strictly later than the latest wave-1 launch.
  expect_gt(min(t[3:4]), max(t[1:2]))
})

test_that("maxTime fires mid-pool with nRuns > nCore returns valid result (PAR-001 regression)", {
  # Pre-fix: crashed with "missing value where TRUE/FALSE needed" because
  # .BuildResult dereferences fields on initial-state stubs for unlaunched
  # pool slots. Post-fix: .RunParallelRuns shrinks the result list to only
  # runs that produced output.
  skip_if_not_installed("callr")
  skip_if_not(.is_mkprime_installed(),
              "MkPrime not installed — callr workers need installed package")
  library("ape")

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)

  # nRuns = 6, nCore = 2 → 3 waves. maxTime = 3L fires before the queue
  # drains, leaving some pool slots unlaunched.
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(
      nRuns        = 6L,
      nCore        = 2L,
      nIter        = Inf,
      maxWarmup    = 1000L,
      minWarmup    = 1000L,
      autoTune     = FALSE,
      pollInterval = 1L,
      maxTime      = 3L
    ))

  expect_s3_class(result, "MkPosterior")
  # The primary regression check is that the code does not crash when maxTime
  # fires with unlaunched pool slots pending. nRuns may be 0 if maxTime fires
  # before any worker finishes (e.g. slow callr startup on CI).
  expect_gte(result$nRuns, 0L)
  expect_lte(result$nRuns, 6L)
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

# --- PAR-012: cancelGrace validation tests ---

test_that("MkPrimeMCMC() validates cancelGrace", {
  expect_error(MkPrimeMCMC(cancelGrace = 0L),   "cancelGrace")
  expect_error(MkPrimeMCMC(cancelGrace = -1L),  "cancelGrace")
  expect_error(suppressWarnings(MkPrimeMCMC(cancelGrace = "foo")), "cancelGrace")
})

test_that("MkPrimeMCMC() accepts cancelGrace as positive integer", {
  mcmc <- MkPrimeMCMC(cancelGrace = 5L)
  expect_equal(mcmc$cancelGrace, 5L)
})

test_that("MkPrimeMCMC() accepts cancelGrace = Inf", {
  mcmc <- MkPrimeMCMC(cancelGrace = Inf)
  expect_true(is.infinite(mcmc$cancelGrace))
})

# --- PAR-008/009: dropped_runs populated on early-stop with short cancelGrace ---

test_that("dropped_runs populated when workers are killed by short cancelGrace (PAR-008/009)", {
  # nRuns = 4, nCore = 2, maxTime = 2s, cancelGrace = 1s.
  # Workers in-flight when maxTime fires won't finish their batch in 1s,
  # so they get hard-killed → dropped_runs should be non-empty.
  # We also check requested_nRuns is preserved.
  skip_if_not_installed("callr")
  skip_if_not(.is_mkprime_installed(),
              "MkPrime not installed — callr workers need installed package")
  library("ape")

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)

  # Use a long checkEvery so workers are mid-batch when cancelled.
  result <- suppressWarnings(RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(
      nRuns        = 4L,
      nCore        = 2L,
      nIter        = Inf,
      maxWarmup    = 50000L,
      minWarmup    = 50000L,
      autoTune     = FALSE,
      pollInterval = 1L,
      maxTime      = 2L,
      cancelGrace  = 1L,
      checkEvery   = 10000L
    )))

  expect_s3_class(result, "MkPosterior")
  # requested_nRuns must equal the configured nRuns (PAR-009)
  expect_equal(result$requested_nRuns, 4L)
  # dropped_runs must be a data.frame with the expected columns
  expect_s3_class(result$dropped_runs, "data.frame")
  expect_true(all(c("run", "reason", "message", "wait_s") %in%
                    names(result$dropped_runs)))
  # Given maxTime=2s, cancelGrace=1s, maxWarmup=50000, at least some runs must
  # be dropped (killed mid-batch or unlaunched before the pool drains).
  expect_gt(nrow(result$dropped_runs), 0L)
  # Conservation invariant: completed + dropped == requested (PAR-009)
  #
  # The fallback must be 1L, not 0L. `.BuildResult()` sets `$nRuns` and
  # `$per_run` only when more than one run completed -- per-run diagnostics
  # need at least two -- and sets `$nRuns <- 0L` explicitly when none did. So
  # an absent `$nRuns` means exactly one, which is the `%||% 1L` convention
  # every other consumer uses (`R/burnin.R`, `R/Convergence.R`,
  # `R/MkPosterior.R`) and which `test-tempering.R` pins directly.
  #
  # With `%||% 0L` this counted a lone surviving run as zero, so whenever the
  # timing left 1 of 4 runs alive the invariant read 0 + 3 != 4. That is not
  # hypothetical: it is why the Windows `Code coverage` job failed on `main`
  # and on every open PR. Instrumentation slows the workers enough to change
  # how many survive the 1s cancel grace, and nothing else in CI runs this
  # test slowly enough to reach the one-survivor case.
  nCompleted <- if (!is.null(result$per_run)) {
    length(result$per_run)
  } else {
    result$nRuns %||% 1L
  }
  expect_equal(result$requested_nRuns, nCompleted + nrow(result$dropped_runs))
  # Each run is accounted for exactly once.
  expect_equal(anyDuplicated(result$dropped_runs$run), 0L)
  expect_true(all(result$dropped_runs$run %in% seq_len(4L)))
  # print() must not error when dropped_runs is populated (covers the new
  # print.MkPosterior branch for PAR-009 display)
  expect_no_error(capture.output(print(result)))
})
