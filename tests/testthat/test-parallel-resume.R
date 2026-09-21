# Tests for parallel-mode resumability after SIGKILL / abnormal parent exit
#
# These tests verify that parallel runs are recoverable when the parent
# orchestrator dies without running R-level cleanup -- the case SLURM
# walltime overrun produces on Hamilton.  Workers write per-run .ckp files
# every checkEvery iterations; ResumeMkPrime() synthesises the master
# checkpoint from those files when the parent did not exit cleanly.

.is_mkprime_installed <- function() {
  if (!requireNamespace("pkgload", quietly = TRUE)) return(TRUE)
  !isTRUE(pkgload::is_dev_package("MkPrime"))
}

# Small fixture re-used across tests
.tiny_fixture <- function() {
  list(
    tree = ape::read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);"),
    pd   = MatrixToPhyDat(matrix(
      c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
      dimnames = list(paste0("t", 1:4), NULL)
    ))
  )
}


# --- .CkpFilePaths helper ---

test_that(".CkpFilePaths derives per-run paths like .TreeFilePaths", {
  # NULL passes through
  expect_null(MkPrime:::.CkpFilePaths(NULL, 4L))
  # Single run: path unchanged
  expect_identical(MkPrime:::.CkpFilePaths("/tmp/run.ckp", 1L),
                   "/tmp/run.ckp")
  # Multi-run: base_1.ckp, base_2.ckp, ...
  expect_identical(
    MkPrime:::.CkpFilePaths("/tmp/run.ckp", 3L),
    c("/tmp/run_1.ckp", "/tmp/run_2.ckp", "/tmp/run_3.ckp")
  )
  # No extension: still numbered suffix
  expect_identical(
    MkPrime:::.CkpFilePaths("/tmp/run", 2L),
    c("/tmp/run_1", "/tmp/run_2")
  )
})


# --- .SaveCheckpoint atomic write ---

test_that(".SaveCheckpoint writes atomically (tmp+rename)", {
  ckp <- tempfile(fileext = ".ckp")
  on.exit(unlink(c(ckp, paste0(ckp, ".tmp"))), add = TRUE)

  # Minimal valid payload via a real (tiny) RunMkPrime invocation
  fx <- .tiny_fixture()
  set.seed(101)
  allow_warning(
    RunMkPrime(fx$pd, fx$tree,
      mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 200L, thin = 5L,
                          maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                          checkEvery = 100L, checkpointFile = ckp)),
    "without stabilisation"
  )

  expect_true(file.exists(ckp))
  # The .tmp file should not exist after a successful write
  expect_false(file.exists(paste0(ckp, ".tmp")))
  # The .ckp should be a valid RDS
  cp <- readRDS(ckp)
  expect_true(is.list(cp$runs))
  expect_equal(cp$version, 3L)
})


# --- .SynthesiseMasterFromPerRun ---

test_that(".SynthesiseMasterFromPerRun rebuilds master from per-run ckps", {
  skip_if_not_installed("callr")
  skip_if_not(.is_mkprime_installed(),
              "MkPrime not installed -- callr workers need installed package")

  fx <- .tiny_fixture()
  td <- tempfile("mkp_synth_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ckp <- file.path(td, "run.ckp")
  log <- file.path(td, "run.log")

  # Run a brief parallel run to completion to populate per-run ckps + master.
  suppressWarnings(RunMkPrime(fx$pd, fx$tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nCore = 2L, nIter = 800L, thin = 5L,
                        maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        checkEvery = 100L,
                        checkpointFile = ckp, logFile = log)))

  expect_true(file.exists(ckp))
  expect_true(file.exists(file.path(td, "run_1.ckp")))
  expect_true(file.exists(file.path(td, "run_2.ckp")))

  # Touch a per-run ckp to make it newer than master -> synthesis fires
  Sys.sleep(1.1)  # ensure mtime resolution distinguishes
  perRun1 <- file.path(td, "run_1.ckp")
  Sys.setFileTime(perRun1, Sys.time())

  origMtime <- file.mtime(ckp)
  synthIter <- MkPrime:::.SynthesiseMasterFromPerRun(ckp, 2L)
  expect_true(synthIter > 0L)
  # Master rewritten
  expect_true(file.mtime(ckp) > origMtime)
})


test_that(".SynthesiseMasterFromPerRun rebuilds when master is missing", {
  skip_if_not_installed("callr")
  skip_if_not(.is_mkprime_installed(),
              "MkPrime not installed -- callr workers need installed package")

  fx <- .tiny_fixture()
  td <- tempfile("mkp_synth2_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ckp <- file.path(td, "run.ckp")
  log <- file.path(td, "run.log")
  tre <- file.path(td, "run.nwk")

  suppressWarnings(RunMkPrime(fx$pd, fx$tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nCore = 2L, nIter = 800L, thin = 5L,
                        maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        checkEvery = 100L,
                        checkpointFile = ckp, logFile = log, treeFile = tre)))

  expect_true(file.exists(file.path(td, "run_1.ckp")))
  unlink(ckp)
  expect_false(file.exists(ckp))

  synthIter <- MkPrime:::.SynthesiseMasterFromPerRun(ckp, 2L)
  expect_true(synthIter > 0L)
  expect_true(file.exists(ckp))

  cp <- readRDS(ckp)
  expect_equal(length(cp$runs), 2L)
  # Critical: synthesised master must carry per-run treeFilePaths, not the
  # master tree path -- otherwise ResumeMkPrime's tree-truncation logic
  # would target the wrong file and trees would drift out of sync.
  expect_equal(length(cp$treeFilePaths), 2L)
  expect_true(grepl("_1\\.nwk$", cp$treeFilePaths[1L]))
  expect_true(grepl("_2\\.nwk$", cp$treeFilePaths[2L]))
})


# --- End-to-end: kill mid-flight, then resume ---

test_that("parallel run resumes after orchestrator process killed (SIGKILL)", {
  skip_on_cran()
  skip_on_os("solaris")
  skip_if_not_installed("callr")
  skip_if_not(.is_mkprime_installed(),
              "MkPrime not installed -- callr workers need installed package")

  fx <- .tiny_fixture()
  td <- tempfile("mkp_kill_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ckp <- file.path(td, "run.ckp")
  log <- file.path(td, "run.log")

  # Spawn the orchestrator in a child R process so we can kill it abruptly
  # without crashing the test session.  supervise=TRUE propagates kill to
  # the callr grandchildren (the actual MkPrime workers).
  proc <- callr::r_bg(
    func = function(pd, tree, ckp, log) {
      MkPrime::RunMkPrime(pd, tree,
        mcmc = MkPrime::MkPrimeMCMC(
          nRuns = 2L, nCore = 2L, nIter = 20000L, thin = 5L,
          maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
          checkEvery = 50L,
          checkpointFile = ckp, logFile = log))
    },
    args = list(pd = fx$pd, tree = fx$tree, ckp = ckp, log = log),
    supervise = TRUE
  )

  # Wait for per-run ckps to materialise.  Poll up to 15s; bail with a
  # diagnostic if workers never produced a checkpoint (probably a
  # platform-specific worker-launch failure rather than the feature under
  # test).
  startTime <- Sys.time()
  perRun1 <- file.path(td, "run_1.ckp")
  perRun2 <- file.path(td, "run_2.ckp")
  while (Sys.time() - startTime < 15) {
    if (file.exists(perRun1) && file.exists(perRun2)) break
    Sys.sleep(0.25)
  }

  # Kill the orchestrator without giving it a chance to clean up.
  if (proc$is_alive()) proc$kill()
  # Give callr's supervisor a moment to reap grandchildren.
  Sys.sleep(0.5)

  if (!file.exists(perRun1) || !file.exists(perRun2)) {
    skip("Per-run ckps never appeared -- worker launch may have failed.")
  }

  # Read iter from each per-run ckp -- both should have made some progress.
  cp1 <- readRDS(perRun1)
  cp2 <- readRDS(perRun2)
  expect_true(cp1$iter > 0L)
  expect_true(cp2$iter > 0L)

  # Resume in the test process (small enough not to recurse into callr).
  resumed <- suppressWarnings(ResumeMkPrime(ckp, fx$pd, fx$tree))
  expect_s3_class(resumed, "MkPosterior")
  # When logFile is user-supplied, samples remain on disk and resumed$samples
  # is NULL.  Check nSamples (always populated) and verify the log files
  # grew.  Both indicate the resume actually produced output.
  expect_gt(resumed$nSamples, 0L)
  expect_true(file.exists(file.path(td, "run_1.log")))
  logSize <- file.info(file.path(td, "run_1.log"))$size
  expect_gt(logSize, 0L)
})


test_that("ResumeMkPrime synthesises master when master.ckp is missing", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not(.is_mkprime_installed(),
              "MkPrime not installed -- callr workers need installed package")

  fx <- .tiny_fixture()
  td <- tempfile("mkp_nomaster_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ckp <- file.path(td, "run.ckp")
  log <- file.path(td, "run.log")

  # Brief parallel run to completion -- populates per-run + master ckps.
  suppressWarnings(RunMkPrime(fx$pd, fx$tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nCore = 2L, nIter = 800L, thin = 5L,
                        maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        checkEvery = 100L,
                        checkpointFile = ckp, logFile = log)))

  expect_true(file.exists(ckp))
  unlink(ckp)
  expect_false(file.exists(ckp))
  # Per-run ckps still present
  expect_true(file.exists(file.path(td, "run_1.ckp")))

  # Resume: should synthesise master from per-run ckps and continue.
  resumed <- suppressWarnings(ResumeMkPrime(ckp, fx$pd, fx$tree))
  expect_s3_class(resumed, "MkPosterior")
  expect_true(file.exists(ckp))  # master synthesised
  expect_gt(resumed$nSamples, 0L)
})
