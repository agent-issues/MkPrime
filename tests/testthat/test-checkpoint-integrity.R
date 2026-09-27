# Checkpoint and resume integrity: per-move counts (#74), adapted thin and
# iteration numbers (#91), run indices and per-run merges (#92), and per-run
# checkpoint cleanup (#107).

.IntegrityFixture <- function() {
  list(
    tree = ape::read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);"),
    pd = MatrixToPhyDat(matrix(
      c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L), 4, 3,
      dimnames = list(paste0("t", 1:4), NULL)
    ))
  )
}

.LogIters <- function(logFile) {
  as.integer(rownames(ReadMkLog(logFile)))
}

# A checkpoint written before `fixTopology` travelled with it resumes on the
# free-topology schedule, which has more moves than the run's counts cover.
.LegacyFixTopologyCheckpoint <- function(mcmc) {
  fx <- .IntegrityFixture()
  set.seed(7401)
  allow_warning(RunMkPrime(fx$pd, fx$tree, mcmc = mcmc, fixTopology = TRUE),
                "without stabilisation")
  ck <- readRDS(mcmc$checkpointFile)
  ck$mcmc$fixTopology <- NULL
  saveRDS(ck, mcmc$checkpointFile)
  ck
}


test_that(".AlignMoveCounts matches counts by name (#74)", {
  moveNames <- c("a", "b", "c")
  expect_identical(.AlignMoveCounts(c(c = 3L, a = 1L, b = 2L), moveNames),
                   c(a = 1L, b = 2L, c = 3L))
  expect_identical(.AlignMoveCounts(1:3, moveNames), c(a = 1L, b = 2L, c = 3L))
  expect_error(.AlignMoveCounts(c(a = 1L, b = 2L), moveNames),
               "do not match the rebuilt move schedule")
  expect_error(.AlignMoveCounts(1:2, moveNames), "do not match")
})

test_that("a post-warmup resume on another schedule aborts, not recycles (#74)", {
  ckp <- tempfile(fileext = ".ckp")
  log <- sub("ckp$", "log", ckp)
  on.exit(unlink(c(ckp, log)), add = TRUE)
  ck <- .LegacyFixTopologyCheckpoint(MkPrimeMCMC(
    nRuns = 1L, nIter = 200L, thin = 5L, minWarmup = 100L, maxWarmup = 100L,
    autoTune = FALSE, checkEvery = 100L, maxTime = 60,
    checkpointFile = ckp, logFile = log))
  expect_identical(ck$runs[[1]]$phase, "Sample")

  expect_error(ResumeMkPrime(ckp, .IntegrityFixture()$pd,
                             mcmc = list(nIter = 400L)),
               "move set differs")
})

test_that("a mid-warmup resume on another schedule aborts clearly (#74)", {
  ckp <- tempfile(fileext = ".ckp")
  log <- sub("ckp$", "log", ckp)
  on.exit(unlink(c(ckp, log)), add = TRUE)
  ck <- .LegacyFixTopologyCheckpoint(MkPrimeMCMC(
    nRuns = 1L, nIter = 200L, thin = 5L, minWarmup = 1000L, maxWarmup = 1000L,
    checkEvery = 100L, maxTime = 60, checkpointFile = ckp, logFile = log))
  expect_identical(ck$phase, "Warmup")
  # Without stored weights only the counts themselves reveal the mismatch.
  ck$moveWeights <- NULL
  saveRDS(ck, ckp)

  expect_error(ResumeMkPrime(ckp, .IntegrityFixture()$pd,
                             mcmc = list(nIter = 400L)),
               "move counts do not match")
})

test_that("Sample records each row's true iteration once thin adapts (#91)", {
  fx <- .IntegrityFixture()
  log <- tempfile(fileext = ".log")
  on.exit(unlink(c(log, sub("log$", "ckp", log))), add = TRUE)
  set.seed(9101)
  res <- allow_warning(RunMkPrime(fx$pd, fx$tree, mcmc = MkPrimeMCMC(
    nRuns = 1L, nIter = 6000L, minWarmup = 200L, maxWarmup = 200L,
    autoTune = FALSE, checkEvery = 2000L, maxTime = 60, logFile = log)),
    "without stabilisation")
  iters <- .LogIters(log)
  gaps <- diff(iters)
  # This chain's autocorrelation always outgrows the initial thin.
  expect_gt(max(gaps), min(gaps))
  expect_true(all(gaps > 0L))
  expect_lte(max(iters), res$actual_iter)
  expect_gt(max(iters) + max(gaps), res$actual_iter)
})

test_that("every checkpoint records each run's adapted thin (#91)", {
  fx <- .IntegrityFixture()
  ckp <- tempfile(fileext = ".ckp")
  log <- sub("ckp$", "log", ckp)
  nwk <- sub("ckp$", "nwk", ckp)
  on.exit(unlink(c(ckp, .LogFilePaths(log, 2L), .TreeFilePaths(nwk, 2L))),
          add = TRUE)

  # The multiple-of-thin rule needs the thin that "auto" resolves to.
  set.seed(9102)
  probe <- allow_warning(RunMkPrime(fx$pd, fx$tree, mcmc = MkPrimeMCMC(
    nRuns = 1L, nIter = 300L, minWarmup = 200L, maxWarmup = 200L,
    autoTune = FALSE, maxTime = 60)), "without stabilisation")
  autoThin <- probe$mcmc$thin

  mcmc <- MkPrimeMCMC(nRuns = 2L, nIter = 6000L, minWarmup = 200L,
                      maxWarmup = 200L, autoTune = FALSE, checkEvery = 2000L,
                      treeThin = 2L * autoThin, maxTime = 60,
                      checkpointFile = ckp, logFile = log, treeFile = nwk)
  set.seed(9103)
  allow_warning(RunMkPrime(fx$pd, fx$tree, mcmc = mcmc),
                "without stabilisation")

  ck <- readRDS(ckp)
  runThin <- vapply(ck$runs, function(r) r$thin %||% NA_integer_, numeric(1))
  expect_true(all(runThin > autoThin))

  # Resuming at the stale thin would space the new rows differently, and
  # the tree rewind would count trees at the stale tree interval.
  resumed <- ResumeMkPrime(ckp, fx$pd, mcmc = list(nIter = 9000L))
  expect_identical(resumed$stop_reason, "max_iter")
  for (run in 1:2) {
    iters <- .LogIters(.LogFilePaths(log, 2L)[run])
    expect_true(all(diff(iters) > 0L))
    expect_identical(unique(diff(iters[iters > 6000L])),
                     as.integer(runThin[run]))
  }
})

test_that("checkpoints look up files by each run's original index (#92)", {
  run <- list(chains = list(), saved_idx = 0L)
  runs <- lapply(c(1L, 2L, 4L), function(i) c(run, run_index = i))
  ckp <- tempfile(fileext = ".ckp")
  on.exit(unlink(ckp), add = TRUE)
  .SaveCheckpoint(runs, list(nRuns = 4L, logFile = "out/base.log",
                             treeFile = "out/base.nwk"),
                  10L, "log_post", ckp)
  ck <- readRDS(ckp)
  expect_identical(ck$logFilePaths,
                   c("out/base_1.log", "out/base_2.log", "out/base_4.log"))
  expect_identical(ck$treeFilePaths,
                   c("out/base_1.nwk", "out/base_2.nwk", "out/base_4.nwk"))
})


# A finished two-run serial job: its master checkpoint and the per-run files
# a parallel worker would have written for the same state.
.TwoRunCheckpoint <- function() {
  fx <- .IntegrityFixture()
  dir <- tempfile("ckp_integrity_")
  dir.create(dir)
  withr::defer(unlink(dir, recursive = TRUE), envir = parent.frame())
  ckp <- file.path(dir, "job.ckp")
  set.seed(9201)
  allow_warning(RunMkPrime(fx$pd, fx$tree, mcmc = MkPrimeMCMC(
    nRuns = 2L, nIter = 400L, thin = 5L, minWarmup = 200L, maxWarmup = 200L,
    autoTune = FALSE, maxTime = 60, checkpointFile = ckp,
    logFile = file.path(dir, "job.log"))), "without stabilisation")
  master <- readRDS(ckp)
  perRun <- .CkpFilePaths(ckp, 2L)
  for (i in 1:2) {
    saveRDS(list(runs = master$runs[i], iter = 400L, mcmc = master$mcmc,
                 paramNames = master$paramNames, version = 3L,
                 timestamp = master$timestamp + 1), perRun[i])
  }
  list(ckp = ckp, dir = dir, perRun = perRun, master = master, pd = fx$pd)
}

test_that("a checkpoint missing a run fails with a clear message (#92)", {
  job <- .TwoRunCheckpoint()
  unlink(job$perRun)
  ck <- job$master
  ck$runs <- ck$runs[1]
  saveRDS(ck, job$ckp)
  expect_error(ResumeMkPrime(job$ckp, job$pd, mcmc = list(nIter = 600L)),
               "state for 1 of 2 runs")

  ck$mcmc$nRuns <- 3L
  saveRDS(ck, job$ckp)
  expect_error(ResumeMkPrime(job$ckp, job$pd, mcmc = list(nIter = 600L)),
               "No state for runs 2 and 3")
})

test_that("a parallel master keeps the runs its workers dropped (#92)", {
  launched <- lapply(1:3, function(i) list(run_index = i, phase = NULL))
  returned <- list(list(run_index = 1L, phase = "Sample"),
                   list(run_index = 3L, phase = "Sample"))
  kept <- .WithDroppedRuns(returned, launched)
  expect_identical(vapply(kept, `[[`, 0L, "run_index"), 1:3)
  expect_identical(kept[[2]], launched[[2]])
  expect_identical(kept[c(1, 3)], returned)
})

test_that("per-run checkpoints that disagree are not merged (#92)", {
  job <- .TwoRunCheckpoint()

  expect_identical(.SynthesiseMasterFromPerRun(job$ckp, 2L), 400L)

  job <- .TwoRunCheckpoint()
  stray <- readRDS(job$perRun[2])
  stray$paramNames <- rev(stray$paramNames)
  saveRDS(stray, job$perRun[2])
  expect_error(.SynthesiseMasterFromPerRun(job$ckp, 2L), "paramNames")

  job <- .TwoRunCheckpoint()
  stray <- readRDS(job$perRun[2])
  stray$runs[[1]]$chains <- stray$runs[[1]]$chains[1]
  stray$runs[[1]]$chains[[2]] <- stray$runs[[1]]$chains[[1]]
  stray$runs[[1]]$chains[[3]] <- stray$runs[[1]]$chains[[1]]
  saveRDS(stray, job$perRun[2])
  expect_error(.SynthesiseMasterFromPerRun(job$ckp, 2L), "nChains")

  # Merged into another run's slot.
  job <- .TwoRunCheckpoint()
  file.copy(job$perRun[1], job$perRun[2], overwrite = TRUE)
  expect_error(.SynthesiseMasterFromPerRun(job$ckp, 2L), "run index")
})

test_that("per-run checkpoints older than the master are refused (#107)", {
  job <- .TwoRunCheckpoint()
  before <- readRDS(job$ckp)
  Backdate <- function(path) {
    ck <- readRDS(path)
    ck$runs[[1]]$saved_idx <- 999L
    ck$timestamp <- before$timestamp - 3600
    saveRDS(ck, path)
  }
  Backdate(job$perRun[2])
  .SynthesiseMasterFromPerRun(job$ckp, 2L)
  expect_identical(readRDS(job$ckp)$runs[[2]]$saved_idx,
                   before$runs[[2]]$saved_idx)

  job <- .TwoRunCheckpoint()
  before <- readRDS(job$ckp)
  lapply(job$perRun, Backdate)
  expect_identical(.SynthesiseMasterFromPerRun(job$ckp, 2L), 0L)
})

test_that("overwrite = TRUE discards the previous run's per-run checkpoints (#107)", {
  fx <- .IntegrityFixture()
  dir <- tempfile("ckp107_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  ckp <- file.path(dir, "job.ckp")
  stale <- file.path(dir, c("job_1.ckp", "job_3.ckp", "job_2.ckp.tmp",
                            "job.ckp.tmp"))
  keep <- file.path(dir, c("job_4.ckp", "job_notes.ckp", "other_1.ckp"))
  saveRDS(list(mcmc = list(nRuns = 3L)), ckp)
  file.create(c(stale, keep))

  set.seed(10701)
  allow_warning(RunMkPrime(fx$pd, fx$tree, overwrite = TRUE,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 300L, thin = 5L, minWarmup = 200L,
                       maxWarmup = 200L, autoTune = FALSE, maxTime = 60,
                       checkpointFile = ckp,
                       logFile = file.path(dir, "job.log"))),
    "without stabilisation")
  expect_false(any(file.exists(stale)))
  expect_true(all(file.exists(keep)))
})

test_that("a temp-log parallel run leaves no per-run checkpoints (#107)", {
  skip_if_not_installed("callr")
  skip_if(requireNamespace("pkgload", quietly = TRUE) &&
            isTRUE(pkgload::is_dev_package("MkPrime")),
          "callr workers need the installed package")
  fx <- .IntegrityFixture()
  CkpFiles <- function() {
    list.files(tempdir(), pattern = "^mkp_run_.*\\.ckp(\\.tmp)?$")
  }
  before <- CkpFiles()
  set.seed(10702)
  allow_warning(RunMkPrime(fx$pd, fx$tree, mcmc = MkPrimeMCMC(
    nRuns = 2L, nCore = 2L, nIter = 1500L, thin = 5L, minWarmup = 200L,
    maxWarmup = 200L, autoTune = FALSE, checkEvery = 500L, pollInterval = 1L,
    maxTime = 60)), "without stabilisation")
  expect_identical(setdiff(CkpFiles(), before), character(0))
})

test_that("a parallel resume on another schedule stops before launching (#74)", {
  skip_if_not_installed("callr")
  skip_if(requireNamespace("pkgload", quietly = TRUE) &&
            isTRUE(pkgload::is_dev_package("MkPrime")),
          "callr workers need the installed package")
  ckp <- tempfile(fileext = ".ckp")
  log <- sub("ckp$", "log", ckp)
  on.exit(unlink(c(ckp, .CkpFilePaths(ckp, 2L), .LogFilePaths(log, 2L))),
          add = TRUE)
  ck <- .LegacyFixTopologyCheckpoint(MkPrimeMCMC(
    nRuns = 2L, nCore = 2L, nIter = 600L, thin = 5L, minWarmup = 200L,
    maxWarmup = 200L, autoTune = FALSE, checkEvery = 200L, pollInterval = 1L,
    maxTime = 60, checkpointFile = ckp, logFile = log))
  ck$moveWeights <- NULL
  saveRDS(ck, ckp)

  expect_error(ResumeMkPrime(ckp, .IntegrityFixture()$pd,
                             mcmc = list(nIter = 800L)),
               "move counts do not match")
  expect_length(readRDS(ckp)$runs, 2L)
})
