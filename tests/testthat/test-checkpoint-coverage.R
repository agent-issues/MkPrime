# Checkpoints taken part way through a job (#14), their size (#106), and the
# paths they record (#31).
skip_under_memcheck()

CkpData <- function() {
  set.seed(1)
  tree <- ape::rtree(8)
  tree$tip.label <- paste0("t", 1:8)
  mat <- matrix(sample(0:1, 8 * 20, TRUE), 8, 20,
                dimnames = list(tree$tip.label, NULL))
  list(tree = tree, pd = MatrixToPhyDat(mat))
}

LogRows <- function(logFile) {
  length(grep("^[0-9]", readLines(logFile)))
}

LogIters <- function(logFile) {
  rows <- grep("^[0-9]", readLines(logFile), value = TRUE)
  as.integer(vapply(strsplit(rows, "\t"), `[`, "", 1L))
}

# A progress callback that fails the job as a kill would, once `Stop(info)`.
Killer <- function(Stop) {
  function(info) if (Stop(info)) stop("killed")
}

test_that("a kill in run 2 keeps run 1's samples (#14)", {
  d <- CkpData()
  for (maxRhat in list(NULL, 1.01)) {
    dir <- withr::local_tempdir()
    log <- file.path(dir, "run.log")
    ckp <- file.path(dir, "run.ckp")
    seen <- new.env()
    seen$warmupEnds <- 0L
    # Run 2's first warmup batch ends the second time iteration 500 is seen.
    Stop <- function(info) {
      if (info$iter == 500L) seen$warmupEnds <- seen$warmupEnds + 1L
      seen$warmupEnds == 2L
    }
    set.seed(2)
    expect_error(allow_warning(RunMkPrime(
      d$pd, d$tree,
      mcmc = MkPrimeMCMC(
        nRuns = 2L, nIter = 1000L, thin = 5L, minWarmup = 500L,
        maxWarmup = 500L, autoTune = FALSE, checkEvery = 500L,
        bufferSize = 50L, logFile = log, maxRhat = maxRhat, maxTime = 60,
        plotEvery = 500L, progressFn = Killer(Stop)
      )
    ), "maxWarmup"), "killed")

    ck <- readRDS(ckp)
    expect_gt(ck$iter, 0)
    expect_length(ck$runs, 2L)
    run1Log <- file.path(dir, "run_1.log")
    run1Rows <- LogRows(run1Log)
    expect_gt(run1Rows, 0L)
    expect_identical(ck$runs[[1]]$saved_idx, run1Rows)

    resumed <- allow_warning(ResumeMkPrime(
      ckp, d$pd, d$tree,
      mcmc = list(progressFn = function(info) NULL, maxTime = 60)
    ), "maxWarmup")
    expect_identical(LogRows(run1Log), run1Rows)
    expect_gt(LogRows(file.path(dir, "run_2.log")), 0L)
  }
})

test_that("a checkpoint resumes a run from the batch it was saved at", {
  d <- CkpData()
  dir <- withr::local_tempdir()
  log <- file.path(dir, "run.log")
  set.seed(3)
  expect_error(allow_warning(RunMkPrime(
    d$pd, d$tree,
    mcmc = MkPrimeMCMC(
      nRuns = 1L, nIter = 3000L, thin = 5L, minWarmup = 2000L,
      maxWarmup = 2000L, autoTune = FALSE, checkEvery = 500L,
      logFile = log, maxTime = 60, plotEvery = 500L,
      progressFn = Killer(function(info) info$iter >= 1000L)
    )
  ), "maxWarmup"), "killed")
  ck <- readRDS(file.path(dir, "run.ckp"))
  expect_equal(ck$iter, 1000)
  expect_equal(ck$runs[[1]]$actual_iter, 1000)
})

test_that("a kill mid-epoch checkpoints the epoch, not the inner run (#14)", {
  d <- CkpData()
  dir <- withr::local_tempdir()
  log <- file.path(dir, "run.log")
  ckp <- file.path(dir, "run.ckp")
  seen <- new.env()
  # The first callback after phase 2 begins ends run 1's first epoch batch.
  Stop <- function(info) {
    phase <- tryCatch(readRDS(ckp)$serialPhase, error = function(e) NULL)
    seen$iter <- info$iter
    identical(phase, 2L)
  }
  set.seed(4)
  expect_error(allow_warning(RunMkPrime(
    d$pd, d$tree,
    mcmc = MkPrimeMCMC(
      nRuns = 2L, nIter = 20000L, thin = 5L, minWarmup = 500L,
      maxWarmup = 500L, autoTune = FALSE, checkEvery = 1000L,
      bufferSize = 50L, logFile = log, maxRhat = 1.00001, maxTime = 60,
      plotEvery = 500L, progressFn = Killer(Stop)
    )
  ), "maxWarmup"), "killed")

  ck <- readRDS(ckp)
  expect_identical(ck$serialPhase, 2L)
  expect_equal(ck$runs[[1]]$actual_iter, seen$iter)
  # The inner run's settings would reroute the resume.
  expect_equal(ck$mcmc$maxRhat, 1.00001)
  expect_equal(ck$mcmc$nIter, 20000L)
})

test_that("checkpoints omit trees already streamed to a tree file (#106)", {
  run <- list(
    chains = list(), saved_idx = 1L, tree_saved_idx = 1L,
    tree_samples = list(ape::rtree(4))
  )
  ckp <- withr::local_tempfile(fileext = ".ckp")
  mcmc <- list(nRuns = 1L, logFile = "run.log")

  .SaveCheckpoint(list(run), mcmc, 1L, "log_post", ckp)
  expect_length(readRDS(ckp)$runs[[1]]$tree_samples, 1L)

  mcmc$treeFile <- "run.nwk"
  .SaveCheckpoint(list(run), mcmc, 1L, "log_post", ckp)
  expect_null(readRDS(ckp)$runs[[1]]$tree_samples)
})

test_that("resume reads streamed trees back from the tree file (#106)", {
  d <- CkpData()
  dir <- withr::local_tempdir()
  log <- file.path(dir, "run.log")
  ckp <- file.path(dir, "run.ckp")
  set.seed(5)
  expect_error(allow_warning(RunMkPrime(
    d$pd, d$tree,
    mcmc = MkPrimeMCMC(
      nRuns = 1L, nIter = 2000L, thin = 5L, minWarmup = 500L,
      maxWarmup = 500L, autoTune = FALSE, checkEvery = 500L,
      bufferSize = 50L, logFile = log, maxTime = 60, plotEvery = 500L,
      progressFn = Killer(function(info) info$iter >= 2000L)
    )
  ), "maxWarmup"), "killed")
  nBefore <- readRDS(ckp)$runs[[1]]$tree_saved_idx
  expect_gt(nBefore, 0L)

  resumed <- ResumeMkPrime(
    ckp, d$pd, d$tree,
    mcmc = list(progressFn = function(info) NULL, nIter = 2500L)
  )
  nTrees <- length(readLines(file.path(dir, "run_trees.nwk")))
  expect_gt(nTrees, nBefore)
  expect_length(resumed$trees, nTrees)
  expect_false(any(vapply(resumed$trees, is.null, logical(1))))
  expect_identical(resumed$trees[[1]]$tip.label, rownames(resumed$data$matrix))
})

test_that("a run started with relative paths resumes from elsewhere (#31)", {
  d <- CkpData()
  runDir <- withr::local_tempdir()
  elsewhere <- withr::local_tempdir()
  withr::with_dir(runDir, {
    set.seed(6)
    allow_warning(RunMkPrime(
      d$pd, d$tree,
      mcmc = MkPrimeMCMC(
        nRuns = 1L, nIter = 1000L, thin = 5L, minWarmup = 500L,
        maxWarmup = 500L, autoTune = FALSE, logFile = "run.log", maxTime = 60
      )
    ), "maxWarmup")
  })
  ckp <- file.path(runDir, "run.ckp")
  log <- file.path(runDir, "run.log")
  expect_true(.IsAbsolutePath(readRDS(ckp)$logFilePaths))
  rows <- LogRows(log)

  withr::with_dir(elsewhere, {
    allow_warning(ResumeMkPrime(ckp, d$pd, d$tree, mcmc = list(nIter = 1500L)),
                  "maxWarmup")
  })
  expect_gt(LogRows(log), rows)
  expect_true(all(diff(LogIters(log)) > 0L))
  expect_length(list.files(elsewhere), 0L)
})

test_that("relative paths in older checkpoints resolve to the run's directory (#31)", {
  d <- CkpData()
  runDir <- withr::local_tempdir()
  dir.create(file.path(runDir, "out"))
  elsewhere <- withr::local_tempdir()
  withr::with_dir(runDir, {
    set.seed(7)
    allow_warning(RunMkPrime(
      d$pd, d$tree,
      mcmc = MkPrimeMCMC(
        nRuns = 1L, nIter = 1000L, thin = 5L, minWarmup = 500L,
        maxWarmup = 500L, autoTune = FALSE, logFile = "out/run.log",
        maxTime = 60
      )
    ), "maxWarmup")
  })
  ckp <- file.path(runDir, "out", "run.ckp")
  log <- file.path(runDir, "out", "run.log")
  ck <- readRDS(ckp)
  ck$logFilePaths <- "out/run.log"
  ck$treeFilePaths <- "out/run_trees.nwk"
  ck$mcmc$logFile <- "out/run.log"
  ck$mcmc$treeFile <- "out/run_trees.nwk"
  ck$mcmc$checkpointFile <- "out/run.ckp"
  saveRDS(ck, ckp)
  rows <- LogRows(log)

  withr::with_dir(elsewhere, {
    allow_warning(ResumeMkPrime(ckp, d$pd, d$tree, mcmc = list(nIter = 1500L)),
                  "maxWarmup")
  })
  expect_gt(LogRows(log), rows)
  expect_length(list.files(elsewhere), 0L)

  # With the logs gone, a relative path must not pass for a cleaned temp log.
  saveRDS(ck, ckp)
  file.rename(log, file.path(runDir, "moved.log"))
  withr::with_dir(elsewhere, {
    expect_error(ResumeMkPrime(ckp, d$pd, d$tree), "relative paths")
  })
})

test_that("resuming with the relative paths a run started with is silent", {
  withr::with_dir(tempdir(), {
    stored <- list(logFile = file.path(getwd(), "run.log"), nIter = 10L)
    expect_no_warning(
      expect_identical(.ResumeMcmc(stored, list(logFile = "run.log")), stored)
    )
    expect_warning(.ResumeMcmc(stored, list(logFile = "other.log")),
                   "keeps the checkpoint")
  })
})

test_that(".RunDirectory recovers the directory a run started in", {
  expect_identical(.RunDirectory("out/run.ckp", "/a/b/out/run.ckp"), "/a/b")
  expect_identical(.RunDirectory("./run.ckp", "/a/b/run.ckp"), "/a/b")
  withr::with_dir(tempdir(), {
    here <- getwd()
    expect_identical(.RunDirectory(NULL, "/a/run.ckp"), here)
    expect_identical(.RunDirectory("/x/run.ckp", "/a/run.ckp"), here)
    expect_identical(.RunDirectory("../run.ckp", "/a/run.ckp"), here)
    expect_identical(.RunDirectory("out/run.ckp", "/a/renamed/run.ckp"), here)
  })
  expect_identical(.AbsolutePath(c("/a/x", "y"), "/b"), c("/a/x", "/b/y"))
  expect_identical(.AbsolutePath(character(0)), character(0))
  expect_null(.AbsolutePath(NULL))
})

test_that("a checkpoint without a log keeps its samples past the session (#31)", {
  d <- CkpData()
  dir <- withr::local_tempdir()
  ckp <- file.path(dir, "run.ckp")
  set.seed(8)
  expect_error(allow_warning(RunMkPrime(
    d$pd, d$tree,
    mcmc = MkPrimeMCMC(
      nRuns = 1L, nIter = 3000L, thin = 5L, minWarmup = 500L,
      maxWarmup = 500L, autoTune = FALSE, bufferSize = 50L,
      checkpointFile = ckp, maxTime = 60, plotEvery = 500L,
      progressFn = Killer(function(info) info$iter >= 3000L)
    )
  ), "maxWarmup"), "killed")

  log <- readRDS(ckp)$logFilePaths
  expect_identical(dirname(log), dirname(.AbsolutePath(ckp)))
  expect_true(file.exists(log))
  rows <- LogRows(log)
  expect_gt(rows, 0L)

  resumed <- ResumeMkPrime(
    ckp, d$pd, d$tree,
    mcmc = list(progressFn = function(info) NULL, nIter = 3500L)
  )
  expect_gt(LogRows(log), rows)
  expect_true(all(diff(LogIters(log)) > 0L))
})

test_that("a completed run deletes the log it kept beside the checkpoint", {
  d <- CkpData()
  dir <- withr::local_tempdir()
  ckp <- file.path(dir, "run.ckp")
  set.seed(9)
  res <- allow_warning(RunMkPrime(
    d$pd, d$tree,
    mcmc = MkPrimeMCMC(
      nRuns = 1L, nIter = 1000L, thin = 5L, minWarmup = 500L,
      maxWarmup = 500L, autoTune = FALSE, checkpointFile = ckp, maxTime = 60
    )
  ), "maxWarmup")
  expect_gt(nrow(res$samples), 0L)
  expect_null(res$logFile)
  expect_identical(list.files(dir), "run.ckp")
})

test_that("MkPrimeRecover() keeps a log the user named", {
  d <- CkpData()
  dir <- withr::local_tempdir()
  log <- file.path(dir, "run.log")
  set.seed(10)
  expect_warning(allow_warning(RunMkPrime(
    d$pd, d$tree,
    mcmc = MkPrimeMCMC(
      nRuns = 1L, nIter = 3000L, thin = 5L, minWarmup = 500L,
      maxWarmup = 500L, autoTune = FALSE, bufferSize = 50L, logFile = log,
      maxTime = 60, plotEvery = 500L,
      progressFn = function(info) if (info$iter >= 3000L) rlang::interrupt()
    )
  ), "maxWarmup"), "interrupted")
  recovered <- MkPrimeRecover()
  expect_gt(nrow(recovered$samples), 0L)
  expect_true(file.exists(log))
})

test_that("a run paused on maxTime keeps the log beside its checkpoint", {
  d <- CkpData()
  dir <- withr::local_tempdir()
  ckp <- file.path(dir, "run.ckp")
  set.seed(11)
  res <- allow_warning(RunMkPrime(
    d$pd, d$tree,
    mcmc = MkPrimeMCMC(
      nRuns = 1L, nIter = Inf, thin = 5L, minWarmup = 500L, maxWarmup = 500L,
      autoTune = FALSE, bufferSize = 50L, checkpointFile = ckp, maxTime = 0.5
    )
  ), "maxWarmup")
  expect_identical(res$stop_reason, "max_time")
  log <- readRDS(ckp)$logFilePaths
  expect_true(file.exists(log))
  rows <- LogRows(log)

  allow_warning(ResumeMkPrime(ckp, d$pd, d$tree, mcmc = list(maxTime = 0.5)),
                "maxWarmup")
  expect_gt(LogRows(log), rows)
})
