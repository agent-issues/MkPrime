# Resume and interrupt recovery: thin read by exact name (#314), a dropped
# worker's per-run checkpoint (#313), nIter at or below a checkpoint and the
# final save (#318), and the interrupt handlers' save and flush (#317).

.RecoveryFixture <- function() {
  list(
    tree = ape::read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);"),
    pd = MatrixToPhyDat(matrix(
      c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L), 4, 3,
      dimnames = list(paste0("t", 1:4), NULL)
    ))
  )
}

.RecoveryJob <- function(dir, nRuns = 1L, nIter = 400L, ...) {
  MkPrimeMCMC(
    nRuns = nRuns, nIter = nIter, thin = 5L, minWarmup = 200L,
    maxWarmup = 200L, autoTune = FALSE, maxTime = 60,
    checkpointFile = file.path(dir, "job.ckp"),
    logFile = file.path(dir, "job.log"), ...
  )
}

.DataRows <- function(logFile) {
  lines <- readLines(logFile, warn = FALSE)[-1L]
  sum(!grepl("^(#|\\s*$)", lines))
}

test_that("resume reads a run's thin by exact name (#314)", {
  fx <- .RecoveryFixture()
  dir <- withr::local_tempdir()
  mcmc <- .RecoveryJob(dir, treeThin = 10L, checkEvery = 100L)
  set.seed(3141)
  allow_warning(RunMkPrime(fx$pd, fx$tree, mcmc = mcmc),
                "without stabilisation")
  ck <- readRDS(mcmc$checkpointFile)
  # Thin adaptation fired but left thin as it was.
  ck$runs[[1]]$thinAdapted <- TRUE
  ck$runs[[1]][["thin"]] <- NULL
  ck$runs[[1]][["treeThin"]] <- NULL
  saveRDS(ck, mcmc$checkpointFile)

  allow_warning(ResumeMkPrime(mcmc$checkpointFile, fx$pd,
                              mcmc = list(nIter = 600L)),
                "without stabilisation")
  iters <- as.integer(rownames(ReadMkLog(mcmc$logFile)))
  expect_identical(unique(diff(iters)), 5L)
  expect_identical(max(iters), 600L)
})

test_that("a dropped worker keeps its own last checkpoint (#313)", {
  dir <- withr::local_tempdir()
  ckp <- file.path(dir, "job.ckp")
  perRun <- .CkpFilePaths(ckp, 3L)
  launched <- lapply(1:3, function(i) list(run_index = i, saved_idx = 0L))
  returned <- list(list(run_index = 1L, actual_iter = 16000L))
  launchTimes <- as.numeric(Sys.time()) - c(100, 100, 100)
  Write <- function(i, iter, when) {
    saveRDS(list(runs = list(list(run_index = i, actual_iter = iter,
                                  saved_idx = iter %/% 10L)),
                 iter = iter, timestamp = when), perRun[i])
  }
  # Run 2 was killed after checkpointing; run 3's file predates its launch.
  Write(2L, 11000L, Sys.time() - 50)
  Write(3L, 11000L, Sys.time() - 200)

  kept <- .RecoverDroppedRuns(launched, returned, ckp, launchTimes)
  expect_identical(kept[[2]]$actual_iter, 11000L)
  expect_identical(kept[[3]], launched[[3]])
  expect_identical(kept[[1]], launched[[1]])

  # The master, written after the worker's checkpoint, carries it through.
  master <- .WithDroppedRuns(returned, kept)
  expect_identical(master[[2]]$saved_idx, 1100L)
  expect_identical(master[[1]], returned[[1]])

  # A per-run state no further on than the launch state is not taken.
  launched[[2]]$actual_iter <- 12000L
  kept <- .RecoverDroppedRuns(launched, returned, ckp, launchTimes)
  expect_identical(kept[[2]], launched[[2]])
})

# A progress callback that raises an interrupt once `When(info)` holds, so the
# interrupt handlers can be driven without sending a signal.
.InterruptWhen <- function(When) {
  function(info) {
    if (When(info)) {
      signalCondition(structure(class = c("interrupt", "condition"),
                                list(message = "", call = NULL)))
    }
  }
}

test_that("resume to an nIter at or below the checkpoint runs nothing (#318)", {
  fx <- .RecoveryFixture()
  dir <- withr::local_tempdir()
  mcmc <- .RecoveryJob(dir, checkEvery = 100L)
  set.seed(3142)
  allow_warning(RunMkPrime(fx$pd, fx$tree, mcmc = mcmc),
                "without stabilisation")
  ck <- readRDS(mcmc$checkpointFile)
  expect_identical(as.integer(ck$runs[[1]]$actual_iter), 400L)
  nRows <- length(readLines(mcmc$logFile))
  # A checkpointed run keeps its buffered-row count but not the buffer.
  ck$runs[[1]]$flush_idx <- 3L
  saveRDS(ck, mcmc$checkpointFile)

  for (nIter in c(300L, 399L)) {
    res <- ResumeMkPrime(mcmc$checkpointFile, fx$pd,
                         mcmc = list(nIter = nIter))
    expect_identical(as.integer(res[["actual_iter"]]), 400L)
    expect_identical(length(readLines(mcmc$logFile)), nRows)
  }
})

test_that("a run that ends at nIter checkpoints its final rows (#318)", {
  fx <- .RecoveryFixture()
  dir <- withr::local_tempdir()
  # 40 samples: neither a full buffer nor a checkEvery crossing.
  mcmc <- .RecoveryJob(dir, checkEvery = 1000L)
  set.seed(3143)
  allow_warning(RunMkPrime(fx$pd, fx$tree, mcmc = mcmc),
                "without stabilisation")
  ck <- readRDS(mcmc$checkpointFile)
  expect_identical(as.integer(ck$runs[[1]]$actual_iter), 400L)
  expect_identical(ck$runs[[1]]$saved_idx, .DataRows(mcmc$logFile))

  allow_warning(ResumeMkPrime(mcmc$checkpointFile, fx$pd,
                              mcmc = list(nIter = 600L)),
                "without stabilisation")
  iters <- as.integer(rownames(ReadMkLog(mcmc$logFile)))
  expect_identical(iters, seq.int(205L, 600L, by = 5L))
})

test_that("an interrupt flushes before it checkpoints (#317)", {
  fx <- .RecoveryFixture()
  dir <- withr::local_tempdir()
  # The first Sample batch runs 201-5200.
  mcmc <- .RecoveryJob(dir, nIter = 8000L, checkEvery = 10000L,
                       bufferSize = 2000L, plotEvery = 100L,
                       progressFn = .InterruptWhen(function(info) {
                         info$iter >= 5000L
                       }))
  set.seed(3144)
  caught <- new.env()
  withCallingHandlers(
    RunMkPrime(fx$pd, fx$tree, mcmc = mcmc),
    warning = function(w) {
      if (grepl("interrupted", conditionMessage(w))) {
        caught$msg <- conditionMessage(w)
      }
      invokeRestart("muffleWarning")
    }
  )
  ck <- readRDS(mcmc$checkpointFile)
  run <- ck$runs[[1]]
  expect_identical(as.integer(run$actual_iter), 5200L)
  expect_identical(run$flush_idx, 0L)
  expect_identical(run$saved_idx, .DataRows(mcmc$logFile))
  # Comment lines are not samples.
  expect_match(caught$msg, paste0(" ", .DataRows(mcmc$logFile), " samples saved"))
})

test_that("a resume interrupted mid-run checkpoints the run in progress (#317)", {
  fx <- .RecoveryFixture()
  dir <- withr::local_tempdir()
  mcmc <- .RecoveryJob(dir, nRuns = 2L, checkEvery = 10000L,
                       bufferSize = 2000L)
  set.seed(3145)
  allow_warning(RunMkPrime(fx$pd, fx$tree, mcmc = mcmc),
                "without stabilisation")

  # Run 1 finishes at nIter; run 2 is interrupted part way.
  seen <- new.env()
  seen$secondRun <- FALSE
  seen$last <- 0
  Stop <- .InterruptWhen(function(info) {
    if (info$iter < seen$last) seen$secondRun <- TRUE
    seen$last <- info$iter
    seen$secondRun && info$iter >= 5000L
  })
  allow_warning(
    ResumeMkPrime(mcmc$checkpointFile, fx$pd,
                  mcmc = list(nIter = 8000L, plotEvery = 100L,
                              progressFn = Stop)),
    "interrupted|without stabilisation|Resuming"
  )
  ck <- readRDS(mcmc$checkpointFile)
  # Run 2's first batch runs 401-5400.
  expect_identical(as.integer(ck$runs[[1]]$actual_iter), 8000L)
  expect_identical(as.integer(ck$runs[[2]]$actual_iter), 5400L)
  expect_identical(ck$runs[[2]]$saved_idx,
                   .DataRows(sub("\\.log$", "_2.log", mcmc$logFile)))
})
