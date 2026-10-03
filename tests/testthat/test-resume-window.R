# A resumed run carries on its convergence window and skips a run that has
# already converged (#316).

.WindowFixture <- function() {
  list(
    tree = ape::read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);"),
    pd = MatrixToPhyDat(matrix(
      c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L), 4, 3,
      dimnames = list(paste0("t", 1:4), NULL)
    ))
  )
}

.WindowJob <- function(dir, nIter, ...) {
  MkPrimeMCMC(nRuns = 1L, nIter = nIter, thin = 10L, minWarmup = 100L,
              maxWarmup = 100L, autoTune = FALSE, maxTime = 60,
              checkpointFile = file.path(dir, "job.ckp"),
              logFile = file.path(dir, "job.log"), ...)
}

# The window a resume to `nIter` holds, as its progress callback sees it.
.ResumedWindow <- function(dir, pd, nIter, ...) {
  seen <- new.env()
  Capture <- function(info) seen$rows <- info$runSamples[[1]]
  ResumeMkPrime(file.path(dir, "job.ckp"), pd,
                mcmc = list(nIter = nIter, plotEvery = nIter,
                            progressFn = Capture, ...))
  seen$rows
}

test_that("a log's tail seeds the window of a run resumed from it", {
  params <- c("log_posterior", "tree_length")
  log <- withr::local_tempfile(fileext = ".log")
  rows <- matrix(as.numeric(seq_len(60)), 30, 2,
                 dimnames = list(seq(10, 300, 10), params))
  writeLines(c(paste(c("Sample", params), collapse = "\t"),
               paste(rownames(rows), rows[, 1], rows[, 2], sep = "\t")), log)
  r <- list(saved_idx = 30L)

  tail20 <- .PriorWindowRows(r, log, params, 20L)
  expect_equal(unname(tail20), unname(rows[11:30, ]))
  full <- .SeedConvWindow(c(r, .InitStreamBuffers(2L, params, 10L, 20L)),
                          tail20)
  expect_true(full$conv_filled)
  expect_equal(unname(.ConvWindowRows(full)), unname(rows[11:30, ]))
  # The oldest row is the first to go.
  full <- .AddToStreamBuffer(full, c(-1, -2), 310L, log, 10L)
  expect_equal(unname(.ConvWindowRows(full)),
               unname(rbind(rows[12:30, ], c(-1, -2))))

  part <- .SeedConvWindow(c(r, .InitStreamBuffers(2L, params, 10L, 50L)),
                          .PriorWindowRows(r, log, params, 50L))
  expect_false(part$conv_filled)
  expect_identical(part$conv_head, 30L)
  expect_equal(unname(.ConvWindowRows(part)), unname(rows))

  # A run carried on in memory keeps its own window.
  expect_identical(.PriorWindowRows(part, NULL, params, 50L),
                   .ConvWindowRows(part))
  # Nothing to seed from a fresh run, or a log of other columns.
  expect_null(.PriorWindowRows(list(saved_idx = 0L), log, params, 20L))
  expect_null(.PriorWindowRows(r, log, rev(params), 20L))
})

test_that("a resume keeps the samples already in the window", {
  fx <- .WindowFixture()
  dir <- withr::local_tempdir()
  set.seed(5772)
  allow_warning(RunMkPrime(fx$pd, fx$tree, fixTopology = TRUE,
                           mcmc = .WindowJob(dir, 600L)),
                "without stabilisation")
  before <- ReadMkLog(file.path(dir, "job.log"))
  expect_identical(nrow(before), 50L)

  window <- .ResumedWindow(dir, fx$pd, 700L)
  expect_identical(nrow(window), 60L)
  expect_equal(unname(window[1:50, ]), unname(before))
})

test_that("a resume keeps a window that grew", {
  fx <- .WindowFixture()
  dir <- withr::local_tempdir()
  # A 40-row window, short of an unreachable minEss, grows at a check.
  job <- .WindowJob(dir, 1000L, bufferSize = 20L, checkEvery = 100L,
                    minEss = 1e6)
  set.seed(1414)
  allow_warning(RunMkPrime(fx$pd, fx$tree, fixTopology = TRUE, mcmc = job),
                "without stabilisation|minEss")
  expect_gt(readRDS(job$checkpointFile)$runs[[1]]$conv_size, 40L)

  window <- .ResumedWindow(dir, fx$pd, 1100L)
  expect_identical(nrow(window), 100L)
})

test_that("a converged run is not run again", {
  fx <- .WindowFixture()
  dir <- withr::local_tempdir()
  ckp <- file.path(dir, "job.ckp")
  set.seed(2236)
  allow_warning(RunMkPrime(fx$pd, fx$tree, fixTopology = TRUE,
                           mcmc = .WindowJob(dir, 300L)),
                "without stabilisation")
  ck <- readRDS(ckp)
  ck$runs[[1]]$stop_reason <- "converged"
  ck$runs[[1]]$conv_criteria <- list(minEss = 50, minTreeEss = NULL)
  saveRDS(ck, ckp)

  same <- ResumeMkPrime(ckp, fx$pd, mcmc = list(nIter = 600L, minEss = 50))
  expect_identical(same$stop_reason, "converged")
  expect_equal(same$actual_iter, 300)
  expect_identical(nrow(ReadMkLog(file.path(dir, "job.log"))), 20L)

  # A stricter minEss is a reason to run on.
  raised <- ResumeMkPrime(ckp, fx$pd, mcmc = list(nIter = 400L, minEss = 1e6))
  expect_identical(raised$stop_reason, "max_iter")
  expect_equal(raised$actual_iter, 400)
})

test_that(".StillConverged() compares the criteria a run converged under", {
  r <- list(stop_reason = "converged",
            conv_criteria = list(minEss = 100, minTreeEss = NULL))
  expect_true(.StillConverged(r, list(minEss = 100)))
  expect_true(.StillConverged(r, list(minEss = 50)))
  expect_false(.StillConverged(r, list()))
  expect_false(.StillConverged(r, list(minEss = 200)))
  expect_false(.StillConverged(r, list(minEss = 100, minTreeEss = 10)))
  both <- list(stop_reason = "converged",
               conv_criteria = list(minEss = 100, minTreeEss = 50))
  expect_true(.StillConverged(both, list(minTreeEss = 50)))
  expect_false(.StillConverged(list(stop_reason = "converged"),
                               list(minEss = 100)))
  expect_false(.StillConverged(list(stop_reason = "max_iter",
                                    conv_criteria = list()), list()))
})

test_that("serial runs skip a converged run and reset streaks per epoch", {
  calls <- new.env()
  calls$log <- list()
  local_mocked_bindings(
    .RunMkPrimeSingleRun = function(mkd, model, mcmc, initialState, moves,
                                    tipLabels, runIdx, ..., startIter = 1L) {
      calls$log <- c(calls$log, list(list(
        run = runIdx, streak = initialState$convStreak %||% 0L
      )))
      initialState$actual_iter <- min(mcmc$nIter, startIter + 999L)
      initialState$stop_reason <- "max_iter"
      initialState
    },
    .CheckConvergenceFromLogs = function(...) {
      list(converged = FALSE, maxRhat = 2)
    }
  )
  mcmc <- list(nIter = 3000L, checkEvery = 1000L, maxRhat = 1.01,
               minEss = 50)
  runs <- list(
    list(stop_reason = "converged", actual_iter = 1000L, convStreak = 2L,
         conv_criteria = list(minEss = 50)),
    list(actual_iter = 0L, convStreak = 0L)
  )
  out <- .RunSerialRuns(NULL, NULL, mcmc, runs, NULL, NULL, NULL, 0L, 0L,
                        logFilePaths = c("a", "b"), convWindowSize = 0L,
                        treeFilePaths = NULL, startIters = c(1001L, 1L))
  ran <- vapply(calls$log, `[[`, 0L, "run")
  streaks <- vapply(calls$log, `[[`, 0L, "streak")
  # Phase 1 runs run 2 alone; each epoch then runs both from a clean streak.
  expect_identical(ran, c(2L, 1L, 2L, 1L, 2L))
  expect_identical(streaks[-1], c(0L, 0L, 0L, 0L))
  expect_identical(out$stopReason, "max_iter")
})
