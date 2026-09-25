# What each tuning window is scored on, seen through the calls to .MinEssRate()

.windowRecord <- new.env()
.tuningTree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
.tuningData <- MatrixToPhyDat(matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                                     dimnames = list(paste0("t", 1:4), NULL)))

# Tuning on a 4-tip tree, 20000 iterations after a 500-iteration warmup, so a
# budget of 9750. The progress callback sleeps once per 2000-iteration window
# during tuning, so that time charged to a window shows in its `sec`. The pause
# dwarfs a window's own sampling time (~0.1 s on Linux; more on Windows CI).
.tuningPause <- 1

.TuningWindows <- function(pause = .tuningPause) {
  if (!is.null(.windowRecord$windows)) return(.windowRecord$windows)
  windows <- new.env()
  windows$list <- list()
  RealRate <- MkPrime:::.MinEssRate
  local_mocked_bindings(
    .MinEssRate = function(sampleMatrix, wallTimeSec, tuningTrees = NULL,
                           fixedCols = NULL) {
      windows$list[[length(windows$list) + 1L]] <- list(
        n = nrow(sampleMatrix), sec = wallTimeSec,
        trees = !is.null(tuningTrees)
      )
      RealRate(sampleMatrix, wallTimeSec, tuningTrees, fixedCols)
    },
    .package = "MkPrime"
  )
  Pause <- function(info) if (identical(info$phase, "Tuning")) Sys.sleep(pause)
  set.seed(7820)
  allow_warning(
    RunMkPrime(.tuningData, .tuningTree,
               mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 20000L, thin = 17L,
                                  maxWarmup = 500L, minWarmup = 500L,
                                  autoTune = TRUE, progressFn = Pause,
                                  plotEvery = 2000L, maxTime = 60)),
    "maxWarmup"
  )
  .windowRecord$windows <- windows$list
  .windowRecord$windows
}

test_that("a tuning window holds ~100 samples, not one batch (#78)", {
  windows <- .TuningWindows()
  expect_gte(length(windows), 4L)
  expect_true(all(vapply(windows, `[[`, integer(1), "n") >= 100L))
})

test_that("tuning starts no round that would overrun its budget", {
  windows <- .TuningWindows()
  # Samples per window times thin: the iterations tuning spent.
  expect_lte(sum(vapply(windows, `[[`, integer(1), "n")) * 17L, 9750L)
})

test_that("a tuning window is not charged for the progress callback (#220)", {
  windows <- .TuningWindows()
  expect_lt(max(vapply(windows, `[[`, numeric(1), "sec")), .tuningPause)
})

test_that("tree ESS enters the tuning score without minTreeEss (#222)", {
  skip_if_not_installed("TreeDist")
  windows <- .TuningWindows()
  expect_true(all(vapply(windows, `[[`, logical(1), "trees")))
})

test_that("a mid-round resume keeps the round's best so far (#221)", {
  cpFile <- tempfile(fileext = ".ckp")
  snapFile <- tempfile(fileext = ".ckp")
  on.exit(unlink(c(cpFile, snapFile)), add = TRUE)

  # Keep the last checkpoint written during tuning: mid-way through the last
  # candidate's window, with the incumbent and earlier candidates scored.
  RealSave <- MkPrime:::.SaveCheckpoint
  local_mocked_bindings(
    .SaveCheckpoint = function(runs, mcmc, iter, paramNames, file, ...) {
      RealSave(runs, mcmc, iter, paramNames, file, ...)
      if (identical(runs[[1]]$phase, "Tuning")) {
        file.copy(file, snapFile, overwrite = TRUE)
      }
    },
    .package = "MkPrime"
  )
  set.seed(2211)
  allow_warning(
    RunMkPrime(.tuningData, .tuningTree,
               mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 20000L, thin = 17L,
                                  maxWarmup = 500L, minWarmup = 500L,
                                  autoTune = TRUE, tuningRounds = 1L,
                                  checkEvery = 500L, checkpointFile = cpFile,
                                  maxTime = 60)),
    "maxWarmup"
  )
  expect_true(file.exists(snapFile))
  saved <- readRDS(snapFile)$runs[[1]]$tuningRound
  expect_gt(saved$candIdx, 0L)
  expect_true(is.finite(saved$bestRate))

  seen <- new.env()
  local_mocked_bindings(
    .BeatsIncumbent = function(candRate, candEss, bestRate, bestEss, ...) {
      if (is.null(seen$bestRate)) seen$bestRate <- bestRate
      FALSE
    },
    .package = "MkPrime"
  )
  ResumeMkPrime(snapFile, .tuningData, .tuningTree)
  # The candidate on trial is judged against the best found before the
  # checkpoint, not against nothing.
  expect_identical(seen$bestRate, saved$bestRate)
})
