# What each tuning window is scored on, seen through the calls to .MinEssRate()

.windowRecord <- new.env()

# One tuning round on a 4-tip tree. The progress callback sleeps during
# tuning, so that time charged to a window shows in its `sec`.
.TuningWindows <- function(pause = 0.2) {
  if (!is.null(.windowRecord$windows)) return(.windowRecord$windows)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
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
    RunMkPrime(MatrixToPhyDat(mat), tree,
               mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 20000L, thin = 17L,
                                  maxWarmup = 500L, minWarmup = 500L,
                                  autoTune = TRUE, tuningRounds = 1L,
                                  progressFn = Pause, plotEvery = 1000L,
                                  maxTime = 60)),
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

test_that("a tuning window is not charged for the progress callback (#220)", {
  windows <- .TuningWindows()
  expect_true(all(vapply(windows, `[[`, numeric(1), "sec") < 0.2))
})

test_that("tree ESS enters the tuning score without minTreeEss (#222)", {
  skip_if_not_installed("TreeDist")
  windows <- .TuningWindows()
  expect_true(all(vapply(windows, `[[`, logical(1), "trees")))
})

test_that("a mid-round resume restarts the round from the best schedule (#221)", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)
  cpFile <- tempfile(fileext = ".ckp")
  cpFile2 <- tempfile(fileext = ".ckp")
  on.exit(unlink(c(cpFile, cpFile2)), add = TRUE)

  set.seed(2211)
  allow_warning(
    RunMkPrime(pd, tree,
               mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 3000L, thin = 5L,
                                  maxWarmup = 500L, minWarmup = 500L,
                                  autoTune = FALSE, checkEvery = 300L,
                                  checkpointFile = cpFile, maxTime = 60)),
    "maxWarmup"
  )
  cp <- readRDS(cpFile)
  r <- cp$runs[[1]]
  best <- r$moveWeights
  # The checkpoint caught a candidate on trial, not yet compared.
  onTrial <- MkPrime:::.PerturbMoveWeights(best, NULL, names(best))[[1]]
  r$phase <- "Tuning"
  r$tuningIterUsed <- 0L
  r$tuningRoundsDone <- 0L
  r$moveWeights <- onTrial
  r$tuningBestWeights <- best
  cp$runs[[1]] <- r
  cp$moveWeights <- onTrial
  cp$iter <- 500L
  saveRDS(cp, cpFile2)

  seen <- new.env()
  local_mocked_bindings(
    .PerturbMoveWeights = function(currentWeights, ...) {
      if (is.null(seen$first)) seen$first <- currentWeights
      list()
    },
    .package = "MkPrime"
  )
  ResumeMkPrime(cpFile2, pd, tree)
  expect_equal(seen$first, best)
})
