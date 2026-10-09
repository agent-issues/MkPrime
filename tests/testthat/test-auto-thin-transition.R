# Auto-thin adapts at the first check with enough Sample rows (#422).

ThinData <- function() {
  set.seed(1)
  tree <- ape::rtree(8)
  tree$tip.label <- paste0("t", 1:8)
  mat <- matrix(sample(0:1, 8 * 20, TRUE), 8, 20,
                dimnames = list(tree$tip.label, NULL))
  list(tree = tree, pd = MatrixToPhyDat(mat))
}

LoggedIters <- function(logFile) {
  rows <- grep("^[0-9]", readLines(logFile), value = TRUE)
  as.integer(vapply(strsplit(rows, "\t"), `[`, "", 1L))
}

# Whatever the chain's autocorrelation, an adapting check doubles thin.
DoubleThin <- function(env = parent.frame()) {
  local_mocked_bindings(
    .AdaptThinning = function(sampleMatrix, currentThin, nMoves) {
      2L * currentThin
    },
    .package = "MkPrime", .env = env
  )
}

ThinMcmc <- function(dir, nIter, maxWarmup, ...) {
  MkPrimeMCMC(nRuns = 1L, nIter = nIter, minWarmup = maxWarmup,
              maxWarmup = maxWarmup, autoTune = FALSE, checkEvery = 1000L,
              logFile = file.path(dir, "run.log"),
              checkpointFile = file.path(dir, "run.ckp"), bufferSize = 1000L,
              maxTime = 60, ...)
}

test_that("thin adapts when Warmup ends on a checkEvery multiple", {
  skip_under_memcheck()
  d <- ThinData()
  dir <- withr::local_tempdir()
  DoubleThin()
  set.seed(2)
  res <- allow_warning(
    RunMkPrime(d$pd, d$tree, mcmc = ThinMcmc(dir, 3000L, 1000L)),
    "maxWarmup"
  )
  startThin <- diff(LoggedIters(file.path(dir, "run.log")))[[1]]
  expect_equal(res$thin, 2L * startThin)
})

test_that("the check that adapts thin checkpoints the new value", {
  skip_under_memcheck()
  d <- ThinData()
  dir <- withr::local_tempdir()
  DoubleThin()
  set.seed(3)
  # Fails the job, as a kill would, after the check at 2000 adapts thin.
  Kill <- function(info) if (info$iter == 2100L) stop("killed")
  expect_error(allow_warning(
    RunMkPrime(d$pd, d$tree,
               mcmc = ThinMcmc(dir, 3000L, 1000L, plotEvery = 100L,
                               progressFn = Kill)),
    "maxWarmup"
  ), "killed")
  ckp <- readRDS(file.path(dir, "run.ckp"))
  startThin <- diff(LoggedIters(file.path(dir, "run.log")))[[1]]
  expect_equal(ckp$iter, 2000)
  expect_equal(ckp$runs[[1]][["thin"]], 2L * startThin)
})

test_that("a resume adapts thin for a run that has not yet adapted it", {
  skip_under_memcheck()
  d <- ThinData()
  dir <- withr::local_tempdir()
  set.seed(4)
  # Stops in Sample before its first check.
  allow_warning(
    RunMkPrime(d$pd, d$tree, mcmc = ThinMcmc(dir, 900L, 500L)),
    "maxWarmup"
  )
  startThin <- diff(LoggedIters(file.path(dir, "run.log")))[[1]]

  DoubleThin()
  res <- ResumeMkPrime(file.path(dir, "run.ckp"), d$pd, d$tree,
                       mcmc = list(nIter = 3000L))
  expect_equal(res$thin, 2L * startThin)
})
