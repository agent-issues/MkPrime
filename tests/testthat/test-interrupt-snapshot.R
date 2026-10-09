# An interrupt must not write a batch to the log twice (#419).

SnapData <- function() {
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

# Ctrl-C while the first convergence check runs: after it has flushed the
# buffer, before the run's next batch.
InterruptAtCheck <- function(env = parent.frame()) {
  local_mocked_bindings(
    .CheckConvergence = function(...) rlang::interrupt(),
    .package = "MkPrime", .env = env
  )
}

SnapMcmc <- function(dir, nIter) {
  MkPrimeMCMC(
    nRuns = 1L, nIter = nIter, thin = 5L, minWarmup = 500L,
    maxWarmup = 500L, autoTune = FALSE, checkEvery = 1000L,
    bufferSize = 1000L, logFile = file.path(dir, "run.log"),
    checkpointFile = file.path(dir, "run.ckp"), maxTime = 60
  )
}

test_that("an interrupt after a check's flush logs each sample once", {
  d <- SnapData()
  dir <- withr::local_tempdir()
  InterruptAtCheck()
  set.seed(2)
  expect_warning(allow_warning(
    RunMkPrime(d$pd, d$tree, mcmc = SnapMcmc(dir, 3000L)), "maxWarmup"
  ), "interrupted")
  iters <- LoggedIters(file.path(dir, "run.log"))
  expect_gt(length(iters), 0L)
  expect_false(anyDuplicated(iters) > 0L)
  expect_true(all(diff(iters) > 0L))
})

test_that("an interrupted resume logs each sample once", {
  d <- SnapData()
  dir <- withr::local_tempdir()
  set.seed(3)
  allow_warning(RunMkPrime(d$pd, d$tree, mcmc = SnapMcmc(dir, 1500L)),
                "maxWarmup")
  before <- LoggedIters(file.path(dir, "run.log"))

  InterruptAtCheck()
  expect_warning(
    ResumeMkPrime(file.path(dir, "run.ckp"), d$pd, d$tree,
                  mcmc = list(nIter = 3000L)),
    "interrupted"
  )
  iters <- LoggedIters(file.path(dir, "run.log"))
  expect_gt(length(iters), length(before))
  expect_false(anyDuplicated(iters) > 0L)
  expect_true(all(diff(iters) > 0L))
})
