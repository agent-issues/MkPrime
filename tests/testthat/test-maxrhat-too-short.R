# The maxRhat orchestrator's Phase 2 (#420).

RhatData <- function() {
  set.seed(1)
  tree <- ape::rtree(8)
  tree$tip.label <- paste0("t", 1:8)
  mat <- matrix(sample(0:1, 8 * 20, TRUE), 8, 20,
                dimnames = list(tree$tip.label, NULL))
  list(tree = tree, pd = MatrixToPhyDat(mat))
}

test_that("a run too short to tune ends the job rather than each epoch", {
  skip_under_memcheck()
  d <- RhatData()
  dir <- withr::local_tempdir()
  ckp <- file.path(dir, "run.ckp")
  setTimeLimit(elapsed = 180, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  set.seed(2)
  res <- allow_warning(RunMkPrime(
    d$pd, d$tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nCore = 1L, nIter = 1500L,
                       minWarmup = 500L, maxWarmup = 500L, maxRhat = 1.01,
                       maxTime = 20, checkpointFile = ckp)
  ), "Tuning needs|without drawing|maxWarmup")
  expect_identical(res$stop_reason, "too_short")

  # Given room, a resume tunes across Phase 2's shorter epochs and samples.
  res <- allow_warning(
    ResumeMkPrime(ckp, d$pd, d$tree,
                  mcmc = list(nIter = 8000L, maxTime = 60)),
    "maxWarmup"
  )
  expect_false(res$stop_reason == "too_short")
  expect_gt(res$nSamples, 0L)
})

test_that("Phase 2 runs every run to nIter", {
  skip_under_memcheck()
  d <- RhatData()
  dir <- withr::local_tempdir()
  ckp <- file.path(dir, "run.ckp")
  setTimeLimit(elapsed = 240, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  set.seed(3)
  allow_warning(RunMkPrime(
    d$pd, d$tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nCore = 1L, nIter = 600L, thin = 5L,
                       minWarmup = 200L, maxWarmup = 200L, autoTune = FALSE,
                       maxRhat = 1 + 1e-12, maxTime = 60, checkEvery = 200L,
                       checkpointFile = ckp, logFile = file.path(dir, "run.log"))
  ), "maxWarmup")

  # Interrupt run 2 part way through the second epoch, so that it lags run 1.
  Interrupter <- function(info) {
    if (info$iter == 1700L) {
      seen$n <- seen$n + 1L
      if (seen$n == 2L) rlang::interrupt()
    }
  }
  environment(Interrupter) <- list2env(
    list(seen = list2env(list(n = 0L))), parent = baseenv()
  )
  allow_warning(
    ResumeMkPrime(ckp, d$pd, d$tree,
                  mcmc = list(nIter = 6000L, plotEvery = 100L,
                              progressFn = Interrupter)),
    "interrupted"
  )
  lagged <- vapply(readRDS(ckp)$runs, `[[`, 0, "actual_iter")
  expect_gt(lagged[[1]], lagged[[2]])

  nIter <- lagged[[1]] + 1000L
  res <- ResumeMkPrime(ckp, d$pd, d$tree,
                       mcmc = list(nIter = nIter, progressFn = NULL))
  expect_identical(res$stop_reason, "max_iter")
  expect_equal(vapply(readRDS(ckp)$runs, `[[`, 0, "actual_iter"),
               c(nIter, nIter))
})
