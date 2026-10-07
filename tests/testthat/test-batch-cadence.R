# Callbacks, checks and the cancel file are polled between C++ calls, so each
# call must end on a plotEvery or checkEvery boundary (#327).

.CadenceFixture <- function() {
  list(
    tree = ape::read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);"),
    pd = MatrixToPhyDat(matrix(
      c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L), 4, 3,
      dimnames = list(paste0("t", 1:4), NULL)
    ))
  )
}

test_that(".NextMultiple() finds the nearest boundary of any period", {
  expect_equal(.NextMultiple(1L, 1000L), 1000)
  expect_equal(.NextMultiple(1000L, 1000L), 1000)
  expect_equal(.NextMultiple(1001L, c(1000L, 300L)), 1200)
  expect_equal(.NextMultiple(5L, integer(0)), Inf)
})

test_that("the progress callback fires every plotEvery, warmup included", {
  fx <- .CadenceFixture()
  seen <- new.env()
  seen$iter <- integer(0)
  seen$inWarmup <- logical(0)
  Track <- function(info) {
    seen$iter <- c(seen$iter, info$iter)
    seen$inWarmup <- c(seen$inWarmup, info$inWarmup)
  }
  mcmc <- MkPrimeMCMC(nRuns = 1L, nIter = 1000L, thin = 10L,
                      minWarmup = 200L, maxWarmup = 200L, autoTune = FALSE,
                      maxTime = 60, plotEvery = 100L, progressFn = Track)
  set.seed(2718)
  allow_warning(RunMkPrime(fx$pd, fx$tree, mcmc = mcmc, fixTopology = TRUE),
                "without stabilisation")
  expect_equal(seen$iter, seq(100, 1000, by = 100))
  expect_equal(seen$inWarmup, seen$iter <= 200)
})

test_that("a cancel file is seen at the first checkEvery boundary", {
  fx <- .CadenceFixture()
  dir <- withr::local_tempdir()
  cancel <- file.path(dir, "stop")
  file.create(cancel)
  Job <- function(...) {
    MkPrimeMCMC(nRuns = 1L, thin = 10L, minWarmup = 1000L,
                maxWarmup = 1000L, autoTune = FALSE, maxTime = 60,
                checkEvery = 300L,
                checkpointFile = file.path(dir, "job.ckp"),
                logFile = file.path(dir, "job.log"), ...)
  }

  set.seed(1618)
  warm <- RunMkPrime(fx$pd, fx$tree, fixTopology = TRUE,
                     mcmc = Job(nIter = 2000L, cancelFile = cancel))
  expect_identical(warm$stop_reason, "cancelled")
  expect_equal(warm$actual_iter, 300)

  # In the Sample phase, too, rather than at the 5000-iteration batch end.
  unlink(cancel)
  set.seed(1618)
  allow_warning(
    RunMkPrime(fx$pd, fx$tree, fixTopology = TRUE, overwrite = TRUE,
               mcmc = Job(nIter = 1000L)),
    "without stabilisation"
  )
  file.create(cancel)
  resumed <- ResumeMkPrime(file.path(dir, "job.ckp"), fx$pd,
                           mcmc = list(nIter = 3000L, cancelFile = cancel))
  expect_identical(resumed$stop_reason, "cancelled")
  expect_equal(resumed$actual_iter, 1200)
})
