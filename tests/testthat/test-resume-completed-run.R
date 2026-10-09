# Extending a finished run that kept only a checkpoint (#423).
skip_under_memcheck()

DoneData <- function() {
  set.seed(1)
  tree <- ape::rtree(8)
  tree$tip.label <- paste0("t", 1:8)
  mat <- matrix(sample(0:1, 8 * 20, TRUE), 8, 20,
                dimnames = list(tree$tip.label, NULL))
  list(tree = tree, pd = MatrixToPhyDat(mat))
}

test_that("a resume extends the log beside a finished run's checkpoint", {
  d <- DoneData()
  dir <- withr::local_tempdir()
  ckp <- file.path(dir, "job.ckp")
  set.seed(2)
  r1 <- allow_warning(RunMkPrime(
    d$pd, d$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000L, thin = 10L,
                       minWarmup = 500L, maxWarmup = 500L, autoTune = FALSE,
                       checkpointFile = ckp, maxTime = 60)
  ), "maxWarmup")
  expect_equal(nrow(r1$samples), 50L)
  expect_true(file.exists(file.path(dir, "job_mkp_run.log")))

  r2 <- ResumeMkPrime(ckp, d$pd, d$tree, mcmc = list(nIter = 1500L))
  expect_equal(r2$nSamples, 100L)
  expect_equal(nrow(r2$samples), 100L)
  expect_equal(r2$samples[seq_len(50), ], r1$samples, ignore_attr = TRUE)
  expect_length(r2$trees, 100L)
  expect_null(r2$logFile)
})

test_that("a resume whose log is lost says so", {
  d <- DoneData()
  dir <- withr::local_tempdir()
  ckp <- file.path(dir, "job.ckp")
  set.seed(3)
  allow_warning(RunMkPrime(
    d$pd, d$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000L, thin = 10L,
                       minWarmup = 500L, maxWarmup = 500L, autoTune = FALSE,
                       checkpointFile = ckp, maxTime = 60)
  ), "maxWarmup")
  unlink(file.path(dir, "job_mkp_run.log"))

  expect_warning(
    r2 <- ResumeMkPrime(ckp, d$pd, d$tree, mcmc = list(nIter = 1500L)),
    "previous samples unavailable"
  )
  expect_equal(nrow(r2$samples), r2$nSamples)
})
