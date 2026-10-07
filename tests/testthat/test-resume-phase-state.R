# A resume reproduces the phase-transient state an uninterrupted run carries
# between batches (#404).

.PhaseStateData <- function() {
  set.seed(404)
  tree <- ape::rtree(10)
  mat <- replicate(25, as.character(
    ape::rTraitDisc(tree, k = 3, rate = 2, states = c("0", "1", "2"))
  ))
  rownames(mat) <- tree$tip.label
  list(tree = tree, pd = MatrixToPhyDat(mat))
}

# Runs the job until the progress hook reaches `stopAt`, and returns the
# checkpoint it leaves.
.StopAt <- function(data, stopAt, ckpFile, nIter = 4000L) {
  cancelFile <- tempfile()
  on.exit(unlink(cancelFile))
  Stopper <- function(info) if (info$iter >= stopAt) file.create(cancel)
  environment(Stopper) <- list2env(list(cancel = cancelFile, stopAt = stopAt),
                                   parent = baseenv())
  allow_warning(
    RunMkPrime(data$pd, data$tree, mcmc = MkPrimeMCMC(
      nIter = nIter, minWarmup = 2000L, maxWarmup = 2000L, nRuns = 1L,
      nChains = 1L, nCore = 1L, thin = 10L, maxTime = 100,
      checkpointFile = ckpFile, cancelFile = cancelFile,
      progressFn = Stopper, plotEvery = 100L
    )),
    "stabilis|never started sampling"
  )
  readRDS(ckpFile)
}

# Resumes `cp` until its first stop: the end of one C++ call, which a progress
# hook every `plotEvery` iterations cuts short.
.ResumeOneCall <- function(cp, data, ckpFile, plotEvery = NULL) {
  cancelFile <- tempfile()
  on.exit(unlink(cancelFile))
  file.create(cancelFile)
  saveRDS(cp, ckpFile)
  allow_warning(
    ResumeMkPrime(ckpFile, data$pd, mcmc = list(
      cancelFile = cancelFile, plotEvery = plotEvery,
      progressFn = if (!is.null(plotEvery)) function(info) NULL
    )),
    "stabilis|never started sampling"
  )
  readRDS(ckpFile)
}

.WithTimeLimit <- function(env = parent.frame()) {
  setTimeLimit(elapsed = 240, transient = TRUE)
  withr::defer(setTimeLimit(elapsed = Inf, transient = TRUE), envir = env)
}

test_that("a mid-batch warmup resume keeps the partial batch (#404)", {
  skip_on_cran()
  data <- .PhaseStateData()
  ckpFile <- tempfile(fileext = ".ckp")
  on.exit(unlink(ckpFile), add = TRUE)
  .WithTimeLimit()

  set.seed(4041)
  cp <- .StopAt(data, 300L, ckpFile)
  run <- cp$runs[[1]]
  expect_identical(cp$phase, "Warmup")
  expect_equal(run$actual_iter, 300L)
  expect_equal(run$warmupBatchEnd, 500L)
  expect_gt(sum(run$warmupCounts$propose_counts), 0)

  after <- .ResumeOneCall(cp, data, ckpFile)$runs[[1]]
  # The batch ends at 500, as it would have without the resume, not at 800.
  expect_equal(after$actual_iter, 500L)
  expect_length(after$logPostHistory, 1L)
})

test_that("a mid-warmup resume keeps gibbs_kPrime capped (#404)", {
  skip_on_cran()
  data <- .PhaseStateData()
  ckpFile <- tempfile(fileext = ".ckp")
  on.exit(unlink(ckpFile), add = TRUE)
  .WithTimeLimit()

  set.seed(4042)
  cp <- .StopAt(data, 700L, ckpFile)
  run <- cp$runs[[1]]
  capped <- run$moveWeights[["gibbs_kPrime"]]
  expect_lt(capped, run$autoPins[["gibbs_kPrime"]])

  after <- .ResumeOneCall(cp, data, ckpFile, plotEvery = 100L)$runs[[1]]
  expect_equal(after$actual_iter, 800L)
  expect_equal(after$moveWeights[["gibbs_kPrime"]], capped)

  # A checkpoint that predates `autoPins` recovers the full pin.
  cp$runs[[1]]$autoPins <- NULL
  legacy <- .ResumeOneCall(cp, data, ckpFile, plotEvery = 100L)$runs[[1]]
  expect_equal(legacy$autoPins[["gibbs_kPrime"]],
               run$autoPins[["gibbs_kPrime"]])
  expect_equal(legacy$moveWeights[["gibbs_kPrime"]], capped)
})

test_that("resumes shorter than a tuning window still advance Tuning (#404)", {
  skip_on_cran()
  data <- .PhaseStateData()
  ckpFile <- tempfile(fileext = ".ckp")
  on.exit(unlink(ckpFile), add = TRUE)
  .WithTimeLimit()

  set.seed(4043)
  # A tuning window at thin = 10 spans two batches; stop after the first.
  cp <- .StopAt(data, 2500L, ckpFile, nIter = 12000L)
  expect_identical(cp$phase, "Tuning")
  expect_equal(cp$runs[[1]]$tuningIterUsed, 500L)
  expect_equal(cp$runs[[1]]$tuningRound$windowStart, 0L)

  # Each resume runs one batch, half a window.
  for (i in 1:3) cp <- .ResumeOneCall(cp, data, ckpFile)
  run <- cp$runs[[1]]
  expect_equal(run$tuningIterUsed, 2000L)
  expect_equal(run$tuningRound$windowStart, 2000L)
})
