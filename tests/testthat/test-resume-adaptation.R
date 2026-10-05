# Resume rebuilds the adaptation state that lives only for a phase (#302), and
# tuning is priced against the right target and spend (#201).

.ResumeAdaptData <- function() {
  set.seed(302)
  tree <- ape::rtree(10)
  mat <- replicate(25, as.character(
    ape::rTraitDisc(tree, k = 3, rate = 2, states = c("0", "1", "2"))
  ))
  rownames(mat) <- tree$tip.label
  list(tree = tree, pd = MatrixToPhyDat(mat))
}

.ResumeAdaptMcmc <- function(...) {
  do.call(MkPrimeMCMC, utils::modifyList(list(
    nIter = 4000L, minWarmup = 2000L, maxWarmup = 2000L, nRuns = 1L,
    nChains = 1L, nCore = 1L, thin = 10L, maxTime = 100
  ), list(...)))
}

# Runs the job, cancelling it once the progress hook reaches `stopAt`, and
# returns the checkpoint with the hook removed.
.InterruptAt <- function(data, stopAt, ckpFile, ...) {
  cancelFile <- tempfile()
  on.exit(unlink(cancelFile))
  Stopper <- function(info) if (info$iter >= stopAt) file.create(cancel)
  environment(Stopper) <- list2env(list(cancel = cancelFile, stopAt = stopAt),
                                   parent = baseenv())
  allow_warning(
    RunMkPrime(data$pd, data$tree, mcmc = .ResumeAdaptMcmc(
      checkpointFile = ckpFile, cancelFile = cancelFile,
      progressFn = Stopper, plotEvery = 500L, ...
    )),
    "stabilis|never started sampling"
  )
  cp <- readRDS(ckpFile)
  cp$mcmc$progressFn <- NULL
  cp$mcmc$plotEvery <- NULL
  cp$mcmc$cancelFile <- NULL
  cp
}

# Resumes `cp` for one batch: a cancel file already in place stops the run as
# soon as that batch is checkpointed.
.ResumeOneBatch <- function(cp, data, ckpFile) {
  cancelFile <- tempfile()
  on.exit(unlink(cancelFile))
  file.create(cancelFile)
  cp$mcmc$cancelFile <- cancelFile
  saveRDS(cp, ckpFile)
  allow_warning(ResumeMkPrime(ckpFile, data$pd), "stabilis")
  readRDS(ckpFile)
}

test_that(".SchedulePins keeps stored auto-pins over capped weights", {
  weights <- c(nni = 0.4, gibbs_kPrime = 0.1, slice_tl = 0.1, spr = 0.4)
  types <- c("nni", "gibbs_kprime_sweep", "slice", "spr")
  capped <- MkPrime:::.WarmupGibbsCap(
    weights, MkPrime:::.SchedulePins(weights, types, NULL), 2L
  )
  base <- MkPrime:::.SchedulePins(weights, types, NULL)
  expect_equal(MkPrime:::.SchedulePins(capped, types, NULL, basePins = base),
               base)
  expect_lt(MkPrime:::.SchedulePins(capped, types, NULL)[["gibbs_kPrime"]],
            base[["gibbs_kPrime"]])
  # A user pin still overrides the stored base.
  expect_equal(
    MkPrime:::.SchedulePins(capped, types, c(gibbs_kPrime = 0.05),
                            basePins = base)[["gibbs_kPrime"]],
    0.05
  )
})

test_that("a mid-warmup resume freezes the same gibbs_kPrime weight (#302)", {
  skip_under_memcheck()
  skip_on_cran()
  data <- .ResumeAdaptData()
  ckpFile <- tempfile(fileext = ".ckp")
  on.exit(unlink(ckpFile), add = TRUE)
  setTimeLimit(elapsed = 240, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)

  set.seed(3021)
  full <- allow_warning(
    RunMkPrime(data$pd, data$tree, mcmc = .ResumeAdaptMcmc()), "stabilis"
  )
  cp <- .InterruptAt(data, 1000L, ckpFile)
  expect_identical(cp$phase, "Warmup")
  # Two warmup batches, two rho snapshots, checkpointed with the run.
  expect_equal(nrow(cp$runs[[1]]$rhoSampleBuf), 2L)

  saveRDS(cp, ckpFile)
  resumed <- allow_warning(ResumeMkPrime(ckpFile, data$pd), "stabilis")
  # Before #302 each resume cut the pin by the warmup cap: 0.0776 v. 0.2329.
  expect_equal(resumed$moveWeights[["gibbs_kPrime"]],
               full$moveWeights[["gibbs_kPrime"]])
})

test_that("a mid-warmup resume keeps adapted rhos and the rho buffer (#302)", {
  skip_on_cran()
  data <- .ResumeAdaptData()
  ckpFile <- tempfile(fileext = ".ckp")
  on.exit(unlink(ckpFile), add = TRUE)
  setTimeLimit(elapsed = 240, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)

  set.seed(3022)
  cp <- .InterruptAt(data, 1000L, ckpFile)
  adapted <- list(rho_tl_rls = 0.4, rho_tl_rl = 0, rho_tl_rn = 0)
  cp$runs[[1]]$chain_rhos <- list(adapted)
  after <- .ResumeOneBatch(cp, data, ckpFile)

  expect_equal(after$iter, 1500)
  expect_equal(nrow(after$runs[[1]]$rhoSampleBuf), 3L)
  expect_equal(after$runs[[1]]$chain_rhos[[1]], adapted)
})

test_that("a mid-Tuning resume does not charge the interrupted window (#302)", {
  skip_on_cran()
  data <- .ResumeAdaptData()
  ckpFile <- tempfile(fileext = ".ckp")
  on.exit(unlink(ckpFile), add = TRUE)
  setTimeLimit(elapsed = 240, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)

  set.seed(3023)
  # A tuning window at thin = 10 spans two batches; stop after the first.
  cp <- .InterruptAt(data, 2500L, ckpFile, nIter = 12000L)
  expect_identical(cp$phase, "Tuning")
  expect_equal(cp$runs[[1]]$tuningIterUsed, 500L)
  after <- .ResumeOneBatch(cp, data, ckpFile)
  # The window restarts, so only the batch that re-runs it is charged.
  expect_equal(after$runs[[1]]$tuningIterUsed, 500L)
})

test_that("tuning is priced on the full minEss and all spend (#201)", {
  skip_under_memcheck()
  skip_on_cran()
  data <- .ResumeAdaptData()
  ckpFile <- tempfile(fileext = ".ckp")
  on.exit(unlink(ckpFile), add = TRUE)
  setTimeLimit(elapsed = 240, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)

  seen <- new.env()
  seen$nRuns <- integer(0)
  seen$spent <- numeric(0)
  payback <- MkPrime:::.TuningPayback
  local_mocked_bindings(
    .TuningPayback = function(streak, spentSec, bestRate, targetEss,
                              nRuns = 1L, frozen = FALSE) {
      seen$nRuns <- c(seen$nRuns, nRuns)
      seen$spent <- c(seen$spent, spentSec)
      payback(streak, spentSec, bestRate, targetEss, nRuns, frozen)
    }
  )
  set.seed(2011)
  cp <- .InterruptAt(data, 2500L, ckpFile, nIter = 12000L, nRuns = 2L,
                     tuningBudget = 4000L, minEss = 1e6)
  spentBefore <- cp$runs[[1]]$adaptSec
  expect_gt(spentBefore, 0)
  saveRDS(cp, ckpFile)
  seen$nRuns <- integer(0)
  seen$spent <- numeric(0)
  allow_warning(ResumeMkPrime(ckpFile, data$pd), "stabilis")

  # Each serial run stops on its own minEss, so none prices a per-run share.
  expect_gt(length(seen$nRuns), 0L)
  expect_true(all(seen$nRuns == 1L))
  # The resumed run is charged for what it spent before the interruption.
  expect_gt(seen$spent[[1]], spentBefore)
})

test_that(".SampleElapsed averages sampling time over runs still sampling", {
  # Run 1 sampling since 10, run 2 sampled from 20 to 50, run 3 not yet.
  expect_equal(MkPrime:::.SampleElapsed(c(10, 20, NA), c(NA, 50, NA), 100),
               (90 + 30) / 2)
  expect_equal(MkPrime:::.SampleElapsed(c(NA, NA), c(NA, NA), 100), 0)
})
