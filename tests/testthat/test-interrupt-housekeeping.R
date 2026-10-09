# Interrupt, cancel and resume bookkeeping (#424, #425).

.HousekeepingData <- function() {
  set.seed(424)
  tree <- ape::rtree(8)
  tree$tip.label <- paste0("t", 1:8)
  mat <- matrix(sample(0:1, 8 * 20, TRUE), 8, 20,
                dimnames = list(tree$tip.label, NULL))
  list(tree = tree, pd = MatrixToPhyDat(mat))
}

.HousekeepingMcmc <- function(...) {
  do.call(MkPrimeMCMC, utils::modifyList(list(
    nRuns = 1L, nIter = 3000L, thin = 5L, minWarmup = 500L, maxWarmup = 500L,
    autoTune = FALSE, checkEvery = 1000L, maxTime = 60
  ), list(...)))
}

# Ctrl-C while the first convergence check runs.
.InterruptAtCheck <- function(env = parent.frame()) {
  local_mocked_bindings(
    .CheckConvergence = function(...) rlang::interrupt(),
    .env = env
  )
}

test_that("an interrupted temporary run points at its checkpoint (#424)", {
  skip_on_cran()
  d <- .HousekeepingData()
  mkpEnv <- environment(RunMkPrime)$.mkp_env
  on.exit(MkPrime:::.CleanupStaleTempLogs(), add = TRUE)
  .InterruptAtCheck()
  set.seed(4241)
  msg <- paste(capture_warnings(
    allow_warning(RunMkPrime(d$pd, d$tree, mcmc = .HousekeepingMcmc()),
                  "maxWarmup")
  ), collapse = "\n")
  ckpFile <- mkpEnv$recovery$mcmc$checkpointFile
  expect_true(file.exists(ckpFile))
  # A re-run would start afresh in a new temporary file.
  expect_no_match(msg, "Re-run the same")
  expect_match(msg, "ResumeMkPrime", fixed = TRUE)
  expect_match(msg, basename(ckpFile), fixed = TRUE)
})

test_that("MkPrimeRecover() finds an interrupted resume's samples (#424)", {
  skip_on_cran()
  d <- .HousekeepingData()
  dir <- withr::local_tempdir()
  logFile <- file.path(dir, "run.log")
  ckpFile <- file.path(dir, "run.ckp")
  mkpEnv <- environment(RunMkPrime)$.mkp_env
  set.seed(4242)
  allow_warning(
    RunMkPrime(d$pd, d$tree, mcmc = .HousekeepingMcmc(
      nIter = 1500L, logFile = logFile, checkpointFile = ckpFile
    )),
    "maxWarmup"
  )
  mkpEnv$recovery <- NULL

  .InterruptAtCheck()
  expect_warning(
    ResumeMkPrime(ckpFile, d$pd, d$tree, mcmc = list(nIter = 3000L)),
    "ResumeMkPrime"
  )
  local_mkp_verbosity()
  expect_message(rec <- MkPrimeRecover(), "Recovered")
  expect_s3_class(rec, "MkPosterior")
  expect_gt(nrow(rec$samples), 0L)
  # The logs belong to the checkpoint, which can still resume.
  expect_true(file.exists(logFile))
})

test_that("a resume clears the cancel file its run stopped for (#424)", {
  skip_on_cran()
  d <- .HousekeepingData()
  dir <- withr::local_tempdir()
  cancelFile <- file.path(dir, "cancel")
  Job <- function(nIter) {
    .HousekeepingMcmc(nIter = nIter, cancelFile = cancelFile,
                      logFile = file.path(dir, "run.log"))
  }
  file.create(cancelFile)
  set.seed(4243)
  first <- allow_warning(RunMkPrime(d$pd, d$tree, mcmc = Job(3000L)),
                         "maxWarmup|without drawing|before the first batch")
  expect_identical(first$stop_reason, "cancelled")
  expect_lt(first$actual_iter, 3000)
  expect_true(file.exists(cancelFile))

  # Re-running the same call resumes, and runs on to nIter.
  local_mkp_verbosity()
  expect_message(
    again <- allow_warning(RunMkPrime(d$pd, d$tree, mcmc = Job(3000L)),
                           "maxWarmup"),
    "Removed cancel file"
  )
  expect_false(file.exists(cancelFile))
  expect_identical(again$stop_reason, "max_iter")
  expect_equal(again$actual_iter, 3000)

  # A cancel file written since is a new request, and is obeyed.
  ckpFile <- file.path(dir, "run.ckp")
  cp <- readRDS(ckpFile)
  file.create(cancelFile)
  expect_false(identical(cp$runs[[1]]$cancelSeen$mtime,
                         file.mtime(cancelFile)))
  later <- ResumeMkPrime(ckpFile, d$pd, mcmc = list(nIter = 6000L))
  expect_identical(later$stop_reason, "cancelled")
  expect_lt(later$actual_iter, 6000)
})

test_that("a per-run log finds the checkpoint it was written beside (#424)", {
  dir <- withr::local_tempdir()
  ckp <- file.path(dir, "x.ckp")
  file.create(ckp)
  expect_equal(MkPrime:::.DerivedCkpFile(file.path(dir, "x_mkp_run_1.log")),
               ckp)
  expect_equal(MkPrime:::.DerivedCkpFile(file.path(dir, "x_mkp_run.log")), ckp)
  expect_equal(MkPrime:::.DerivedCkpFile(file.path(dir, "y_2.log")),
               file.path(dir, "y.ckp"))
  # No checkpoint beside it: the log's own name, as before.
  expect_equal(MkPrime:::.DerivedCkpFile(file.path(dir, "z_mkp_run_1.log")),
               file.path(dir, "z_mkp_run.ckp"))
})

test_that("parallel runs leave no cancel files behind (#424)", {
  skip_on_cran()
  skip_under_memcheck()
  skip_if_not_installed("callr")
  d <- .HousekeepingData()
  dir <- withr::local_tempdir()
  cancelFile <- file.path(dir, "cancel")
  file.create(cancelFile)
  # The per-run cancel files are bare tempfile()s.
  Bare <- function() list.files(tempdir(), "^file[0-9a-f]+$", full.names = TRUE)
  before <- Bare()
  set.seed(4245)
  allow_warning(
    RunMkPrime(d$pd, d$tree, mcmc = .HousekeepingMcmc(
      nRuns = 2L, nCore = 2L, cancelFile = cancelFile,
      logFile = file.path(dir, "run.log"), pollInterval = 1L
    )),
    "maxWarmup|never launched|cancel|without drawing|before the first batch"
  )
  expect_identical(setdiff(Bare(), before), character(0))
})

test_that("every master checkpoint names the move set (#425)", {
  mkdir <- withr::local_tempdir()
  file <- file.path(mkdir, "master.ckp")
  weights <- c(nni = 0.6, slice_tl = 0.4)
  runs <- list(list(moveWeights = weights, phase = "Sample"), NULL)
  MkPrime:::.SaveCheckpoint(runs, list(), 10L, "tree_length", file)
  payload <- readRDS(file)
  expect_equal(payload$moveWeights, weights)
  expect_identical(payload$phase, "Sample")

  skip_on_cran()
  d <- .HousekeepingData()
  dir <- withr::local_tempdir()
  ckpFile <- file.path(dir, "run.ckp")
  set.seed(4251)
  allow_warning(
    RunMkPrime(d$pd, d$tree, mcmc = .HousekeepingMcmc(
      nRuns = 2L, nIter = 1000L, checkpointFile = ckpFile,
      logFile = file.path(dir, "run.log")
    )),
    "maxWarmup"
  )
  payload <- readRDS(ckpFile)
  expect_setequal(names(payload$moveWeights),
                  names(payload$runs[[2]]$moveWeights))
  expect_identical(payload$phase, "Sample")
})

test_that("a resume is converged only if every run is (#425)", {
  skip_on_cran()
  d <- .HousekeepingData()
  dir <- withr::local_tempdir()
  ckpFile <- file.path(dir, "run.ckp")
  set.seed(4252)
  allow_warning(
    RunMkPrime(d$pd, d$tree, mcmc = .HousekeepingMcmc(
      nRuns = 2L, nIter = 1000L, checkpointFile = ckpFile,
      logFile = file.path(dir, "run.log")
    )),
    "maxWarmup"
  )
  cp <- readRDS(ckpFile)
  Converged <- function(r) {
    r$stop_reason <- "converged"
    r$conv_criteria <- list(minEss = 10, minTreeEss = NULL)
    r$actual_iter <- 900L
    r
  }
  cp$runs[[2]] <- Converged(cp$runs[[2]])

  # Run 1 reached nIter without converging; run 2, skipped, did converge.
  saveRDS(cp, ckpFile)
  resumed <- ResumeMkPrime(ckpFile, d$pd, mcmc = list(minEss = 10))
  expect_identical(resumed$stop_reason, "max_iter")

  # Both converged and are skipped: the job converged.
  cp$runs[[1]] <- Converged(cp$runs[[1]])
  saveRDS(cp, ckpFile)
  again <- ResumeMkPrime(ckpFile, d$pd, mcmc = list(minEss = 10))
  expect_identical(again$stop_reason, "converged")
})
