# Job-level stop decisions: a parallel job honours minTreeEss and its runs'
# own stops (#355), and maxTime bounds the whole job on the sequential paths
# (#356).

.IsMkPrimeInstalled <- function() {
  if (!requireNamespace("pkgload", quietly = TRUE)) return(TRUE)
  !isTRUE(pkgload::is_dev_package("MkPrime"))
}

# Valgrind inflates the warmup and set-up that precede the clock's first
# check, so a job that honours maxTime still overruns any margin that would
# tell one budget from two.
.ExpectWithinBudget <- function(elapsed, budget) {
  if (identical(Sys.getenv("MKPRIME_MEMCHECK"), "true")) {
    return(invisible())
  }
  expect_lt(elapsed, 1.6 * budget)
}

# A hang guard, not an assertion: valgrind runners vary enough in speed to
# push a healthy sequential job past the native limit.
.SetHangGuard <- function(seconds) {
  if (identical(Sys.getenv("MKPRIME_MEMCHECK"), "true")) {
    seconds <- 4 * seconds
  }
  setTimeLimit(elapsed = seconds, transient = TRUE)
}

test_that(".WorkersStopReason() takes the job's reason from its runs (#355)", {
  expect_identical(.WorkersStopReason(c("converged", "converged"), 2L),
                   "converged")
  # A run that returned nothing did not converge.
  expect_identical(.WorkersStopReason("converged", 2L), "max_iter")
  expect_identical(.WorkersStopReason(c("converged", "max_iter"), 2L),
                   "max_iter")
  expect_identical(.WorkersStopReason(c("max_time", "converged"), 2L),
                   "max_time")
  expect_identical(.WorkersStopReason(c("max_time", "cancelled"), 2L),
                   "cancelled")
})

test_that("a run's minTreeEss verdict reaches the parallel parent (#355)", {
  files <- replicate(3L, tempfile())
  on.exit(unlink(files), add = TRUE)
  mcmc <- list(minTreeEss = 100)
  .WriteTreeVerdict(files[1], list(treeEss = 150, treeEssStatus = "ok"),
                    mcmc, 1000L)
  .WriteTreeVerdict(files[2], list(treeEss = 50, treeEssStatus = "ok"),
                    mcmc, 1000L)
  # Not yet computed, as when the run's scalars are far off.
  .WriteTreeVerdict(files[3], list(treeEss = NA_real_,
                                   treeEssStatus = NA_character_),
                    mcmc, 1000L)
  expect_identical(.ReadTreeVerdicts(files), c(TRUE, FALSE, FALSE))
  unlink(files[3])
  expect_identical(.ReadTreeVerdicts(files), c(TRUE, FALSE, FALSE))
  expect_identical(.ReadTreeVerdicts(NULL), logical(0))
})

test_that("a parallel run computes its tree ESS on its share of minEss (#355)", {
  skip_if_not_installed("TreeDist")
  set.seed(5)
  n <- 200L
  run <- list(
    samples = matrix(rnorm(2L * n), n, 2L,
                     dimnames = list(NULL, c("tree_length", "rate_log_sd"))),
    saved_idx = n,
    tree_samples = lapply(seq_len(n), function(i) ape::rtree(6L)),
    tree_saved_idx = n
  )
  mcmc <- list(minEss = 1000, minTreeEss = 50)
  # About 200 per run: far off the full 1000, but half of a quarter-share.
  alone <- .CheckConvergence(list(run), colnames(run$samples), mcmc)
  expect_identical(alone$treeEssPrecision, "skip")
  shared <- .CheckConvergence(list(run), colnames(run$samples), mcmc,
                              essRuns = 4L)
  expect_false(identical(shared$treeEssPrecision, "skip"))
  expect_true(is.finite(shared$treeEss))
  # The share decides only whether to look, not what the run needs.
  expect_false(shared$converged)
})

test_that("maxTime bounds a sequential job without maxRhat (#356)", {
  tree <- ape::read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  budget <- 6
  set.seed(2718)
  .SetHangGuard(120)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)

  t0 <- proc.time()[["elapsed"]]
  warns <- testthat::capture_warnings(res <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nCore = 1L, nIter = 1e7, thin = 5L,
                       maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                       checkEvery = 500L, maxTime = budget,
                       # Unreachable, so only the clock stops the job.
                       minEss = 1e6)))
  elapsed <- proc.time()[["elapsed"]] - t0

  expect_identical(res$stop_reason, "max_time")
  expect_identical(res$requested_nRuns, 2L)
  expect_true(any(grepl("maxTime.+run 1 of 2", warns)))
  # A fresh budget per run took about 2 * budget.
  .ExpectWithinBudget(elapsed, budget)
})

test_that("maxTime bounds a sequential resume (#356)", {
  tree <- ape::read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  ckpFile <- file.path(dir, "job.ckp")

  set.seed(1618)
  .SetHangGuard(120)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  allow_warning(RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nCore = 1L, nIter = 400L, thin = 5L,
                       maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                       checkEvery = 100L, checkpointFile = ckpFile,
                       logFile = file.path(dir, "job.log"), maxTime = 30)),
                "maxWarmup")
  expect_true(file.exists(ckpFile))

  budget <- 6
  t0 <- proc.time()[["elapsed"]]
  warns <- testthat::capture_warnings(res <- ResumeMkPrime(ckpFile, pd,
    mcmc = list(nIter = 1e7, minEss = 1e6, maxTime = budget)))
  elapsed <- proc.time()[["elapsed"]] - t0

  expect_identical(res$stop_reason, "max_time")
  expect_true(any(grepl("maxTime.+run 1 of 2", warns)))
  .ExpectWithinBudget(elapsed, budget)
})

test_that("a parallel job does not stop converged short of minTreeEss (#355)", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("TreeDist")
  skip_if_not(.IsMkPrimeInstalled(),
              "MkPrime not installed; callr workers need the installed package")

  set.seed(42)
  mat <- matrix(sample(0:1, 10 * 30, TRUE), 10, 30,
                dimnames = list(paste0("t", 1:10), NULL))
  pd <- MatrixToPhyDat(mat)

  budget <- 25
  .SetHangGuard(180)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  set.seed(1)
  # The pooled scalar check passes within the budget, but no run can sample
  # enough trees in it to reach this tree ESS.
  res <- RunMkPrime(pd,
    mcmc = MkPrimeMCMC(nRuns = 2L, nCore = 2L, nIter = 1e7, minEss = 20,
                       minTreeEss = 5000, checkEvery = 1000L,
                       pollInterval = 1L, maxTime = budget))

  expect_identical(res$stop_reason, "max_time")
})

test_that("a parallel job whose runs all converge reports it (#355)", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("TreeDist")
  skip_if_not(.IsMkPrimeInstalled(),
              "MkPrime not installed; callr workers need the installed package")

  set.seed(42)
  mat <- matrix(sample(0:1, 10 * 30, TRUE), 10, 30,
                dimnames = list(paste0("t", 1:10), NULL))
  pd <- MatrixToPhyDat(mat)

  .SetHangGuard(240)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  set.seed(1)
  # minTreeEss alone gives the parent nothing to judge, so each run stops
  # itself.
  res <- RunMkPrime(pd,
    mcmc = MkPrimeMCMC(nRuns = 2L, nCore = 2L, nIter = 1e7, minTreeEss = 30,
                       checkEvery = 1000L, pollInterval = 1L,
                       maxTime = 120))

  expect_identical(res$stop_reason, "converged")
  expect_identical(vapply(res$per_run, `[[`, "", "stop_reason"),
                   c("converged", "converged"))
})
