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

# --- #195: runs stuck on one topology are compared across runs ---

# 6 tips, 160 characters each congruent with one tree: every run settles on it.
.CongruentData <- function() {
  tree <- BalancedTree(6)
  splits <- as.logical(as.Splits(tree))
  mat <- t(splits[rep_len(seq_len(nrow(splits)), 160L), ]) * 1L
  rownames(mat) <- TipLabels(tree)
  MatrixToPhyDat(mat)
}

.StuckRun195 <- function(tree, reason = "single_topology", n = 50L) {
  list(stop_reason = reason, tree_samples = rep(list(tree), n),
       tree_saved_idx = n, actual_iter = 1000L)
}

test_that(".TopologyFingerprint() ignores rooting and tip order (#195)", {
  tree <- as.phylo(0, 8)
  key <- .TopologyFingerprint(tree)
  expect_identical(.TopologyFingerprint(RootTree(tree, "t5")), key)
  expect_identical(.TopologyFingerprint(UnrootTree(tree)), key)
  expect_identical(.TopologyFingerprint(RenumberTips(tree, rev(TipLabels(tree)))), key)
  expect_false(identical(.TopologyFingerprint(as.phylo(1, 8)), key))
})

test_that("runs stuck on one shared topology meet minTreeEss (#195)", {
  t1 <- .TopologyFingerprint(as.phylo(0, 8))
  t2 <- .TopologyFingerprint(as.phylo(1, 8))
  Stuck <- function(key) list(met = FALSE, stuck = TRUE, topology = key)
  moving <- list(met = FALSE, stuck = FALSE, topology = NA_character_)

  expect_identical(.TreeVerdictsMet(list(Stuck(t1), Stuck(t1))), c(TRUE, TRUE))
  expect_identical(.TreeVerdictsMet(list(Stuck(t1), Stuck(t2))),
                   c(FALSE, FALSE))
  expect_identical(.TreeVerdictsMet(list(Stuck(t1), moving)), c(FALSE, FALSE))
  expect_identical(.TreeVerdictsMet(list(Stuck(t1), NULL)), c(FALSE, FALSE))
  # One run alone cannot agree with anything.
  expect_false(.TreeVerdictsMet(list(Stuck(t1))))
  # A verdict from before #195 carries no topology.
  expect_identical(.TreeVerdictsMet(list(list(met = TRUE), list(met = FALSE))),
                   c(TRUE, FALSE))
})

test_that("a stuck worker's topology reaches the parallel parent (#195)", {
  skip_if_not_installed("TreeDist")
  files <- replicate(2L, tempfile())
  on.exit(unlink(files), add = TRUE)
  mcmc <- list(minTreeEss = 100)
  run <- list(saved_idx = 50L, tree_saved_idx = 50L,
              tree_samples = rep(list(as.phylo(0, 8)), 50L),
              samples = matrix(rnorm(50), 50, 1,
                               dimnames = list(NULL, "log_posterior")))
  diag <- .CheckConvergence(list(run), "log_posterior", mcmc)
  expect_identical(diag$treeEssStatus, "stuck")
  expect_identical(diag$treeTopology, .TopologyFingerprint(as.phylo(0, 8)))
  # Workers never hand over: the parent decides.
  expect_false(diag$singleTopology)

  .WriteTreeVerdict(files[1], diag, mcmc, 1000L)
  .WriteTreeVerdict(files[2], diag, mcmc, 1000L)
  expect_identical(.ReadTreeVerdicts(files), c(TRUE, TRUE))
  run$tree_samples <- rep(list(as.phylo(1, 8)), 50L)
  .WriteTreeVerdict(files[2], .CheckConvergence(list(run), "log_posterior",
                                                mcmc), mcmc, 2000L)
  expect_identical(.ReadTreeVerdicts(files), c(FALSE, FALSE))
})

test_that("a serial run hands over only once its scalars are met (#195)", {
  skip_if_not_installed("TreeDist")
  set.seed(195)
  run <- list(saved_idx = 50L, tree_saved_idx = 50L,
              tree_samples = rep(list(as.phylo(0, 8)), 50L),
              samples = matrix(rnorm(50), 50, 1,
                               dimnames = list(NULL, "log_posterior")))
  Check <- function(...) {
    .CheckConvergence(list(run), "log_posterior",
                      list(minTreeEss = 100, stuckHandOff = TRUE, ...))
  }
  expect_true(Check()$singleTopology)
  expect_true(Check(minEss = 20)$singleTopology)
  expect_false(Check(minEss = 1000)$singleTopology)
  expect_false(Check(minEss = 20)$converged)
})

test_that(".SingleTopologyStop() compares the runs (#195)", {
  skip_if_not_installed("TreeDist")
  t1 <- as.phylo(0, 8)
  t2 <- as.phylo(1, 8)
  mcmc <- list(nRuns = 2L, minTreeEss = 100)

  none <- .SingleTopologyStop(list(.StuckRun195(t1, "max_iter")), mcmc)
  expect_null(none$stopReason)

  agree <- .SingleTopologyStop(list(.StuckRun195(t1), .StuckRun195(t1)), mcmc)
  expect_identical(agree$stopReason, "converged")
  expect_identical(vapply(agree$runs, `[[`, "", "stop_reason"),
                   c("converged", "converged"))

  expect_warning(
    apart <- .SingleTopologyStop(list(.StuckRun195(t1), .StuckRun195(t2)),
                                 mcmc),
    "settled"
  )
  expect_identical(apart$stopReason, "single_topology")
  expect_identical(apart$runs[[1]]$stop_reason, "single_topology")

  expect_warning(
    mixed <- .SingleTopologyStop(
      list(.StuckRun195(t1), .StuckRun195(t1, "converged")), mcmc
    ),
    "1 of 2"
  )
  expect_identical(mixed$stopReason, "single_topology")

  expect_warning(
    alone <- .SingleTopologyStop(list(.StuckRun195(t1)),
                                 list(nRuns = 1L, minTreeEss = 100)),
    "nRuns"
  )
  expect_identical(alone$stopReason, "single_topology")
  expect_identical(.WorkersStopReason(c("single_topology", "converged"), 2L),
                   "single_topology")
})

test_that("serial runs on one shared topology converge (#195)", {
  skip_if_not_installed("TreeDist")
  pd <- .CongruentData()
  .SetHangGuard(180)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  Run <- function(...) {
    set.seed(195)
    RunMkPrime(pd, mcmc = MkPrimeMCMC(
      nIter = 1e7, thin = 1L, minWarmup = 200L, maxWarmup = 200L,
      autoTune = FALSE, checkEvery = 200L, minTreeEss = 100, maxTime = 40,
      nCore = 1L, ...
    ))
  }

  warns <- testthat::capture_warnings(alone <- Run(nRuns = 1L))
  expect_identical(alone$stop_reason, "single_topology")
  expect_true(any(grepl("nRuns", warns)))

  allow_warning(both <- Run(nRuns = 2L), "single tree topology|maxWarmup")
  expect_identical(both$stop_reason, "converged")
  expect_identical(vapply(both$per_run, `[[`, "", "stop_reason"),
                   c("converged", "converged"))

  # Loose: this checks that agreeing runs reach the R-hat phase, not mixing.
  allow_warning(rhat <- Run(nRuns = 2L, maxRhat = 5),
                "single tree topology|maxWarmup")
  expect_identical(rhat$stop_reason, "converged")
})

test_that("parallel runs on one shared topology converge (#195)", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("TreeDist")
  skip_if_not(.IsMkPrimeInstalled(),
              "MkPrime not installed; callr workers need the installed package")
  pd <- .CongruentData()
  .SetHangGuard(180)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  set.seed(195)
  allow_warning(
    res <- RunMkPrime(pd, mcmc = MkPrimeMCMC(
      nRuns = 2L, nCore = 2L, nIter = 1e7, thin = 1L, minWarmup = 200L,
      maxWarmup = 200L, autoTune = FALSE, checkEvery = 200L,
      minTreeEss = 100, pollInterval = 1L, maxTime = 90
    )),
    "single tree topology|maxWarmup"
  )
  expect_identical(res$stop_reason, "converged")
})

test_that("serial runs reach the R-hat phase only if they agree (#195)", {
  skip_if_not_installed("TreeDist")
  Serial <- function(trees) {
    phase2 <- new.env()
    phase2$calls <- 0L
    local_mocked_bindings(
      .RunMkPrimeSingleRun = function(mkd, model, mcmc, initialState, moves,
                                      tipLabels, runIdx, ..., startIter = 1L) {
        run <- .StuckRun195(trees[[runIdx]])
        run$actual_iter <- min(mcmc$nIter, startIter + 999L)
        run
      },
      .CheckConvergenceFromLogs = function(...) {
        phase2$calls <- phase2$calls + 1L
        list(converged = TRUE, maxRhat = 1)
      }
    )
    mcmc <- list(nRuns = 2L, nIter = 5000L, checkEvery = 1000L,
                 maxRhat = 1.01, minTreeEss = 100)
    runs <- list(list(actual_iter = 0L), list(actual_iter = 0L))
    out <- .RunSerialRuns(NULL, NULL, mcmc, runs, NULL, NULL, NULL, 0L, 0L,
                          logFilePaths = c("a", "b"), convWindowSize = 0L,
                          treeFilePaths = NULL)
    c(out, phase2Calls = phase2$calls)
  }
  t1 <- as.phylo(0, 8)

  agree <- Serial(list(t1, t1))
  expect_identical(agree$stopReason, "converged")
  expect_gt(agree$phase2Calls, 0L)

  expect_warning(apart <- Serial(list(t1, as.phylo(1, 8))), "settled")
  expect_identical(apart$stopReason, "single_topology")
  expect_identical(apart$phase2Calls, 0L)
})
