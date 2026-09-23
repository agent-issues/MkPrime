# #8 and #23: a resumed run must reconstruct the same layout and the same
# sampler as the run it resumes.

test_that(".BrColStart selects exactly the br_ columns (#8)", {
  # The invariant that broke: the window `brColStart:(brColStart + nEdge - 1)`
  # must land on br_1..br_nEdge. ResumeMkPrime's hand-written formula was off
  # by one whenever the data had no neomorphic characters, which shifted every
  # reconstructed edge length by one column and read one past the end --
  # silently, because the tree still parsed.
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  transOnly <- MkPrimeData(MatrixToPhyDat(mat))
  expect_false(any(transOnly$type == "neomorphic"))

  for (nEdge in c(5L, 6L)) {
    for (prior in c("geometric", "beta_geometric", "logseries")) {
      for (qHet in c(FALSE, TRUE)) {
        nms <- .ParamNames(transOnly, nEdge, kPrimePrior = prior,
                           qHeterogeneity = qHet)
        start <- .BrColStart(nms)
        expect_identical(nms[start:(start + nEdge - 1L)],
                         paste0("br_", seq_len(nEdge)))
      }
    }
  }
})

test_that(".BrColStart also holds when neomorphic characters are present", {
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  mkd <- MkPrimeData(MatrixToPhyDat(mat), neomorphic = 1L)
  skip_if_not(any(mkd$type == "neomorphic"),
              "Could not construct a neomorphic character for this check.")

  nEdge <- 5L
  nms <- .ParamNames(mkd, nEdge)
  start <- .BrColStart(nms)
  expect_identical(nms[start:(start + nEdge - 1L)],
                   paste0("br_", seq_len(nEdge)))
})



test_that(".RemainingBudget spends one budget across the job (#13)", {
  # No budget set: nothing to apportion.
  expect_null(.RemainingBudget(NULL, proc.time()["elapsed"]))
  expect_identical(.RemainingBudget(Inf, proc.time()["elapsed"]), Inf)

  # A budget started "10 seconds ago" has ~10s less left, and never goes
  # negative once spent.
  now <- proc.time()["elapsed"]
  expect_lt(.RemainingBudget(60, now - 10), 60)
  expect_gt(.RemainingBudget(60, now - 10), 40)
  expect_identical(.RemainingBudget(5, now - 600), 0)

  # The value is unnamed, so it can be assigned into mcmc$maxTime cleanly.
  expect_null(names(.RemainingBudget(60, now - 10)))
})

test_that(".BuildResult returns a partial posterior when a run never started (#93)", {
  # Cancelling during run 1 of a multi-run job leaves runs 2..n as the bare
  # `.InitRun` structure: chain counters, but no `saved_idx`, `flush_idx` or
  # `tree_samples`.  Every consumer below dereferences those, so `NULL > 0L`
  # reached `if ()` as `logical(0)` and the job died instead of returning the
  # partial posterior that cancellation exists to produce.
  paramNames <- c("log_posterior", "log_likelihood", "tree_length")
  started <- list(
    tree_samples = vector("list", 4L),
    saved_idx = 4L, tree_saved_idx = 4L, flush_idx = 0L,
    chain_accept = list(c(nni = 1)), chain_propose = list(c(nni = 2)),
    chain_tuning = list(list()), logPostHistory = numeric(0),
    moveWeights = c(nni = 0.5, spr = 0.5)
  )
  neverStarted <- started[c("chain_accept", "chain_propose", "chain_tuning")]

  logs <- vapply(1:2, function(i) tempfile(fileext = ".log"), character(1L))
  on.exit(unlink(logs), add = TRUE)
  for (f in logs) {
    writeLines(paste(c("Sample", paramNames), collapse = "\t"), f)
  }

  mcmc <- list(nChains = 1L, nRuns = 2L, warmup = 0L, thin = 1L, treeThin = 1L,
               logFile = logs[[1]])

  expect_warning(
    result <- .BuildResult(list(started, neverStarted), model = NULL,
                          mkd = NULL, mcmc = mcmc, paramNames = paramNames,
                          logFilePaths = logs, actualIter = 4L,
                          stopReason = "cancelled"),
    "never started"
  )
  expect_s3_class(result, "MkPosterior")
  expect_identical(result$nSamples, 4L)
  expect_identical(result$requested_nRuns, 2L)
  expect_identical(result$dropped_runs$run, 2L)

  # More than one unstarted run: the report names them all.
  mcmc$nRuns <- 3L
  logs3 <- c(logs, logs[[2]])
  expect_warning(
    result <- .BuildResult(list(started, neverStarted, neverStarted),
                          model = NULL, mkd = NULL, mcmc = mcmc,
                          paramNames = paramNames,
                          logFilePaths = logs3, actualIter = 4L,
                          stopReason = "cancelled"),
    "never started"
  )
  expect_identical(result$dropped_runs$run, c(2L, 3L))

  # An unstarted run in the middle shifts every later run's position, so
  # `logFile` must be filtered in step: callers pair `per_run[[i]]` with
  # `logFile[i]`, and against the unfiltered vector run 3 would be read from
  # run 2's never-written log.
  expect_warning(
    result <- .BuildResult(list(started, neverStarted, started),
                          model = NULL, mkd = NULL, mcmc = mcmc,
                          paramNames = paramNames,
                          logFilePaths = logs3, actualIter = 4L,
                          stopReason = "cancelled"),
    "never started"
  )
  expect_identical(result$logFile, logs3[c(1L, 3L)])
  expect_length(result$per_run, 2L)
})


test_that("a run with no recorded position resumes from iteration 1 (#94)", {
  # `checkpoint$iter` is a job-wide maximum.  Handing it to a run that never
  # launched starts it at the fastest run's iteration -- past `warmup`, so it
  # samples with none.  Only a run's own position may be inherited.
  runs <- list(
    list(actual_iter = 60000L, phase = "Sample"),   # finished
    list(actual_iter =  5000L, phase = "Warmup"),   # stopped mid-warmup
    list()                                          # never launched
  )
  expect_identical(.ResumeStartIters(runs, fallback = 60001L),
                   c(60001L, 5001L, 1L))
})


skip_slow_tests()

test_that("a fixTopology run resumes with the same move set (#23)", {
  # A fixTopology run must not gain nni/spr/tbr/pspr on resume: the resumed
  # segment would explore topology space the user asked to hold fixed.
  tree <- ape::read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  ckpFile <- tempfile(fileext = ".rds")
  on.exit(unlink(ckpFile), add = TRUE)

  set.seed(4242)
  setTimeLimit(elapsed = 240, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)

  RunMkPrime(pd, tree, fixTopology = TRUE,
             mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000L, thin = 5L,
                                maxWarmup = 200L, minWarmup = 200L,
                                autoTune = FALSE, maxTime = 120,
                                checkEvery = 300L, checkpointFile = ckpFile))

  cp <- readRDS(ckpFile)

  # The run's own configuration now travels with the checkpoint.
  expect_true(isTRUE(cp$mcmc$fixTopology))

  # No topology move was ever scheduled, and resuming must not introduce one.
  topoMoves <- c("nni", "spr", "tbr", "pspr")
  expect_length(intersect(names(cp$moveWeights), topoMoves), 0L)

  # ResumeMkPrime warns when the rebuilt move set differs from the stored one,
  # so a silent resume is the assertion.
  expect_no_warning(
    resumed <- ResumeMkPrime(ckpFile, pd),
    message = "move set"
  )
  expect_s3_class(resumed, "MkPosterior")
})

test_that("maxTime bounds the job, not each run (#13)", {
  # maxTime must bound the job: with a fresh budget per run the real ceiling
  # is nRuns * maxTime, overrunning the external wall clock it was sized for.
  tree <- ape::read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  budget <- 15
  set.seed(31415)
  setTimeLimit(elapsed = 240, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)

  t0 <- proc.time()["elapsed"]
  warns <- testthat::capture_warnings(res <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 1e7, thin = 5L,
                       maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                       checkEvery = 500L, maxTime = budget,
                       # Unreachable, so the clock is what stops the job.
                       maxRhat = 1 + 1e-12)))
  elapsed <- unname(proc.time()["elapsed"] - t0)

  expect_identical(res$stop_reason, "max_time")
  # The budget stopped the job, so the second run never started; the result
  # records that two were requested.
  expect_identical(res$requested_nRuns, 2L)
  expect_true(any(grepl("maxTime", warns)))
  # The discriminating assertion: under the old behaviour this took about
  # 2 * budget. The margin is loose because the inner clock is only consulted
  # once per sampling batch.
  expect_lt(elapsed, 1.8 * budget)
})

test_that("resumed trees have finite edge lengths on neo-free data (#64)", {
  # A wrong brColStart shifts every edge length by one column and reads past
  # the end, writing `t_last:NA` into each post-resume tree; the newick still
  # parses, so nothing downstream fails. Only neo-free data is affected, and
  # the streaming fixture is all-transformational -- the existing resume tests
  # walk this path and assert nothing about it.
  tree <- ape::read.tree(
    text = "(((t1:0.1,t2:0.2):0.1,t3:0.15):0.1,(t4:0.1,t5:0.2):0.1,t6:0.3);"
  )
  mat <- matrix(c(0, 0, 1, 1, 0, 1,
                  0, 1, 0, 1, 1, 0,
                  1, 0, 0, 1, 0, 1,
                  1, 1, 0, 0, 1, 0,
                  0, 1, 1, 0, 0, 1,
                  1, 0, 1, 0, 1, 0),
                nrow = 6, byrow = TRUE,
                dimnames = list(paste0("t", 1:6), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  expect_false(any(mkd$type == "neomorphic"))

  ckpFile  <- tempfile(fileext = ".ckp")
  logFile <- tempfile(fileext = ".log")
  treeFile <- tempfile(fileext = ".trees")
  on.exit(unlink(c(ckpFile, logFile, treeFile)), add = TRUE)

  set.seed(606)
  setTimeLimit(elapsed = 300, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)

  # Stop on the clock so there is sampling left to resume into.
  first <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 200000L, thin = 5L,
                       maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                       maxTime = 8, checkEvery = 500L,
                       checkpointFile = ckpFile, logFile = logFile,
                       treeFile = treeFile))
  expect_false(anyNA(unlist(lapply(first$trees, `[[`, "edge.length"))))

  resumed <- ResumeMkPrime(ckpFile, pd)
  postEdges <- unlist(lapply(resumed$trees, `[[`, "edge.length"))
  expect_gt(length(postEdges), 0L)
  expect_false(anyNA(postEdges))

  # And the newick actually written to disk, which is what a user reads back.
  written <- ape::read.tree(treeFile)
  if (inherits(written, "phylo")) written <- list(written)
  expect_false(anyNA(unlist(lapply(written, `[[`, "edge.length"))))
})


test_that("a never-started run resumes into Warmup, not mid-Sample (#94)", {
  tree <- .mkp_test_tree()
  pd   <- .mkp_test_pd()

  ckpFile  <- tempfile(fileext = ".ckp")
  logFile  <- tempfile(fileext = ".log")
  logPaths <- .LogFilePaths(logFile, 2L)
  on.exit(unlink(c(ckpFile, logPaths)), add = TRUE)

  set.seed(94)
  setTimeLimit(elapsed = 300, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)

  warmup <- 200L
  thin   <- 5L
  allow_warning(
    RunMkPrime(pd, tree,
      mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 600L, thin = thin,
                         maxWarmup = warmup, minWarmup = warmup,
                         autoTune = FALSE, maxTime = 120, checkEvery = 200L,
                         checkpointFile = ckpFile, logFile = logFile)),
    "stabilis"
  )

  # Recast run 2 as one that never launched: exactly the field set `.InitRun`
  # returns, with no position, phase or sample count.  Its chain state is run
  # 2's own final state, which is immaterial -- what is under test is the
  # iteration the resume hands it.
  cp <- readRDS(ckpFile)
  cp$runs[[2]] <- cp$runs[[2]][c(
    "chains", "betas", "chain_accept", "chain_propose", "chain_time_ns",
    "chain_slice_exp", "chain_tuning", "chain_rhos", "swap_accept",
    "swap_propose"
  )]
  cp$mcmc$nIter <- 1200L
  saveRDS(cp, ckpFile)

  allow_warning(ResumeMkPrime(ckpFile, pd), "stabilis")

  nRows <- vapply(logPaths, function(f) length(readLines(f)) - 1L,
                  integer(1L), USE.NAMES = FALSE)
  # Run 2 warmed from scratch, so it covers the whole post-warmup budget.
  # Inheriting run 1's position put it at iteration 601 in the Sample phase
  # with `cppWarmup = 0`: it covered only the tail, and every draw it did take
  # was un-warmed, pooled into the posterior and into the cross-run R-hat.
  expect_gt(nRows[[1]], 0L)
  expect_gte(nRows[[2]], (1200L - warmup) %/% thin)
})


test_that("interrupting a resumed run does not rewind the checkpoint (#96)", {
  # `.RunSerialRuns` writes an epoch checkpoint carrying real `actual_iter`
  # values, but `ResumeMkPrime`'s own `runs` is refreshed only on normal
  # return.  The interrupt handler therefore held the pre-resume state, and
  # wrote it back over the master -- destroying exactly the work it exists to
  # preserve, and silently, since the next resume reports the log rewind as
  # routine post-checkpoint cleanup.
  tree <- .mkp_test_tree()
  pd   <- .mkp_test_pd()

  ckpFile  <- tempfile(fileext = ".ckp")
  logFile  <- tempfile(fileext = ".log")
  logPaths <- .LogFilePaths(logFile, 2L)
  on.exit(unlink(c(ckpFile, logPaths)), add = TRUE)

  set.seed(96)
  setTimeLimit(elapsed = 600, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)

  # An unreachable maxRhat routes both segments through `.RunSerialRuns` and
  # keeps its Phase 2 epoch loop running.
  RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 600L, thin = 5L,
                       maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                       maxRhat = 1 + 1e-12, maxTime = 240, checkEvery = 200L,
                       checkpointFile = ckpFile, logFile = logFile))
  before <- readRDS(ckpFile)$iter
  expect_equal(before, 600)

  # Simulated Ctrl-C.  `progressFn` is the only user hook inside the batch
  # loop, and it fires after the batch has published its state.  Injecting it
  # into the checkpoint leaves the first segment unhooked; its environment is
  # rooted at `baseenv()` so the checkpoint carries a closure, not the test
  # frame.  Phase 2 epochs are 1000 iterations, so the first epoch (both runs
  # to 1600) completes and is checkpointed before the second one trips this.
  Interrupter <- function(info) {
    if (info$iter > stopAt) {
      stop(structure(class = c("interrupt", "condition"),
                     list(message = "simulated Ctrl-C", call = NULL)))
    }
  }
  environment(Interrupter) <- list2env(list(stopAt = 1600), parent = baseenv())

  cp <- readRDS(ckpFile)
  cp$mcmc$nIter      <- 6000L
  cp$mcmc$plotEvery  <- 100L
  cp$mcmc$progressFn <- Interrupter
  saveRDS(cp, ckpFile)

  allow_warning(ResumeMkPrime(ckpFile, pd), "interrupted")
  after <- readRDS(ckpFile)
  expect_gt(after$iter, before)
  # `checkpoint$moveWeights` is what the move-set-drift check reads, so an
  # interrupt that writes a checkpoint without it silences that warning for
  # every later resume of the file.
  expect_false(is.null(after$moveWeights))
})
