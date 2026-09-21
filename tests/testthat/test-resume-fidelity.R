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
