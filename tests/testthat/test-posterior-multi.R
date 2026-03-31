# Tests for MkPosterior multi-run updates (Phase 5: M-036)
skip_slow_tests()

test_that("print.MkPosterior works for single run", {
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(5194)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE))

  expect_no_error(print(result))
  expect_no_error(summary(result))
  expect_no_error(plot(result))
})


test_that("print.MkPosterior works for multi-run", {
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(2204)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 1000L, thin = 5L, maxWarmup = 500L, minWarmup = 500L, autoTune = FALSE))

  expect_no_error(print(result))
  # Multi-run should have nRuns and per_run
  expect_equal(result$nRuns, 2L)
  expect_true(!is.null(result$per_run))
})


test_that("print.MkPosterior shows tempering info", {
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(6842)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nChains = 2L, heat = 0.3,
                        nIter = 500L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE))

  # Swap rates should be present in result
  expect_true(!is.null(result$swap_rates))
  expect_no_error(print(result))
})


test_that("summary.MkPosterior includes ESS and R-hat for multi-run", {
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(4487)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 1000L, thin = 5L, maxWarmup = 500L, minWarmup = 500L, autoTune = FALSE))

  s <- summary(result)
  expect_true(is.data.frame(s))
  expect_true("ESS" %in% names(s))
  expect_true("Rhat" %in% names(s))
  # ESS should be positive (or NA for constant params)
  expect_true(all(is.na(s$ESS) | s$ESS > 0))
})


test_that("summary.MkPosterior works for single run (no R-hat)", {
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(1498)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE))

  s <- summary(result)
  expect_true("ESS" %in% names(s))
  expect_false("Rhat" %in% names(s))
})


test_that("plot.MkPosterior shows multi-run traces", {
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(7218)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 500L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE))

  expect_no_error(plot(result))
})


test_that("Result includes stop_reason and actual_iter", {
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(3382)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 100000L, thin = 5L,
                        maxWarmup = 100L, minWarmup = 100L,
                        autoTune = FALSE, maxTime = 0.5))

  expect_equal(result$stop_reason, "max_time")
  expect_lt(result$actual_iter, 100000L)
  expect_no_error(print(result))
})
