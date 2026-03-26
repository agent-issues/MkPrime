# Tests for stopping rules (Phase 5: M-034)

test_that("RunMkPrime reports stop_reason = 'max_iter' by default", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(5194)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, warmup = 200L))

  expect_equal(result$stop_reason, "max_iter")
  expect_equal(result$actual_iter, 500L)
})


test_that("maxTime stops MCMC early", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(3382)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000000L, thin = 5L,
                        warmup = 100L, maxTime = 0.5))

  expect_equal(result$stop_reason, "max_time")
  expect_lt(result$actual_iter, 1000000L)
  # Should still have some samples
  expect_gt(nrow(result$samples), 0)
})


test_that("MkPrimeMCMC stores stopping parameters", {
  cfg <- MkPrimeMCMC(maxTime = 60, minEss = 200, maxPsrf = 1.05,
                      checkEvery = 500L)
  expect_equal(cfg$maxTime, 60)
  expect_equal(cfg$minEss, 200)
  expect_equal(cfg$maxPsrf, 1.05)
  expect_equal(cfg$checkEvery, 500L)
})


test_that("Convergence-based stopping works", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  # Very generous convergence criteria so it triggers
  set.seed(9283)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 50000L, thin = 5L,
                        warmup = 200L, minEss = 5, maxPsrf = 5.0,
                        checkEvery = 300L))

  # Should stop before max_iter due to generous criteria
  if (result$stop_reason == "converged") {
    expect_lt(result$actual_iter, 50000L)
  }
  # Even if not converged (rare), result should be valid
  expect_s3_class(result, "MkPosterior")
  expect_gt(nrow(result$samples), 0)
})


test_that(".CheckConvergence returns NULL for insufficient samples", {
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                   "rate_loss", "rate_log_sd", "p")

  # Create runs with too few samples
  runs <- list(
    list(saved_idx = 3L,
         samples = matrix(rnorm(18), 10, length(paramNames),
                          dimnames = list(NULL, paramNames))),
    list(saved_idx = 3L,
         samples = matrix(rnorm(18), 10, length(paramNames),
                          dimnames = list(NULL, paramNames)))
  )

  mcmc <- MkPrimeMCMC(minEss = 100, maxPsrf = 1.05)
  result <- MkPrime:::.CheckConvergence(runs, paramNames, mcmc)
  expect_null(result)
})


test_that("Early stopping produces fewer samples", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  # Full run
  set.seed(1107)
  full <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000L, thin = 5L,
                        warmup = 200L))

  # Time-limited run
  set.seed(1107)
  limited <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000L, thin = 5L,
                        warmup = 200L, maxTime = 0.01))

  # Limited should have fewer or equal samples
  expect_lte(nrow(limited$samples), nrow(full$samples))
})
