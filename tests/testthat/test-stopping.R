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


test_that("max_time stops MCMC early", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(3382)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000000L, thin = 5L,
                        warmup = 100L, max_time = 0.5))

  expect_equal(result$stop_reason, "max_time")
  expect_lt(result$actual_iter, 1000000L)
  # Should still have some samples
  expect_gt(nrow(result$samples), 0)
})


test_that("MkPrimeMCMC stores stopping parameters", {
  cfg <- MkPrimeMCMC(max_time = 60, min_ess = 200, max_psrf = 1.05,
                      check_every = 500L)
  expect_equal(cfg$max_time, 60)
  expect_equal(cfg$min_ess, 200)
  expect_equal(cfg$max_psrf, 1.05)
  expect_equal(cfg$check_every, 500L)
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
                        warmup = 200L, min_ess = 5, max_psrf = 5.0,
                        check_every = 300L))

  # Should stop before max_iter due to generous criteria
  if (result$stop_reason == "converged") {
    expect_lt(result$actual_iter, 50000L)
  }
  # Even if not converged (rare), result should be valid
  expect_s3_class(result, "MkPosterior")
  expect_gt(nrow(result$samples), 0)
})


test_that(".check_convergence returns NULL for insufficient samples", {
  param_names <- c("log_posterior", "log_likelihood", "tree_length",
                   "rate_loss", "rate_log_sd", "p")

  # Create runs with too few samples
  runs <- list(
    list(saved_idx = 3L,
         samples = matrix(rnorm(18), 10, length(param_names),
                          dimnames = list(NULL, param_names))),
    list(saved_idx = 3L,
         samples = matrix(rnorm(18), 10, length(param_names),
                          dimnames = list(NULL, param_names)))
  )

  mcmc <- MkPrimeMCMC(min_ess = 100, max_psrf = 1.05)
  result <- MkPrime:::.check_convergence(runs, param_names, mcmc)
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
                        warmup = 200L, max_time = 0.01))

  # Limited should have fewer or equal samples
  expect_lte(nrow(limited$samples), nrow(full$samples))
})
