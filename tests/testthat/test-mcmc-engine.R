# Smoke tests for MCMC engine (M-016, M-017, M-018, M-019, M-020)

test_that("RunMkPrime runs on simple transformational data", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(5194)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nIter = 500L, thin = 5L, warmup = 200L))

  expect_s3_class(result, "MkPosterior")
  expect_equal(nrow(result$samples), 60L)
  expect_equal(length(result$trees), 60L)
  expect_true(all(is.finite(result$samples[, "log_posterior"])))
})


test_that("RunMkPrime handles neomorphic characters", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(4781)
  result <- RunMkPrime(pd, tree, neomorphic = 1L,
    mcmc = MkPrimeMCMC(nIter = 500L, thin = 5L, warmup = 200L))

  expect_s3_class(result, "MkPosterior")
  expect_true("rate_loss" %in% names(result$acceptance))
  # rate_loss should vary (neomorphic chars present)
  rate_loss_vals <- result$samples[, "rate_loss"]
  expect_true(sd(rate_loss_vals) > 0)
})


test_that("RunMkPrime handles known state-space characters", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(2837)
  result <- RunMkPrime(pd, tree, known_states = c("1" = 4L),
    mcmc = MkPrimeMCMC(nIter = 500L, thin = 5L, warmup = 200L))

  expect_s3_class(result, "MkPosterior")
  # No kPrime columns (no transformational chars)
  expect_false(any(grepl("kPrime", colnames(result$samples))))
})


test_that("RunMkPrime handles mixed character types", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1,
                  0, 1, 2, 0,
                  0, 1, 0, 1), 4, 3,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(9201)
  result <- RunMkPrime(pd, tree, neomorphic = 1L,
    known_states = c("3" = 3L),
    mcmc = MkPrimeMCMC(nIter = 500L, thin = 5L, warmup = 200L))

  expect_s3_class(result, "MkPosterior")
  # Should have kPrime_2 (char 2 is transformational)
  expect_true(any(grepl("kPrime", colnames(result$samples))))
})


test_that("MkPosterior print, summary, plot methods work", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat <- matrix(c(0, 1, 0), 3, 1,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(7712)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nIter = 300L, thin = 3L, warmup = 150L))

  expect_no_error(print(result))
  s <- summary(result)
  expect_true(is.data.frame(s))
  expect_true(all(c("parameter", "mean", "median") %in% names(s)))

  # plot should not error
  expect_no_error(plot(result))
})


test_that("Acceptance rates are non-degenerate", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(8371)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nIter = 2000L, thin = 10L, warmup = 1000L))

  # No move type should have 0% or 100% acceptance
  for (nm in names(result$acceptance)) {
    expect_gt(result$acceptance[nm], 0, label = paste(nm, "acceptance > 0"))
    expect_lt(result$acceptance[nm], 1, label = paste(nm, "acceptance < 1"))
  }
})
