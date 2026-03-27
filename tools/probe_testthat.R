#!/usr/bin/env Rscript
# Mimics tests exactly using test_that() to see if testthat env teardown is the trigger.
# Run from the mkp package root directory.

cat("PROBE_TT: starting\n")
suppressMessages(devtools::load_all(quiet = TRUE))
library(ape)
library(testthat)

tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
               dimnames = list(paste0("t", 1:4), NULL))
pd   <- TreeTools::MatrixToPhyDat(mat)

cat("PROBE_TT: test 1 — single run\n")
test_that("test1", {
  set.seed(5194)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, warmup = 200L))
  expect_no_error(print(result))
  expect_no_error(summary(result))
  expect_no_error(plot(result))
})
cat("PROBE_TT: test 1 done\n")

cat("PROBE_TT: test 2 — multi-run\n")
test_that("test2", {
  set.seed(2204)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 1000L, thin = 5L, warmup = 500L))
  expect_no_error(print(result))
  expect_equal(result$nRuns, 2L)
  expect_true(!is.null(result$per_run))
})
cat("PROBE_TT: test 2 done\n")

cat("PROBE_TT: RESULT=OK\n")
