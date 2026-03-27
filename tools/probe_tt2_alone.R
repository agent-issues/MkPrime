#!/usr/bin/env Rscript
# Test 2 (multi-run, nIter=1000) alone in a test_that() block.
# Run from the mkp package root directory.

cat("PROBE_TT2: starting\n")
suppressMessages(devtools::load_all(quiet = TRUE))
library(ape)
library(testthat)

tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
               dimnames = list(paste0("t", 1:4), NULL))
pd   <- TreeTools::MatrixToPhyDat(mat)

cat("PROBE_TT2: running test\n")
test_that("test2_alone", {
  set.seed(2204)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 1000L, thin = 5L, warmup = 500L))
  cat("PROBE_TT2: RunMkPrime returned\n")
  expect_no_error(print(result))
  cat("PROBE_TT2: expect_no_error(print) passed\n")
  expect_equal(result$nRuns, 2L)
  cat("PROBE_TT2: expect_equal passed\n")
  expect_true(!is.null(result$per_run))
  cat("PROBE_TT2: expect_true passed — about to exit test_that block\n")
})
cat("PROBE_TT2: test_that returned — RESULT=OK\n")
