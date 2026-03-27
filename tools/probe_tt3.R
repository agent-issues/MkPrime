#!/usr/bin/env Rscript
# Narrow down: does print(result) crash inside test_that, without expect_no_error?
# Run from the mkp package root directory.

cat("PROBE_TT3: starting\n")
suppressMessages(devtools::load_all(quiet = TRUE))
library(ape)
library(testthat)

tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
               dimnames = list(paste0("t", 1:4), NULL))
pd   <- TreeTools::MatrixToPhyDat(mat)

# ── Variant A: print() without expect_no_error ──
cat("PROBE_TT3: variant A — bare print() in test_that\n")
test_that("variantA", {
  set.seed(2204)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 1000L, thin = 5L, warmup = 500L))
  cat("PROBE_TT3: RunMkPrime returned\n")
  print(result)
  cat("PROBE_TT3: print returned\n")
})
cat("PROBE_TT3: variant A test_that returned\n")

cat("PROBE_TT3: RESULT=OK\n")
