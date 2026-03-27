#!/usr/bin/env Rscript
# Mimics the test sequence in test-posterior-multi.R.
# Run from the mkp package root directory.

cat("PROBE_SEQ: starting\n")
suppressMessages(devtools::load_all(quiet = TRUE))
library(ape)

tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
               dimnames = list(paste0("t", 1:4), NULL))
pd   <- TreeTools::MatrixToPhyDat(mat)

# ── Test 1: single run, nIter=500 ──
cat("PROBE_SEQ: test 1 start\n")
local({
  set.seed(5194)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, warmup = 200L))
  cat("PROBE_SEQ: test 1 RunMkPrime returned\n")
})
cat("PROBE_SEQ: test 1 scope exited — gc()\n")
gc(verbose = FALSE)
cat("PROBE_SEQ: test 1 gc done\n")

# ── Test 2: multi-run, nIter=1000 ──
cat("PROBE_SEQ: test 2 start\n")
local({
  set.seed(2204)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 1000L, thin = 5L, warmup = 500L))
  cat("PROBE_SEQ: test 2 RunMkPrime returned\n")
})
cat("PROBE_SEQ: test 2 scope exited — gc()\n")
gc(verbose = FALSE)
cat("PROBE_SEQ: RESULT=OK\n")
