#!/usr/bin/env Rscript
# Probe WITHOUT checkEvery=NULL — uses the default, matching the real test.
# Usage: Rscript tools/probe_crash2.R <nIter>
# Run from the mkp package root directory.

args <- commandArgs(trailingOnly = TRUE)
nIter <- as.integer(args[[1]])

cat(sprintf("CRASH_PROBE2: starting nIter=%d\n", nIter))

suppressMessages(devtools::load_all(quiet = TRUE))
library(ape)

tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
               dimnames = list(paste0("t", 1:4), NULL))
pd   <- TreeTools::MatrixToPhyDat(mat)

set.seed(2204)
result <- RunMkPrime(pd, tree,
  mcmc = MkPrimeMCMC(
    nRuns  = 2L,
    nIter  = nIter,
    thin   = 5L,
    warmup = as.integer(nIter * 0.5)
    # checkEvery = default
  )
)

cat(sprintf("CRASH_PROBE2: nIter=%d RESULT=OK\n", nIter))
