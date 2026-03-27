#!/usr/bin/env Rscript
# Test whether forced GC after RunMkPrime triggers the crash.
# Run from the mkp package root directory.

cat("PROBE_GC: starting\n")
suppressMessages(devtools::load_all(quiet = TRUE))
library(ape)

tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
               dimnames = list(paste0("t", 1:4), NULL))
pd   <- TreeTools::MatrixToPhyDat(mat)

cat("PROBE_GC: before RunMkPrime\n")

# Run in a local scope so result goes out of scope
local({
  set.seed(2204)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(
      nRuns  = 2L,
      nIter  = 1000L,
      thin   = 5L,
      warmup = 500L
    )
  )
  cat("PROBE_GC: RunMkPrime returned\n")
  # result goes out of scope here
})

cat("PROBE_GC: after local() — forcing gc()\n")
gc()
cat("PROBE_GC: gc() returned — RESULT=OK\n")
