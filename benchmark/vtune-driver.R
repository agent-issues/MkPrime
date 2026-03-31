# VTune driver script for S-PROF round 4
# Exercises the MCMC hot path on Sun2018 (54 taxa, 225 chars)
# Target: ~30s of CPU time in the C++ inner loop

library(MkPrime, lib.loc = ".vtune-lib")

# Load Sun2018 dataset from TreeSearch's bundled nexus file
nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
phyDat <- TreeTools::ReadAsPhyDat(nexFile)
mkd <- MkPrimeData(phyDat)

cat("Dataset:", nrow(mkd$matrix), "taxa,", ncol(mkd$matrix), "characters\n")
cat("Types:", table(mkd$type), "\n")

# Use a fixed NJ starting tree for reproducibility
tree <- TreeTools::NJTree(phyDat, edgeLengths = TRUE)
tree$edge.length[tree$edge.length <= 0] <- 1e-8

# Run MCMC: fixed nIter, no convergence criteria, no tempering overhead.
# 15000 iterations at ~2ms/iter ≈ 30s of hot-path CPU time.
set.seed(4619)
posterior <- RunMkPrime(
  mkd, tree,
  nIter    = 15000L,
  warmup   = 1000L,
  nRuns    = 1L,
  nChains  = 1L,
  thin     = 50L,
  maxTime  = 60,
  plotEvery = 0L
)

cat("Completed", posterior$nIter, "iterations,", posterior$nSamples, "samples\n")
cat("Stop reason:", posterior$stopReason, "\n")
