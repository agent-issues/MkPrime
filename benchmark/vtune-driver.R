# VTune driver script for M-166 (S-PROF round 5)
# Exercises the MCMC hot path on Sun2018 (54 taxa, 225 chars)
# Target: ~60-90s of CPU time in the C++ inner loop

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
# Short warmup, no tuning, straight to sampling for clean VTune profile.
set.seed(4619)
posterior <- RunMkPrime(
  mkd, tree,
  nIter      = 5000L,
  maxWarmup  = 500L,
  minWarmup  = 200L,
  autoTune   = FALSE,
  nRuns      = 1L,
  nChains    = 1L,
  thin       = 50L,
  maxTime    = 180,
  plotEvery  = 0L
)

cat("Completed", posterior$nIter, "iterations,", posterior$nSamples, "samples\n")
cat("Stop reason:", posterior$stopReason, "\n")
