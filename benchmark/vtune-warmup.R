# VTune driver: profile the warmup phase of RunMkPrime on Sun2018
# Target: ~60-90s of CPU time in warmup iterations
.libPaths(c(".vtune-lib", .libPaths()))
library(MkPrime)
library(TreeTools)

# Load Sun2018 NEXUS from TreeSearch
nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
pd <- TreeTools::ReadAsPhyDat(nexFile)

# Classify: 2-state = neomorphic, >2-state = transformational
nStates <- vapply(pd, function(x) length(attr(x, "levels")), integer(1))
neomorphicChars <- which(nStates == 2)
transformationalChars <- which(nStates > 2)

mkd <- MkPrimeData(pd, neomorphic = neomorphicChars)

model <- MkPrimeModel(nCat = 6)

# Warmup-focused config: nIter = 5000 forces exit after 5000 iters.
# minWarmup = 5000 ensures we stay in warmup the whole time.
# maxTime = 120 as safety net.
mcmc <- MkPrimeMCMC(
  nIter = 5000L,
  minWarmup = 5000L,
  maxWarmup = 5000L,
  nRuns = 1L,
  nChains = 1L,
  maxTime = 120,
  plotEvery = NULL,
  logFile = NULL,
  treeFile = NULL,
  checkpointFile = NULL
)

cat("Starting warmup profiling run...\n")
t0 <- proc.time()
posterior <- RunMkPrime(mkd, model = model, mcmc = mcmc)
elapsed <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("Completed in %.1f seconds\n", elapsed))
