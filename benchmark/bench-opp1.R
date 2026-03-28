# A/B benchmark for OPP-1 (JC O(k) product in flat pruning)
library(MkPrime)

nex <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
pd  <- TreeTools::ReadAsPhyDat(nex)
mkd <- MkPrimeData(pd)
mdl <- MkPrimeModel()

set.seed(7382)
tree <- ape::rtree(length(pd))
tree$tip.label <- names(pd)

cfg <- MkPrimeMCMC(nIter = 4000L, nRuns = 1L, nChains = 1L,
                   warmup = 0L, thin = 200L,
                   gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE,
                   tbr = FALSE)

# Warm up
invisible(RunMkPrime(mkd, tree, model = mdl,
                     mcmc = MkPrimeMCMC(nIter = 100L, nRuns = 1L,
                                        nChains = 1L, warmup = 0L,
                                        thin = 100L, gibbsSpr = TRUE,
                                        gibbsSubtreeSwap = TRUE,
                                        tbr = FALSE)))

# 3 timed runs
times <- numeric(3)
for (i in 1:3) {
  t <- system.time({
    r <- RunMkPrime(mkd, tree, model = mdl, mcmc = cfg)
  })
  times[i] <- t["elapsed"]
  cat(sprintf("Run %d: %.2f sec (%.1f iter/s)\n", i, times[i], 4000/times[i]))
}
cat(sprintf("\nMedian: %.2f sec (%.1f iter/s)\n",
            median(times), 4000/median(times)))
