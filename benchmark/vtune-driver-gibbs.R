# VTune hotspot driver for Gibbs/Weighted moves
#
# Exercises specific Gibbs/weighted move configurations for ~30s CPU each.
# Complements vtune-driver.R (which profiles the standard move schedule).
#
# Usage:
#   vtune -collect hotspots -result-dir vtune-gibbs-b \
#         -- Rscript benchmark/vtune-driver-gibbs.R b
#
#   vtune -collect hotspots -result-dir vtune-gibbs-c \
#         -- Rscript benchmark/vtune-driver-gibbs.R c
#
# Configs:
#   b — +GibbsSPR + GibbsSubtreeSwap (default ON; isolates clone + pruning)
#   c — +WeightedSPR (bin marginalisation overhead)
#   d — (b) + BlockGibbsBranch (full sweep cost)
#   e — +WeightedBranchScale (bin-based branch scaling)
#
# Build for profiling first (see PROFILING.md):
#   1. Add -g -fno-omit-frame-pointer to src/Makevars.win PKG_CXXFLAGS
#   2. Install into .vtune-lib with DLLFLAGS override
#   3. Remove src/Makevars.win after profiling

args <- commandArgs(trailingOnly = TRUE)
config <- if (length(args) >= 1) tolower(args[1]) else "b"
if (!config %in% c("b", "c", "d", "e")) {
  stop("Config must be one of: b, c, d, e. Got: ", config)
}

vtune_lib <- file.path(getwd(), ".vtune-lib")
.libPaths(c(vtune_lib, .libPaths()))
library(MkPrime)
library(ape)

# ---- Load Sun2018 hyolith data from TreeSearch (all 3 partition types) ------
nex_file <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
if (!nzchar(nex_file)) {
  stop("TreeSearch package required for Sun2018.nex dataset")
}

pd  <- TreeTools::ReadAsPhyDat(nex_file)
mkd <- MkPrimeData(pd)
mdl <- MkPrimeModel()

# Config-specific MCMC settings.
# Gibbs/weighted moves are ~7-60x more expensive per iteration (M-102).
# Iteration counts calibrated for ~30s wall time at -O2 on i7-10700:
#   b: 142 iter/s → 4000 iter;  c: 41 iter/s → 1200 iter
#   d:  13 iter/s →  400 iter;  e: 476 iter/s → 14000 iter
cfg <- switch(config,
  b = MkPrimeMCMC(nIter = 4000L, nRuns = 1L, nChains = 1L,
                  warmup = 0L, thin = 100L,
                  gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE),
  c = MkPrimeMCMC(nIter = 1200L, nRuns = 1L, nChains = 1L,
                  warmup = 0L, thin = 50L,
                  gibbsSpr = FALSE, gibbsSubtreeSwap = FALSE,
                  weightedSpr = TRUE),
  d = MkPrimeMCMC(nIter = 400L, nRuns = 1L, nChains = 1L,
                  warmup = 0L, thin = 20L,
                  gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE,
                  blockGibbsBranch = TRUE),
  e = MkPrimeMCMC(nIter = 14000L, nRuns = 1L, nChains = 1L,
                  warmup = 0L, thin = 200L,
                  gibbsSpr = FALSE, gibbsSubtreeSwap = FALSE,
                  weightedBranchScale = TRUE)
)

start_tree <- ape::rtree(mkd$nTip, tip.label = mkd$taxon_names,
                         rooted = FALSE)

cat(sprintf("Config %s: running %s iterations...\n", config,
            format(cfg$nIter, big.mark = ",")))
post <- RunMkPrime(mkd, start_tree, model = mdl, mcmc = cfg)
cat("Done.\n")
