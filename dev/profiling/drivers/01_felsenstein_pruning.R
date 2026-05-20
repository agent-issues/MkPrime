# Driver: area #1 — Felsenstein pruning / CL accumulation
#
# Representative workload exercising the C++ pruning hot path
# (`pruning_jc_acrv_flat`, `pruning_jc_flat`, `cpp_partition_log_likelihood`)
# through `run_mcmc_batch_cpp`.  Deterministic seed; short bare run.
#
# Per PROFILING.md (2026-03-28 VTune) the pre-M-145 hot path was:
#   pruning_jc_acrv  42.0 %  +  exp() 15.2 %  +  std::fill 6.3 %
# That baseline is now 6 weeks stale (M-145..M-172 + EG/BG priors landed
# in between); refresh before any new optimisation decision.

suppressPackageStartupMessages({
  library(MkPrime, lib.loc = "dev/profiling/.vtune-lib-LATEST")
  library(TreeTools)
})

set.seed(5813)

# Sun2018 hyoliths: the standard mkp profiling workload (54 taxa, 225 chars).
nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
if (!nzchar(nexFile) || !file.exists(nexFile)) {
  stop("Sun2018.nex not found via TreeSearch::datasets path.")
}
pd  <- TreeTools::ReadAsPhyDat(nexFile)
mkd <- MkPrimeData(pd)
cat(sprintf("Dataset: %d taxa x %d chars\n",
            nrow(mkd$matrix), ncol(mkd$matrix)))
cat(sprintf("Types:   %s\n",
            paste(names(table(mkd$type)), table(mkd$type),
                  sep = ":", collapse = "  ")))

# Deterministic NJ start tree.
tree <- TreeTools::NJTree(pd, edgeLengths = TRUE)
tree$edge.length[tree$edge.length <= 0] <- 1e-8

# Short bare run — target ~3-5 s wall on Rscript.
# nIter = 800 is enough for the C++ hot loop to dominate dispatch overhead
# but short enough to keep VTune software-sampling collection time tractable.
nIter <- 800L
profPath <- "dev/profiling/drivers/01_felsenstein_pruning.Rprof"

t0 <- proc.time()
Rprof(profPath, interval = 0.02, line.profiling = FALSE,
      memory.profiling = FALSE)
res <- RunMkPrime(
  mkd, tree,
  model = MkPrimeModel(coding = "variable",
                        kPrimePrior = "empirical_geometric"),
  mcmc  = MkPrimeMCMC(nIter      = nIter,
                       thin       = 50L,
                       maxWarmup  = 200L,
                       minWarmup  = 200L,
                       autoTune   = FALSE,
                       nRuns      = 1L,
                       nChains    = 1L,
                       progressFn = function(...) invisible())
)
Rprof(NULL)
elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("Bare wall: %.2f s for %d iter\n", elapsed, nIter))

# Rprof summary — confirms the path is dominated by .Call -> run_mcmc_batch_cpp.
prof <- summaryRprof(profPath)
cat("\n--- Rprof by.self (top 5) ---\n")
print(head(prof$by.self, 5))
cat("\n--- Rprof by.total (top 5) ---\n")
print(head(prof$by.total, 5))

cpp_share <- prof$by.self[".Call", "self.pct"]
if (!is.na(cpp_share)) {
  cat(sprintf("\nC++ self.pct: %.2f %% (Rprof attribution to .Call)\n",
              cpp_share))
}
cat("\nNote: Rprof cannot resolve below the .Call boundary; use VTune\n")
cat("(Windows) or `perf` (Linux) for C++-internal hotspots.\n")
