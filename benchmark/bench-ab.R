# A/B benchmark: compare two builds of mkp head-to-head.
#
# Avoids DLL-lock issues (Windows) and environmental drift (cross-session)
# by loading both package copies in the same subprocess.
#
# Usage (from mkp/ root):
#   Rscript benchmark/bench-ab.R
#
# Prerequisites:
#   1. Install the reference build:
#        ref_lib <- ".bench-lib/ref"
#        stale <- Sys.glob(file.path("src", c("*.o", "*.dll", "*.so")))
#        if (length(stale)) file.remove(stale)
#        install.packages(".", lib = ref_lib, repos = NULL, type = "source",
#                         INSTALL_opts = "--no-multiarch")
#
#   2. Make your changes to src/, then install the dev build:
#        dev_lib <- ".bench-lib/dev"
#        stale <- Sys.glob(file.path("src", c("*.o", "*.dll", "*.so")))
#        if (length(stale)) file.remove(stale)
#        install.packages(".", lib = dev_lib, repos = NULL, type = "source",
#                         INSTALL_opts = "--no-multiarch")
#
#   3. Rename one copy (so both can be loaded simultaneously) — see
#      PROFILING.md §A/B renaming recipe.
#
# Because both mkp variants have the same package name, you must rename one.
# The rename_and_install() helper in PROFILING.md handles this.

ref_lib <- ".bench-lib/ref"
dev_lib <- ".bench-lib/dev"
.libPaths(c(dev_lib, ref_lib, .libPaths()))

library(mkpRef)   # reference build (renamed)
library(mkpDev)   # dev build (renamed)
library(ape)
library(bench)

# ---- Shared test data -------------------------------------------------------
nex_file <- file.path("vignettes", "hyoliths.nex")
if (!file.exists(nex_file)) stop("Run from mkp/ root")

pd  <- ape::read.nexus.data(nex_file)

# Reference mkd/mdl/cfg
mkd_ref <- mkpRef::MkPrimeData(pd, outgroup = "Triplicatella_dentifera")
mdl_ref <- mkpRef::MkPrimeModel(mkd_ref, rateLogSd = 0.5)
cfg_ref  <- mkpRef::MkPrimeMCMC(nIter = 1L, batchSize = 5000L, nChains = 1L,
                                  warmup = 0L, thin = 500L)
start_tree <- ape::rtree(length(mkd_ref$tip_labels),
                          tip.label = mkd_ref$tip_labels, rooted = FALSE)

# Dev mkd/mdl/cfg
mkd_dev <- mkpDev::MkPrimeData(pd, outgroup = "Triplicatella_dentifera")
mdl_dev <- mkpDev::MkPrimeModel(mkd_dev, rateLogSd = 0.5)
cfg_dev  <- mkpDev::MkPrimeMCMC(nIter = 1L, batchSize = 5000L, nChains = 1L,
                                  warmup = 0L, thin = 500L)

# ---- Canary: something that should NOT change --------------------------------
# (If canary moves > 3%, discard the run — environmental drift)
canary <- bench::mark(
  ref = { mkpRef::RunMkPrime(mkd_ref, mdl_ref, cfg_ref, startTree = start_tree) },
  dev = { mkpDev::RunMkPrime(mkd_dev, mdl_dev, cfg_dev, startTree = start_tree) },
  min_iterations = 3L,
  check = FALSE
)
cat("=== Canary (same code, should be <3% difference) ===\n")
print(canary[, c("expression", "min", "median", "n_itr")])

# ---- Primary benchmark: full MCMC batch -------------------------------------
cfg_long_ref <- mkpRef::MkPrimeMCMC(nIter = 1L, batchSize = 20000L,
                                      nChains = 1L, warmup = 0L, thin = 1000L)
cfg_long_dev <- mkpDev::MkPrimeMCMC(nIter = 1L, batchSize = 20000L,
                                      nChains = 1L, warmup = 0L, thin = 1000L)

b <- bench::mark(
  ref = { mkpRef::RunMkPrime(mkd_ref, mdl_ref, cfg_long_ref, startTree = start_tree) },
  dev = { mkpDev::RunMkPrime(mkd_dev, mdl_dev, cfg_long_dev, startTree = start_tree) },
  min_iterations = 5L,
  check = FALSE
)
cat("\n=== Primary benchmark: 20 000 MCMC iterations ===\n")
print(b[, c("expression", "min", "median", "n_itr")])
cat(sprintf("Speedup (ref/dev median): %.2f×\n",
            as.numeric(b$median[1]) / as.numeric(b$median[2])))
