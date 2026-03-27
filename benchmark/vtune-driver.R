# VTune hotspot driver for MkPrime
#
# Exercises the C++ MCMC hot path for ~30 seconds of CPU time.
# Run via:
#
#   vtune -collect hotspots -result-dir vtune-out -- Rscript benchmark/vtune-driver.R
#
# Or for software sampling (no driver needed):
#
#   vtune -collect hotspots -knob sampling-mode=sw -result-dir vtune-out -- Rscript benchmark/vtune-driver.R
#
# Filter to package DLL:
#   vtune -report hotspots -result-dir vtune-out -filter "module=mkp.dll"
#
# Build for profiling first (see PROFILING.md):
#   1. Add -g -fno-omit-frame-pointer to src/Makevars.win PKG_CXXFLAGS
#   2. Run the install block below with MAKEFLAGS override
#   3. Remove src/Makevars.win after profiling

# ---- Install profiling build (run once, separately) -------------------------
# vtune_lib <- file.path(getwd(), ".vtune-lib")
# old_mf <- Sys.getenv("MAKEFLAGS")
# Sys.setenv(MAKEFLAGS = "DLLFLAGS=-static-libgcc")
# install.packages(".", lib = vtune_lib, repos = NULL, type = "source",
#                  INSTALL_opts = "--no-multiarch")
# Sys.setenv(MAKEFLAGS = old_mf)
# -----------------------------------------------------------------------------

vtune_lib <- file.path(getwd(), ".vtune-lib")
.libPaths(c(vtune_lib, .libPaths()))
library(mkp)
library(ape)

# ---- Load hyoliths data (exercises all three partition types) ----------------
nex_file <- file.path("vignettes", "hyoliths.nex")
if (!file.exists(nex_file)) {
  stop("Run from mkp/ root; vignettes/hyoliths.nex not found")
}

pd  <- ape::read.nexus.data(nex_file)
mkd <- MkPrimeData(pd, outgroup = "Triplicatella_dentifera")
mdl <- MkPrimeModel(mkd, rateLogSd = 0.5)
cfg <- MkPrimeMCMC(nIter = 1L, batchSize = 50000L, nChains = 1L,
                   warmup = 0L, thin = 500L)

# Use a fixed start tree
start_tree <- ape::rtree(length(mkd$tip_labels), tip.label = mkd$tip_labels,
                         rooted = FALSE)

cat("Warming up...\n")
post <- RunMkPrime(mkd, mdl, cfg, startTree = start_tree)

cat("Profiling hot path (50 000 iterations)...\n")
# Re-run with fresh state so VTune captures steady-state behaviour
cfg2 <- MkPrimeMCMC(nIter = 1L, batchSize = 200000L, nChains = 1L,
                    warmup = 0L, thin = 1000L)
post2 <- RunMkPrime(mkd, mdl, cfg2, startTree = start_tree)
cat("Done.\n")
