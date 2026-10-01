#!/usr/bin/env Rscript
# M-131 warmup stabilisation validation: single dataset x seed run.
#
# Usage:
#   DATASET=Sun2018 SEED=2847 OUTDIR=results Rscript m131-warmup-validation.R
#
# Runs MkPrime for NITER (default 200k) iterations, all of them warmup
# (minWarmup = maxWarmup = nIter, autoTune off), so the stabilisation detector
# cannot end warmup early. The logP snapshot trace is kept for post-hoc offline
# replay of the detector.

.start <- Sys.time()

library(MkPrime)
library(TreeTools)

dataset_name <- Sys.getenv("DATASET")
seed         <- as.integer(Sys.getenv("SEED"))
outdir       <- Sys.getenv("OUTDIR", "results")

if (nchar(dataset_name) == 0L || is.na(seed)) {
  stop("Set DATASET and SEED environment variables.")
}

nex_file <- file.path("data-raw", paste0(dataset_name, ".nex"))
if (!file.exists(nex_file)) {
  stop("Nexus file not found: ", nex_file)
}

cat(sprintf("M-131 validation: %s, seed %d\n", dataset_name, seed))

dat <- ReadAsPhyDat(nex_file)
mkd <- MkPrimeData(dat)

nIter <- as.integer(Sys.getenv("NITER", "200000"))

# A temp log file avoids accumulating large streaming files; only the
# warmup_trace on the returned posterior is needed.
mcmc_val <- MkPrimeMCMC(
  nIter       = nIter,
  thin        = 50L,
  maxWarmup   = nIter,
  minWarmup   = nIter,
  autoTune    = FALSE,
  nRuns       = 1L,
  nChains     = 1L,
  logFile     = tempfile()
)

set.seed(seed)

posterior <- RunMkPrime(mkd, mcmc = mcmc_val)

nTip  <- length(mkd$taxon_names)
nEdge <- 2L * nTip - 3L

result <- list(
  dataset       = dataset_name,
  seed          = seed,
  nTip          = nTip,
  nEdge         = nEdge,
  nChar         = mkd$nChar,
  warmup_trace  = posterior$warmup_trace[[1]],
  wall_time_sec = as.numeric(difftime(Sys.time(), .start, units = "secs"))
)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(outdir, sprintf("m131_%s_seed%d.rds",
                                       dataset_name, seed))
saveRDS(result, out_file)
cat(sprintf("Saved: %s (%.1f s)\n", out_file, result$wall_time_sec))
