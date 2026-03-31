#!/usr/bin/env Rscript
# M-131 warmup stabilisation validation: single dataset x seed run.
#
# Usage:
#   DATASET=Sun2018 SEED=2847 OUTDIR=results Rscript m131-warmup-validation.R
#
# Runs MkPrime with warmup disabled (minWarmup = maxWarmup = 200000)
# for 200k iterations, capturing the logP snapshot trace for post-hoc
# offline replay of the stabilisation detector.

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

dat <- ReadCharacters(nex_file)
mkd <- MkPrimeData(dat)

mcmc_val <- MkPrimeMCMC(
  nIter       = 200000L,
  thin        = 50L,
  maxWarmup   = 200000L,
  minWarmup   = 200000L,
  autoTune    = FALSE,
  nRuns       = 1L,
  nChains     = 1L
)

set.seed(seed)

# Use a temp log file to avoid accumulating large streaming files.
# We only need the warmup_trace from the posterior object.
posterior <- RunMkPrime(mkd, mcmc = mcmc_val, logFile = tempfile())

nTip  <- length(mkd$tipLabels)
nEdge <- 2L * nTip - 3L

result <- list(
  dataset       = dataset_name,
  seed          = seed,
  nTip          = nTip,
  nEdge         = nEdge,
  nChar         = sum(lengths(mkd$partitions)),
  warmup_trace  = posterior$warmup_trace[[1]],
  wall_time_sec = as.numeric(difftime(Sys.time(), .start, units = "secs"))
)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(outdir, sprintf("m131_%s_seed%d.rds",
                                       dataset_name, seed))
saveRDS(result, out_file)
cat(sprintf("Saved: %s (%.1f s)\n", out_file, result$wall_time_sec))
