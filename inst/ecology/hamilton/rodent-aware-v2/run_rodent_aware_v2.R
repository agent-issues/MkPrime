# run_rodent_aware_v2.R -------------------------------------------------------
# Ecology-AWARE Mk' MCMC on the rodent dataset — Hamilton HPC version.
#
# Mirrors inst/ecology/scripts/rodent-ecology-v2.R exactly in:
#   - data loading, taxon filtering, character classification
#   - parsimony start tree (same set.seed(20260512))
#   - model spec: ecologyAware=TRUE, magnitudeMode="global",
#     kPrimePrior="geometric", coding="variable", expSteps=10
#   - ecology-aware hyperpriors: rho0Alpha=7, rho0Beta=3,
#     thetaAlpha=2, thetaBeta=2, sigmaPhi=0.5
#   - MCMC config: thin = nIter %/% 1000L, treeThin = nIter %/% 1000L,
#     minWarmup = 10000L, maxWarmup = 50000L
#
# Deviations from local v2 (intentional):
#   - nIter = 1000000L (1M; local default was 500k, still warming at 127k)
#   - treeThin = nIter %/% 1000L (task spec: ~1000 trees; local used /500)
#   - Output paths: /nobackup/pjjg18/mkp-rodent-aware-v2/ (Hamilton convention)
#   - Nexus accepted via command-line arg (reuses blind-v2 data dir)
#   - checkpointFile enabled under outRoot (no checkpoint in local v2)
#   - Library from MKP_LIB (not devtools::load_all)
#   - Post-processing (zPost, heatmap) omitted; runs locally after pulling
#
# Usage (via SBATCH script):
#   Rscript run_rodent_aware_v2.R <nexFile> [nIter]
#
suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKP_LIB",
    file.path("/nobackup", Sys.getenv("USER"),
              "mkp-sim3-multirep-v3/lib"))
  .libPaths(c(libPath, .libPaths()))
  library("MkPrime")
  library("TreeTools")
})

args    <- commandArgs(trailingOnly = TRUE)
nexFile <- if (length(args) >= 1) args[1] else
             "~/downloads/mbank_X24848_2026-5-9-1135.nex"
nIter   <- if (length(args) >= 2) as.integer(args[2]) else 1000000L

stopifnot(file.exists(nexFile))

outRoot <- file.path("/nobackup", Sys.getenv("USER"), "mkp-rodent-aware-v2")
dir.create(file.path(outRoot, "results"), showWarnings = FALSE, recursive = TRUE)

cat("=== Rodent AWARE v2 MCMC ===\n")
cat("nexFile:", nexFile, "\n")
cat("nIter:  ", nIter, "\n")
cat("outRoot:", outRoot, "\n")

# ---------------------------------------------------------------------------
# Data loading — identical to rodent-ecology-v2.R
# ---------------------------------------------------------------------------
mat <- TreeTools::ReadCharacters(nexFile)
cat("Raw matrix:", nrow(mat), "tips x", ncol(mat), "chars\n")

ecologyCol <- 220L
extantCol  <- 221L
ecoVec  <- mat[, ecologyCol]
extant  <- mat[, extantCol]
keepTaxa <- which(extant == "0")  # 0 = extant, 1 = extinct
cat("Extant taxa:", length(keepTaxa), "of", nrow(mat), "\n")

matKeep <- mat[keepTaxa, -extantCol, drop = FALSE]  # drop col 221 only
ecoKeep <- ecoVec[keepTaxa]
# Polymorphic ecology codings like (01), (03) — recode as the first
# state listed (treating polymorphism as the prevailing ecology).
poly <- grepl("^\\(", ecoKeep)
if (any(poly)) {
  cat("Recoding", sum(poly), "polymorphic ecology entries to first listed state:\n")
  recoded <- substr(sub("^\\(", "", ecoKeep[poly]), 1, 1)
  print(data.frame(tip = rownames(matKeep)[poly], orig = ecoKeep[poly],
                   new = recoded))
  ecoKeep[poly] <- recoded
  matKeep[poly, ecologyCol] <- recoded
}
# Drop any remaining missing-ecology tips (MkPrimeData ecology requires
# a complete vector; could relax by passing NA, but cleaner to subset).
hasEco <- !is.na(ecoKeep) & ecoKeep != "?" & !is.na(suppressWarnings(as.integer(ecoKeep)))
if (any(!hasEco)) {
  cat("Dropping", sum(!hasEco), "tips with missing/unparseable ecology:\n")
  print(rownames(matKeep)[!hasEco])
  matKeep <- matKeep[hasEco, , drop = FALSE]
  ecoKeep <- ecoKeep[hasEco]
}
ecoStates <- sort(unique(ecoKeep))
cat("Ecology states:", paste(ecoStates, collapse = ", "), "\n")
cat("Ecology tip counts:\n"); print(table(ecoKeep))

# ---------------------------------------------------------------------------
# Character classification — identical to rodent-ecology-v2.R.
# Build MkPrimeData with ecology col (now col 220 in the reduced matrix).
# ---------------------------------------------------------------------------
pdForDetect <- MatrixToPhyDat(matKeep[, -ecologyCol, drop = FALSE])
neoIdx <- AutoDetectNeomorphic(pdForDetect)
cat("Neomorphic chars:", length(neoIdx), " | transformational:",
    ncol(matKeep) - 1L - length(neoIdx), "\n")

mkd <- MkPrimeData(pdForDetect,
                   ecology = setNames(as.integer(ecoKeep), rownames(matKeep)),
                   neomorphic = neoIdx)
cat("MkPrimeData built:",
    " nTip=", mkd$nTip, " nChar=", mkd$nChar,
    " kEcology=", mkd$kEcology,
    " refEcology=", mkd$refEcology, "\n", sep = "")
cat("Ecology tip-fraction:\n"); print(round(mkd$ecoTipFraction, 3))

# ---------------------------------------------------------------------------
# Parsimony start tree — same seed as rodent-ecology-v2.R
# ---------------------------------------------------------------------------
set.seed(20260512)
randTree  <- Preorder(ape::rtree(mkd$nTip, tip.label = mkd$taxon_names))
psTrees   <- TreeSearch::MaximizeParsimony(mkd$phyDat, tree = randTree,
                                           verbosity = 0)
startTree <- Preorder(psTrees[[1]])
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

# ---------------------------------------------------------------------------
# Model — ecology-aware, matching rodent-ecology-v2.R exactly.
# ---------------------------------------------------------------------------
modelV2 <- MkPrimeModel(
  ecologyAware  = TRUE,
  magnitudeMode = "global",
  kPrimePrior   = "geometric",
  coding        = "variable",
  expSteps      = 10,
  rho0Alpha     = 7,
  rho0Beta      = 3,
  thetaAlpha    = 2,
  thetaBeta     = 2,
  sigmaPhi      = 0.5
)

# ---------------------------------------------------------------------------
# MCMC config — thin/treeThin scale with nIter per task spec.
# treeThin = nIter %/% 1000L (not /500 as local v2) → ~1000 trees saved.
# ---------------------------------------------------------------------------
logFile <- file.path(outRoot, "rodent-aware-v2.log")
ckpFile <- file.path(outRoot, "rodent-aware-v2.ckp")

mcmc <- MkPrimeMCMC(
  nIter          = nIter,
  nChains        = 1L,
  nRuns          = 1L,
  thin           = max(1L, nIter %/% 1000L),
  treeThin       = max(1L, nIter %/% 1000L),
  minWarmup      = 10000L,
  maxWarmup      = 50000L,
  logFile        = logFile,
  checkpointFile = ckpFile
)

cat(sprintf("MCMC: nIter=%d  thin=%d  treeThin=%d  warmup=[%d,%d]\n",
            nIter, mcmc$thin, mcmc$treeThin,
            mcmc$minWarmup, mcmc$maxWarmup))

# Resume from checkpoint if present; otherwise start fresh.
if (file.exists(ckpFile)) {
  cat("Resuming from checkpoint:", ckpFile, "\n")
  t0  <- Sys.time()
  res <- ResumeMkPrime(checkpointFile = ckpFile, data = mkd)
} else {
  # Clean stale log if it exists without a checkpoint.
  if (file.exists(logFile)) file.remove(logFile)
  t0  <- Sys.time()
  res <- RunMkPrime(mkd, tree = startTree, model = modelV2, mcmc = mcmc)
}
cat("Elapsed:", format(Sys.time() - t0), "\n")

# ---------------------------------------------------------------------------
# Save result.
# ---------------------------------------------------------------------------
discard <- function(idx) idx[seq.int(ceiling(length(idx) / 4) + 1L,
                                      length(idx))]
samples <- ReadMkLog(logFile)
samples <- samples[discard(seq_len(nrow(samples))), , drop = FALSE]

saveRDS(
  list(mkd   = mkd,
       model  = modelV2,
       mcmc   = mcmc,
       res    = res,
       samples = samples),
  file.path(outRoot, "results", "rodent-aware-v2-result.rds")
)
cat("Saved result to", file.path(outRoot, "results",
                                  "rodent-aware-v2-result.rds"), "\n")
cat("\nDone.\n")
