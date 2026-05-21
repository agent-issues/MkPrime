# run_rodent_blind_v2.R -------------------------------------------------------
# BLIND Mk' MCMC on the rodent dataset — paired with rodent-ecology-v2.R
# (the ecology-AWARE v2 run).
#
# Mirrors rodent-ecology-v2.R exactly in:
#   - data loading, taxon filtering, character classification
#   - parsimony start tree (same set.seed(20260512))
#   - MCMC config: thin = nIter %/% 1000L, treeThin = nIter %/% 500L,
#     minWarmup = 10000L, maxWarmup = 50000L
#   - kPrimePrior = "geometric", coding = "variable", expSteps = 10
#
# Differences from aware v2 (intentional):
#   - ecologyAware = FALSE  (ecology-blind model)
#   - magnitudeMode dropped (N/A for blind)
#   - rho0Alpha/Beta, thetaAlpha/Beta, sigmaPhi dropped (ecology-only params)
#   - nIter = 200000L  (pilot run; continue if ESS < 200)
#   - Output paths: /nobackup/pjjg18/mkp-rodent-blind-v2/ (Hamilton convention)
#   - Nexus path accepted via command-line arg (no local ~/downloads assumption)
#
# Usage (via SBATCH script):
#   Rscript run_rodent_blind_v2.R <nexFile> [nIter]
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
nIter   <- if (length(args) >= 2) as.integer(args[2]) else 200000L

stopifnot(file.exists(nexFile))

outRoot <- file.path("/nobackup", Sys.getenv("USER"), "mkp-rodent-blind-v2")
dir.create(file.path(outRoot, "results"), showWarnings = FALSE, recursive = TRUE)

cat("=== Rodent BLIND v2 MCMC ===\n")
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
# Drop any remaining missing-ecology tips.
hasEco <- !is.na(ecoKeep) & ecoKeep != "?" &
          !is.na(suppressWarnings(as.integer(ecoKeep)))
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
# Character classification — identical to aware v2.
# Build MkPrimeData WITHOUT ecology (blind model does not use it).
# ---------------------------------------------------------------------------
pdForDetect <- MatrixToPhyDat(matKeep[, -ecologyCol, drop = FALSE])
neoIdx <- AutoDetectNeomorphic(pdForDetect)
cat("Neomorphic chars:", length(neoIdx), " | transformational:",
    ncol(matKeep) - 1L - length(neoIdx), "\n")

# Blind: pass no ecology argument — same phyDat, same neomorphic classification.
mkd <- MkPrimeData(pdForDetect, neomorphic = neoIdx)
cat("MkPrimeData built:",
    " nTip=", mkd$nTip, " nChar=", mkd$nChar, "\n", sep = "")

# ---------------------------------------------------------------------------
# Parsimony start tree — same seed as aware v2.
# ---------------------------------------------------------------------------
set.seed(20260512)
randTree  <- Preorder(ape::rtree(mkd$nTip, tip.label = mkd$taxon_names))
startTree <- Preorder(TreeSearch::AdditionTree(mkd$phyDat))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

# ---------------------------------------------------------------------------
# Model — blind, matching non-ecology parameters of aware v2.
# ---------------------------------------------------------------------------
modelBlind <- MkPrimeModel(
  ecologyAware = FALSE,
  kPrimePrior  = "geometric",
  coding       = "variable",
  expSteps     = 10
)

# ---------------------------------------------------------------------------
# MCMC config — same formulas as aware v2, thin/treeThin scale with nIter.
# ---------------------------------------------------------------------------
logFile <- file.path(outRoot, "rodent-blind-v2.log")
ckpFile <- file.path(outRoot, "rodent-blind-v2.ckp")

mcmc <- MkPrimeMCMC(
  nIter          = nIter,
  nChains        = 1L,
  nRuns          = 1L,
  thin           = max(1L, nIter %/% 1000L),
  treeThin       = max(1L, nIter %/% 500L),
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
  res <- ResumeMkPrime(mkd, mcmc = mcmc, model = modelBlind)
} else {
  # Clean stale log if it exists without a checkpoint.
  if (file.exists(logFile)) file.remove(logFile)
  t0  <- Sys.time()
  res <- RunMkPrime(mkd, tree = startTree, model = modelBlind, mcmc = mcmc)
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
  list(mkd     = mkd,
       model   = modelBlind,
       mcmc    = mcmc,
       res     = res,
       samples = samples),
  file.path(outRoot, "results", "rodent-blind-v2-result.rds")
)
cat("Saved result to", file.path(outRoot, "results",
                                  "rodent-blind-v2-result.rds"), "\n")
cat("\nDone.\n")
