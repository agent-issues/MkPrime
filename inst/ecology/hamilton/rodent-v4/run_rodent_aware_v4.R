# run_rodent_aware_v4.R -----------------------------------------------------
# Ecology-AWARE Mk' MCMC on the rodent dataset — v4 PRODUCTION.
#
# Configuration follows the v4 validation result (job 17245693):
#   - Streaming fix (commit 326bb13) confirmed working end-to-end.
#   - Cold-chain throughput ~7.5k iter/h with nChains=4 PT on rodent.
#
# This script: nRuns=4 (Rhat), nChains=4 (PT mode-trap resilience),
# nCore=4 (parallel runs via callr::r_bg). Per-run checkpoints support
# resume if walltime is hit (feat/parallel-checkpointing).
#
# Output: /nobackup/pjjg18/mkp-rodent-v4/aware/.
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
nIter   <- if (length(args) >= 2) as.integer(args[2]) else 100000L

stopifnot(file.exists(nexFile))

outRoot <- file.path("/nobackup", Sys.getenv("USER"),
                     "mkp-rodent-v4", "aware")
dir.create(file.path(outRoot, "results"),
           showWarnings = FALSE, recursive = TRUE)

cat("=== Rodent AWARE v4 PRODUCTION (nRuns=4, nChains=4 PT, nCore=4) ===\n")
cat("nexFile:", nexFile, "\n")
cat("nIter:  ", nIter, "(per run)\n")
cat("outRoot:", outRoot, "\n")

# ---------------------------------------------------------------------------
# Data loading — identical to v3.
# ---------------------------------------------------------------------------
mat <- TreeTools::ReadCharacters(nexFile)
cat("Raw matrix:", nrow(mat), "tips x", ncol(mat), "chars\n")

ecologyCol <- 220L
extantCol  <- 221L
ecoVec  <- mat[, ecologyCol]
extant  <- mat[, extantCol]
keepTaxa <- which(extant == "0")
cat("Extant taxa:", length(keepTaxa), "of", nrow(mat), "\n")

matKeep <- mat[keepTaxa, -extantCol, drop = FALSE]
ecoKeep <- ecoVec[keepTaxa]
poly <- grepl("^\\(", ecoKeep)
if (any(poly)) {
  cat("Recoding", sum(poly), "polymorphic ecology entries to first listed state:\n")
  recoded <- substr(sub("^\\(", "", ecoKeep[poly]), 1, 1)
  print(data.frame(tip = rownames(matKeep)[poly], orig = ecoKeep[poly],
                   new = recoded))
  ecoKeep[poly] <- recoded
  matKeep[poly, ecologyCol] <- recoded
}
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
# Character classification.
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

# ---------------------------------------------------------------------------
# Parsimony start tree.
# ---------------------------------------------------------------------------
set.seed(20260512)
randTree  <- Preorder(ape::rtree(mkd$nTip, tip.label = mkd$taxon_names))
psTrees   <- TreeSearch::MaximizeParsimony(mkd$phyDat, tree = randTree,
                                           verbosity = 0)
startTree <- Preorder(psTrees[[1]])
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

# ---------------------------------------------------------------------------
# Model — ecology-aware, auto-expSteps.
# ---------------------------------------------------------------------------
modelV4 <- MkPrimeModel(
  ecologyAware  = TRUE,
  magnitudeMode = "global",
  kPrimePrior   = "geometric",
  coding        = "variable",
  rho0Alpha     = 7,
  rho0Beta      = 3,
  thetaAlpha    = 2,
  thetaBeta     = 2,
  sigmaPhi      = 0.5
)

# ---------------------------------------------------------------------------
# MCMC config — PRODUCTION: 4 parallel runs, each with nChains=4 PT.
# Resumable: per-run checkpoints under feat/parallel-checkpointing.
# ---------------------------------------------------------------------------
logFile <- file.path(outRoot, "rodent-aware-v4.log")
ckpFile <- file.path(outRoot, "rodent-aware-v4.ckp")

mcmc <- MkPrimeMCMC(
  nIter          = nIter,
  nChains        = 4L,
  nRuns          = 4L,
  nCore          = 4L,
  heat           = 0.2,
  thin           = max(1L, nIter %/% 1000L),
  treeThin       = max(1L, nIter %/% 1000L),
  minWarmup      = 10000L,
  maxWarmup      = 30000L,
  logFile        = logFile,
  checkpointFile = ckpFile
)

cat(sprintf("MCMC: nIter=%d  nRuns=%d  nChains=%d (PT)  nCore=%d  thin=%d  warmup=[%d,%d]\n",
            nIter, mcmc$nRuns, mcmc$nChains, mcmc$nCore, mcmc$thin,
            mcmc$minWarmup, mcmc$maxWarmup))

if (file.exists(ckpFile) || length(Sys.glob(paste0(tools::file_path_sans_ext(ckpFile), "_*.ckp")))) {
  cat("Resuming from checkpoint(s): ", ckpFile, " (and/or per-run files)\n")
  t0  <- Sys.time()
  res <- ResumeMkPrime(checkpointFile = ckpFile, data = mkd)
} else {
  t0  <- Sys.time()
  res <- RunMkPrime(mkd, tree = startTree, model = modelV4, mcmc = mcmc)
}
cat("Elapsed:", format(Sys.time() - t0), "\n")

# ---------------------------------------------------------------------------
# Save result.
# ---------------------------------------------------------------------------
saveRDS(
  list(mkd     = mkd,
       model   = modelV4,
       mcmc    = mcmc,
       res     = res),
  file.path(outRoot, "results", "rodent-aware-v4-result.rds")
)
cat("Saved result.\n\nDone.\n")
