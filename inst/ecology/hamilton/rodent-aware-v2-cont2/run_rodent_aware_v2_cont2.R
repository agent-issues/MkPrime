# run_rodent_aware_v2_cont2.R --------------------------------------------------
# Second continuation of the rodent AWARE Mk' MCMC v2 (1M iter reached,
# minESS=88, 160 trees).  Extends to 2.3M total iterations (~1.3M additional)
# to reach minESS >= 200.
#
# Strategy: the checkpoint has mcmc$nIter = 1000000L and ckp$iter == 1000000L
# (chain is exhausted).  We must PATCH mcmc$nIter to 2300000L before calling
# ResumeMkPrime so that it runs the additional ~1.3M iterations.
# thin and treeThin are left unchanged — the chain accumulates more samples
# naturally as the iteration count grows.
#
# Checkpoint: /nobackup/pjjg18/mkp-rodent-aware-v2/rodent-aware-v2.ckp
#
# Usage (via SBATCH script):
#   Rscript run_rodent_aware_v2_cont2.R <nexFile>
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

stopifnot(file.exists(nexFile))

outRoot <- file.path("/nobackup", Sys.getenv("USER"), "mkp-rodent-aware-v2")
ckpFile <- file.path(outRoot, "rodent-aware-v2.ckp")
logFile <- file.path(outRoot, "rodent-aware-v2.log")

if (!file.exists(ckpFile)) {
  stop("Checkpoint not found: ", ckpFile,
       "\nRun run_rodent_aware_v2.R (and cont) first.")
}

nIterTotal <- 2300000L  # 1M done + 1.3M continuation = 2.3M total

cat("=== Rodent AWARE v2 CONTINUATION 2 ===\n")
cat("nexFile:   ", nexFile, "\n")
cat("ckpFile:   ", ckpFile, "\n")
cat("nIterTotal:", nIterTotal, "\n")
cat("outRoot:   ", outRoot, "\n")

# ---------------------------------------------------------------------------
# Patch checkpoint: extend nIter to 2.3M so ResumeMkPrime runs the remainder.
# ResumeMkPrime reads mcmc$nIter from the checkpoint; we update it in-place.
# thin and treeThin are preserved from the checkpoint — do NOT recalculate them
# here as that would change the sampling density relative to what was used in
# the first 1M iterations.
# ---------------------------------------------------------------------------
cat("Patching checkpoint nIter to", nIterTotal, "...\n")
ckp <- readRDS(ckpFile)
oldNIter <- ckp$mcmc$nIter
cat("  Previous nIter:", oldNIter, "-> new nIter:", nIterTotal, "\n")
cat("  iter reached:  ", ckp$iter, "\n")
cat("  thin:          ", ckp$mcmc$thin, "(unchanged)\n")
cat("  treeThin:      ", ckp$mcmc$treeThin, "(unchanged)\n")
ckp$mcmc$nIter <- nIterTotal
saveRDS(ckp, ckpFile)
cat("Checkpoint patched. Remaining iterations:", nIterTotal - ckp$iter, "\n")

# ---------------------------------------------------------------------------
# Rebuild MkPrimeData — identical to run_rodent_aware_v2.R so character
# indices and taxon order align with the checkpoint state.
# ---------------------------------------------------------------------------
mat <- TreeTools::ReadCharacters(nexFile)
cat("Raw matrix:", nrow(mat), "tips x", ncol(mat), "chars\n")

ecologyCol <- 220L
extantCol  <- 221L
ecoVec  <- mat[, ecologyCol]
extant  <- mat[, extantCol]
keepTaxa <- which(extant == "0")  # 0 = extant, 1 = extinct
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
# Resume from patched checkpoint.
# ResumeMkPrime uses mcmc (incl. nIter = 2.3M) from the checkpoint.
# ---------------------------------------------------------------------------
cat("Resuming from checkpoint...\n")
t0  <- Sys.time()
res <- ResumeMkPrime(checkpointFile = ckpFile, data = mkd)
cat("Elapsed:", format(Sys.time() - t0), "\n")

# ---------------------------------------------------------------------------
# Save combined result (all samples from the log, post-burnin).
# ---------------------------------------------------------------------------
discard <- function(idx) idx[seq.int(ceiling(length(idx) / 4) + 1L,
                                      length(idx))]
samples <- ReadMkLog(logFile)
samples <- samples[discard(seq_len(nrow(samples))), , drop = FALSE]

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

saveRDS(
  list(mkd     = mkd,
       model   = modelV2,
       mcmc    = ckp$mcmc,
       res     = res,
       samples = samples),
  file.path(outRoot, "results", "rodent-aware-v2-result.rds")
)
cat("Saved result to",
    file.path(outRoot, "results", "rodent-aware-v2-result.rds"), "\n")
cat("\nDone.\n")
