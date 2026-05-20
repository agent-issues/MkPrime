# run_rodent_aware_v3_cont.R -------------------------------------------------
# Continuation of the rodent AWARE Mk' MCMC v3 (auto-expSteps).
#
# Previous run timed out at iter ~346k / 1M with minESS=27. ResumeMkPrime
# reads mcmc$nIter=1M from the checkpoint, so will run the remaining
# ~650k iterations automatically. Multiple back-to-back continuations
# may be required; this script is safe to run repeatedly — each invocation
# resumes from the most recent checkpoint.
#
# Data rebuilt identically to run_rodent_aware_v3.R so character indices
# and taxon order align with the checkpoint state.
#
# Usage (via SBATCH script):
#   Rscript run_rodent_aware_v3_cont.R <nexFile>
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

outRoot <- file.path("/nobackup", Sys.getenv("USER"), "mkp-rodent-v3", "aware")
ckpFile <- file.path(outRoot, "rodent-aware-v3.ckp")
logFile <- file.path(outRoot, "rodent-aware-v3.log")

if (!file.exists(ckpFile)) {
  stop("Checkpoint not found: ", ckpFile,
       "\nRun the original run_rodent_aware_v3.R first.")
}

cat("=== Rodent AWARE v3 CONTINUATION ===\n")
cat("nexFile:", nexFile, "\n")
cat("ckpFile:", ckpFile, "\n")
cat("logFile:", logFile, "\n")
cat("outRoot:", outRoot, "\n")

ckp <- readRDS(ckpFile)
cat("Checkpoint nIter:", ckp$mcmc$nIter, "\n")
cat("Checkpoint iter reached:", ckp$iter, "\n")
cat("Remaining iterations:   ", ckp$mcmc$nIter - ckp$iter, "\n")

# ---------------------------------------------------------------------------
# Rebuild MkPrimeData identically to run_rodent_aware_v3.R.
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
  recoded <- substr(sub("^\\(", "", ecoKeep[poly]), 1, 1)
  ecoKeep[poly] <- recoded
  matKeep[poly, ecologyCol] <- recoded
}
hasEco <- !is.na(ecoKeep) & ecoKeep != "?" &
          !is.na(suppressWarnings(as.integer(ecoKeep)))
if (any(!hasEco)) {
  matKeep <- matKeep[hasEco, , drop = FALSE]
  ecoKeep <- ecoKeep[hasEco]
}

pdForDetect <- MatrixToPhyDat(matKeep[, -ecologyCol, drop = FALSE])
neoIdx <- AutoDetectNeomorphic(pdForDetect)

mkd <- MkPrimeData(pdForDetect,
                   ecology = setNames(as.integer(ecoKeep), rownames(matKeep)),
                   neomorphic = neoIdx)
cat("MkPrimeData built:",
    " nTip=", mkd$nTip, " nChar=", mkd$nChar,
    " kEcology=", mkd$kEcology,
    " refEcology=", mkd$refEcology, "\n", sep = "")

# ---------------------------------------------------------------------------
# Resume.
# ---------------------------------------------------------------------------
cat("Resuming from checkpoint...\n")
t0  <- Sys.time()
res <- ResumeMkPrime(checkpointFile = ckpFile, data = mkd)
cat("Elapsed:", format(Sys.time() - t0), "\n")

# ---------------------------------------------------------------------------
# Save combined result (best-effort; the chain may not have completed).
# ---------------------------------------------------------------------------
discard <- function(idx) idx[seq.int(ceiling(length(idx) / 4) + 1L,
                                      length(idx))]
samples <- tryCatch(
  ReadMkLog(logFile),
  error = function(e) {
    cat("ReadMkLog error (continuation may still be incomplete):",
        conditionMessage(e), "\n")
    NULL
  })
if (!is.null(samples) && nrow(samples) > 0L) {
  samples <- samples[discard(seq_len(nrow(samples))), , drop = FALSE]
}

modelV3 <- MkPrimeModel(
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

dir.create(file.path(outRoot, "results"), showWarnings = FALSE,
           recursive = TRUE)
saveRDS(
  list(mkd     = mkd,
       model   = modelV3,
       mcmc    = ckp$mcmc,
       res     = res,
       samples = samples),
  file.path(outRoot, "results", "rodent-aware-v3-result.rds")
)
cat("Saved result to",
    file.path(outRoot, "results", "rodent-aware-v3-result.rds"), "\n")
cat("\nDone.\n")
