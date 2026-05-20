# run_rodent_aware_v4_validate.R --------------------------------------------
# Ecology-AWARE Mk' MCMC on the rodent dataset — v4 VALIDATION run.
#
# Purpose: first end-to-end run on rodent under the post-merge code, after
# the streaming bug fixes in commit 326bb13 (atomic flush, z_samples trim,
# RelabelEcology diagnostics). Serial, single PT run, short — confirms no
# torn rows / z_samples drift / RelabelEcology abort. RESUMABLE if it
# overruns.
#
# Config: nRuns=1, nChains=4 (PT sequential), nIter=20000.
# Output: /nobackup/pjjg18/mkp-rodent-v4/aware-validate/.
#
# After this passes cleanly, the production parallel script
# (run_rodent_aware_v4.R) runs nRuns=4, nChains=4, nCore=4.
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
nIter   <- if (length(args) >= 2) as.integer(args[2]) else 20000L

stopifnot(file.exists(nexFile))

outRoot <- file.path("/nobackup", Sys.getenv("USER"),
                     "mkp-rodent-v4", "aware-validate")
dir.create(file.path(outRoot, "results"),
           showWarnings = FALSE, recursive = TRUE)

cat("=== Rodent AWARE v4 VALIDATION (nChains=4 PT, serial) ===\n")
cat("nexFile:", nexFile, "\n")
cat("nIter:  ", nIter, "\n")
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
# Character classification — identical to v3.
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
# MCMC config — VALIDATION: nRuns=1, nChains=4 (PT), serial.
# Resumable: serial mode (nCore=1) supports ResumeMkPrime() if walltime hits.
# ---------------------------------------------------------------------------
logFile <- file.path(outRoot, "rodent-aware-v4-validate.log")
ckpFile <- file.path(outRoot, "rodent-aware-v4-validate.ckp")

mcmc <- MkPrimeMCMC(
  nIter          = nIter,
  nChains        = 4L,
  nRuns          = 1L,
  nCore          = 1L,
  heat           = 0.2,
  thin           = max(1L, nIter %/% 1000L),
  treeThin       = max(1L, nIter %/% 1000L),
  minWarmup      = 5000L,
  maxWarmup      = 20000L,
  logFile        = logFile,
  checkpointFile = ckpFile
)

cat(sprintf("MCMC: nIter=%d  nChains=%d (PT)  thin=%d  warmup=[%d,%d]\n",
            nIter, mcmc$nChains, mcmc$thin,
            mcmc$minWarmup, mcmc$maxWarmup))

if (file.exists(ckpFile)) {
  cat("Resuming from checkpoint:", ckpFile, "\n")
  t0  <- Sys.time()
  res <- ResumeMkPrime(checkpointFile = ckpFile, data = mkd)
} else {
  if (file.exists(logFile)) file.remove(logFile)
  t0  <- Sys.time()
  res <- RunMkPrime(mkd, tree = startTree, model = modelV4, mcmc = mcmc)
}
cat("Elapsed:", format(Sys.time() - t0), "\n")

# ---------------------------------------------------------------------------
# Validation checks for the streaming fix.
# ---------------------------------------------------------------------------
cat("\n=== Streaming-fix validation ===\n")
cat("nrow(samples): ", nrow(res$samples), "\n", sep = "")
if (!is.null(res$z_samples)) {
  cat("length(z_samples): ", length(res$z_samples), "\n", sep = "")
  cat("z_samples aligned with samples: ",
      length(res$z_samples) == nrow(res$samples), "\n", sep = "")
}

relabRes <- tryCatch(
  RelabelEcology(res),
  error = function(e) {
    cat("RelabelEcology FAILED:\n"); print(e); NULL
  }
)
if (!is.null(relabRes)) {
  cat("RelabelEcology() succeeded.\n")
}

# ---------------------------------------------------------------------------
# Save result.
# ---------------------------------------------------------------------------
discard <- function(idx) idx[seq.int(ceiling(length(idx) / 4) + 1L,
                                      length(idx))]
samples <- ReadMkLog(logFile)
samples <- samples[discard(seq_len(nrow(samples))), , drop = FALSE]

saveRDS(
  list(mkd     = mkd,
       model   = modelV4,
       mcmc    = mcmc,
       res     = res,
       relab   = relabRes,
       samples = samples),
  file.path(outRoot, "results", "rodent-aware-v4-validate-result.rds")
)
cat("Saved result.\n\nDone.\n")
