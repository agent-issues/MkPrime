# run_rodent_blind_v3.R -------------------------------------------------------
# BLIND Mk' MCMC on the rodent dataset — v3 (auto-expSteps from parsimony).
#
# Same as run_rodent_blind_v2.R except:
#   - REMOVED `expSteps = 10` from MkPrimeModel() call so the new
#     parsimony-based default (commit bc24d52) kicks in.
#   - nIter default raised to 1,000,000 (was 200,000); resume if ESS < 200.
#   - Output paths under /nobackup/pjjg18/mkp-rodent-v3/blind/.
#
# Usage (via SBATCH script):
#   Rscript run_rodent_blind_v3.R <nexFile> [nIter]
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

outRoot <- file.path("/nobackup", Sys.getenv("USER"), "mkp-rodent-v3", "blind")
dir.create(file.path(outRoot, "results"), showWarnings = FALSE, recursive = TRUE)

cat("=== Rodent BLIND v3 MCMC (auto-expSteps) ===\n")
cat("nexFile:", nexFile, "\n")
cat("nIter:  ", nIter, "\n")
cat("outRoot:", outRoot, "\n")

# ---------------------------------------------------------------------------
# Data loading — identical to v2.
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

# ---------------------------------------------------------------------------
# Character classification — identical to v2.
# ---------------------------------------------------------------------------
pdForDetect <- MatrixToPhyDat(matKeep[, -ecologyCol, drop = FALSE])
neoIdx <- AutoDetectNeomorphic(pdForDetect)
cat("Neomorphic chars:", length(neoIdx), " | transformational:",
    ncol(matKeep) - 1L - length(neoIdx), "\n")

mkd <- MkPrimeData(pdForDetect, neomorphic = neoIdx)
cat("MkPrimeData built:",
    " nTip=", mkd$nTip, " nChar=", mkd$nChar, "\n", sep = "")

# ---------------------------------------------------------------------------
# Parsimony start tree — same seed as v2.
# ---------------------------------------------------------------------------
set.seed(20260512)
randTree  <- Preorder(ape::rtree(mkd$nTip, tip.label = mkd$taxon_names))
psTrees   <- TreeSearch::MaximizeParsimony(mkd$phyDat, tree = randTree,
                                           verbosity = 0)
startTree <- Preorder(psTrees[[1]])
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

# ---------------------------------------------------------------------------
# Model — blind, NO expSteps (parsimony-based default applies).
# ---------------------------------------------------------------------------
cli::cli_alert_info("Using new auto-expSteps default (parsimony-based prior)")
modelBlind <- MkPrimeModel(
  ecologyAware = FALSE,
  kPrimePrior  = "geometric",
  coding       = "variable"
  # expSteps intentionally omitted -> auto from parsimony score x 1.05
)

# ---------------------------------------------------------------------------
# MCMC config.
# ---------------------------------------------------------------------------
logFile <- file.path(outRoot, "rodent-blind-v3.log")
ckpFile <- file.path(outRoot, "rodent-blind-v3.ckp")

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
  res <- ResumeMkPrime(mkd, mcmc = mcmc, model = modelBlind)
} else {
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
  file.path(outRoot, "results", "rodent-blind-v3-result.rds")
)
cat("Saved result to", file.path(outRoot, "results",
                                  "rodent-blind-v3-result.rds"), "\n")
cat("\nDone.\n")
