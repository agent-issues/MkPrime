# run_rodent_MkNT.R ---------------------------------------------------------
# Rodent inference under the **validated MkNT model**:
#   - Neomorphic chars: F81-type (asymmetric binary, rate_loss + rate_neo)
#   - Transformational chars: standard Mk with k = kObs (no k' inference)
#
# Toggled aware/blind by `MODE` env var or 2nd CLI arg ("aware" | "blind").
# Same MkNT likelihood under both; only ecology layer differs.
#
# Mk' machinery (kPrimePrior, kprime_alpha/beta, etc.) is bypassed because
# every non-neomorphic character is registered as "known" with known_k = kObs
# via MkPrimeData(knownStates = ...). The result has no kPrime_* log columns.
#
# Usage: Rscript run_rodent_MkNT.R <nexFile> <mode> [nIter]

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
mode    <- if (length(args) >= 2) args[2] else Sys.getenv("MODE", "aware")
nIter   <- if (length(args) >= 3) as.integer(args[3]) else 100000L
stopifnot(file.exists(nexFile), mode %in% c("aware", "blind"))

outRoot <- file.path("/nobackup", Sys.getenv("USER"),
                     "mkp-rodent-MkNT-v1", mode)
dir.create(file.path(outRoot, "results"),
           showWarnings = FALSE, recursive = TRUE)

cat(sprintf("=== Rodent MkNT %s (nRuns=4 nChains=4 nCore=4) ===\n", toupper(mode)))
cat("nexFile:", nexFile, "\n")
cat("nIter:  ", nIter, "(per run)\n")
cat("outRoot:", outRoot, "\n")

# ---------------------------------------------------------------------------
# Data loading.
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

# ---------------------------------------------------------------------------
# Character classification + kObs pinning to invoke MkNT.
# ---------------------------------------------------------------------------
charMat <- matKeep[, -ecologyCol, drop = FALSE]
pdForDetect <- MatrixToPhyDat(charMat)
neoIdx <- AutoDetectNeomorphic(pdForDetect)
cat("Neomorphic:", length(neoIdx), " | non-neomorphic:",
    ncol(charMat) - length(neoIdx), "\n")

# Compute kObs per character directly from the raw character matrix in
# ORIGINAL coordinates (i.e. before MkPrimeData drops invariants), so the
# knownStates names match the columns MkPrimeData receives. MkPrimeData
# remaps both neomorphic and knownStates indices internally after dropping
# invariants (see MkPrimeData.R:170-193).
.kObsCol <- function(col) {
  vals <- col[!(col %in% c("?", "-", NA))]
  # split polymorphic codes "(01)" → individual states
  vals <- unlist(strsplit(gsub("[()]", "", vals), ""))
  length(unique(vals))
}
kObsRaw <- vapply(seq_len(ncol(charMat)),
                  function(j) .kObsCol(charMat[, j]),
                  integer(1L))
nonNeoOriginal <- setdiff(seq_len(ncol(charMat)), neoIdx)
# Only pin the variable ones (knownStates must be >= kObs >= 2);
# invariants will be dropped by MkPrimeData anyway.
nonNeoVariable <- nonNeoOriginal[kObsRaw[nonNeoOriginal] >= 2L]
cat("Non-neomorphic, variable: ", length(nonNeoVariable),
    " (invariants to be dropped: ",
    length(nonNeoOriginal) - length(nonNeoVariable), ")\n", sep = "")
cat("Non-neomorphic kObs distribution (variable chars):\n")
print(table(kObsRaw[nonNeoVariable]))

knownK <- setNames(as.integer(kObsRaw[nonNeoVariable]),
                   as.character(nonNeoVariable))
mkd <- MkPrimeData(pdForDetect,
                   neomorphic  = neoIdx,
                   knownStates = knownK,
                   ecology     = setNames(as.integer(ecoKeep),
                                          rownames(matKeep)))
cat("MkPrimeData rebuilt as MkNT: nTip=", mkd$nTip, " nChar=", mkd$nChar,
    " kEcology=", mkd$kEcology, "\n", sep = "")
stopifnot(all(mkd$type %in% c("neomorphic", "known")))
cat("type counts:\n"); print(table(mkd$type))

# ---------------------------------------------------------------------------
# Parsimony start tree.
# ---------------------------------------------------------------------------
set.seed(20260521)
randTree  <- Preorder(ape::rtree(mkd$nTip, tip.label = mkd$taxon_names))
psTrees   <- TreeSearch::MaximizeParsimony(mkd$phyDat, tree = randTree,
                                           ratchIter = 5, verbosity = 0)
startTree <- Preorder(psTrees[[1]])
startTree$edge.length <- rep(0.1, nrow(startTree$edge))
cat("Parsimony start score:", sum(TreeSearch::CharacterLength(startTree, mkd$phyDat)), "\n")

# ---------------------------------------------------------------------------
# Model — toggle aware/blind.
# ---------------------------------------------------------------------------
modelArgs <- list(
  # kPrimePrior is irrelevant when no transformational chars remain
  # (every non-neomorphic char is "known"); keep default.
  rateLossMeanlog = 0, rateLossSdlog = 2,
  rateLogSdShape  = 1, rateLogSdRate = 1,
  rateNeoMeanlog  = 0, rateNeoSdlog  = 1
)
if (mode == "aware") {
  modelArgs <- c(modelArgs, list(
    ecologyAware  = TRUE,
    magnitudeMode = "global",
    rho0Alpha     = 7, rho0Beta = 3,    # mode 0.7 sparsity
    thetaAlpha    = 2, thetaBeta = 2,
    sigmaPhi      = 1.5
  ))
} else {
  modelArgs <- c(modelArgs, list(ecologyAware = FALSE))
}
modelMkNT <- do.call(MkPrimeModel, modelArgs)

# ---------------------------------------------------------------------------
# MCMC config — 4 parallel runs, each with nChains=4 PT.
# Resumable via per-run checkpoints (feat/parallel-checkpointing).
# ---------------------------------------------------------------------------
logFile <- file.path(outRoot, sprintf("rodent-MkNT-%s.log", mode))
ckpFile <- file.path(outRoot, sprintf("rodent-MkNT-%s.ckp", mode))

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
cat(sprintf("MCMC: nIter=%d nRuns=%d nChains=%d (PT) nCore=%d thin=%d\n",
            nIter, mcmc$nRuns, mcmc$nChains, mcmc$nCore, mcmc$thin))

resume <- file.exists(ckpFile) ||
          length(Sys.glob(paste0(tools::file_path_sans_ext(ckpFile), "_*.ckp")))
t0 <- Sys.time()
if (resume) {
  cat("Resuming from existing checkpoint(s).\n")
  res <- ResumeMkPrime(checkpointFile = ckpFile, data = mkd)
} else {
  res <- RunMkPrime(mkd, tree = startTree, model = modelMkNT, mcmc = mcmc)
}
cat("Elapsed:", format(Sys.time() - t0), "\n")

saveRDS(list(mkd = mkd, model = modelMkNT, mcmc = mcmc, res = res),
        file.path(outRoot, "results",
                  sprintf("rodent-MkNT-%s-result.rds", mode)))
cat("Saved result.\n\nDone.\n")
