# 11_aware_vs_blind_rodent.R ------------------------------------------------
# Profile the ecology-aware path on the rodent matrix vs the blind path.
# Aware reportedly runs ~4-5x slower for the same iteration count.
#
# Driver runs both modes for a short iteration count (deterministic), records
# wall time, and dumps a profvis HTML for the aware run so we can confirm
# where the cost goes. Loads the installed binary, not load_all (per skill).
#
# bare target: ~30-60 s for aware (small nIter), ~10 s for blind.

suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKPRIME_LIBPATH", "dev/profiling/.vtune-lib-prof")
  library("MkPrime", lib.loc = libPath)
  library("TreeTools")
  library("profvis")
})

set.seed(20260521)

nexFile <- "C:/Users/pjjg18/downloads/mbank_X24848_2026-5-9-1135.nex"
stopifnot(file.exists(nexFile))

mat <- TreeTools::ReadCharacters(nexFile)
ecologyCol <- 220L; extantCol <- 221L
ecoVec <- mat[, ecologyCol]; extant <- mat[, extantCol]
keepTaxa <- which(extant == "0")
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
matKeep <- matKeep[hasEco, , drop = FALSE]; ecoKeep <- ecoKeep[hasEco]
charMat <- matKeep[, -ecologyCol, drop = FALSE]
pdForDetect <- MatrixToPhyDat(charMat)
neoIdx <- AutoDetectNeomorphic(pdForDetect)
.kObsCol <- function(col) {
  vals <- col[!(col %in% c("?", "-", NA))]
  vals <- unlist(strsplit(gsub("[()]", "", vals), ""))
  length(unique(vals))
}
kObsRaw <- vapply(seq_len(ncol(charMat)),
                  function(j) .kObsCol(charMat[, j]), integer(1L))
nonNeoOriginal <- setdiff(seq_len(ncol(charMat)), neoIdx)
nonNeoVariable <- nonNeoOriginal[kObsRaw[nonNeoOriginal] >= 2L]
knownK <- setNames(as.integer(kObsRaw[nonNeoVariable]),
                   as.character(nonNeoVariable))
mkd <- MkPrimeData(pdForDetect,
                   neomorphic  = neoIdx,
                   knownStates = knownK,
                   ecology     = setNames(as.integer(ecoKeep),
                                          rownames(matKeep)))
cat(sprintf("MkPrimeData: nTip=%d nChar=%d kEcology=%d (neo=%d known=%d)\n",
            mkd$nTip, mkd$nChar, mkd$kEcology,
            sum(mkd$type == "neomorphic"), sum(mkd$type == "known")))

# Greedy stepwise-addition parsimony start tree — instant, and gives a
# methodologically correct (parsimony-based) starting point for the chain.
# ratchIter parsimony was burning ~30 s per run via the slow Morphy
# delegate; AdditionTree is essentially free.
startTree <- Preorder(TreeSearch::AdditionTree(mkd$phyDat))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

baseModelArgs <- list(
  rateLossMeanlog = 0, rateLossSdlog = 2,
  rateLogSdShape  = 1, rateLogSdRate = 1,
  rateNeoMeanlog  = 0, rateNeoSdlog  = 1
)
modelBlind <- do.call(MkPrimeModel, c(baseModelArgs, list(ecologyAware = FALSE)))
modelAware <- do.call(MkPrimeModel, c(baseModelArgs, list(
  ecologyAware  = TRUE, magnitudeMode = "global",
  rho0Alpha = 7, rho0Beta = 3, thetaAlpha = 2, thetaBeta = 2, sigmaPhi = 1.5
)))

# Small iter count — just enough to see per-iter cost ratio.
nIter <- 200L
mcmcArgs <- list(nIter = nIter, nChains = 1L, nRuns = 1L, nCore = 1L,
                 thin = nIter, treeThin = nIter,
                 minWarmup = 10L, maxWarmup = 20L)
mcmc <- do.call(MkPrimeMCMC, mcmcArgs)

cat("\n=== BLIND run ===\n")
tB <- proc.time()
resB <- RunMkPrime(mkd, tree = startTree, model = modelBlind, mcmc = mcmc)
blindElapsed <- (proc.time() - tB)["elapsed"]
cat(sprintf("Blind elapsed: %.2f s for %d iter (%.2f iter/s)\n",
            blindElapsed, nIter, nIter / blindElapsed))

cat("\n=== AWARE run (profvis) ===\n")
tA <- proc.time()
p <- profvis::profvis({
  resA <- RunMkPrime(mkd, tree = startTree, model = modelAware, mcmc = mcmc)
})
awareElapsed <- (proc.time() - tA)["elapsed"]
cat(sprintf("Aware elapsed: %.2f s for %d iter (%.2f iter/s)\n",
            awareElapsed, nIter, nIter / awareElapsed))
cat(sprintf("\n>>> Ratio aware/blind = %.2fx <<<\n",
            awareElapsed / blindElapsed))

htmlOut <- "dev/profiling/drivers/11_aware_rodent-profvis.html"
htmlwidgets::saveWidget(p, htmlOut, selfcontained = TRUE)
cat("profvis HTML written to:", htmlOut, "\n")
