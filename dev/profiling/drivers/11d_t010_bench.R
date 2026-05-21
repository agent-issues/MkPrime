# 11d_t010_bench.R -----------------------------------------------------------
# T-010 verification: run aware + blind rodent runs against the T-010 build,
# report wall time, ratio, and the eco drift counter from a 1000-iter chain.
suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKPRIME_LIBPATH", "dev/profiling/.vtune-lib-t010")
  library("MkPrime", lib.loc = libPath)
  library("TreeTools")
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
cat(sprintf("MkPrimeData: nTip=%d nChar=%d kEcology=%d\n",
            mkd$nTip, mkd$nChar, mkd$kEcology))

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

nIter <- 200L
mcmcArgs <- list(nIter = nIter, nChains = 1L, nRuns = 1L, nCore = 1L,
                 thin = nIter, treeThin = nIter,
                 minWarmup = 10L, maxWarmup = 20L)
mcmc <- do.call(MkPrimeMCMC, mcmcArgs)

cat("\n=== BLIND run (200 iter) ===\n")
tB <- proc.time()
resB <- RunMkPrime(mkd, tree = startTree, model = modelBlind, mcmc = mcmc)
blindElapsed <- (proc.time() - tB)["elapsed"]
cat(sprintf("Blind elapsed: %.2f s (%.2f iter/s)\n",
            blindElapsed, nIter / blindElapsed))

cat("\n=== AWARE run (200 iter) ===\n")
tA <- proc.time()
resA <- RunMkPrime(mkd, tree = startTree, model = modelAware, mcmc = mcmc)
awareElapsed <- (proc.time() - tA)["elapsed"]
cat(sprintf("Aware elapsed: %.2f s (%.2f iter/s)\n",
            awareElapsed, nIter / awareElapsed))
cat(sprintf(">>> Ratio aware/blind = %.2fx\n",
            awareElapsed / blindElapsed))

cat("\n=== AWARE drift check (1000 iter) ===\n")
mcmcDrift <- do.call(MkPrimeMCMC, list(
  nIter = 1000L, nChains = 1L, nRuns = 1L, nCore = 1L,
  thin = 1000L, treeThin = 1000L,
  minWarmup = 10L, maxWarmup = 20L))
tD <- proc.time()
resD <- RunMkPrime(mkd, tree = startTree, model = modelAware, mcmc = mcmcDrift)
driftElapsed <- (proc.time() - tD)["elapsed"]
cat(sprintf("Aware 1000-iter elapsed: %.2f s\n", driftElapsed))
# Attempt to extract the diagnostic counter from result list. Names vary by
# version; print whatever drift information is available.
if (!is.null(resD$diag$drift)) {
  cat(sprintf("diagDriftCount = %d (target: 0)\n", resD$diag$drift))
} else if (!is.null(resD$drift)) {
  cat(sprintf("drift = %d (target: 0)\n", resD$drift))
} else {
  cat("drift counter not in result; check stderr above for [eco-resync] lines\n")
}
