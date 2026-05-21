# rodent_aware_timing.R -----------------------------------------------------
# Quick wall-clock timing of MkNT aware vs blind MCMC on the rodent dataset
# (post-T-013 local patches). Goal: catch the catastrophic slowdown that
# stalled Hamilton job 17254914 (~17 iter/h vs blind ~500k iter/h).
#
# Single-chain, serial (nCore=1, nRuns=1) — so we just measure raw per-iter
# cost without parallel-runs overhead.

suppressPackageStartupMessages({
  library(MkPrime)
  library(TreeTools)
})

nIterShort <- 500L     # tiny — just want iter/sec
nexFile    <- "inst/ecology/data/rodent-X24848.nex"
stopifnot(file.exists(nexFile))

# --- Data prep mirroring run_rodent_MkNT.R ---------------------------------
mat <- TreeTools::ReadCharacters(nexFile)
ecoVec  <- mat[, 220]
extant  <- mat[, 221]
keepTaxa <- which(extant == "0")
matKeep <- mat[keepTaxa, -221, drop = FALSE]
ecoKeep <- ecoVec[keepTaxa]
poly <- grepl("^\\(", ecoKeep)
if (any(poly)) {
  recoded <- substr(sub("^\\(", "", ecoKeep[poly]), 1, 1)
  ecoKeep[poly] <- recoded
  matKeep[poly, 220] <- recoded
}
hasEco <- !is.na(ecoKeep) & ecoKeep != "?" &
            !is.na(suppressWarnings(as.integer(ecoKeep)))
matKeep <- matKeep[hasEco, , drop = FALSE]
ecoKeep <- ecoKeep[hasEco]

charMat <- matKeep[, -220, drop = FALSE]
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
stopifnot(all(mkd$type %in% c("neomorphic", "known")))
cat("MkNT data prep OK. nTip=", mkd$nTip, " nChar=", mkd$nChar,
    " kEcology=", mkd$kEcology, "\n", sep = "")

set.seed(20260521)
startTree <- Preorder(TreeSearch::AdditionTree(mkd$phyDat))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

# --- Blind timing -----------------------------------------------------------
modelBlind <- MkPrimeModel(ecologyAware = FALSE,
                            rateLossMeanlog = 0, rateLossSdlog = 2,
                            rateLogSdShape  = 1, rateLogSdRate = 1,
                            rateNeoMeanlog  = 0, rateNeoSdlog  = 1)
mcmcShort <- MkPrimeMCMC(nIter = nIterShort, nChains = 1L, nRuns = 1L,
                          nCore = 1L,
                          thin = 1L, treeThin = 1L,
                          minWarmup = 100L, maxWarmup = 200L,
                          logFile = tempfile(fileext = ".log"),
                          checkpointFile = tempfile(fileext = ".ckp"))

cat("\n=== BLIND timing (nIter=", nIterShort, ", serial, nChains=1) ===\n", sep="")
tB <- system.time(
  resB <- RunMkPrime(mkd, tree = startTree, model = modelBlind, mcmc = mcmcShort)
)
cat("blind elapsed:", round(tB["elapsed"], 2), "s; rate:",
    round(nIterShort / tB["elapsed"], 1), "iter/s\n")

# --- Aware timing -----------------------------------------------------------
modelAware <- MkPrimeModel(ecologyAware  = TRUE,
                            magnitudeMode = "global",
                            rho0Alpha = 7, rho0Beta = 3,
                            thetaAlpha = 2, thetaBeta = 2,
                            sigmaPhi   = 1.5,
                            rateLossMeanlog = 0, rateLossSdlog = 2,
                            rateLogSdShape  = 1, rateLogSdRate = 1,
                            rateNeoMeanlog  = 0, rateNeoSdlog  = 1)
mcmcShort$logFile <- tempfile(fileext = ".log")
mcmcShort$checkpointFile <- tempfile(fileext = ".ckp")

cat("\n=== AWARE timing (nIter=", nIterShort, ", serial, nChains=1) ===\n", sep="")
tA <- system.time(
  resA <- RunMkPrime(mkd, tree = startTree, model = modelAware, mcmc = mcmcShort)
)
cat("aware elapsed:", round(tA["elapsed"], 2), "s; rate:",
    round(nIterShort / tA["elapsed"], 1), "iter/s\n")

cat("\n=== Ratio aware/blind:",
    round(tA["elapsed"] / tB["elapsed"], 2), "x ===\n")
cat("(Hamilton stalled run was ~30,000x; <50x is acceptable.)\n")
