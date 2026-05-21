# 11e_t011_bench.R — bench T-011 against T-010 baseline.
# Mirrors 11_aware_vs_blind_rodent.R sans profvis (faster + scriptable).
suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKPRIME_LIBPATH", "dev/profiling/.vtune-lib-t011")
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
mkd <- MkPrimeData(pdForDetect, neomorphic = neoIdx, knownStates = knownK,
                   ecology = setNames(as.integer(ecoKeep),
                                      rownames(matKeep)))
startTree <- Preorder(TreeSearch::AdditionTree(mkd$phyDat))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

baseModelArgs <- list(
  rateLossMeanlog = 0, rateLossSdlog = 2,
  rateLogSdShape  = 1, rateLogSdRate = 1,
  rateNeoMeanlog  = 0, rateNeoSdlog  = 1
)
modelBlind <- do.call(MkPrimeModel,
                      c(baseModelArgs, list(ecologyAware = FALSE)))
modelAware <- do.call(MkPrimeModel, c(baseModelArgs, list(
  ecologyAware = TRUE, magnitudeMode = "global",
  rho0Alpha = 7, rho0Beta = 3, thetaAlpha = 2, thetaBeta = 2, sigmaPhi = 1.5
)))
nIter <- 200L
mcmcArgs <- list(nIter = nIter, nChains = 1L, nRuns = 1L, nCore = 1L,
                 thin = nIter, treeThin = nIter,
                 minWarmup = 10L, maxWarmup = 20L)
mcmc <- do.call(MkPrimeMCMC, mcmcArgs)
tB <- proc.time()
resB <- RunMkPrime(mkd, tree = startTree, model = modelBlind, mcmc = mcmc)
blindE <- (proc.time() - tB)["elapsed"]
tA <- proc.time()
resA <- RunMkPrime(mkd, tree = startTree, model = modelAware, mcmc = mcmc)
awareE <- (proc.time() - tA)["elapsed"]
cat(sprintf("LIB=%s  BENCH:  blind %.2fs  aware %.2fs  ratio %.2fx\n",
            libPath, blindE, awareE, awareE / blindE))
