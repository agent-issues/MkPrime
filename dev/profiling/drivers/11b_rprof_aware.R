suppressPackageStartupMessages({
  library("MkPrime", lib.loc = "dev/profiling/.vtune-lib-prof")
  library("TreeTools")
})
set.seed(20260521)
nexFile <- "C:/Users/pjjg18/downloads/mbank_X24848_2026-5-9-1135.nex"
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
# Greedy stepwise-addition parsimony start tree — instant, parsimony-based.
# Parsimony+ratchet was burning ~30 s per run via the slow Morphy delegate;
# AdditionTree is essentially free.
startTree <- Preorder(TreeSearch::AdditionTree(mkd$phyDat))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))
modelAware <- MkPrimeModel(
  rateLossMeanlog = 0, rateLossSdlog = 2,
  rateLogSdShape = 1, rateLogSdRate = 1,
  rateNeoMeanlog = 0, rateNeoSdlog = 1,
  ecologyAware = TRUE, magnitudeMode = "global",
  rho0Alpha = 7, rho0Beta = 3, thetaAlpha = 2, thetaBeta = 2, sigmaPhi = 1.5)
mcmc <- MkPrimeMCMC(nIter = 200L, nChains = 1L, nRuns = 1L, nCore = 1L,
  thin = 200L, treeThin = 200L, minWarmup = 10L, maxWarmup = 20L)

Rprof("dev/profiling/aware_rprof.out", interval = 0.02)
res <- RunMkPrime(mkd, tree = startTree, model = modelAware, mcmc = mcmc)
Rprof(NULL)
summ <- summaryRprof("dev/profiling/aware_rprof.out")
cat("\n=== by.self (top 15) ===\n")
print(head(summ$by.self, 15))
cat("\n=== by.total (top 15) ===\n")
print(head(summ$by.total, 15))
cat("\nsampling.time:", summ$sampling.time, "\n")
