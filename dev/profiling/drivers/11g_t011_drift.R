# 11g_t011_drift.R — 5000-iter aware MCMC drift check.
# Required by the T-011 spec: assert zero `[eco-resync]` warnings.
suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKPRIME_LIBPATH", "dev/profiling/.vtune-lib-t011")
  library("MkPrime", lib.loc = libPath)
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
startTree <- Preorder(TreeSearch::AdditionTree(mkd$phyDat))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))
modelAware <- MkPrimeModel(
  rateLossMeanlog = 0, rateLossSdlog = 2,
  rateLogSdShape = 1, rateLogSdRate = 1,
  rateNeoMeanlog = 0, rateNeoSdlog = 1,
  ecologyAware = TRUE, magnitudeMode = "global",
  rho0Alpha = 7, rho0Beta = 3, thetaAlpha = 2, thetaBeta = 2, sigmaPhi = 1.5
)
mcmc <- MkPrimeMCMC(nIter = 5000L, nChains = 1L, nRuns = 1L, nCore = 1L,
                    thin = 5000L, treeThin = 5000L,
                    minWarmup = 10L, maxWarmup = 20L)
tA <- proc.time()
# Capture stderr so we can count [eco-resync] warnings emitted by the
# periodic resync check inside run_mcmc_batch_cpp.
log <- tempfile(fileext = ".log")
con <- file(log, "wt")
sink(con, type = "message")
resA <- RunMkPrime(mkd, tree = startTree, model = modelAware, mcmc = mcmc)
sink(type = "message")
close(con)
awareE <- (proc.time() - tA)["elapsed"]
cat(sprintf("AWARE 5000-iter wall: %.2f s\n", awareE))
logLines <- readLines(log)
resync <- grep("eco-resync", logLines, value = TRUE)
cat(sprintf("[eco-resync] warnings: %d\n", length(resync)))
if (length(resync) > 0) {
  cat("First few:\n")
  cat(head(resync, 10), sep = "\n")
}
