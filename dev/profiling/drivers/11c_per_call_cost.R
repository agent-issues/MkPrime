# Isolated micro-bench: time one full `cpp_log_likelihood_ecology` call vs
# one `cpp_log_likelihood` (blind) call on the rodent matrix. Tells us the
# per-call cost gap independent of MCMC move scheduling.
suppressPackageStartupMessages({
  library("MkPrime", lib.loc = Sys.getenv("MKPRIME_LIBPATH",
                                          "dev/profiling/.vtune-lib-prof"))
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
# Parsimony+ratchet was burning ~30 s; AdditionTree is essentially free.
startTree <- Preorder(TreeSearch::AdditionTree(mkd$phyDat))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

nTip <- mkd$nTip
parent <- as.integer(startTree$edge[, 1])
child  <- as.integer(startTree$edge[, 2])
edgeLen <- as.numeric(startTree$edge.length)
kPrime <- as.integer(mkd$kObs)  # MkNT: pin k=kObs
rateLoss <- 1; rateLogSd <- 0; rateNeo <- 1
phi <- c(1)
kEco <- mkd$kEcology
zMatrix <- matrix(0L, nrow = mkd$nChar, ncol = kEco - 1)
pi0 <- 0.7
theta <- rep(0.5, kEco - 1)
nCat <- 1L

# Warm-up
invisible(MkPrime:::.MkpEcologyLogLikelihood(
  startTree, mkd, kPrime, rateLoss, rateLogSd, nCat, rateNeo, FALSE,
  phi, zMatrix, "global", "none", refEcology = NULL, theta = theta, pi0 = pi0))
invisible(MkPrime:::.MkpLogLikelihood(
  startTree, mkd, kPrime, rateLoss, rateLogSd, nCat, "none", rateNeo, FALSE))

nRep <- 30L
tA <- system.time(for (i in seq_len(nRep)) {
  MkPrime:::.MkpEcologyLogLikelihood(
    startTree, mkd, kPrime, rateLoss, rateLogSd, nCat, rateNeo, FALSE,
    phi, zMatrix, "global", "none",
    refEcology = NULL, theta = theta, pi0 = pi0)
})["elapsed"]
tB <- system.time(for (i in seq_len(nRep)) {
  MkPrime:::.MkpLogLikelihood(
    startTree, mkd, kPrime, rateLoss, rateLogSd, nCat, "none", rateNeo, FALSE)
})["elapsed"]
cat(sprintf("Per-call cost over %d reps:\n", nRep))
cat(sprintf("  Aware orchestrator: %.2f ms/call\n", 1000 * tA / nRep))
cat(sprintf("  Blind likelihood:   %.2f ms/call\n", 1000 * tB / nRep))
cat(sprintf("  Per-call ratio:     %.2fx\n\n", tA / tB))

# Per-character ecology cost (used by gibbs_z sweep)
zCols <- kEco - 1L
zRow <- integer(zCols)
nCallsZSweep <- mkd$nChar * zCols * 3L
cat(sprintf("Gibbs z sweep calls per_char_log_lik_ecology: %d × %d × 3 = %d\n",
            mkd$nChar, zCols, nCallsZSweep))
