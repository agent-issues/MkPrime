# t019_pt_clamp_test.R - sanity test for tighter heat clamp + round-trip count
suppressPackageStartupMessages({
  library(MkPrime); library(TreeTools)
})
nIter <- 500L
nexFile <- "inst/ecology/data/rodent-X24848.nex"
mat <- TreeTools::ReadCharacters(nexFile)
ecoVec  <- mat[, 220]; extant <- mat[, 221]
keepTaxa <- which(extant == "0")
matKeep <- mat[keepTaxa, -221, drop = FALSE]
ecoKeep <- ecoVec[keepTaxa]
poly <- grepl("^\\(", ecoKeep)
ecoKeep[poly] <- substr(sub("^\\(", "", ecoKeep[poly]), 1, 1)
matKeep[poly, 220] <- ecoKeep[poly]
hasEco <- !is.na(ecoKeep) & ecoKeep != "?" &
          !is.na(suppressWarnings(as.integer(ecoKeep)))
matKeep <- matKeep[hasEco, , drop = FALSE]; ecoKeep <- ecoKeep[hasEco]
charMat <- matKeep[, -220, drop = FALSE]
pdF <- MatrixToPhyDat(charMat); neoIdx <- AutoDetectNeomorphic(pdF)
.kObsCol <- function(col) {
  v <- col[!(col %in% c("?","-",NA))]
  v <- unlist(strsplit(gsub("[()]","",v),""))
  length(unique(v))
}
kObsRaw <- vapply(seq_len(ncol(charMat)),
                  function(j) .kObsCol(charMat[,j]), integer(1L))
nonNeoVar <- setdiff(seq_len(ncol(charMat)), neoIdx)
nonNeoVar <- nonNeoVar[kObsRaw[nonNeoVar] >= 2L]
mkd <- MkPrimeData(pdF, neomorphic = neoIdx,
  knownStates = setNames(as.integer(kObsRaw[nonNeoVar]),
                         as.character(nonNeoVar)),
  ecology = setNames(as.integer(ecoKeep), rownames(matKeep)))
set.seed(20260521)
startTree <- Preorder(TreeSearch::AdditionTree(mkd$phyDat))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))
modelA <- MkPrimeModel(ecologyAware = TRUE, magnitudeMode = "global",
  rho0Alpha = 7, rho0Beta = 3, thetaAlpha = 2, thetaBeta = 2, sigmaPhi = 1.5)
mcmc <- MkPrimeMCMC(nIter = nIter, nChains = 4L, nRuns = 1L, nCore = 1L,
  thin = 1L, treeThin = 1L,
  minWarmup = 100L, maxWarmup = 200L,
  logFile = tempfile(fileext = ".log"),
  checkpointFile = tempfile(fileext = ".ckp"))
t0 <- Sys.time()
res <- RunMkPrime(mkd, tree = startTree, model = modelA, mcmc = mcmc)
cat("Elapsed:", round(as.numeric(Sys.time() - t0, units = "secs"), 2), "s\n")
cat("Final betas:", paste(round(res$betas, 4), collapse = ", "), "\n")
cat("Swap rates:", paste(round(res$swap_rates, 3), collapse = ", "), "\n")
cat("Round trips:", res$round_trip_count, "\n")
cat("ladder_warning attr:", attr(res$betas, "ladder_warning"), "\n")
