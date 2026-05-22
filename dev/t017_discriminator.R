suppressPackageStartupMessages({
  library(MkPrime); library(TreeTools)
})
mat <- TreeTools::ReadCharacters("inst/ecology/data/rodent-X24848.nex")
ecoVec  <- mat[, 220]; extant <- mat[, 221]
keepTaxa <- which(extant == "0")
matKeep <- mat[keepTaxa, -221, drop = FALSE]
ecoKeep <- ecoVec[keepTaxa]
poly <- grepl("^\\(", ecoKeep)
if (any(poly)) {
  recoded <- substr(sub("^\\(", "", ecoKeep[poly]), 1, 1)
  ecoKeep[poly] <- recoded; matKeep[poly, 220] <- recoded
}
hasEco <- !is.na(ecoKeep) & ecoKeep != "?" & !is.na(suppressWarnings(as.integer(ecoKeep)))
matKeep <- matKeep[hasEco, , drop = FALSE]; ecoKeep <- ecoKeep[hasEco]
charMat <- matKeep[, -220, drop = FALSE]
pdForDetect <- MatrixToPhyDat(charMat)
neoIdx <- AutoDetectNeomorphic(pdForDetect)
.kObsCol <- function(col) {
  vals <- col[!(col %in% c("?", "-", NA))]
  vals <- unlist(strsplit(gsub("[()]", "", vals), ""))
  length(unique(vals))
}
kObsRaw <- vapply(seq_len(ncol(charMat)), function(j) .kObsCol(charMat[, j]), integer(1L))
nonNeoOriginal <- setdiff(seq_len(ncol(charMat)), neoIdx)
nonNeoVariable <- nonNeoOriginal[kObsRaw[nonNeoOriginal] >= 2L]
knownK <- setNames(as.integer(kObsRaw[nonNeoVariable]), as.character(nonNeoVariable))
mkd <- MkPrimeData(pdForDetect, neomorphic = neoIdx, knownStates = knownK,
                   ecology = setNames(as.integer(ecoKeep), rownames(matKeep)))
cat("nChar=", mkd$nChar, "  nTip=", mkd$nTip, "\n")
cat("type counts:\n"); print(table(mkd$type))
cat("\nDistinct partition types:", length(unique(mkd$type)), "\n")
trans_idx <- which(mkd$type %in% c("transformational", "known"))
cat("\nTrans/known-type kObs distribution (kPrime initialised to kObs):\n")
print(table(mkd$kObs[trans_idx]))
cat("\nNumber of distinct kObs in trans/known:",
    length(unique(mkd$kObs[trans_idx])), "\n")
