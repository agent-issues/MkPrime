#!/usr/bin/env Rscript
# 12c_rodent_dirty_fraction.R — diagnostic: how many internal nodes does the
# NNI partial-eval mark dirty on the rodent matrix?  If most of the tree
# becomes dirty after every NNI, partial-eval cannot beat full-eval and the
# T-013 plan is dead.

suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKPRIME_LIBPATH", "dev/profiling/.vtune-lib-t013")
  library("MkPrime", lib.loc = libPath)
  library("TreeTools")
})

`%||%` <- function(a, b) if (is.null(a)) b else a

make_eco_dataptr <- function(mkd, model) {
  parts <- lapply(mkd$partitions, function(p) {
    list(type = p$type, k = p$k, kObs = p$kObs,
         char_indices = p$char_indices,
         tip_states = p$tip_states,
         unique_tip_states = p$unique_tip_states,
         pattern_index = p$pattern_index)
  })
  ecoTip <- as.integer(mkd$ecology)
  ecoTip[is.na(ecoTip)] <- -1L
  MkPrime:::prepare_mcmc_data(
    parts, as.integer(mkd$kObs), mkd$type,
    any(mkd$type == "neomorphic"),
    model$nCat, model$coding, model$relabel,
    model$treeLengthShape %||% 2.0, model$treeLengthRate %||% 1.0,
    model$rateLossMeanlog, model$rateLossSdlog,
    model$rateLogSdShape, model$rateLogSdRate,
    model$rateNeoMeanlog, model$rateNeoSdlog,
    model$kprimeHyperA, model$kprimeHyperB,
    identical(model$kPrimePrior, "logseries"),
    model$kprimeLogseriesC %||% 0.7,
    FALSE, FALSE, 4L, 1.0, 1.0,
    FALSE, numeric(0), 1L, 0L, 0, -1e308,
    TRUE,
    ecoTip, as.integer(mkd$kEcology),
    model$magnitudeMode %||% "global",
    model$rho0Alpha %||% 75.0, model$rho0Beta %||% 25.0,
    model$sigmaPhi %||% 1.0, as.integer(model$gibbsZEvery %||% 50L),
    model$thetaAlpha %||% 1.0, model$thetaBeta %||% 1.0
  )
}

# ------------------------ Rodent fixture (matches 11e/11g) ------------------
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

model <- MkPrimeModel(
  rateLossMeanlog = 0, rateLossSdlog = 2,
  rateLogSdShape = 1, rateLogSdRate = 1,
  rateNeoMeanlog = 0, rateNeoSdlog = 1,
  ecologyAware = TRUE, magnitudeMode = "global",
  rho0Alpha = 7, rho0Beta = 3,
  thetaAlpha = 2, thetaBeta = 2,
  sigmaPhi = 1.5, coding = "none"
)
dataPtr <- make_eco_dataptr(mkd, model)

# ---------------- Enumerate in-place NNIs on the addition tree --------------
parent <- as.integer(startTree$edge[, 1])
child  <- as.integer(startTree$edge[, 2])
edgeLen <- as.numeric(startTree$edge.length)
nTip <- as.integer(length(startTree$tip.label))

enumerate_inplace_nni <- function(parent, child, nTip) {
  nEdge <- length(parent)
  out <- list()
  for (edgeRow in seq_len(nEdge)) {
    if (!(parent[edgeRow] > nTip && child[edgeRow] > nTip)) next
    u <- parent[edgeRow]; v <- child[edgeRow]
    vCh <- which(parent == v)
    uSib <- which(parent == u & child != v)
    if (length(vCh) == 0 || length(uSib) == 0) next
    for (cRow in vCh) for (wRow in uSib) {
      if (wRow > edgeRow) {
        out[[length(out) + 1L]] <- list(
          edgeRow = edgeRow, u = u, v = v,
          cRow = cRow, wRow = wRow,
          cNode = child[cRow], wNode = child[wRow]
        )
      }
    }
  }
  out
}
apply_nni <- function(parent, swap) {
  pNew <- parent
  pNew[swap$cRow] <- swap$u
  pNew[swap$wRow] <- swap$v
  pNew
}

swaps <- enumerate_inplace_nni(parent, child, nTip)
cat(sprintf("Available in-place NNI swaps: %d\n", length(swaps)))
stopifnot(length(swaps) > 0)

zCols <- as.integer(mkd$kEcology - 1L)
zMat <- matrix(sample(0:2, mkd$nChar * zCols, replace = TRUE),
               nrow = mkd$nChar, ncol = zCols)
storage.mode(zMat) <- "integer"
kPrime <- as.integer(mkd$kObs)
phi <- 1.6; pi0 <- 0.5; theta <- rep(0.5, zCols)

# Sample N NNI moves; record dirty_count, nWedgeDirty, fallback (with threshold 0.99).
N <- min(30L, length(swaps))
idx <- sample.int(length(swaps), N)
res <- lapply(idx, function(i) {
  s <- swaps[[i]]
  pNew <- apply_nni(parent, s)
  r <- MkPrime:::.CppPartialEvalEcologyNNI(
    dataPtr, parent, child, pNew, child, edgeLen,
    kPrime, 1.0, 0.5, 1.0, phi, zMat, pi0, theta,
    s$v, s$u, 0.99
  )
  list(dirty = r$dirty_count, wEdge = r$n_wedge_dirty_edges,
       n_internal = r$n_internal, fallback = r$fallback)
})

dirty <- vapply(res, function(x) x$dirty, integer(1))
wEdge <- vapply(res, function(x) x$wEdge, integer(1))
nInt  <- res[[1]]$n_internal
cat(sprintf("nTip = %d, nInternal = %d, nEdge = %d\n", nTip, nInt, length(parent)))
cat(sprintf("dirty count  : mean = %.1f, median = %.0f, max = %d (fraction of nInternal: mean %.3f, max %.3f)\n",
            mean(dirty), median(dirty), max(dirty),
            mean(dirty)/nInt, max(dirty)/nInt))
cat(sprintf("wEdge dirty  : mean = %.1f, median = %.0f, max = %d\n",
            mean(wEdge), median(wEdge), max(wEdge)))
cat(sprintf("fallbacks @ 0.99 threshold: %d / %d\n",
            sum(vapply(res, function(x) isTRUE(x$fallback), logical(1))), N))

# Print expected savings: if dirty_count = X, naive expectation is 1 - X/nInternal
# of orchestrator pruning work is saved.
cat(sprintf("\nExpected naive partial-eval saving = 1 - dirty/nInternal = %.1f%% (median)\n",
            100 * (1 - median(dirty)/nInt)))
