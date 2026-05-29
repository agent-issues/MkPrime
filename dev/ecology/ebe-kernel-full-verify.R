# EBE kernel full verification (post all edits): cache-path parity (entry 7),
# partial-NNI parity (entry 8, R3/T5), and R8 (transformational z-invariance).
# Compiles in THIS subprocess via load_all. Standalone driver, kernel agent.

suppressMessages({
  pkgbuild::compile_dll(".")
  devtools::load_all(".", quiet = TRUE)
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

cat("=== EBE full verification (mixed neo + trans data) ===\n")

set.seed(42)
nTip <- 12L
tips <- paste0("t", seq_len(nTip))
mat <- cbind(
  sample(0:2, nTip, replace = TRUE),   # trans (char 1)
  sample(0:2, nTip, replace = TRUE),   # trans
  sample(0:3, nTip, replace = TRUE),   # trans (k=4)
  sample(0:1, nTip, replace = TRUE),   # neo (char 4)
  sample(0:1, nTip, replace = TRUE)    # neo (char 5)
)
rownames(mat) <- tips
pd <- TreeTools::MatrixToPhyDat(mat)
ecoVec <- setNames(sample(0:2, nTip, replace = TRUE), tips)
mkd <- MkPrimeData(pd, neomorphic = c(4L, 5L), ecology = ecoVec)
model <- MkPrimeModel(
  rateLossMeanlog = 0, rateLossSdlog = 2,
  rateLogSdShape = 1, rateLogSdRate = 1,
  rateNeoMeanlog = 0, rateNeoSdlog = 1,
  ecologyAware = TRUE, magnitudeMode = "global",
  rho0Alpha = 7, rho0Beta = 3, thetaAlpha = 2, thetaBeta = 2,
  sigmaPhi = 1.5, coding = "none"
)
dataPtr <- make_eco_dataptr(mkd, model)

tree <- TreeTools::Preorder(ape::rtree(nTip, tip.label = tips))
parent <- as.integer(tree$edge[, 1]); child <- as.integer(tree$edge[, 2])
edgeLen <- as.numeric(tree$edge.length)
kPrime <- as.integer(mkd$kObs)
zCols <- as.integer(mkd$kEcology - 1L)
phi <- 1.6; pi0 <- 0.5; theta <- rep(0.5, zCols)
rateLoss <- 1.1; rateLogSd <- 0.5; rateNeo <- 0.9

zZero <- matrix(0L, nrow = mkd$nChar, ncol = zCols)
set.seed(7)
zMix <- matrix(as.integer(sample(0:2, mkd$nChar * zCols, replace = TRUE)),
               nrow = mkd$nChar, ncol = zCols)

orch <- function(p, c, z) MkPrime:::.CppLogLikelihoodEcology(
  dataPtr, p, c, edgeLen, kPrime, rateLoss, rateLogSd, rateNeo,
  phi = phi, zMatrix = z, pi0 = pi0, theta = theta)
cached <- function(p, c, z) MkPrime:::.CppLogLikelihoodEcologyCached(
  dataPtr, p, c, edgeLen, kPrime, rateLoss, rateLogSd, rateNeo,
  phi = phi, zMatrix = z, pi0 = pi0, theta = theta)

# ---- Entry-7 cache parity: cached-full == orchestrator (z=0 AND z!=0) ----
for (lab in c("z=0", "z!=0")) {
  z <- if (lab == "z=0") zZero else zMix
  lo <- orch(parent, child, z); lc <- cached(parent, child, z)
  cat(sprintf("[cache parity %-5s] orch=%.10f cached=%.10f  |diff|=%.3g  %s\n",
              lab, lo, lc, abs(lo - lc),
              if (abs(lo - lc) < 1e-10) "PASS" else "FAIL"))
  stopifnot(abs(lo - lc) < 1e-10)
}

# ---- R8: transformational/known carry NO ecology effect (z ignored) ----
# Build a TRANS-ONLY dataset; orchestrator logLik must be identical for any z.
matT <- mat[, 1:3, drop = FALSE]
pdT <- TreeTools::MatrixToPhyDat(matT)
mkdT <- MkPrimeData(pdT, ecology = ecoVec)   # no neomorphic => all trans
dataPtrT <- make_eco_dataptr(mkdT, model)
kPrimeT <- as.integer(mkdT$kObs)
zCT <- as.integer(mkdT$kEcology - 1L)
zT0 <- matrix(0L, nrow = mkdT$nChar, ncol = zCT)
set.seed(99)
zT1 <- matrix(as.integer(sample(0:2, mkdT$nChar * zCT, replace = TRUE)),
              nrow = mkdT$nChar, ncol = zCT)
llT0 <- MkPrime:::.CppLogLikelihoodEcology(
  dataPtrT, parent, child, edgeLen, kPrimeT, rateLoss, rateLogSd, rateNeo,
  phi = phi, zMatrix = zT0, pi0 = pi0, theta = theta)
llT1 <- MkPrime:::.CppLogLikelihoodEcology(
  dataPtrT, parent, child, edgeLen, kPrimeT, rateLoss, rateLogSd, rateNeo,
  phi = phi, zMatrix = zT1, pi0 = pi0, theta = theta)
cat(sprintf("[R8 trans z-invariance] z0=%.10f zMix=%.10f  |diff|=%.3g  %s\n",
            llT0, llT1, abs(llT0 - llT1),
            if (abs(llT0 - llT1) < 1e-12) "PASS (z ignored)" else "FAIL"))
stopifnot(abs(llT0 - llT1) < 1e-12)
# And R8 with phi varied (z!=0): still invariant => phi has no trans effect.
llT1b <- MkPrime:::.CppLogLikelihoodEcology(
  dataPtrT, parent, child, edgeLen, kPrimeT, rateLoss, rateLogSd, rateNeo,
  phi = 3.3, zMatrix = zT1, pi0 = pi0, theta = theta)
cat(sprintf("[R8 trans phi-invariance] phi1.6=%.10f phi3.3=%.10f |diff|=%.3g %s\n",
            llT1, llT1b, abs(llT1 - llT1b),
            if (abs(llT1 - llT1b) < 1e-12) "PASS" else "FAIL"))
stopifnot(abs(llT1 - llT1b) < 1e-12)

# ---- R3 / T5: partial-NNI == full recompute (z!=0), partial path ENGAGED ----
enumerate_inplace_nni <- function(parent, child, nTip) {
  nEdge <- length(parent); out <- list()
  for (edgeRow in seq_len(nEdge)) {
    if (!(parent[edgeRow] > nTip && child[edgeRow] > nTip)) next
    u <- parent[edgeRow]; v <- child[edgeRow]
    vCh <- which(parent == v); uSib <- which(parent == u & child != v)
    if (length(vCh) == 0 || length(uSib) == 0) next
    for (cRow in vCh) for (wRow in uSib) if (wRow > edgeRow)
      out[[length(out) + 1L]] <- list(edgeRow = edgeRow, u = u, v = v,
        cRow = cRow, wRow = wRow, cNode = child[cRow], wNode = child[wRow])
  }
  out
}
apply_nni <- function(parent, swap) {
  pNew <- parent; pNew[swap$cRow] <- swap$u; pNew[swap$wRow] <- swap$v; pNew
}

maxDP <- 0; maxDR <- 0; done <- 0L; attempt <- 0L; engaged <- 0L
while (done < 15L && attempt < 200L) {
  attempt <- attempt + 1L
  tr <- TreeTools::Preorder(ape::rtree(nTip, tip.label = tips))
  p <- as.integer(tr$edge[, 1]); c <- as.integer(tr$edge[, 2])
  el <- as.numeric(tr$edge.length)
  swaps <- enumerate_inplace_nni(p, c, nTip)
  if (length(swaps) == 0) next
  swap <- swaps[[sample.int(length(swaps), 1)]]
  pNew <- apply_nni(p, swap)
  # NOTE: fallback threshold 2.0 forces partial-eval to ALWAYS run (no fallback),
  # so this genuinely exercises the partial path under z != 0.
  pr <- MkPrime:::.CppPartialEvalEcologyNNI(
    dataPtr, p, c, pNew, c, el, kPrime, rateLoss, rateLogSd, rateNeo,
    phi, zMix, pi0, theta, swap$v, swap$u, 2.0)
  if (isTRUE(pr$fallback)) next
  engaged <- engaged + 1L
  freshNew <- MkPrime:::.CppLogLikelihoodEcology(
    dataPtr, pNew, c, el, kPrime, rateLoss, rateLogSd, rateNeo,
    phi = phi, zMatrix = zMix, pi0 = pi0, theta = theta)
  freshOld <- MkPrime:::.CppLogLikelihoodEcology(
    dataPtr, p, c, el, kPrime, rateLoss, rateLogSd, rateNeo,
    phi = phi, zMatrix = zMix, pi0 = pi0, theta = theta)
  dP <- abs(pr$partial_ll - freshNew); dR <- abs(pr$restored_ll - freshOld)
  maxDP <- max(maxDP, dP); maxDR <- max(maxDR, dR)
  done <- done + 1L
}
cat(sprintf("[R3/T5 partial-NNI z!=0] samples=%d engaged=%d  max|partial-full|=%.3g  max|restored-full|=%.3g  %s\n",
            done, engaged, maxDP, maxDR,
            if (done >= 10L && maxDP < 1e-10 && maxDR < 1e-10) "PASS" else "FAIL"))
stopifnot(done >= 10L, maxDP < 1e-10, maxDR < 1e-10)

# ---------------------------------------------------------------------------
# coding = "variable" (ascertainment) — the spec's emphasized correctness gate
# (§8 T2-with-coding=1).  Re-run cache-parity + one partial-NNI pass at z!=0
# so the orchestrator vs cached const-site paths (which have different guards)
# are exercised on the EBE kernel.
# ---------------------------------------------------------------------------
cat("\n--- coding = 'variable' (ascertainment) ---\n")
modelV <- MkPrimeModel(
  rateLossMeanlog = 0, rateLossSdlog = 2,
  rateLogSdShape = 1, rateLogSdRate = 1,
  rateNeoMeanlog = 0, rateNeoSdlog = 1,
  ecologyAware = TRUE, magnitudeMode = "global",
  rho0Alpha = 7, rho0Beta = 3, thetaAlpha = 2, thetaBeta = 2,
  sigmaPhi = 1.5, coding = "variable"
)
dataPtrV <- make_eco_dataptr(mkd, modelV)
orchV <- function(p, c, z) MkPrime:::.CppLogLikelihoodEcology(
  dataPtrV, p, c, edgeLen, kPrime, rateLoss, rateLogSd, rateNeo,
  phi = phi, zMatrix = z, pi0 = pi0, theta = theta)
cachedV <- function(p, c, z) MkPrime:::.CppLogLikelihoodEcologyCached(
  dataPtrV, p, c, edgeLen, kPrime, rateLoss, rateLogSd, rateNeo,
  phi = phi, zMatrix = z, pi0 = pi0, theta = theta)
for (lab in c("z=0", "z!=0")) {
  z <- if (lab == "z=0") zZero else zMix
  lo <- orchV(parent, child, z); lc <- cachedV(parent, child, z)
  cat(sprintf("[cache parity coding=1 %-5s] orch=%.10f cached=%.10f  |diff|=%.3g  %s\n",
              lab, lo, lc, abs(lo - lc),
              if (abs(lo - lc) < 1e-10) "PASS" else "FAIL"))
  stopifnot(abs(lo - lc) < 1e-10)
}
# One partial-NNI pass under coding=1, z!=0, partial engaged.
maxDPv <- 0; doneV <- 0L; attemptV <- 0L; engagedV <- 0L
while (doneV < 8L && attemptV < 200L) {
  attemptV <- attemptV + 1L
  tr <- TreeTools::Preorder(ape::rtree(nTip, tip.label = tips))
  p <- as.integer(tr$edge[, 1]); c <- as.integer(tr$edge[, 2])
  el <- as.numeric(tr$edge.length)
  swaps <- enumerate_inplace_nni(p, c, nTip)
  if (length(swaps) == 0) next
  swap <- swaps[[sample.int(length(swaps), 1)]]
  pNew <- apply_nni(p, swap)
  pr <- MkPrime:::.CppPartialEvalEcologyNNI(
    dataPtrV, p, c, pNew, c, el, kPrime, rateLoss, rateLogSd, rateNeo,
    phi, zMix, pi0, theta, swap$v, swap$u, 2.0)
  if (isTRUE(pr$fallback)) next
  engagedV <- engagedV + 1L
  freshNew <- MkPrime:::.CppLogLikelihoodEcology(
    dataPtrV, pNew, c, el, kPrime, rateLoss, rateLogSd, rateNeo,
    phi = phi, zMatrix = zMix, pi0 = pi0, theta = theta)
  maxDPv <- max(maxDPv, abs(pr$partial_ll - freshNew))
  doneV <- doneV + 1L
}
cat(sprintf("[R3/T5 coding=1 z!=0] samples=%d engaged=%d  max|partial-full|=%.3g  %s\n",
            doneV, engagedV, maxDPv,
            if (doneV >= 5L && maxDPv < 1e-10) "PASS" else "FAIL"))
stopifnot(doneV >= 5L, maxDPv < 1e-10)

cat("=== ALL FULL-VERIFY CHECKS PASS ===\n")
