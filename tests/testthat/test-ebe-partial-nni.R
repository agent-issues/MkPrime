# EBE Phase 2b / R3 / T5: the partial-NNI cache path (eco_cache_partial_eval_nni)
# must reconcile with a full recompute under the NON-STATIONARY EBE kernel.
# Promoted from the 1A dev driver dev/ecology/ebe-kernel-full-verify.R into the
# suite as a persisted regression guard (the path was 1A-verified at 0-7e-15 but
# not previously auto-run).  Forces the partial path (dirtyFracThreshold = 2.0,
# which never falls back) at z != 0, for coding "none" and "variable".
#
# This is the NNI leg of T5; the phi leg (R2: scale_phi cache invalidation) is
# pinned in test-ecology-plumbing.R ("scale_phi move keeps logLik in sync ...").

library("TreeTools")

`%||%` <- function(a, b) if (is.null(a)) b else a

# Self-contained ecology McmcData XPtr builder (mirrors the 1A driver and
# test-ecology-plumbing.R's .MakeEcoDataPtr; kept local — testthat does not
# share top-level defs across test files).
.EboNniDataPtr <- function(mkd, model) {
  parts <- lapply(mkd$partitions, function(p) {
    list(type = p$type, k = p$k, kObs = p$kObs,
         char_indices = p$char_indices,
         tip_states = p$tip_states,
         unique_tip_states = p$unique_tip_states,
         pattern_index = p$pattern_index)
  })
  ecoTip <- as.integer(mkd$ecology); ecoTip[is.na(ecoTip)] <- -1L
  MkPrime:::prepare_mcmc_data(
    parts, as.integer(mkd$kObs), mkd$type, any(mkd$type == "neomorphic"),
    model$nCat, model$coding, model$relabel,
    model$treeLengthShape %||% 2.0, model$treeLengthRate %||% 1.0,
    model$rateLossMeanlog, model$rateLossSdlog,
    model$rateLogSdShape, model$rateLogSdRate,
    model$rateNeoMeanlog, model$rateNeoSdlog,
    model$kprimeHyperA, model$kprimeHyperB,
    identical(model$kPrimePrior, "logseries"), model$kprimeLogseriesC %||% 0.7,
    FALSE, FALSE, 4L, 1.0, 1.0, FALSE, numeric(0), 1L, 0L, 0, -1e308, TRUE,
    ecoTip, as.integer(mkd$kEcology), model$magnitudeMode %||% "global",
    model$rho0Alpha %||% 75.0, model$rho0Beta %||% 25.0,
    model$sigmaPhi %||% 1.0, as.integer(model$gibbsZEvery %||% 50L),
    model$thetaAlpha %||% 1.0, model$thetaBeta %||% 1.0)
}

# Enumerate in-place NNIs (both endpoints internal) and apply one.
.EboEnumNni <- function(parent, child, nTip) {
  nEdge <- length(parent); out <- list()
  for (edgeRow in seq_len(nEdge)) {
    if (!(parent[edgeRow] > nTip && child[edgeRow] > nTip)) next
    u <- parent[edgeRow]; v <- child[edgeRow]
    vCh <- which(parent == v); uSib <- which(parent == u & child != v)
    if (length(vCh) == 0 || length(uSib) == 0) next
    for (cRow in vCh) for (wRow in uSib) if (wRow > edgeRow)
      out[[length(out) + 1L]] <- list(cRow = cRow, wRow = wRow, u = u, v = v)
  }
  out
}
.EboApplyNni <- function(parent, swap) {
  pNew <- parent; pNew[swap$cRow] <- swap$u; pNew[swap$wRow] <- swap$v; pNew
}

.EboPartialNniFixture <- function(coding) {
  set.seed(42)
  nTip <- 12L
  tips <- paste0("t", seq_len(nTip))
  mat <- cbind(
    sample(0:2, nTip, replace = TRUE), sample(0:2, nTip, replace = TRUE),
    sample(0:3, nTip, replace = TRUE),
    sample(0:1, nTip, replace = TRUE), sample(0:1, nTip, replace = TRUE))
  rownames(mat) <- tips
  pd <- TreeTools::MatrixToPhyDat(mat)
  ecoVec <- setNames(sample(0:2, nTip, replace = TRUE), tips)
  mkd <- MkPrimeData(pd, neomorphic = c(4L, 5L), ecology = ecoVec)
  model <- MkPrimeModel(
    rateLossMeanlog = 0, rateLossSdlog = 2, rateLogSdShape = 1, rateLogSdRate = 1,
    rateNeoMeanlog = 0, rateNeoSdlog = 1, ecologyAware = TRUE,
    magnitudeMode = "global", rho0Alpha = 7, rho0Beta = 3,
    thetaAlpha = 2, thetaBeta = 2, sigmaPhi = 1.5, coding = coding)
  list(mkd = mkd, model = model, tips = tips, nTip = nTip)
}

.EboCheckPartialNni <- function(coding) {
  f <- .EboPartialNniFixture(coding)
  dataPtr <- .EboNniDataPtr(f$mkd, f$model)
  kPrime <- as.integer(f$mkd$kObs)
  zCols  <- as.integer(f$mkd$kEcology - 1L)
  phi <- 1.6; pi0 <- 0.5; theta <- rep(0.5, zCols)
  rateLoss <- 1.1; rateLogSd <- 0.5; rateNeo <- 0.9
  set.seed(7L)
  zMix <- matrix(as.integer(sample(0:2, f$mkd$nChar * zCols, replace = TRUE)),
                 nrow = f$mkd$nChar, ncol = zCols)
  full <- function(p, c, el) MkPrime:::.CppLogLikelihoodEcology(
    dataPtr, p, c, el, kPrime, rateLoss, rateLogSd, rateNeo,
    phi = phi, zMatrix = zMix, pi0 = pi0, theta = theta)
  maxDP <- 0; maxDR <- 0; done <- 0L; attempt <- 0L
  set.seed(202L)
  while (done < 8L && attempt < 250L) {
    attempt <- attempt + 1L
    tr <- TreeTools::Preorder(ape::rtree(f$nTip, tip.label = f$tips))
    p <- as.integer(tr$edge[, 1]); c <- as.integer(tr$edge[, 2])
    el <- as.numeric(tr$edge.length)
    swaps <- .EboEnumNni(p, c, f$nTip)
    if (length(swaps) == 0) next
    swap <- swaps[[sample.int(length(swaps), 1L)]]
    pNew <- .EboApplyNni(p, swap)
    # threshold 2.0 -> dirty fraction (<=1) never exceeds it -> partial ALWAYS runs.
    pr <- MkPrime:::.CppPartialEvalEcologyNNI(
      dataPtr, p, c, pNew, c, el, kPrime, rateLoss, rateLogSd, rateNeo,
      phi, zMix, pi0, theta, swap$v, swap$u, 2.0)
    if (isTRUE(pr$fallback)) next
    maxDP <- max(maxDP, abs(pr$partial_ll  - full(pNew, c, el)))
    maxDR <- max(maxDR, abs(pr$restored_ll - full(p,    c, el)))
    done <- done + 1L
  }
  list(done = done, maxDP = maxDP, maxDR = maxDR)
}

test_that("EBE partial-NNI eval == full recompute at z != 0, coding = none (R3/T5)", {
  r <- .EboCheckPartialNni("none")
  expect_gte(r$done, 5L)
  expect_lt(r$maxDP, 1e-10)   # partial-eval (new tree) == full recompute
  expect_lt(r$maxDR, 1e-10)   # restored-eval (old tree) == full recompute
})

test_that("EBE partial-NNI eval == full recompute at z != 0, coding = variable (R3/T5)", {
  r <- .EboCheckPartialNni("variable")   # ascertainment correction engaged
  expect_gte(r$done, 5L)
  expect_lt(r$maxDP, 1e-10)
})
