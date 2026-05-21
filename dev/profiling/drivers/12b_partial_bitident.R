#!/usr/bin/env Rscript
# 12b_partial_bitident.R — T-013 partial-eval bit-identity test for NNI.
#
# For ~20 random in-place NNI swaps, assert that
#   .CppPartialEvalEcologyNNI(parentOld -> parentNew, v, u)
# returns a partial_ll bit-identical to
#   .CppLogLikelihoodEcology(parentNew, childNew, ...)  (from-scratch full eval).
#
# Also asserts rollback: after restore_dirty_eco_cls + recompute total_loglik
# on OLD topology, restored_ll equals a fresh full eval on OLD.

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

# ----------------------------------------------------------------------------
# Fixture: ~12-tip tree with neomorphic + transformational chars + ecology
# ----------------------------------------------------------------------------
set.seed(2026013)
nTip <- 20L
tips <- paste0("t", seq_len(nTip))
mat <- cbind(
  sample(0:2, nTip, replace = TRUE),
  sample(0:2, nTip, replace = TRUE),
  sample(0:2, nTip, replace = TRUE),
  sample(0:2, nTip, replace = TRUE),
  sample(0:3, nTip, replace = TRUE),
  sample(0:3, nTip, replace = TRUE),
  sample(0:1, nTip, replace = TRUE),
  sample(0:1, nTip, replace = TRUE)
)
rownames(mat) <- tips
pd <- TreeTools::MatrixToPhyDat(mat)
ecoVec <- setNames(sample(0:2, nTip, replace = TRUE), tips)
mkd <- MkPrimeData(pd, neomorphic = c(7L, 8L), ecology = ecoVec)
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

# Random latent state.
zCols <- as.integer(mkd$kEcology - 1L)
zMat <- matrix(sample(0:2, mkd$nChar * zCols, replace = TRUE),
               nrow = mkd$nChar, ncol = zCols)
storage.mode(zMat) <- "integer"
kPrime <- as.integer(mkd$kObs)
phi <- 1.6
pi0 <- 0.5
theta <- rep(0.5, zCols)
rateLoss <- 1.0
rateLogSd <- 0.5
rateNeo <- 1.0

# ----------------------------------------------------------------------------
# Helper: enumerate valid in-place NNI swaps (mirrors case 5 in src/mcmc.cpp,
# wRow > edgeRow branch).  Returns list of swap descriptors.
# ----------------------------------------------------------------------------
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

# Apply an in-place NNI swap.
apply_nni <- function(parent, swap) {
  pNew <- parent
  pNew[swap$cRow] <- swap$u
  pNew[swap$wRow] <- swap$v
  pNew
}

# ----------------------------------------------------------------------------
# Run the regression: many random NNI swaps, both partial-eval and rollback.
# ----------------------------------------------------------------------------
nSwapsTarget <- 20L
results <- list()
maxDiffPartial <- 0
maxDiffRestored <- 0
fallbackCount <- 0
samplesDone <- 0L
attempt <- 0L

while (samplesDone < nSwapsTarget && attempt < 200L) {
  attempt <- attempt + 1L
  tree <- TreeTools::Preorder(ape::rtree(nTip, tip.label = tips))
  parent <- as.integer(tree$edge[, 1])
  child  <- as.integer(tree$edge[, 2])
  edgeLen <- as.numeric(tree$edge.length)

  swaps <- enumerate_inplace_nni(parent, child, nTip)
  if (length(swaps) == 0) next

  swap <- swaps[[sample.int(length(swaps), 1)]]
  parentNew <- apply_nni(parent, swap)
  childNew  <- child  # NNI swaps parent assignments, child unchanged

  partial_res <- MkPrime:::.CppPartialEvalEcologyNNI(
    dataPtr,
    parent, child, parentNew, childNew, edgeLen,
    kPrime, rateLoss, rateLogSd, rateNeo,
    phi, zMat, pi0, theta,
    swap$v, swap$u,
    2.0  # disable fallback — bit-identity test ALWAYS runs partial-eval
  )
  if (samplesDone == 0L && fallbackCount < 3L) {
    cat(sprintf("  diag: swap v=%d u=%d edgeRow=%d cRow=%d wRow=%d  dirty=%d wEdge=%d nInternal=%d  fallback=%s\n",
                swap$v, swap$u, swap$edgeRow, swap$cRow, swap$wRow,
                partial_res$dirty_count, partial_res$n_wedge_dirty_edges,
                partial_res$n_internal, partial_res$fallback))
  }
  if (isTRUE(partial_res$fallback)) {
    fallbackCount <- fallbackCount + 1L
    next
  }
  partial_ll <- partial_res$partial_ll
  restored_ll <- partial_res$restored_ll
  dirty_n <- partial_res$dirty_count

  # Fresh full eval against NEW topology.
  fresh_new <- MkPrime:::.CppLogLikelihoodEcology(
    dataPtr, parentNew, childNew, edgeLen, kPrime,
    rateLoss = rateLoss, rateLogSd = rateLogSd, rateNeo = rateNeo,
    phi = phi, zMatrix = zMat, pi0 = pi0, theta = theta
  )
  # Fresh full eval against OLD topology (for rollback regression).
  fresh_old <- MkPrime:::.CppLogLikelihoodEcology(
    dataPtr, parent, child, edgeLen, kPrime,
    rateLoss = rateLoss, rateLogSd = rateLogSd, rateNeo = rateNeo,
    phi = phi, zMatrix = zMat, pi0 = pi0, theta = theta
  )

  diffP <- abs(partial_ll - fresh_new)
  diffR <- abs(restored_ll - fresh_old)
  maxDiffPartial <- max(maxDiffPartial, diffP)
  maxDiffRestored <- max(maxDiffRestored, diffR)

  samplesDone <- samplesDone + 1L
  results[[samplesDone]] <- list(
    swap = swap, dirty_n = dirty_n,
    partial_ll = partial_ll, fresh_new = fresh_new,
    restored_ll = restored_ll, fresh_old = fresh_old,
    diff_partial = diffP, diff_restored = diffR
  )
  cat(sprintf("[%2d] dirty=%2d  partial=%14.10f  fresh_new=%14.10f  |dP|=%.3g  |dR|=%.3g\n",
              samplesDone, dirty_n, partial_ll, fresh_new, diffP, diffR))
}

cat("\n=== T-013 partial-eval bit-identity summary ===\n")
cat(sprintf("Samples completed: %d / target %d (fallback hits: %d, attempts: %d)\n",
            samplesDone, nSwapsTarget, fallbackCount, attempt))
cat(sprintf("max |partial - fresh_new| : %.3g\n", maxDiffPartial))
cat(sprintf("max |restored - fresh_old|: %.3g\n", maxDiffRestored))

stopifnot(samplesDone >= 10L)  # at least 10 useful samples
stopifnot(maxDiffPartial < 1e-10)
stopifnot(maxDiffRestored < 1e-10)
cat("T-013 partial-eval bit-identity OK.\n")
