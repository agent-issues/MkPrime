#!/usr/bin/env Rscript
# 12a_bitident.R — T-012 cached-vs-legacy bit-identity check.
# Builds an MkPrimeData + ecology model, then calls both .CppLogLikelihoodEcology
# (legacy full pruner) and .CppLogLikelihoodEcologyCached (T-012 cache path) on
# the same inputs and asserts |diff| < 1e-10.

suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKPRIME_LIBPATH", "dev/profiling/.vtune-lib-t012")
  library("MkPrime", lib.loc = libPath)
  library("TreeTools")
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# Mirror the test helper from tests/testthat/test-ecology-plumbing.R.
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

set.seed(42)
nTip <- 10L
tips <- paste0("t", seq_len(nTip))
mat <- cbind(
  sample(0:2, nTip, replace = TRUE),
  sample(0:2, nTip, replace = TRUE),
  sample(0:2, nTip, replace = TRUE),
  sample(0:2, nTip, replace = TRUE),
  sample(0:1, nTip, replace = TRUE),
  sample(0:1, nTip, replace = TRUE)
)
rownames(mat) <- tips
pd <- TreeTools::MatrixToPhyDat(mat)

ecoVec <- setNames(sample(0:2, nTip, replace = TRUE), tips)
mkd <- MkPrimeData(pd, neomorphic = c(5L, 6L), ecology = ecoVec)
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

tree <- TreeTools::Preorder(ape::rtree(nTip, tip.label = tips))
parent <- as.integer(tree$edge[, 1])
child  <- as.integer(tree$edge[, 2])
edgeLen <- as.numeric(tree$edge.length)
kPrime <- as.integer(mkd$kObs)
zCols <- as.integer(mkd$kEcology - 1L)
zMat <- matrix(sample(0:2, mkd$nChar * zCols, replace = TRUE),
               nrow = mkd$nChar, ncol = zCols)
storage.mode(zMat) <- "integer"
phi <- 1.6
pi0 <- 0.5
theta <- rep(0.5, zCols)

ll_legacy <- MkPrime:::.CppLogLikelihoodEcology(
  dataPtr, parent, child, edgeLen, kPrime,
  rateLoss = 1.0, rateLogSd = 0.5, rateNeo = 1.0,
  phi = phi, zMatrix = zMat, pi0 = pi0, theta = theta
)
ll_cached <- MkPrime:::.CppLogLikelihoodEcologyCached(
  dataPtr, parent, child, edgeLen, kPrime,
  rateLoss = 1.0, rateLogSd = 0.5, rateNeo = 1.0,
  phi = phi, zMatrix = zMat, pi0 = pi0, theta = theta
)
diff <- abs(ll_legacy - ll_cached)
cat(sprintf("legacy  = %.10f\n", ll_legacy))
cat(sprintf("cached  = %.10f\n", ll_cached))
cat(sprintf("|diff|  = %.3g\n", diff))
stopifnot(diff < 1e-10)
cat("T-012 bit-identity OK.\n")
