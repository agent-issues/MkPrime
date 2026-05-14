# truth-init-deeper.R — peel back warmup behaviour + recompute truth logL.
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/simulations/ecology/sim3-helpers.R")
source("inst/simulations/ecology/sim3-simulate.R")

res <- readRDS("inst/simulations/ecology/truth-init-results/truth-init-result.rds")
cat("Result names:", paste(names(res), collapse = ", "), "\n")
cat("nSamples:", res$nSamples, "  warmup:", res$warmup, "\n")

if (!is.null(res$warmup_trace)) {
  wt <- res$warmup_trace
  cat("\n=== warmup_trace dim ===\n"); print(dim(wt))
  cat("=== warmup_trace columns ===\n"); print(colnames(wt))
  cat("=== first 5 rows ===\n"); print(head(wt, 5))
  cat("=== last 5 rows ===\n"); print(tail(wt, 5))
}

cat("\n=== Stop reason ===\n"); print(res$stop_reason)
cat("actual_iter:", res$actual_iter, "\n")

# Recompute log_lik + log_prior at truth (with theta=0.999 nudge)
nEco <- 120L; nBase <- 360L; phi <- 4
stemBr <- 0.30; rootBr <- 0.15; tipBr <- 0.5
tree <- .BuildConvergentTree(tipBranch = tipBr,
                              stemBranch = stemBr, rootBranch = rootBr)
eco  <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)
set.seed(20260602L)
zFull <- matrix(0L, nrow = nEco + nBase, ncol = 2)
zFull[seq_len(nEco), 2] <- 1L
cFull <- c(rep("neomorphic", nEco), rep("transformational", nBase))
datSim <- .SimulateMkPrimeEcology(tree, edgeEc, zFull, phi = phi,
                                  type = cFull, baseRate = 1.0,
                                  normalize = TRUE,
                                  pi0 = 0.75, theta = 1.0,
                                  refEcology = 0L, rateLoss = 1)
pdSim <- MatrixToPhyDat(datSim)
mkdBlind <- MkPrimeData(pdSim, neomorphic = seq_len(nEco))
ecoTipVec <- eco[mkdBlind$taxon_names]
mkdAware <- MkPrimeData(pdSim, neomorphic = seq_len(nEco), ecology = ecoTipVec)

trueTree  <- Preorder(tree)
trueZ <- matrix(0L, nrow = nEco + nBase, ncol = 1L)
trueZ[seq_len(nEco), 1L] <- 1L

modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)

logLikTruth <- .MkpEcologyLogLikelihood(
  trueTree, mkdAware,
  kPrime = as.integer(mkdAware$kObs),
  rate_loss = 1.0, rate_log_sd = 0.5,
  nCat = modelAware$nCat, rate_neo = 1.0,
  relabel = modelAware$relabel,
  phi = 4.0, zMat = trueZ,
  magnitudeMode = modelAware$magnitudeMode,
  coding = modelAware$coding,
  refEcology = mkdAware$refEcology,
  theta = 0.999, pi0 = 0.75)
cat("\n=== Truth log_lik (theta=0.999):", logLikTruth, "===\n")

stateTruth <- list(
  tree = trueTree,
  tree_length = sum(trueTree$edge.length),
  rel_br_lengths = trueTree$edge.length / sum(trueTree$edge.length),
  rate_loss = 1.0, rate_log_sd = 0.5, rate_neo = 1.0, p = 0.5,
  kPrime = as.integer(mkdAware$kObs),
  phi = 4.0, pi0 = 0.75, theta = 0.999, z = trueZ
)
logPriorTruth <- LogPrior(stateTruth, modelAware, mkdAware)
cat("Truth log_prior:", logPriorTruth, "\n")
cat("Truth log_post:",  logLikTruth + logPriorTruth, "\n")

# Compare to drifted mode (median of last 50%)
cat("\n=== Drifted mode median log_post (from log): -5117 ===\n")
cat("Gap (drifted - truth):", -5117 - (logLikTruth + logPriorTruth), "nats\n")
cat("Interpretation: negative = truth higher posterior (drift is genuinely worse)\n")
