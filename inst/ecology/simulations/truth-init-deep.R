# truth-init-deep.R — dig into why log_posterior = -Inf at truth.
# Loads the result, isolates which prior component is -Inf, and reports
# whether the chain state actually changed (or whether the log columns
# are stale).
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/ecology/simulations/sim3-helpers.R")
source("inst/ecology/simulations/sim3-simulate.R")

res <- readRDS("inst/ecology/simulations/truth-init-results/truth-init-result.rds")
cat("Class:", class(res), "\n")
cat("Names:", paste(names(res), collapse = ", "), "\n")
cat("Samples dim: ", paste(dim(res$samples), collapse = " x "),
    if (!is.null(res$samples)) "" else " (NULL)", "\n")
cat("trees length:", length(res$trees), "\n")

cat("\n=== Last sampled state (chain end) ===\n")
endState <- res$chains[[1]]
if (is.null(endState)) {
  endState <- res$chainStates[[1]]
}
cat("Names of end state:\n"); print(names(endState))
cat("phi:", endState$phi, "\n")
cat("pi0:", endState$pi0, "\n")
cat("theta:", endState$theta, "\n")
cat("rate_loss:", endState$rate_loss, "\n")
cat("rate_neo:", endState$rate_neo, "\n")
cat("tree_length:", endState$tree_length, "\n")
cat("log_lik:", endState$log_lik, "\n")
cat("log_prior:", endState$log_prior, "\n")
cat("log_post:", endState$log_post, "\n")

if (!is.null(endState$z)) {
  zMat <- endState$z
  cat("z table:\n"); print(table(zMat))
}

# Rebuild data + model to invoke LogPrior cleanly
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

modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)

# Construct truth state matching what initOverrides would build
trueTree  <- Preorder(tree)
trueTL    <- sum(trueTree$edge.length)
trueRelBr <- trueTree$edge.length / trueTL
trueZ <- matrix(0L, nrow = nEco + nBase, ncol = 1L)
trueZ[seq_len(nEco), 1L] <- 1L

stateTruth <- list(
  tree = trueTree,
  tree_length = trueTL,
  rel_br_lengths = trueRelBr,
  rate_loss = 1.0,
  rate_log_sd = 0.5,
  rate_neo = 1.0,
  p = 0.5,
  kPrime = as.integer(mkdAware$kObs),
  phi = 4.0,
  pi0 = 0.75,
  theta = 1.0,
  z = trueZ
)

cat("\n=== LogPrior at truth (manual call) ===\n")
lp <- LogPrior(stateTruth, modelAware, mkdAware)
cat("Total log_prior:", lp, "\n")

# Decompose by component if accessible
if (exists("LogPrior")) {
  # Try invoking sub-priors if defined
  cat("\nPrior contributions (try to manually decompose):\n")
  cat("  pi0~Beta(7,3): ", dbeta(0.75, 7, 3, log = TRUE), "\n")
  cat("  phi~LN(0,1):   ", dlnorm(4.0, 0, 1, log = TRUE), "\n")
  cat("  theta~U(0,1):  ", dbeta(1.0, 1, 1, log = TRUE), "\n")
  cat("  z (cells = 480):\n")
  pi0v <- 0.75; thetaV <- 1.0
  pNone <- pi0v
  pEnc  <- (1 - pi0v) * thetaV
  pDisc <- (1 - pi0v) * (1 - thetaV)
  cat("    P(none)=", pNone, ", P(enc)=", pEnc, ", P(disc)=", pDisc, "\n")
  zCounts <- table(factor(trueZ, levels = 0:2))
  cat("    counts (0/1/2):", zCounts, "\n")
  cat("    log_prior_z = nNone*log(pNone) + nEnc*log(pEnc) =",
      zCounts[1] * log(pNone) + zCounts[2] * log(pEnc), "\n")
  cat("    (z=2 would give log(0)=-Inf; but truth has none)\n")
}
