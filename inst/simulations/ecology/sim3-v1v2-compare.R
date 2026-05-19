# sim3-v1v2-compare.R --------------------------------------------------------
# Direct v1 vs v2 comparison driver for the convergent-ecology Sim 3.
#
# Re-runs the same fixed dataset under both the legacy v1 parameterisation
# (un-normalised rate factors, symmetric slab, K columns of z) and the v2
# parameterisation (gamma-normalised, reference ecology, asymmetric slab,
# K-1 columns of z). Reports head-to-head:
#   * Mixing diagnostics: minESS, logL median/IQR, tree_length stability
#   * Posterior on hyperparameters: phi, pi0 (v2 adds theta)
#   * Bipartition support: P(true), P(wrong); CID gap to truth
#
# Headline expectations at the Goldilocks config (nEco=60, nBase=180,
# phi=4, stem=0.10, root=0.15):
#   * v1: aware logL << blind logL (rate-time identifiability ridge),
#     tree_length blows up, posterior diffuse
#   * v2: aware logL >= blind logL (ridge closed), tree_length tight,
#     posterior concentrated near truth or split between truth and wrong
#
# Requires v2 to be working. v1 results are loaded from the existing
# sim3-mcmc-result.rds (already on disk). v2 runs fresh.
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
  library("TreeDist")
})
source("inst/simulations/ecology/sim3-helpers.R")
source("inst/simulations/ecology/sim3-simulate.R")

args <- commandArgs(trailingOnly = TRUE)
nIter <- if (length(args) >= 1) as.integer(args[1]) else 200000L

# Same fixed dataset used by sim3-mcmc.R: seed 20260512, Goldilocks config.
set.seed(20260512)
nEco <- 60L; nBase <- 180L; phi <- 4; stemBr <- 0.10; rootBr <- 0.15
tree   <- .BuildConvergentTree(stemBranch = stemBr, rootBranch = rootBr)
eco    <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)

zFull <- matrix(0L, nrow = nEco + nBase, ncol = 2)
zFull[seq_len(nEco), 2] <- 1L
cFull <- c(rep("neomorphic", nEco), rep("transformational", nBase))
datSim <- .SimulateMkPrimeEcology(tree, edgeEc, zFull, phi = phi,
                                  type = cFull, baseRate = 0.5,
                                  rateLoss = 1)
pdSim <- MatrixToPhyDat(datSim)
mkdBlind <- MkPrimeData(pdSim, neomorphic = seq_len(nEco))
ecoTipVec <- eco[mkdBlind$taxon_names]
mkd <- MkPrimeData(pdSim, neomorphic = seq_len(nEco), ecology = ecoTipVec)

# Parsimony start tree (matches sim3-mcmc.R).
set.seed(1)
randTree <- Preorder(ape::rtree(mkd$nTip, tip.label = mkd$taxon_names))
psTrees  <- TreeSearch::MaximizeParsimony(pdSim, tree = randTree,
                                           verbosity = 0)
startTree <- Preorder(psTrees[[1]])
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

modelBlind <- MkPrimeModel(ecologyAware = FALSE,
                           kPrimePrior = "geometric", coding = "variable")
modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           kPrimePrior = "geometric", coding = "variable")
mcmc <- MkPrimeMCMC(nIter = nIter, nChains = 1L, nRuns = 1L,
                    thin = max(1L, nIter %/% 500L),
                    treeThin = max(1L, nIter %/% 500L),
                    minWarmup = 5000L, maxWarmup = 40000L,
                    logFile = NULL, checkpointFile = NULL)

cat("== Running v2 blind chain ==\n")
m <- mcmc; m$logFile <- "sim3-v2-blind.log"
if (file.exists(m$logFile)) file.remove(m$logFile)
t0 <- Sys.time()
resBlind <- RunMkPrime(mkdBlind, tree = startTree,
                       model = modelBlind, mcmc = m)
cat("Blind elapsed:", format(Sys.time() - t0), "\n")

cat("\n== Running v2 aware chain ==\n")
m <- mcmc; m$logFile <- "sim3-v2-aware.log"
if (file.exists(m$logFile)) file.remove(m$logFile)
t0 <- Sys.time()
resAware <- RunMkPrime(mkd, tree = startTree,
                       model = modelAware, mcmc = m)
cat("Aware elapsed:", format(Sys.time() - t0), "\n")

# Compare to v1 result if present.
v1File <- "inst/simulations/ecology/sim3-mcmc-result.rds"
v1Exists <- file.exists(v1File)
if (v1Exists) {
  v1 <- readRDS(v1File)
} else {
  cat("\nNote: v1 sim3-mcmc-result.rds not found; will report v2 only.\n")
}

# Diagnostic summary.
discard <- function(x) x[seq.int(ceiling(length(x) / 4) + 1L, length(x))]
samples <- function(logFile) {
  s <- ReadMkLog(logFile)
  s[discard(seq_len(nrow(s))), , drop = FALSE]
}
sB2 <- samples("sim3-v2-blind.log")
sA2 <- samples("sim3-v2-aware.log")

cat("\n== Hyperparameter posterior (v2 aware) ==\n")
for (nm in c("phi", "pi0", "theta_1", "log_likelihood", "tree_length")) {
  if (nm %in% colnames(sA2)) {
    x <- sA2[, nm]
    cat(sprintf("  %-15s mean=%.3f  CI=[%.3f, %.3f]\n",
                nm, mean(x),
                stats::quantile(x, 0.025, names = FALSE),
                stats::quantile(x, 0.975, names = FALSE)))
  }
}

cat("\n== logL distribution: aware vs blind (v2) ==\n")
cat(sprintf("  blind  median=%.0f  IQR=[%.0f, %.0f]\n",
            median(sB2[, "log_likelihood"]),
            stats::quantile(sB2[, "log_likelihood"], 0.25, names = FALSE),
            stats::quantile(sB2[, "log_likelihood"], 0.75, names = FALSE)))
cat(sprintf("  aware  median=%.0f  IQR=[%.0f, %.0f]\n",
            median(sA2[, "log_likelihood"]),
            stats::quantile(sA2[, "log_likelihood"], 0.25, names = FALSE),
            stats::quantile(sA2[, "log_likelihood"], 0.75, names = FALSE)))
gap <- median(sA2[, "log_likelihood"]) - median(sB2[, "log_likelihood"])
cat(sprintf("  aware - blind median logL gap: %.1f  (positive = aware better)\n",
            gap))

cat("\n== Tree length stability ==\n")
cat(sprintf("  blind  median=%.2f  IQR=[%.2f, %.2f]  max=%.2f\n",
            median(sB2[, "tree_length"]),
            stats::quantile(sB2[, "tree_length"], 0.25, names = FALSE),
            stats::quantile(sB2[, "tree_length"], 0.75, names = FALSE),
            max(sB2[, "tree_length"])))
cat(sprintf("  aware  median=%.2f  IQR=[%.2f, %.2f]  max=%.2f\n",
            median(sA2[, "tree_length"]),
            stats::quantile(sA2[, "tree_length"], 0.25, names = FALSE),
            stats::quantile(sA2[, "tree_length"], 0.75, names = FALSE),
            max(sA2[, "tree_length"])))

if (v1Exists) {
  cat("\n== v1 reference (for comparison) ==\n")
  v1samp <- ReadMkLog("sim3-aware.log")
  v1samp <- v1samp[discard(seq_len(nrow(v1samp))), , drop = FALSE]
  cat(sprintf("  v1 aware logL median=%.0f  tree_length median=%.2f  max=%.2f\n",
              median(v1samp[, "log_likelihood"]), median(v1samp[, "tree_length"]),
              max(v1samp[, "tree_length"])))
  v1B <- ReadMkLog("sim3-blind.log")
  v1B <- v1B[discard(seq_len(nrow(v1B))), , drop = FALSE]
  v1gap <- median(v1samp[, "log_likelihood"]) - median(v1B[, "log_likelihood"])
  cat(sprintf("  v1 (aware - blind) median logL gap: %.1f\n", v1gap))
  cat(sprintf("  v2 (aware - blind) median logL gap: %.1f\n", gap))
  cat("  v2 should have a non-negative gap (aware nests blind);\n")
  cat("  v1 had aware ~76 nats worse than blind from the rate-time ridge.\n")
}

# Bipartition support.
# Use root-invariant HasBipartSplits() — see sim3-scoring.R.
source("inst/simulations/ecology/sim3-scoring.R")
trueSplit  <- c(paste0("A", 1:4), paste0("C", 1:4))
wrongSplit <- c(paste0("A", 1:4), paste0("B", 1:4))
hasBipart <- HasBipartSplits
for (nm in c("blind", "aware")) {
  trees <- discard(if (nm == "blind") resBlind$trees else resAware$trees)
  pTrue  <- mean(hasBipart(trees, trueSplit))
  pWrong <- mean(hasBipart(trees, wrongSplit))
  cat(sprintf("\nv2 %s : P(true)=%.3f  P(wrong)=%.3f  (n=%d)\n",
              nm, pTrue, pWrong, length(trees)))
}

saveRDS(list(tree = tree, eco = eco, z = zFull, charType = cFull,
             phi = phi, data = datSim, edgeEco = edgeEc,
             resBlind = resBlind, resAware = resAware,
             samplesBlind = sB2, samplesAware = sA2),
        "inst/simulations/ecology/sim3-v2-result.rds")
cat("\nDone. v2 result saved.\n")
