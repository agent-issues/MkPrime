# sim1-null.R -----------------------------------------------------------------
# Sim 1: null simulation.  Simulate transformational + neomorphic characters
# on a random 16-tip tree with NO ecology effect (z = 0 everywhere) and a
# random binary tip-ecology assignment.  Run the aware model on the data.
#
# Pass conditions (identifiability / non-pathology):
#   - pi0 posterior stays near prior mean (0.7), CI overlapping it
#   - phi posterior concentrated near 1
#   - posterior recovers the simulating tree (or includes it with high
#     probability), regardless of ecology
#
# Failure of any of these would indicate the prior is too weak or that
# the ecology layer has a sign / scale problem.
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/ecology/simulations/sim3-simulate.R")

set.seed(20260513)
nTip  <- 16L
nEco  <- 30L
nBase <- 90L
phi   <- 1   # neutral; with z=0 everywhere it makes no difference
# Random tree, ultrametric-ish via rcoal so branch lengths are reasonable.
tree <- Preorder(ape::rcoal(nTip,
                            tip.label = paste0("T", seq_len(nTip))))

# Random binary ecology, ~half / half tips.
ecoTip <- setNames(sample(c(0L, 1L), nTip, replace = TRUE),
                   tree$tip.label)
edgeEc <- .AssignEdgeEcology(tree, ecoTip)
cat("Tip ecology:\n"); print(ecoTip)
cat("Edge ecology counts:\n"); print(table(edgeEc))

zNull <- matrix(0L, nrow = nEco + nBase, ncol = 2)
cFull <- c(rep("neomorphic", nEco), rep("transformational", nBase))
datSim <- .SimulateMkPrimeEcology(tree, edgeEc, zNull, phi = phi,
                                  type = cFull, baseRate = 0.5,
                                  rateLoss = 1)
cat("Simulated", nrow(datSim), "tips x", ncol(datSim), "chars\n")

pdSim <- MatrixToPhyDat(datSim)
mkdAware <- MkPrimeData(pdSim, neomorphic = seq_len(nEco),
                        ecology = ecoTip[rownames(datSim)])
cat("kEcology:", mkdAware$kEcology, " | nChar:", mkdAware$nChar, "\n")

set.seed(1)
randTree  <- Preorder(ape::rtree(mkdAware$nTip,
                                 tip.label = mkdAware$taxon_names))
psTrees   <- TreeSearch::MaximizeParsimony(pdSim, tree = randTree,
                                           verbosity = 0)
startTree <- Preorder(psTrees[[1]])
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)
mcmc <- MkPrimeMCMC(nIter = 100000L, nChains = 1L, nRuns = 1L,
                    thin = 200L, treeThin = 200L,
                    minWarmup = 5000L, maxWarmup = 30000L,
                    logFile = NULL, checkpointFile = NULL)

for (f in c("sim1-null.log", "sim1-null.ckp")) {
  if (file.exists(f)) file.remove(f)
}
mcmc$logFile <- "sim1-null.log"
t0 <- Sys.time()
resAware <- RunMkPrime(mkdAware, tree = startTree,
                       model = modelAware, mcmc = mcmc)
cat("Aware elapsed:", format(Sys.time() - t0), "\n")

# Posterior summary.
samp <- ReadMkLog("sim1-null.log")
samp <- samp[seq.int(ceiling(nrow(samp) / 4) + 1L, nrow(samp)), ,
             drop = FALSE]
cat("\n== Aware posterior (null sim: expect phi~1, pi0~0.7) ==\n")
for (nm in c("phi", "pi0", "log_likelihood", "tree_length")) {
  if (nm %in% colnames(samp)) {
    x <- samp[, nm]
    cat(sprintf("  %-15s mean=%.3f  CI=[%.3f, %.3f]\n",
                nm, mean(x),
                stats::quantile(x, 0.025, names = FALSE),
                stats::quantile(x, 0.975, names = FALSE)))
  }
}

# Recovery of the simulating tree.
discard <- function(x) x[seq.int(ceiling(length(x) / 4) + 1L, length(x))]
trees <- discard(resAware$trees)
class(trees) <- "multiPhylo"
cidTrue <- as.numeric(TreeDist::ClusteringInfoDist(
  trees, tree, normalize = TRUE))
cat(sprintf("\n== Tree recovery (n=%d posterior) ==\n", length(trees)))
cat(sprintf("  CID to TRUE: mean=%.3f  median=%.3f  min=%.3f  max=%.3f\n",
            mean(cidTrue), median(cidTrue), min(cidTrue), max(cidTrue)))
cat(sprintf("  Trees identical to TRUE: %d (%.1f%%)\n",
            sum(cidTrue == 0), 100 * mean(cidTrue == 0)))

saveRDS(list(tree = tree, eco = ecoTip, z = zNull, charType = cFull,
             phi = phi, data = datSim, edgeEco = edgeEc,
             resAware = resAware),
        "inst/ecology/simulations/sim1-null-result.rds")
cat("\nDone. Sim 1 null result saved.\n")
