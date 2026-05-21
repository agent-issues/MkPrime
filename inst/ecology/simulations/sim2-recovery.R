# sim2-recovery.R -------------------------------------------------------------
# Sim 2: parameter recovery.  Simulate on a random 16-tip tree with a
# random binary tip ecology.  Half the characters are ecology-driven
# (z = 1 in ecology 1, z = 0 in ecology 0); the other half are z = 0
# everywhere.  Simulate with phi = 3, then run the aware model.
#
# Pass conditions:
#   - phi posterior covers ~3 (or, given the relabel non-identifiability
#     for neomorphic chars, the |log phi| posterior is well away from 0)
#   - pi0 posterior is well below the prior mean (0.7), reflecting the
#     ~50% of characters with z != 0
#   - posterior recovers the simulating tree
#
# This is a complementary check to Sim 3: there is no convergent
# ecology pattern here, so the *blind* model is also expected to do
# fine on tree recovery; the diagnostic value is on phi and pi0.
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/ecology/simulations/sim3-simulate.R")

set.seed(20260514)
nTip  <- 16L
phi   <- 3
# Use roughly the same total char count as Sim 1, but balanced across
# z = 0 / z = 1.
nEcoDriven   <- 60L  # transformational chars with z = 1 in eco 1
nNeutralTr   <- 60L  # transformational chars with z = 0 everywhere
nNeutralNeo  <- 30L  # neomorphic chars with z = 0 everywhere
nChar <- nEcoDriven + nNeutralTr + nNeutralNeo
tree <- Preorder(ape::rcoal(nTip,
                            tip.label = paste0("T", seq_len(nTip))))

# Random binary ecology, roughly balanced.
repeat {
  ecoTip <- setNames(sample(c(0L, 1L), nTip, replace = TRUE),
                     tree$tip.label)
  if (sum(ecoTip == 1L) >= 5L && sum(ecoTip == 0L) >= 5L) break
}
edgeEc <- .AssignEdgeEcology(tree, ecoTip)
cat("Tip ecology counts:\n"); print(table(ecoTip))
cat("Edge ecology counts:\n"); print(table(edgeEc))

# Build z matrix.  Neomorphic (first nNeutralNeo) entries have z = 0
# everywhere.  Then 60 transformational chars with z[c, eco=1] = 1.
# Then 60 transformational neutral.  Ordering matters because
# MkPrimeData expects neomorphic chars at the leading indices.
zMat <- matrix(0L, nrow = nChar, ncol = 2)
neoIdx <- seq_len(nNeutralNeo)
ecoIdx <- nNeutralNeo + seq_len(nEcoDriven)
neuIdx <- nNeutralNeo + nEcoDriven + seq_len(nNeutralTr)
zMat[ecoIdx, 2] <- 1L  # column 2 = ecology state 1
cType <- character(nChar)
cType[neoIdx] <- "neomorphic"
cType[ecoIdx] <- "transformational"
cType[neuIdx] <- "transformational"

datSim <- .SimulateMkPrimeEcology(tree, edgeEc, zMat, phi = phi,
                                  type = cType, baseRate = 0.5,
                                  rateLoss = 1)
cat("Simulated", nrow(datSim), "tips x", ncol(datSim), "chars\n")

pdSim <- MatrixToPhyDat(datSim)
mkdAware <- MkPrimeData(pdSim, neomorphic = neoIdx,
                        ecology = ecoTip[rownames(datSim)])
cat("kEcology:", mkdAware$kEcology, " | nChar:", mkdAware$nChar, "\n")

set.seed(1)
randTree  <- Preorder(ape::rtree(mkdAware$nTip,
                                 tip.label = mkdAware$taxon_names))
startTree <- Preorder(TreeSearch::AdditionTree(pdSim))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)
mcmc <- MkPrimeMCMC(nIter = 100000L, nChains = 1L, nRuns = 1L,
                    thin = 200L, treeThin = 200L,
                    minWarmup = 5000L, maxWarmup = 30000L,
                    logFile = NULL, checkpointFile = NULL)

for (f in c("sim2-recovery.log", "sim2-recovery.ckp")) {
  if (file.exists(f)) file.remove(f)
}
mcmc$logFile <- "sim2-recovery.log"
t0 <- Sys.time()
resAware <- RunMkPrime(mkdAware, tree = startTree,
                       model = modelAware, mcmc = mcmc)
cat("Aware elapsed:", format(Sys.time() - t0), "\n")

samp <- ReadMkLog("sim2-recovery.log")
samp <- samp[seq.int(ceiling(nrow(samp) / 4) + 1L, nrow(samp)), ,
             drop = FALSE]
cat("\n== Aware posterior (recovery sim: expect |log phi| away from 0, pi0 < 0.7) ==\n")
for (nm in c("phi", "pi0", "log_likelihood", "tree_length")) {
  if (nm %in% colnames(samp)) {
    x <- samp[, nm]
    cat(sprintf("  %-15s mean=%.3f  CI=[%.3f, %.3f]\n",
                nm, mean(x),
                stats::quantile(x, 0.025, names = FALSE),
                stats::quantile(x, 0.975, names = FALSE)))
  }
}
if ("phi" %in% colnames(samp)) {
  lx <- abs(log(samp[, "phi"]))
  cat(sprintf("  |log phi|       mean=%.3f  CI=[%.3f, %.3f]\n",
              mean(lx),
              stats::quantile(lx, 0.025, names = FALSE),
              stats::quantile(lx, 0.975, names = FALSE)))
}

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

saveRDS(list(tree = tree, eco = ecoTip, z = zMat, charType = cType,
             phi = phi, data = datSim, edgeEco = edgeEc,
             neoIdx = neoIdx, ecoIdx = ecoIdx, neuIdx = neuIdx,
             resAware = resAware),
        "inst/ecology/simulations/sim2-recovery-result.rds")
cat("\nDone. Sim 2 recovery result saved.\n")
