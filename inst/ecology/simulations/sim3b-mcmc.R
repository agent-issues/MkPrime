# sim3b-mcmc.R ---------------------------------------------------------------
# Sim 3b: negative control for the convergent-clades study.  Same tree,
# same tip-ecology assignment (A and B both in ecology 1), same character
# counts (60 + 180), same phi=4 *but* z = 0 everywhere — meaning NO
# character actually depends on ecology.  A and B share ecology by
# coincidence, not by mechanism.
#
# Tests: does the ecology layer spuriously break the true (A, C)
# bipartition because A and B share ecology?  Answer should be NO —
# blind and aware should give comparable posteriors, both recovering
# the true topology.  pi0 posterior should stay near prior mean (0.7),
# phi posterior near 1.
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/ecology/simulations/sim3-helpers.R")
source("inst/ecology/simulations/sim3-simulate.R")

set.seed(20260512)
nEco <- 60L; nBase <- 180L; stemBr <- 0.10; rootBr <- 0.15
tree   <- .BuildConvergentTree(stemBranch = stemBr, rootBranch = rootBr)
eco    <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)

# z = 0 everywhere: no character depends on ecology.
zNull <- matrix(0L, nrow = nEco + nBase, ncol = 2)
cFull <- c(rep("neomorphic", nEco), rep("transformational", nBase))
# phi value doesn't matter when z is all zero; pass 4 for parity.
datSim <- .SimulateMkPrimeEcology(tree, edgeEc, zNull, phi = 4,
                                  type = cFull, baseRate = 0.5,
                                  rateLoss = 1)
cat("Sim3b: ", nrow(datSim), "tips x", ncol(datSim), "chars\n")
cat("z is all-zero by design (no ecology effect).\n")

treeAB <- Preorder(ape::read.tree(text =
  paste0("(((A1,A2),(A3,A4)),((B1,B2),(B3,B4)),",
         "(((C1,C2),(C3,C4)),((D1,D2),(D3,D4))));")))
pdSim <- MatrixToPhyDat(datSim)
cat("Parsimony: TRUE =", phangorn::parsimony(tree, pdSim),
    " WRONG =", phangorn::parsimony(treeAB, pdSim), "\n")

mkdBlind <- MkPrimeData(pdSim, neomorphic = seq_len(nEco))
ecoTipVec <- eco[mkdBlind$taxon_names]
mkdAware <- MkPrimeData(pdSim, neomorphic = seq_len(nEco),
                        ecology = ecoTipVec)
cat("kEcology:", mkdAware$kEcology, "\n")

set.seed(1)
randTree <- Preorder(ape::rtree(mkdBlind$nTip,
                                tip.label = mkdBlind$taxon_names))
startTree <- Preorder(TreeSearch::AdditionTree(pdSim))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

modelBlind <- MkPrimeModel(ecologyAware = FALSE,
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)
modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)
mcmc <- MkPrimeMCMC(nIter = 200000L, nChains = 1L, nRuns = 1L,
                    thin = 400L, treeThin = 400L,
                    minWarmup = 5000L, maxWarmup = 40000L,
                    logFile = NULL, checkpointFile = NULL)

for (f in c("sim3b-blind.log", "sim3b-blind.ckp",
            "sim3b-aware.log", "sim3b-aware.ckp")) {
  if (file.exists(f)) file.remove(f)
}

cat("\n== Blind run ==\n")
m <- mcmc; m$logFile <- "sim3b-blind.log"
t0 <- Sys.time()
resBlind <- RunMkPrime(mkdBlind, tree = startTree,
                       model = modelBlind, mcmc = m)
cat("Blind elapsed:", format(Sys.time() - t0), "\n")

cat("\n== Aware run ==\n")
m <- mcmc; m$logFile <- "sim3b-aware.log"
t0 <- Sys.time()
resAware <- RunMkPrime(mkdAware, tree = startTree,
                       model = modelAware, mcmc = m)
cat("Aware elapsed:", format(Sys.time() - t0), "\n")

# Bipartition support.
# Use root-invariant HasBipartSplits() — see sim3-scoring.R.
source("inst/ecology/simulations/sim3-scoring.R")
trueSplit  <- c(paste0("A", 1:4), paste0("C", 1:4))
wrongSplit <- c(paste0("A", 1:4), paste0("B", 1:4))
hasBipart <- HasBipartSplits
discard <- function(x) x[seq.int(ceiling(length(x) / 4) + 1L, length(x))]

cat("\n== Posterior bipartition support ==\n")
for (nm in c("blind", "aware")) {
  trees <- if (nm == "blind") resBlind$trees else resAware$trees
  trees <- discard(trees)
  pTrue  <- mean(hasBipart(trees, trueSplit))
  pWrong <- mean(hasBipart(trees, wrongSplit))
  cat(sprintf("  %-5s : n=%d trees  P(true)=%.3f  P(wrong)=%.3f\n",
              nm, length(trees), pTrue, pWrong))
}

cat("\n== Aware posterior phi / pi0 (should stay near prior under H0) ==\n")
samp <- ReadMkLog("sim3b-aware.log")
samp <- samp[seq.int(ceiling(nrow(samp) / 4) + 1L, nrow(samp)), ,
             drop = FALSE]
for (nm in c("phi", "pi0", "log_likelihood", "tree_length")) {
  if (nm %in% colnames(samp)) {
    x <- samp[, nm]
    cat(sprintf("  %-15s mean=%.3f  CI=[%.3f, %.3f]\n",
                nm, mean(x),
                stats::quantile(x, 0.025, names = FALSE),
                stats::quantile(x, 0.975, names = FALSE)))
  }
}

saveRDS(list(tree = tree, eco = eco, z = zNull, charType = cFull,
             phi = 4, baseRate = 0.5, stem = stemBr, root = rootBr,
             data = datSim, edgeEco = edgeEc,
             resBlind = resBlind, resAware = resAware),
        "inst/ecology/simulations/sim3b-mcmc-result.rds")
cat("\nDone. Sim 3b result saved.\n")
