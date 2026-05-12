# sim3-mcmc.R -----------------------------------------------------------------
# Headline Sim 3 demonstration.  Generates one simulated dataset at the
# Goldilocks cell identified in sim3-pilot.R:
#
#   nEco = 60 neomorphic, nBase = 180 transformational,
#   phi  = 4, stem = 0.10, root = 0.15.
#
# Runs two MCMC chains on the same data:
#   * blind:  MkPrimeModel(ecologyAware = FALSE)
#   * aware:  MkPrimeModel(ecologyAware = TRUE,  magnitudeMode = "global")
#
# Reports, for each model:
#   - posterior probability of the TRUE bipartition (A1..A4, C1..C4)
#   - posterior probability of the WRONG bipartition (A1..A4, B1..B4)
#   - chain log-likelihood summary, phi / pi0 summaries (aware only)
#
# The expectation: blind concentrates posterior mass on the wrong tree;
# aware spreads mass to include the truth (or recovers it outright).
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/simulations/ecology/sim3-helpers.R")
source("inst/simulations/ecology/sim3-simulate.R")

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

cat("Simulated", nrow(datSim), "tips x", ncol(datSim), "chars\n")

# Parsimony sanity check before running MCMC.
treeAB <- Preorder(ape::read.tree(text =
  paste0("(((A1,A2),(A3,A4)),((B1,B2),(B3,B4)),",
         "(((C1,C2),(C3,C4)),((D1,D2),(D3,D4))));")))
pdSim <- MatrixToPhyDat(datSim)
cat("Parsimony: TRUE =", phangorn::parsimony(tree, pdSim),
    " WRONG =", phangorn::parsimony(treeAB, pdSim), "\n")

# Tip ecology as a vector for the ecology-aware model.
mkdBlind <- MkPrimeData(pdSim,
                        neomorphic = seq_len(nEco))
ecoTipVec <- eco[mkdBlind$taxon_names]
mkdAware <- MkPrimeData(pdSim,
                        neomorphic = seq_len(nEco),
                        ecology = ecoTipVec)
cat("Blind nChar:", mkdBlind$nChar,
    " | Aware nChar:", mkdAware$nChar,
    " | kEcology:", mkdAware$kEcology, "\n")

# Parsimony starting tree (much closer to posterior modes than random;
# a random rtree leaves the chain wandering through tree space for many
# tens of thousands of iterations before locating the high-likelihood
# region).
set.seed(1)
randTree  <- Preorder(ape::rtree(mkdBlind$nTip,
                                 tip.label = mkdBlind$taxon_names))
# Use TreeSearch::MaximizeParsimony instead of phangorn::optim.parsimony:
# the latter triggers a reorderRcpp type-mismatch crash with the installed
# phangorn 2.12.1 + ape 5.8.1 combination on Windows.
psTrees   <- TreeSearch::MaximizeParsimony(pdSim, tree = randTree,
                                           verbosity = 0)
startTree <- Preorder(psTrees[[1]])
# MaximizeParsimony returns trees without branch lengths; RunMkPrime requires
# them.  Initialise to a small uniform value -- MCMC will quickly relax.
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

modelBlind <- MkPrimeModel(
  ecologyAware = FALSE,
  kPrimePrior  = "geometric",
  coding       = "variable",
  expSteps     = 10
)
modelAware <- MkPrimeModel(
  ecologyAware  = TRUE,
  magnitudeMode = "global",
  kPrimePrior   = "geometric",
  coding        = "variable",
  expSteps      = 10
)
mcmc <- MkPrimeMCMC(
  nIter          = 50000L,
  nChains        = 1L,
  nRuns          = 1L,
  thin           = 100L,
  treeThin       = 100L,
  minWarmup      = 2000L,
  maxWarmup      = 15000L,
  logFile        = NULL,
  checkpointFile = NULL
)

cat("\n== Blind run ==\n")
for (f in c("sim3-blind.log", "sim3-blind.ckp")) {
  if (file.exists(f)) file.remove(f)
}
mcmcBlind <- mcmc; mcmcBlind$logFile <- "sim3-blind.log"
t0 <- Sys.time()
resBlind <- RunMkPrime(mkdBlind, tree = startTree,
                       model = modelBlind, mcmc = mcmcBlind)
cat("Blind elapsed:", format(Sys.time() - t0), "\n")

cat("\n== Aware run ==\n")
for (f in c("sim3-aware.log", "sim3-aware.ckp")) {
  if (file.exists(f)) file.remove(f)
}
mcmcAware <- mcmc; mcmcAware$logFile <- "sim3-aware.log"
t0 <- Sys.time()
resAware <- RunMkPrime(mkdAware, tree = startTree,
                       model = modelAware, mcmc = mcmcAware)
cat("Aware elapsed:", format(Sys.time() - t0), "\n")

# Score posterior on key bipartitions.
trueSplit  <- c(paste0("A", 1:4), paste0("C", 1:4))
wrongSplit <- c(paste0("A", 1:4), paste0("B", 1:4))
hasBipart <- function(treeList, splitTips) {
  vapply(treeList, function(tr) {
    cl <- ape::prop.part(tr)
    tips <- attr(cl, "labels")
    splitSet <- which(tips %in% splitTips)
    any(vapply(cl, function(p) setequal(p, splitSet), logical(1)))
  }, logical(1))
}
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

# Aware-specific posterior summary.
cat("\n== Aware posterior phi / pi0 ==\n")
samp <- ReadMkLog("sim3-aware.log")
samp <- samp[seq.int(ceiling(nrow(samp) / 4) + 1L, nrow(samp)), , drop = FALSE]
for (nm in c("phi", "pi0", "log_likelihood", "tree_length")) {
  if (nm %in% colnames(samp)) {
    x <- samp[, nm]
    cat(sprintf("  %-15s mean=%.3f  CI=[%.3f, %.3f]\n",
                nm, mean(x),
                stats::quantile(x, 0.025, names = FALSE),
                stats::quantile(x, 0.975, names = FALSE)))
  }
}

saveRDS(list(tree = tree, eco = eco, z = zFull, charType = cFull,
             phi = phi, baseRate = 0.5, stem = stemBr, root = rootBr,
             data = datSim, edgeEco = edgeEc,
             resBlind = resBlind, resAware = resAware),
        "inst/simulations/ecology/sim3-mcmc-result.rds")
cat("\nDone. Result saved.\n")
