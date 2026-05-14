# Sim 3 v3 Bayesian pilot — aware vs blind at top redesign candidate
#
# Candidate: 2x-3x-stem
#   nEco = 120, nBase = 360 (480 chars total — 2x v2 Goldilocks)
#   stem = 0.30 (3x v2), root = 0.15, tipBranch = 0.5
#   phi = 4, tree = ((A,C),(B,D)) convergent
#
# Two chains:
#   - blind: standard MkPrime, no ecology
#   - aware: ecology-aware MkPrime, magnitudeMode = "global"
#
# Each runs from a parsimony start tree (not truth-init — that comes
# later as a separate diagnostic). 100k iter; aligned simulator
# (baseRate=1.0, normalize=TRUE, pi0=0.75, theta=1.0).
#
# Success criteria (paper-relevant):
#   1. blind chain shows P(wrong AB bipartition) > P(true AC) — fooled
#   2. aware chain shows P(true AC) > P(wrong AB) — rescued
#   3. aware CID-to-truth lower than blind
#   4. Posterior tl/phi/pi0 nearer truth on aware than blind
#
# Reality check (mixing):
#   - ESS of phi, pi0, theta, tl should be > 100 (we had ~0 on v2)
#   - Modal topology should hold >5% of samples (was 0.4% on v2)
#
# Run from worktree root:
#   Rscript -e "setwd('C:/Users/pjjg18/GitHub/mkp/.claude/worktrees/ecology-aware'); source('dev/pilots/2026-05-14-sim3-v3-redesign/bayesian-pilot.R')"

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/simulations/ecology/sim3-helpers.R")
source("inst/simulations/ecology/sim3-simulate.R")

PILOT_DIR <- "dev/pilots/2026-05-14-sim3-v3-redesign"

# v3 candidate
nEco <- 120L
nBase <- 360L
phi <- 4
stemBr <- 0.30
rootBr <- 0.15
tipBr  <- 0.5

tree <- .BuildConvergentTree(tipBranch = tipBr,
                             stemBranch = stemBr,
                             rootBranch = rootBr)
eco <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)
truthTL <- sum(tree$edge.length)
cat(sprintf("Sim 3 v3 candidate: 2x-3x-stem\n"))
cat(sprintf("  nChar = %d (neo %d + trans %d)\n", nEco + nBase, nEco, nBase))
cat(sprintf("  stem = %.2f  root = %.2f  tip = %.2f\n", stemBr, rootBr, tipBr))
cat(sprintf("  truth tl = %.2f  phi = %g\n", truthTL, phi))

# Simulate aligned-units data
set.seed(20260514)
zFull <- matrix(0L, nrow = nEco + nBase, ncol = 2)
zFull[seq_len(nEco), 2] <- 1L
cFull <- c(rep("neomorphic", nEco), rep("transformational", nBase))
datSim <- .SimulateMkPrimeEcology(tree, edgeEc, zFull, phi = phi,
                                  type = cFull, baseRate = 1.0,
                                  normalize = TRUE,
                                  pi0 = 0.75, theta = 1.0,
                                  refEcology = 0L,
                                  rateLoss = 1)
pdSim <- MatrixToPhyDat(datSim)
mkdBlind <- MkPrimeData(pdSim, neomorphic = seq_len(nEco))
ecoTipVec <- eco[mkdBlind$taxon_names]
mkdAware <- MkPrimeData(pdSim, neomorphic = seq_len(nEco),
                        ecology = ecoTipVec)

# Parsimony start tree (NOT truth-init)
set.seed(1)
randTree <- Preorder(ape::rtree(mkdBlind$nTip,
                                tip.label = mkdBlind$taxon_names))
psTrees <- TreeSearch::MaximizeParsimony(pdSim, tree = randTree,
                                         verbosity = 0)
startTree <- Preorder(psTrees[[1]])
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

modelBlind <- MkPrimeModel(ecologyAware = FALSE,
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)
modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)

nIter <- 100000L
thin <- 40L
mcmc <- MkPrimeMCMC(nIter = nIter, nChains = 1L, nRuns = 1L,
                    thin = thin, treeThin = thin,
                    minWarmup = 2000L, maxWarmup = 5000L,
                    logFile = NULL, checkpointFile = NULL)

# Clean
for (sub in c("blind", "aware")) {
  for (ext in c(".log", ".ckp")) {
    f <- file.path(PILOT_DIR, paste0(sub, "-chain", ext))
    if (file.exists(f)) file.remove(f)
  }
}

# ---- Blind chain ----
cat("\n== Running BLIND chain ==\n")
mB <- mcmc
mB$logFile <- file.path(PILOT_DIR, "blind-chain.log")
mB$checkpointFile <- file.path(PILOT_DIR, "blind-chain.ckp")
t0 <- Sys.time()
resBlind <- RunMkPrime(mkdBlind, tree = startTree,
                       model = modelBlind, mcmc = mB)
cat("Elapsed (blind):", format(Sys.time() - t0), "\n")

# ---- Aware chain ----
cat("\n== Running AWARE chain ==\n")
mA <- mcmc
mA$logFile <- file.path(PILOT_DIR, "aware-chain.log")
mA$checkpointFile <- file.path(PILOT_DIR, "aware-chain.ckp")
t0 <- Sys.time()
resAware <- RunMkPrime(mkdAware, tree = startTree,
                       model = modelAware, mcmc = mA)
cat("Elapsed (aware):", format(Sys.time() - t0), "\n")

# Apply RelabelEcology to aware
resAwareRl <- RelabelEcology(resAware)

# Save
saveRDS(list(
  config = list(nEco = nEco, nBase = nBase, phi = phi,
                stemBr = stemBr, rootBr = rootBr, tipBr = tipBr,
                nIter = nIter, thin = thin),
  tree = tree, eco = eco, edgeEco = edgeEc,
  data = datSim, z = zFull, charType = cFull,
  truthTL = truthTL,
  resBlind = resBlind, resAware = resAwareRl
), file.path(PILOT_DIR, "result.rds"))
cat(sprintf("\nSaved to %s/result.rds\n", PILOT_DIR))
cat("Done. Run analyse.R for diagnostics.\n")
