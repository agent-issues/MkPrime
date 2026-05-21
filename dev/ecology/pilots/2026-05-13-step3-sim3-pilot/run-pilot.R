# run-pilot.R ------------------------------------------------------------------
# Step 3 pilot: 20k-iteration Sim 3 aware chain with the two bug fixes:
#   Bug A: asymmetric-slab prior (cpp_log_prior + LogPrior)
#   Bug B: theta now logged in .ParamNames / .StateToRow
#
# Config:
#   seed=20260512 (canonical pre-fix seed, direct comparability)
#   kEco=2 (A/B=ecology1, C/D=ecology0), nTheta=1
#   nEco=60 neomorphic, nBase=180 transformational
#   phi=4, stem=0.10, root=0.15 (Goldilocks)
#   nIter=20000, thin=40 → ~500 raw samples, ~375 post-25%-discard
#   Single chain, aware only (diagnostic, not comparative)
#
# Run from the worktree root:
#   Rscript -e "setwd('C:/Users/pjjg18/GitHub/mkp/.claude/worktrees/ecology-aware'); source('dev/pilots/2026-05-13-step3-sim3-pilot/run-pilot.R')"

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/ecology/simulations/sim3-helpers.R")
source("inst/ecology/simulations/sim3-simulate.R")

PILOT_DIR <- "dev/pilots/2026-05-13-step3-sim3-pilot"

# Canonical seed (same as sim3-v1v2-compare.R for direct comparability)
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
cat("kEcology:", mkd$kEcology, "  nTheta:", mkd$kEcology - 1L, "\n")

# Parsimony start tree
set.seed(1)
randTree <- Preorder(ape::rtree(mkd$nTip, tip.label = mkd$taxon_names))
startTree <- Preorder(TreeSearch::AdditionTree(MatrixToPhyDat(datSim)))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

modelAware <- MkPrimeModel(
  ecologyAware  = TRUE,
  magnitudeMode = "global",
  kPrimePrior   = "geometric",
  coding        = "variable",
  expSteps      = 10
)

# 20k iterations, thin=40 -> ~500 rows; 25%-discard leaves ~375
nIter <- 20000L
thin  <- 40L
mcmc  <- MkPrimeMCMC(
  nIter          = nIter,
  nChains        = 1L,
  nRuns          = 1L,
  thin           = thin,
  treeThin       = thin,
  minWarmup      = 2000L,
  maxWarmup      = 5000L,
  logFile        = NULL,
  checkpointFile = NULL
)
logFile  <- file.path(PILOT_DIR, "sim3-pilot.log")
ckpFile  <- file.path(PILOT_DIR, "sim3-pilot.ckp")

# Clean any previous partial run
for (f in c(logFile, ckpFile)) if (file.exists(f)) file.remove(f)

mcmc$logFile        <- logFile
mcmc$checkpointFile <- ckpFile

cat("== Running Step 3 pilot aware chain (nIter=", nIter, ", thin=", thin, ") ==\n")
t0 <- Sys.time()
resAware <- RunMkPrime(mkd, tree = startTree, model = modelAware, mcmc = mcmc)
elapsed <- Sys.time() - t0
cat("Elapsed:", format(elapsed), "\n")

# Save RDS
rdsFile <- file.path(PILOT_DIR, "sim3-pilot-result.rds")
saveRDS(list(
  config = list(seed = 20260512, nEco = nEco, nBase = nBase, phi = phi,
                stemBr = stemBr, rootBr = rootBr, nIter = nIter, thin = thin,
                kEcology = mkd$kEcology),
  tree = tree, eco = eco, z = zFull, charType = cFull,
  phi_truth = phi, data = datSim, edgeEco = edgeEc,
  res = resAware, elapsed = elapsed
), rdsFile)
cat("Saved RDS to", rdsFile, "\n")
cat("Log file:", logFile, "\n")
cat("Done.\n")
