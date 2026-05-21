# run.R -------------------------------------------------------------------------
# Step 3b pilot: 100k-iteration mode-persistence test.
#
# Question: Does a single chain stay in one phi-mode for 100k iterations?
# This determines whether Option 1 (post-hoc relabel) is sufficient or whether
# the (phi <-> 1/phi, theta <-> 1-theta) symmetry must be broken via the prior.
#
# Config:
#   seed=20260513 (one more than the 20k pilot; independent chain)
#   kEco=2 (A/B=ecology1, C/D=ecology0), nTheta=1
#   nEco=60 neomorphic, nBase=180 transformational
#   phi=4, stem=0.10, root=0.15 (Goldilocks — identical to 20k pilot)
#   nIter=100000, thin=40 -> ~2500 raw samples, ~1875 post-25%-discard
#   Initial state: default (whatever .InitState yields)
#   Single chain, simple Bactrians only (no new joint scalers)
#
# Run from the worktree root:
#   Rscript -e "setwd('C:/Users/pjjg18/GitHub/mkp/.claude/worktrees/ecology-aware'); source('dev/pilots/2026-05-13-step3b-mode-persistence-A/run.R')"

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/ecology/simulations/sim3-helpers.R")
source("inst/ecology/simulations/sim3-simulate.R")

PILOT_DIR <- "dev/pilots/2026-05-13-step3b-mode-persistence-A"

# Seed: 20260513 (one more than the 20k pilot — independent chain)
set.seed(20260513)
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

# 100k iterations, thin=40 -> ~2500 rows; 25%-discard leaves ~1875
nIter <- 100000L
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
logFile  <- file.path(PILOT_DIR, "chain.log")
ckpFile  <- file.path(PILOT_DIR, "chain.ckp")

# Clean any previous partial run
for (f in c(logFile, ckpFile)) if (file.exists(f)) file.remove(f)

mcmc$logFile        <- logFile
mcmc$checkpointFile <- ckpFile

cat("== Running Step 3b mode-persistence chain (nIter=", nIter, ", thin=", thin, ") ==\n")
cat("   Seed: 20260513\n")
cat("   Expected runtime: ~22 min (5x the 20k pilot's 4.4 min)\n")
t0 <- Sys.time()
resAware <- tryCatch(
  RunMkPrime(mkd, tree = startTree, model = modelAware, mcmc = mcmc),
  error = function(e) {
    cat("ERROR during RunMkPrime:", conditionMessage(e), "\n")
    NULL
  }
)
elapsed <- Sys.time() - t0
cat("Elapsed:", format(elapsed), "\n")

if (is.null(resAware)) {
  cat("Chain crashed. Saving NULL result and partial log if it exists.\n")
  saveRDS(list(
    config = list(seed = 20260513, nEco = nEco, nBase = nBase, phi = phi,
                  stemBr = stemBr, rootBr = rootBr, nIter = nIter, thin = thin,
                  kEcology = mkd$kEcology),
    tree = tree, eco = eco, z = zFull, charType = cFull,
    phi_truth = phi, data = datSim, edgeEco = edgeEc,
    res = NULL, elapsed = elapsed, crashed = TRUE
  ), file.path(PILOT_DIR, "result.rds"))
  stop("Chain crashed — check log above.")
}

# Save RDS (result object needed for downstream z_samples and relabeling code)
rdsFile <- file.path(PILOT_DIR, "result.rds")
saveRDS(list(
  config = list(seed = 20260513, nEco = nEco, nBase = nBase, phi = phi,
                stemBr = stemBr, rootBr = rootBr, nIter = nIter, thin = thin,
                kEcology = mkd$kEcology),
  tree = tree, eco = eco, z = zFull, charType = cFull,
  phi_truth = phi, data = datSim, edgeEco = edgeEc,
  res = resAware, elapsed = elapsed
), rdsFile)
cat("Saved RDS to", rdsFile, "\n")
cat("Log file:", logFile, "\n")
cat("Done. Now run analyse.R to generate REPORT.md\n")
