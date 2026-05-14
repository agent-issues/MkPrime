# run.R -------------------------------------------------------------------------
# Discriminating test: aligned-units re-sim.
#
# Hypothesis: the ~2x tree_length upward bias in chain A (seed 20260513) is NOT
# a model bug but a unit mismatch between simulator and model.
#
# Previous chain A used .SimulateMkPrimeEcology(..., baseRate=0.5, normalize=FALSE).
# With baseRate=0.5 the simulator's rate-time product is half the model's
# per-edge substitution count, so the model must inflate tree_length to fit.
#
# This run fixes the mismatch by:
#   1. baseRate = 1.0  (model's JC kernel absorbs no separate baseRate)
#   2. normalize = TRUE (divide per-cell rate factor by gammaE so that at truth
#      the simulator's expected substitutions equal the model's expectation)
#   3. pi0 = 0.75, theta = 1.0 EXPLICIT — the simulator's default pi0 would
#      compute mean(z==0) over ALL columns including the reference column
#      (giving 420/480 = 0.875, not 0.75). The model counts only the K-1
#      non-reference columns (180/240 = 0.75). Passing explicit values makes
#      simulator gammaE identical to model gammaE at truth.
#
# Config:
#   seed=20260514 (fresh; same seed as truth-init chain B, but with aligned sim
#   it generates different data — documented clearly)
#   kEco=2 (A/B=ecology1, C/D=ecology0), nTheta=1
#   nEco=60 neomorphic, nBase=180 transformational
#   phi=4, stem=0.10, root=0.15 (Goldilocks — identical to previous chains)
#   nIter=100000, thin=40 -> ~2500 raw samples, ~1875 post-25%-discard
#   Initial state: DEFAULT (no phi=4 force — let chain find its own mode)
#   Single chain, simple Bactrians only
#
# If posterior tl_median collapses to ~truth tl: unit mismatch confirmed.
# If pi0_median ~0.75: pi0 bias was also unit-mismatch artefact.
# If pi0 still ~0.235: pi0 bias is real and separate.
#
# Run from the worktree root:
#   Rscript -e "setwd('C:/Users/pjjg18/GitHub/mkp/.claude/worktrees/ecology-aware'); source('dev/pilots/2026-05-14-aligned-units-resim/run.R')"

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/simulations/ecology/sim3-helpers.R")
source("inst/simulations/ecology/sim3-simulate.R")

PILOT_DIR <- "dev/pilots/2026-05-14-aligned-units-resim"

# ---- Simulation ---------------------------------------------------------------
# Seed for data generation: 20260514.
# NOTE: same integer as truth-init chain B, but that chain used a different
# (misaligned) dataset. This run generates a NEW dataset with aligned params.
set.seed(20260514)
nEco <- 60L; nBase <- 180L; phi <- 4; stemBr <- 0.10; rootBr <- 0.15
tree   <- .BuildConvergentTree(stemBranch = stemBr, rootBranch = rootBr)
eco    <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)

zFull <- matrix(0L, nrow = nEco + nBase, ncol = 2)
zFull[seq_len(nEco), 2] <- 1L  # col 2 = non-ref ecology; z=1 for neo chars
cFull <- c(rep("neomorphic", nEco), rep("transformational", nBase))

# --- ALIGNED SIMULATOR CALL ---
# baseRate = 1.0 : matches model's JC kernel (no separate baseRate scaling)
# normalize = TRUE: divides per-cell factor by gammaE so expected subs = model
# pi0 = 0.75, theta = 1.0: MUST be explicit (see header; default gives 0.875)
# refEcology = 0L: default, matches model; reference-ecology branches factor=1
datSim <- .SimulateMkPrimeEcology(tree, edgeEc, zFull, phi = phi,
                                  type      = cFull,
                                  baseRate  = 1.0,
                                  normalize = TRUE,
                                  pi0       = 0.75,
                                  theta     = 1.0,
                                  refEcology = 0L,
                                  rateLoss  = 1)

# Report simulator-computed truth tree_length (model units)
truthTL <- sum(tree$edge.length)
cat("Truth tree_length (sum edge.length, model units):", truthTL, "\n")
cat("  [stem=", stemBr, " x4 + root=", rootBr, " x2 + within-clade edges]\n", sep="")

pdSim <- MatrixToPhyDat(datSim)

mkdBlind <- MkPrimeData(pdSim, neomorphic = seq_len(nEco))
ecoTipVec <- eco[mkdBlind$taxon_names]
mkd <- MkPrimeData(pdSim, neomorphic = seq_len(nEco), ecology = ecoTipVec)
cat("kEcology:", mkd$kEcology, "  nTheta:", mkd$kEcology - 1L, "\n")

# Parsimony start tree
set.seed(1)
randTree <- Preorder(ape::rtree(mkd$nTip, tip.label = mkd$taxon_names))
psTrees  <- TreeSearch::MaximizeParsimony(MatrixToPhyDat(datSim), tree = randTree,
                                          verbosity = 0)
startTree <- Preorder(psTrees[[1]])
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

# ---- Model --------------------------------------------------------------------
modelAware <- MkPrimeModel(
  ecologyAware  = TRUE,
  magnitudeMode = "global",
  kPrimePrior   = "geometric",
  coding        = "variable",
  expSteps      = 10
)

# ---- MCMC config --------------------------------------------------------------
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

# ---- Run ----------------------------------------------------------------------
cat("== Running aligned-units re-sim (nIter=", nIter, ", thin=", thin, ") ==\n")
cat("   Seed (data): 20260514 — NEW data with aligned sim params\n")
cat("   Seed (MCMC): default init (no phi=4 override)\n")
cat("   baseRate=1.0, normalize=TRUE, pi0=0.75, theta=1.0\n")
cat("   Expected runtime: ~30 min\n")
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

# ---- Save result --------------------------------------------------------------
rdsFile <- file.path(PILOT_DIR, "result.rds")
saveRDS(list(
  config = list(
    seed_data  = 20260514,
    seed_mcmc  = "default",
    nEco       = nEco,
    nBase      = nBase,
    phi        = phi,
    stemBr     = stemBr,
    rootBr     = rootBr,
    nIter      = nIter,
    thin       = thin,
    kEcology   = mkd$kEcology,
    baseRate   = 1.0,
    normalize  = TRUE,
    pi0_sim    = 0.75,
    theta_sim  = 1.0,
    refEcology = 0L,
    note       = "Aligned-units re-sim: unit mismatch test. baseRate=1 + normalize=TRUE + explicit pi0/theta"
  ),
  tree       = tree,
  eco        = eco,
  z          = zFull,
  charType   = cFull,
  phi_truth  = phi,
  truth_tl   = truthTL,
  data       = datSim,
  edgeEco    = edgeEc,
  res        = resAware,
  elapsed    = elapsed,
  crashed    = is.null(resAware)
), rdsFile)
cat("Saved RDS to", rdsFile, "\n")
cat("Log file:", logFile, "\n")

if (is.null(resAware)) stop("Chain crashed — check log above.")

cat("Done. Now run analyse.R to generate REPORT.md\n")
