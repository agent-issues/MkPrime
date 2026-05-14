# run.R -----------------------------------------------------------------------
# Step 3b pilot: truth-init mode-persistence check
#
# Same config as step3 20k pilot EXCEPT:
#   seed     = 20260514  (independent of 20k pilot and sibling 100k chain)
#   nIter    = 100,000
#   phi_init = 4.0   (truth)
#   theta_init = 0.97 (truth-side slab balance)
#
# The phi/theta override is injected via a temporary monkey-patch of
# .InitState (package namespace), which re-computes log_lik/log_prior/log_post
# after the override so the chain starts with a valid cached posterior.
#
# Question: does the truth-init chain stay in phi > 1 or cross to phi < 1?
#
# Run from the worktree root:
#   Rscript -e "setwd('C:/Users/pjjg18/GitHub/mkp/.claude/worktrees/ecology-aware'); source('dev/pilots/2026-05-13-step3b-mode-persistence-B-truthInit/run.R')"

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/simulations/ecology/sim3-helpers.R")
source("inst/simulations/ecology/sim3-simulate.R")

OUT_DIR <- "dev/pilots/2026-05-13-step3b-mode-persistence-B-truthInit"

# ---------------------------------------------------------------------------
# Monkey-patch .InitState to inject truth phi/theta and recompute posteriors
# ---------------------------------------------------------------------------
ns <- "MkPrime"
orig_InitState <- get(".InitState", envir = asNamespace(ns))

patched_InitState <- function(tree, mkd, model) {
  state <- orig_InitState(tree, mkd, model)
  if (isTRUE(model$ecologyAware)) {
    # Override phi and theta to the truth-side mode
    state$phi   <- if (length(state$phi) > 1L) rep(4.0, length(state$phi)) else 4.0
    state$theta <- rep(0.97, length(state$theta))
    # z stays at the zero matrix — Gibbs resamples it within the first few iters

    # Recompute log_lik with the new phi/theta (z=0 so slab contribution is
    # via pi0 only at this moment; phi enters via gamma_e on the edge rates)
    state$log_lik <- MkPrime:::.MkpEcologyLogLikelihood(
      tree, mkd,
      kPrime        = state$kPrime,
      rate_loss     = state$rate_loss,
      rate_log_sd   = state$rate_log_sd,
      nCat          = model$nCat,
      rate_neo      = if (is.null(state$rate_neo)) 1.0 else state$rate_neo,
      relabel       = model$relabel,
      phi           = state$phi,
      zMat          = state$z,
      magnitudeMode = model$magnitudeMode,
      coding        = model$coding,
      refEcology    = mkd$refEcology,
      theta         = state$theta,
      pi0           = state$pi0
    )
    state$log_prior <- MkPrime:::LogPrior(state, model, mkd)
    state$log_post  <- state$log_lik + state$log_prior

    cat(sprintf(
      "[truth-init] phi=%.4f  theta=%.4f  log_lik=%.2f  log_prior=%.2f\n",
      state$phi[[1]], state$theta[[1]], state$log_lik, state$log_prior
    ))
  }
  state
}
assignInNamespace(".InitState", patched_InitState, ns = ns)
cat("[patch] .InitState replaced in MkPrime namespace\n")

# ---------------------------------------------------------------------------
# Simulate data (identical to 20k pilot)
# ---------------------------------------------------------------------------
set.seed(20260512)   # same data seed as 20k pilot for direct comparability
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

# Parsimony start tree (same topology as 20k pilot)
set.seed(1)
randTree <- Preorder(ape::rtree(mkd$nTip, tip.label = mkd$taxon_names))
psTrees  <- TreeSearch::MaximizeParsimony(MatrixToPhyDat(datSim), tree = randTree,
                                          verbosity = 0)
startTree <- Preorder(psTrees[[1]])
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

# ---------------------------------------------------------------------------
# Model & MCMC config
# ---------------------------------------------------------------------------
modelAware <- MkPrimeModel(
  ecologyAware  = TRUE,
  magnitudeMode = "global",
  kPrimePrior   = "geometric",
  coding        = "variable",
  expSteps      = 10
)

set.seed(20260514)   # chain seed independent of both other pilots

nIter <- 100000L
thin  <- 40L
logFile <- file.path(OUT_DIR, "chain.log")
ckpFile <- file.path(OUT_DIR, "chain.ckp")

for (f in c(logFile, ckpFile)) if (file.exists(f)) file.remove(f)

mcmc <- MkPrimeMCMC(
  nIter          = nIter,
  nChains        = 1L,
  nRuns          = 1L,
  thin           = thin,
  treeThin       = thin,
  minWarmup      = 2000L,
  maxWarmup      = 5000L,
  logFile        = logFile,
  checkpointFile = ckpFile
)

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
cat("== Running Step 3b truth-init chain (nIter=", nIter, ", thin=", thin, ") ==\n")
t0 <- Sys.time()
res <- RunMkPrime(mkd, tree = startTree, model = modelAware, mcmc = mcmc)
elapsed <- Sys.time() - t0
cat("Elapsed:", format(elapsed), "\n")

# ---------------------------------------------------------------------------
# Restore original .InitState (good practice, though Rscript exits anyway)
# ---------------------------------------------------------------------------
assignInNamespace(".InitState", orig_InitState, ns = ns)
cat("[patch] .InitState restored\n")

# ---------------------------------------------------------------------------
# Save RDS
# ---------------------------------------------------------------------------
rdsFile <- file.path(OUT_DIR, "result.rds")
saveRDS(list(
  config = list(seed_data = 20260512, seed_chain = 20260514,
                nEco = nEco, nBase = nBase, phi_truth = phi,
                stemBr = stemBr, rootBr = rootBr,
                nIter = nIter, thin = thin,
                phi_init = 4.0, theta_init = 0.97,
                kEcology = mkd$kEcology),
  tree = tree, eco = eco, z = zFull, charType = cFull,
  data = datSim, edgeEco = edgeEc,
  res = res, elapsed = elapsed
), rdsFile)
cat("Saved RDS to", rdsFile, "\n")
cat("Log file:", logFile, "\n")
cat("Done.\n")
