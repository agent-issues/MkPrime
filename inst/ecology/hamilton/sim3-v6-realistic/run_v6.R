# run_v6.R — Sim 3 v6-realistic replicate (blind + aware) on Hamilton.
#
# v6-realistic uses a *realistic* morphological-data tree length
# (truth TL ~ 1.16) at 16 tips x 300 characters, with the multirep-v3
# 8-tip eco-1 topology but distinguished eco vs non-eco clade-stem
# lengths (stemBrEco = 0.15, stemBrClade = 0.035) so the eco-shared
# parallel-substitution signal can compete with ancestry at non-saturated
# branch lengths.
#
# Uses the auto-expSteps default (commit bc24d52): expSteps is set
# from 1.05 * parsimony(start tree) at .FinalizeModel time. We DO NOT
# pass expSteps to MkPrimeModel(); we expect a cli_alert_info line of
# the form
#   "Tree length prior: parsimony score = N; expSteps = 1.05N; ..."
# in the .out log when each chain starts.
#
# Root-invariant scoring via inst/ecology/simulations/sim3-scoring.R
# (HasBipartSplits, ScoreTreesUnrooted).
#
# Usage on Hamilton:
#   Rscript run_v6.R <rep_id> <n_iter> <n_chains> <seed_base>
#
# Outputs to /nobackup/$USER/mkp-sim3-v6-realistic/results/rep%02d/
#   - blind-chain.log, blind-chain.ckp, blind-result.rds
#   - aware-chain.log, aware-chain.ckp, aware-result.rds
#   - summary.rds   (named list of corrected/legacy scores)

suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKP_LIB",
    file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-v6-realistic/lib"))
  .libPaths(c(libPath, .libPaths()))
  library("MkPrime")
  library("TreeTools")
  library("TreeDist")
})

args <- commandArgs(trailingOnly = TRUE)
repId    <- if (length(args) >= 1) as.integer(args[1]) else stop("rep_id required")
nIter    <- if (length(args) >= 2) as.integer(args[2]) else 100000L
nChains  <- if (length(args) >= 3) as.integer(args[3]) else 4L
seedBase <- if (length(args) >= 4) as.integer(args[4]) else 20260701L

outRoot <- file.path("/nobackup", Sys.getenv("USER"),
                     "mkp-sim3-v6-realistic",
                     "results", sprintf("rep%02d", repId))
dir.create(outRoot, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== v6-realistic rep %d / seed %d / nIter %d / nChains %d ===\n",
            repId, seedBase + repId, nIter, nChains))
cat("Output:", outRoot, "\n")

# Helper + simulator scripts -- reuse repo checkout
mkpRoot <- Sys.getenv("MKP_REPO_ROOT",
  file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-v6-realistic/MkPrime"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-helpers.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-simulate.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-scoring.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-v6-params.R"))

cfg <- SIM3V6_PARAMS$v6
nNeo  <- cfg$nNeo
nTrans <- cfg$nTrans
nChar <- nNeo + nTrans
phi   <- cfg$phi

tree <- .BuildConvergentTreeV6(tipBr = cfg$tipBr,
                                stemBrClade = cfg$stemBrClade,
                                stemBrEco = cfg$stemBrEco,
                                rootBr = cfg$rootBr)
eco    <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)
truthTL <- sum(tree$edge.length)
cat(sprintf("Truth TL = %.4f (expected ~1.16 for v6)\n", truthTL))
cat("Edge ecology distribution:\n"); print(table(edgeEc))

# --- Simulate -------------------------------------------------------------
set.seed(seedBase + repId)
zFull <- matrix(0L, nrow = nChar, ncol = 2L)
# Eco column for state 1 gets z = 1 in (1 - pi0) fraction of characters.
ecoActive <- stats::runif(nChar) > cfg$pi0
zFull[ecoActive, 2L] <- 1L
cFull <- c(rep("neomorphic", nNeo), rep("transformational", nTrans))
datSim <- .SimulateMkPrimeEcology(tree, edgeEc, zFull, phi = phi,
                                  type = cFull, baseRate = cfg$baseRate,
                                  normalize = TRUE,
                                  pi0 = cfg$pi0, theta = cfg$theta,
                                  refEcology = 0L, rateLoss = cfg$rateLoss)
cat(sprintf("Eco-active chars: %d / %d  (pi0 = %.2f)\n",
            sum(ecoActive), nChar, cfg$pi0))

pdSim    <- MatrixToPhyDat(datSim)
mkdBlind <- MkPrimeData(pdSim, neomorphic = seq_len(nNeo))
ecoTipVec <- eco[mkdBlind$taxon_names]
mkdAware <- MkPrimeData(pdSim, neomorphic = seq_len(nNeo),
                        ecology = ecoTipVec)

# --- Parsimony start tree -------------------------------------------------
set.seed(1)
randTree <- Preorder(ape::rtree(mkdBlind$nTip,
                                tip.label = mkdBlind$taxon_names))
startTree <- Preorder(TreeSearch::AdditionTree(pdSim))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))
cli::cli_alert_info(
  "Parsimony start tree: TL slot will be 0.1/edge; pre-MCMC parsimony \\
   used as the expSteps anchor by .FinalizeModel."
)

# --- Models: AUTO-expSteps (do NOT pass expSteps) -------------------------
# This is the key change for v6: rely on the bc24d52 auto-default so the
# Gamma TL prior is anchored at 1.05 * parsimony, NOT at the old hard-
# coded expSteps = 10.
modelBlind <- MkPrimeModel(ecologyAware = FALSE,
                           kPrimePrior = "geometric",
                           coding = "variable")
modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable")
cli::cli_alert_info("Models constructed with auto-expSteps (will resolve in .FinalizeModel)")

mcmcCfg <- function(logFile, ckpFile) {
  MkPrimeMCMC(nIter = nIter, nChains = nChains, nRuns = 1L,
              thin = max(1L, nIter %/% 500L),
              treeThin = max(1L, nIter %/% 500L),
              minWarmup = 5000L, maxWarmup = 40000L,
              logFile = logFile, checkpointFile = ckpFile)
}

run_one <- function(tag, mkd, model) {
  logFile <- file.path(outRoot, paste0(tag, "-chain.log"))
  ckpFile <- file.path(outRoot, paste0(tag, "-chain.ckp"))
  mcmc <- mcmcCfg(logFile, ckpFile)

  cat(sprintf("\n--- %s chain (nChains = %d) ---\n", tag, nChains))
  t0 <- Sys.time()
  if (file.exists(ckpFile)) {
    cat("Resuming from existing checkpoint:", ckpFile, "\n")
    res <- ResumeMkPrime(checkpointFile = ckpFile, data = mkd, model = model)
  } else {
    res <- RunMkPrime(mkd, tree = startTree, model = model, mcmc = mcmc)
  }
  cat(sprintf("%s elapsed: %s\n", tag, format(Sys.time() - t0)))
  saveRDS(res, file.path(outRoot, paste0(tag, "-result.rds")))
  res
}

resBlind <- run_one("blind", mkdBlind, modelBlind)
resAware <- run_one("aware", mkdAware, modelAware)
resAware <- RelabelEcology(resAware)
saveRDS(resAware, file.path(outRoot, "aware-result.rds"))

# --- Score using root-invariant ScoreTreesUnrooted ------------------------
biparts <- list(
  AC = c(paste0("A", 1:4), paste0("C", 1:4)),
  AB = c(paste0("A", 1:4), paste0("B", 1:4)),
  cladeA = paste0("A", 1:4),
  cladeB = paste0("B", 1:4),
  cladeC = paste0("C", 1:4),
  cladeD = paste0("D", 1:4)
)

scoreBlind <- ScoreTreesUnrooted(resBlind$trees, biparts, refTree = tree)
scoreAware <- ScoreTreesUnrooted(resAware$trees, biparts, refTree = tree)

summary <- list(
  repId = repId, seed = seedBase + repId,
  nIter = nIter, nChains = nChains,
  config = cfg,
  truthTL = truthTL,
  blind = scoreBlind,
  aware = scoreAware
)
saveRDS(summary, file.path(outRoot, "summary.rds"))

cat("\n=== v6-realistic per-rep scores (corrected / legacy) ===\n")
cat("-- blind --\n")
print(scoreBlind)
cat("-- aware --\n")
print(scoreAware)
cat("\nDone. Saved to", outRoot, "\n")
