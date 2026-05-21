# run_pt.R — Sim 3 multirep-v3 PT run, rep 05 only (blind + aware)
#
# Tests whether Parallel Tempering can overcome the mixing problem that
# leaves single-chain aware diffuse over topologies in the multirep-v3
# setup (single-chain: P(AC) mean 0.012 across 8 reps, 94-167 unique
# topologies). 4 chains, sequential (MkPrime PT is single-threaded).
#
# Usage on Hamilton:
#   Rscript run_pt.R   (no args; rep 05 is hardcoded)
#
# Outputs written to /nobackup/$USER/mkp-sim3-multirep-v3-pt/results/rep05/
#   - blind-chain.log, blind-chain.ckp, blind-result.rds
#   - aware-chain.log, aware-chain.ckp, aware-result.rds
#   - summary.rds (data frame from ScoreTreesUnrooted)

suppressPackageStartupMessages({
  # Library path: reuse the existing v3 library (same built package).
  libPath <- Sys.getenv("MKP_LIB",
    file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/lib"))
  .libPaths(c(libPath, .libPaths()))
  library("MkPrime")
  library("TreeTools")
  library("TreeDist")
})

repId    <- 5L
nIter    <- 100000L
nChains  <- 4L
seedBase <- 20260601L

outRoot <- file.path("/nobackup", Sys.getenv("USER"),
                     "mkp-sim3-multirep-v3-pt",
                     "results", sprintf("rep%02d", repId))
dir.create(outRoot, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== PT Rep %d / seed %d / nIter %d / nChains %d ===\n",
            repId, seedBase + repId, nIter, nChains))
cat("Output:", outRoot, "\n")

# Helper and simulator scripts — reuse existing repo checkout
mkpRoot <- Sys.getenv("MKP_REPO_ROOT",
  file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/MkPrime"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-helpers.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-simulate.R"))

# Root-invariant scoring (required — legacy prop.part is rooting-sensitive)
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-scoring.R"))

# V3 config (identical to run_rep.R)
nEco <- 120L; nBase <- 360L; phi <- 4
stemBr <- 0.30; rootBr <- 0.15; tipBr <- 0.5
tree <- .BuildConvergentTree(tipBranch = tipBr,
                              stemBranch = stemBr, rootBranch = rootBr)
eco  <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)

# Simulate rep-05 data (same seed as multirep-v3 rep 05)
set.seed(seedBase + repId)
zFull <- matrix(0L, nrow = nEco + nBase, ncol = 2)
zFull[seq_len(nEco), 2] <- 1L
cFull <- c(rep("neomorphic", nEco), rep("transformational", nBase))
datSim <- .SimulateMkPrimeEcology(tree, edgeEc, zFull, phi = phi,
                                  type = cFull, baseRate = 1.0,
                                  normalize = TRUE,
                                  pi0 = 0.75, theta = 1.0,
                                  refEcology = 0L, rateLoss = 1)

pdSim <- MatrixToPhyDat(datSim)
mkdBlind <- MkPrimeData(pdSim, neomorphic = seq_len(nEco))
ecoTipVec <- eco[mkdBlind$taxon_names]
mkdAware <- MkPrimeData(pdSim, neomorphic = seq_len(nEco),
                        ecology = ecoTipVec)

# Parsimony start tree (same strategy as run_rep.R)
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

mcmcPT <- function(logFile, ckpFile) {
  MkPrimeMCMC(nIter = nIter, nChains = nChains, nRuns = 1L,
              thin = max(1L, nIter %/% 500L),
              treeThin = max(1L, nIter %/% 500L),
              minWarmup = 5000L, maxWarmup = 40000L,
              logFile = logFile, checkpointFile = ckpFile)
}

run_one <- function(tag, mkd, model) {
  logFile <- file.path(outRoot, paste0(tag, "-chain.log"))
  ckpFile <- file.path(outRoot, paste0(tag, "-chain.ckp"))
  mcmc <- mcmcPT(logFile, ckpFile)

  cat(sprintf("\n--- %s PT chain (%d chains) ---\n", tag, nChains))
  t0 <- Sys.time()
  if (file.exists(ckpFile)) {
    cat("Resuming from existing checkpoint:", ckpFile, "\n")
    res <- ResumeMkPrime(mkd, mcmc = mcmc, model = model)
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

# Score using root-invariant ScoreTreesUnrooted
biparts <- list(
  AC = c(paste0("A", 1:4), paste0("C", 1:4)),
  AB = c(paste0("A", 1:4), paste0("B", 1:4))
)

scoreBlind <- ScoreTreesUnrooted(resBlind$trees, biparts, refTree = tree)
scoreAware <- ScoreTreesUnrooted(resAware$trees, biparts, refTree = tree)

summary <- list(
  repId = repId, seed = seedBase + repId, nIter = nIter, nChains = nChains,
  config = list(nEco = nEco, nBase = nBase, phi = phi,
                stemBr = stemBr, rootBr = rootBr, tipBr = tipBr),
  blind = scoreBlind,
  aware = scoreAware
)
saveRDS(summary, file.path(outRoot, "summary.rds"))

cat("\n=== PT per-rep scores (corrected / legacy) ===\n")
cat("-- blind --\n")
print(scoreBlind)
cat("-- aware --\n")
print(scoreAware)
cat("\nDone. Saved to", outRoot, "\n")
