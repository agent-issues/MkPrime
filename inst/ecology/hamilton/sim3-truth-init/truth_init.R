# truth_init.R — diagnostic: does the aware chain hold the truth peak?
#
# Reuses sim3 v3 config (rep 1, seed 20260602) but inits the chain at the
# true topology + true (phi, pi0, theta, z, TL, rate_neo). If posterior
# stays near these values for 10k iter, the priors/landscape are fine and
# topology mixing is the only remaining bottleneck. If it drifts within
# 1k iter, params can't hold the peak — tighten priors first.
#
# Usage on Hamilton:
#   Rscript truth_init.R [n_iter]
suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKP_LIB",
    file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/lib"))
  .libPaths(c(libPath, .libPaths()))
  library("MkPrime")
  library("TreeTools")
  library("TreeDist")
})

args  <- commandArgs(trailingOnly = TRUE)
nIter <- if (length(args) >= 1) as.integer(args[1]) else 10000L

outRoot <- file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3",
                     "truth-init")
dir.create(outRoot, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== Truth-init aware, nIter %d ===\n", nIter))
cat("Output:", outRoot, "\n")

mkpRoot <- Sys.getenv("MKP_REPO_ROOT",
  file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/MkPrime"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-helpers.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-simulate.R"))

# Reproduce rep 1 of multirep v3 exactly (seed 20260602)
nEco <- 120L; nBase <- 360L; phi <- 4
stemBr <- 0.30; rootBr <- 0.15; tipBr <- 0.5
tree <- .BuildConvergentTree(tipBranch = tipBr,
                              stemBranch = stemBr, rootBranch = rootBr)
eco  <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)

set.seed(20260602L)
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

# True params (Phase 1)
trueTree <- Preorder(tree)  # MCMC expects preorder
trueTL   <- sum(trueTree$edge.length)
trueRelBr <- trueTree$edge.length / trueTL
# zMatrix in model layout: nChar × (kEcology - 1) = nChar × 1.
# Truth: chars 1..120 z=1 (encouraged on non-ref ecology), rest z=0.
trueZ <- matrix(0L, nrow = nEco + nBase, ncol = 1L)
trueZ[seq_len(nEco), 1L] <- 1L

modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)

mcmc <- MkPrimeMCMC(nIter = nIter, nChains = 1L, nRuns = 1L,
                    thin = max(1L, nIter %/% 500L),
                    treeThin = max(1L, nIter %/% 500L),
                    # Short warmup: truth-init should already be at the peak,
                    # so we don't need long convergence checks. Want most of
                    # the budget in the sampling phase to see whether truth
                    # holds or drifts.
                    minWarmup = 200L, maxWarmup = 1000L,
                    logFile = file.path(outRoot, "truth-init.log"),
                    checkpointFile = file.path(outRoot, "truth-init.ckp"))

# Wipe old logs/ckpt so we always start fresh
for (f in c(mcmc$logFile, mcmc$checkpointFile)) {
  if (file.exists(f)) file.remove(f)
}

initOverrides <- list(
  tree = trueTree,
  tree_length = trueTL,
  rel_br_lengths = trueRelBr,
  phi = 4.0,
  pi0 = 0.75,
  # Truth is theta = 1, but the model's prior support is the open interval
  # (0, 1) — both R LogPrior and C++ cpp_log_prior reject theta == 1 exactly.
  # Nudge to 0.999 to land just inside the prior support; difference in
  # likelihood is < 0.5 nats with z = 0 in column.
  theta = 0.999,
  z = trueZ,
  rate_neo = 1.0,
  rate_loss = 1.0,
  rate_log_sd = 0.5
)

cat("\n--- Running truth-init aware chain ---\n")
t0 <- Sys.time()
res <- RunMkPrime(mkdAware, tree = trueTree, model = modelAware, mcmc = mcmc,
                  initOverrides = initOverrides)
cat(sprintf("Elapsed: %s\n", format(Sys.time() - t0)))
saveRDS(res, file.path(outRoot, "truth-init-result.rds"))

# Quick diagnostic: trajectory of phi, pi0, theta, TL, log_lik
samp <- res$samples
cat("\n== First & last 5 samples ==\n")
relCols <- intersect(c("log_likelihood", "log_posterior", "tree_length",
                       "phi", "pi0", "theta_1", "rate_neo", "rate_loss"),
                     colnames(samp))
print(round(rbind(head(samp[, relCols], 5), tail(samp[, relCols], 5)), 3))

cat("\n== Median of last 50% ==\n")
n <- nrow(samp)
keep <- samp[seq.int(ceiling(n / 2), n), relCols, drop = FALSE]
print(round(apply(keep, 2, median), 3))

cat("\nDone. Truth-init artefacts in", outRoot, "\n")
