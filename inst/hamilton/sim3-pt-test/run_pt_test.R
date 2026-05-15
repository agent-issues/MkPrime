# run_pt_test.R — single-rep parallel-tempering test for Sim 3 v3.
#
# Reuses the multirep-v3 simulator + config, but runs the aware chain
# with nChains = 4 (geometric heat ladder) instead of nChains = 1.
# Diagnostic question: does PT close the 100+ nat valley between the
# chain's MAP basin and the true AC bipartition (rep 05-type
# pathology in the 1-chain multirep)?
#
# Seed: 20260606 = seedBase 20260601 + rep 05 (the "high P(wrong)"
# case from the 1-chain multirep where truth was +106 nats above
# chain MAP under MAP nuisance).
#
# Targets to inspect after the run:
#   1. Swap acceptance rates across adjacent betas (aim 20-60%)
#   2. Did any chain visit the true AC bipartition?
#   3. logPost trace per chain — do hot chains explore higher-logPost
#      regions that the cold chain then swapped into?
#   4. Cold-chain P(true AC), P(wrong AB), CID — direct comparison to
#      the 1-chain rep 05 numbers (P(true)=0.000, P(wrong)=0.462).
#
# Usage on Hamilton:
#   Rscript run_pt_test.R

suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKP_LIB",
    file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-pt-test/lib"))
  .libPaths(c(libPath, .libPaths()))
  library("MkPrime")
  library("TreeTools")
  library("TreeDist")
})

repId    <- 5L
nIter    <- 100000L
nChains  <- 4L
heat     <- 0.2
seedBase <- 20260601L

outRoot <- file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-pt-test",
                     "results")
dir.create(outRoot, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== PT test: rep%d seed=%d nIter=%d nChains=%d heat=%.2f ===\n",
            repId, seedBase + repId, nIter, nChains, heat))
cat("Output:", outRoot, "\n")

mkpRoot <- Sys.getenv("MKP_REPO_ROOT",
  file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-pt-test/MkPrime"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3-helpers.R"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3-simulate.R"))

# V3 config (matches sim3-multirep-v3)
nEco <- 120L; nBase <- 360L; phi <- 4
stemBr <- 0.30; rootBr <- 0.15; tipBr <- 0.5
tree <- .BuildConvergentTree(tipBranch = tipBr,
                              stemBranch = stemBr, rootBranch = rootBr)
eco  <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)

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
truthTL <- sum(tree$edge.length)

set.seed(1)
randTree <- Preorder(ape::rtree(mkdAware$nTip,
                                tip.label = mkdAware$taxon_names))
psTrees <- TreeSearch::MaximizeParsimony(pdSim, tree = randTree,
                                         verbosity = 0)
startTree <- Preorder(psTrees[[1]])
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)

logFile <- file.path(outRoot, "pt-aware.log")
ckpFile <- file.path(outRoot, "pt-aware.ckp")
mcmc <- MkPrimeMCMC(nIter = nIter, nChains = nChains, heat = heat,
                    nRuns = 1L,
                    thin = max(1L, nIter %/% 500L),
                    treeThin = max(1L, nIter %/% 500L),
                    minWarmup = 5000L, maxWarmup = 40000L,
                    logFile = logFile, checkpointFile = ckpFile)

cat(sprintf("\n--- AWARE PT chain (nChains=%d) ---\n", nChains))
t0 <- Sys.time()
if (file.exists(ckpFile)) {
  cat("Resuming from existing checkpoint:", ckpFile, "\n")
  res <- ResumeMkPrime(mkdAware, mcmc = mcmc, model = modelAware)
} else {
  res <- RunMkPrime(mkdAware, tree = startTree,
                    model = modelAware, mcmc = mcmc)
}
cat(sprintf("Elapsed: %s\n", format(Sys.time() - t0)))
res <- RelabelEcology(res)
saveRDS(res, file.path(outRoot, "pt-aware-result.rds"))

# Score cold-chain trees
discard <- function(x) x[seq.int(ceiling(length(x) / 4) + 1L, length(x))]
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

tr <- discard(res$trees); class(tr) <- "multiPhylo"
pTrue  <- mean(hasBipart(tr, trueSplit))
pWrong <- mean(hasBipart(tr, wrongSplit))
cidTrue <- mean(as.numeric(TreeDist::ClusteringInfoDistance(
  tr, tree, normalize = TRUE)))

cat("\n=== PT test scores (cold chain) ===\n")
cat(sprintf("P(true AC)=%.3f  P(wrong AB)=%.3f  CID-to-truth=%.3f\n",
            pTrue, pWrong, cidTrue))
cat(sprintf("Compare 1-chain rep 05: P(true)=0.000 P(wrong)=0.462 CID=0.473\n"))

# Swap rates if available
if (!is.null(res$swap_accept)) {
  cat("\nSwap accept rates between adjacent betas:\n")
  rates <- res$swap_accept / pmax(1, res$swap_propose)
  for (i in seq_along(rates)) {
    cat(sprintf("  ladder %d<->%d: %d/%d = %.3f\n",
                i, i + 1L,
                res$swap_accept[i], res$swap_propose[i], rates[i]))
  }
}

cat("\nDone. Result + log + ckp in", outRoot, "\n")
