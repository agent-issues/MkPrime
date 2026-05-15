# run_phi2.R — softened landscape (phi=2 instead of phi=4) single-chain test.
#
# Posture 2 probe: does softening the convergent signal let standard MCMC
# comfortably recover AC under the aware model, while blind still fails?
#
# Same v3 config except phi=2. Single chain (nChains=1) each for blind
# and aware so we get the direct aware-vs-blind comparison on identical
# data.
#
# Seed: 20260602 (= seedBase 20260601 + rep 01). Matches multirep-v3 rep 01
# slot for cross-phi comparison if we extend to a sweep later.

suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKP_LIB",
    file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/lib"))
  .libPaths(c(libPath, .libPaths()))
  library("MkPrime")
  library("TreeTools")
  library("TreeDist")
})

repId    <- 1L
nIter    <- 100000L
nChains  <- 1L
seedBase <- 20260601L
phiSim   <- 2  # SOFTENED (multirep-v3 uses 4)

outRoot <- file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-phi2",
                     "results")
dir.create(outRoot, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== phi=%g test: rep%d seed=%d nIter=%d ===\n",
            phiSim, repId, seedBase + repId, nIter))
cat("Output:", outRoot, "\n")

mkpRoot <- Sys.getenv("MKP_REPO_ROOT",
  file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/MkPrime"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3-helpers.R"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3-simulate.R"))

# V3 config except phi
nEco <- 120L; nBase <- 360L
stemBr <- 0.30; rootBr <- 0.15; tipBr <- 0.5
tree <- .BuildConvergentTree(tipBranch = tipBr,
                              stemBranch = stemBr, rootBranch = rootBr)
eco  <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)

set.seed(seedBase + repId)
zFull <- matrix(0L, nrow = nEco + nBase, ncol = 2)
zFull[seq_len(nEco), 2] <- 1L
cFull <- c(rep("neomorphic", nEco), rep("transformational", nBase))
datSim <- .SimulateMkPrimeEcology(tree, edgeEc, zFull, phi = phiSim,
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

runOne <- function(mkd, model, label) {
  cat(sprintf("\n--- %s chain (phi_sim=%g, nChains=1) ---\n", label, phiSim))
  logFile <- file.path(outRoot, sprintf("%s.log", label))
  ckpFile <- file.path(outRoot, sprintf("%s.ckp", label))
  mcmc <- MkPrimeMCMC(nIter = nIter, nChains = nChains, nRuns = 1L,
                      thin = max(1L, nIter %/% 500L),
                      treeThin = max(1L, nIter %/% 500L),
                      minWarmup = 5000L, maxWarmup = 40000L,
                      logFile = logFile, checkpointFile = ckpFile)
  t0 <- Sys.time()
  if (file.exists(ckpFile)) {
    cat("Resuming from checkpoint\n")
    res <- ResumeMkPrime(mkd, mcmc = mcmc, model = model)
  } else {
    res <- RunMkPrime(mkd, tree = startTree, model = model, mcmc = mcmc)
  }
  cat(sprintf("Elapsed: %s\n", format(Sys.time() - t0)))
  if (isTRUE(model$ecologyAware)) res <- RelabelEcology(res)
  saveRDS(res, file.path(outRoot, sprintf("%s-result.rds", label)))

  tr <- discard(res$trees); class(tr) <- "multiPhylo"
  pTrue  <- mean(hasBipart(tr, trueSplit))
  pWrong <- mean(hasBipart(tr, wrongSplit))
  cidTrue <- mean(as.numeric(TreeDist::ClusteringInfoDistance(
    tr, tree, normalize = TRUE)))
  cat(sprintf("[%s] P(true AC)=%.3f  P(wrong AB)=%.3f  CID=%.3f\n",
              label, pTrue, pWrong, cidTrue))
  invisible(list(pTrue = pTrue, pWrong = pWrong, cid = cidTrue))
}

modelBlind <- MkPrimeModel(ecologyAware = FALSE,
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)
modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)

bScores <- runOne(mkdBlind, modelBlind, "blind")
aScores <- runOne(mkdAware, modelAware, "aware")

cat(sprintf("\n=== phi=%g summary (rep %d, seed %d) ===\n",
            phiSim, repId, seedBase + repId))
cat(sprintf("  BLIND: P(true)=%.3f  P(wrong)=%.3f  CID=%.3f\n",
            bScores$pTrue, bScores$pWrong, bScores$cid))
cat(sprintf("  AWARE: P(true)=%.3f  P(wrong)=%.3f  CID=%.3f\n",
            aScores$pTrue, aScores$pWrong, aScores$cid))
cat(sprintf("Compare phi=4 multirep-v3 mean: aware P(true)=0.002 P(wrong)=mixed\n"))
cat("\nPosture 2 success criterion: aware P(true) substantially > blind P(true).\n")

cat("\nDone. Results in", outRoot, "\n")
