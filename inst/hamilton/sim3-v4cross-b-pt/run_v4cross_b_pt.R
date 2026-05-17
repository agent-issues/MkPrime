# run_v4cross_b_pt.R — Sim 3 v4-cross balanced regime, PT rerun (nChains=4).
#
# Clone of run_v4cross_b.R with nChains=1 → 4 to diagnose the mode-trap
# observed in the single-chain aware run (clade A OK, clade B P=0.121;
# samples 1-83 mode 1, 84+ stuck in mode 2 at lower likelihood).
#
# PT should break the trap if the problem is a mixing barrier.  All other
# parameters are identical to v4cross-b: same seed (20260601+1), same
# nIter=100000, same v4b balanced regime.

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
nChains  <- 4L
seedBase <- 20260601L
PARAM_SET <- "v4b"  # reuse v4b balanced regime

outRoot <- file.path("/nobackup", Sys.getenv("USER"),
                     "mkp-sim3-v4cross-b-pt", "results")
dir.create(outRoot, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== Sim3 v4cross-b-pt (nChains=%d): rep%d seed=%d nIter=%d ===\n",
            nChains, repId, seedBase + repId, nIter))
cat("Output:", outRoot, "\n")

mkpRoot <- Sys.getenv("MKP_REPO_ROOT",
  file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/MkPrime"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3-simulate.R"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3v4-helpers.R"))      # .BuildV4Data
source(file.path(mkpRoot, "inst/simulations/ecology/sim3v4cross-helpers.R")) # new
source(file.path(mkpRoot, "inst/simulations/ecology/sim3v4-params.R"))

params <- SIM3V4_PARAMS[[PARAM_SET]]
cat(sprintf("Params: tipBr=%.3f stemBrEco=%.3f stemBrClade=%.3f rootBr=%.3f\n",
            params$tipBr, params$stemBrEco, params$stemBrClade, params$rootBr))
cat(sprintf("        nNeo=%d nTrans=%d phi=%g pi0=%.2f theta=%.2f\n",
            params$nNeo, params$nTrans, params$phi, params$pi0, params$theta))

# v4-cross tree + ecology
tree   <- .BuildConvergentTreeV4Cross(params$tipBr, params$stemBrEco,
                                       params$stemBrClade, params$rootBr)
eco    <- .ConvergentEcologyV4Cross(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)
biparts <- .SimBipartitionsV4Cross()

set.seed(seedBase + repId)
v4dat <- .BuildV4Data(tree, eco, edgeEc, params,
                      simulator = .SimulateMkPrimeEcology)
mkdBlind <- v4dat$mkdBlind
mkdAware <- v4dat$mkdAware
cat(sprintf("Tree TL=%.3f  nChar=%d (%d neo + %d trans)\n",
            v4dat$truthTL, mkdBlind$nChar,
            sum(mkdBlind$type == "neomorphic"),
            sum(mkdBlind$type == "transformational")))

set.seed(1)
randTree <- Preorder(ape::rtree(mkdBlind$nTip,
                                tip.label = mkdBlind$taxon_names))
psTrees  <- TreeSearch::MaximizeParsimony(v4dat$pdSim, tree = randTree,
                                          verbosity = 0)
startTree <- Preorder(psTrees[[1]])
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

hasBipart <- function(treeList, splitTips) {
  vapply(treeList, function(tr) {
    cl   <- ape::prop.part(tr)
    tips <- attr(cl, "labels")
    splitSet <- which(tips %in% splitTips)
    any(vapply(cl, function(p) setequal(p, splitSet), logical(1)))
  }, logical(1))
}
discard <- function(x) x[seq.int(ceiling(length(x) / 4) + 1L, length(x))]

scoreTrees <- function(treeList, label) {
  tr <- discard(treeList); class(tr) <- "multiPhylo"
  pA   <- mean(hasBipart(tr, biparts$trueClade_A))
  pB   <- mean(hasBipart(tr, biparts$trueClade_B))
  pC   <- mean(hasBipart(tr, biparts$trueClade_C))
  pD   <- mean(hasBipart(tr, biparts$trueClade_D))
  pAC  <- mean(hasBipart(tr, biparts$trueSister_AC))
  pBD  <- mean(hasBipart(tr, biparts$trueSister_BD))
  pFi  <- mean(hasBipart(tr, biparts$falseInner))
  pFc  <- mean(hasBipart(tr, biparts$falseClade_AB))
  cid  <- mean(as.numeric(
    TreeDist::ClusteringInfoDistance(tr, tree, normalize = TRUE)))
  cat(sprintf(
    "[%s] clA=%.3f clB=%.3f clC=%.3f clD=%.3f | AC=%.3f BD=%.3f | falseInner=%.3f falseAB=%.3f CID=%.3f\n",
    label, pA, pB, pC, pD, pAC, pBD, pFi, pFc, cid))
  invisible(list(pA=pA, pB=pB, pC=pC, pD=pD,
                 pAC=pAC, pBD=pBD,
                 pFi=pFi, pFc=pFc, cid=cid))
}

runOne <- function(mkd, model, label) {
  cat(sprintf("\n--- %s chain (nChains=%d) ---\n", label, nChains))
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
  scoreTrees(res$trees, label)
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

cat(sprintf("\n=== v4cross-b-pt summary (rep %d, seed %d, nChains=%d) ===\n",
            repId, seedBase + repId, nChains))
cat(sprintf("  BLIND: clA=%.3f clB=%.3f AC=%.3f BD=%.3f falseInner=%.3f falseAB=%.3f CID=%.3f\n",
            bScores$pA, bScores$pB, bScores$pAC, bScores$pBD,
            bScores$pFi, bScores$pFc, bScores$cid))
cat(sprintf("  AWARE: clA=%.3f clB=%.3f AC=%.3f BD=%.3f falseInner=%.3f falseAB=%.3f CID=%.3f\n",
            aScores$pA, aScores$pB, aScores$pAC, aScores$pBD,
            aScores$pFi, aScores$pFc, aScores$cid))
cat("\nTarget: PT fixes clade B recovery (was P=0.121 under single chain).\n")
cat("Compare v4cross-b single-chain: clB=0.121, mode-trap at iter 84.\n")
cat("\nDone. Results in", outRoot, "\n")
