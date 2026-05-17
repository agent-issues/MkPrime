# run_v4cross_c.R — Sim 3 v4-cross STRESS regime (single-rep blind+aware test).
#
# v4-cross differs from v4 only in ECOLOGY ASSIGNMENT:
#   v4       : eco-1 = {A1, A2, C1, C2}  (inner cherries of SISTER clades A,C)
#   v4-cross : eco-1 = {A1, A2, B1, B2}  (inner cherries of NON-SISTER clades)
#
# True topology ((A,C),(B,D)) unchanged.  In v4-cross the eco signal pulls
# A1,A2 toward B1,B2 ACROSS the true topology divide, so the false
# bipartitions DIRECTLY contradict the truth:
#   falseInner    : {A1,A2,B1,B2}  — breaks clade A and B monophyly
#   falseClade_AB : {A1..A4, B1..B4} — preserves clades, breaks (A,C),(B,D)
#
# Uses SIM3V4_PARAMS$v4c (strong/stress regime: stemBrEco=0.15, stemBrClade=0.05).
# Reuses .BuildV4Data() from sim3v4-helpers.R; only tree-builder and ecology
# vector are v4-cross-specific.
#
# 300 characters (100 neo + 200 trans), phi=4, pi0=0.75, theta=1.

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
PARAM_SET <- "v4c"  # strong/stress regime

outRoot <- file.path("/nobackup", Sys.getenv("USER"),
                     "mkp-sim3-v4cross-c", "results")
dir.create(outRoot, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== Sim3 v4cross-c stress test: rep%d seed=%d nIter=%d ===\n",
            repId, seedBase + repId, nIter))
cat("Output:", outRoot, "\n")

mkpRoot <- Sys.getenv("MKP_REPO_ROOT",
  file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/MkPrime"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3-simulate.R"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3v4-helpers.R"))      # .BuildV4Data
source(file.path(mkpRoot, "inst/simulations/ecology/sim3v4cross-helpers.R")) # v4-cross builders
source(file.path(mkpRoot, "inst/simulations/ecology/sim3v4-params.R"))

params <- SIM3V4_PARAMS[[PARAM_SET]]
cat(sprintf("Params: tipBr=%.3f stemBrEco=%.3f stemBrClade=%.3f rootBr=%.3f\n",
            params$tipBr, params$stemBrEco, params$stemBrClade, params$rootBr))
cat(sprintf("        nNeo=%d nTrans=%d phi=%g pi0=%.2f theta=%.2f\n",
            params$nNeo, params$nTrans, params$phi, params$pi0, params$theta))

# v4-cross tree + ecology (stress: strong eco stem, weak clade stem)
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
  cat(sprintf("\n--- %s chain ---\n", label))
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

cat(sprintf("\n=== v4cross-c summary (rep %d, seed %d) ===\n",
            repId, seedBase + repId))
cat(sprintf("  BLIND: clA=%.3f clB=%.3f AC=%.3f BD=%.3f falseInner=%.3f falseAB=%.3f CID=%.3f\n",
            bScores$pA, bScores$pB, bScores$pAC, bScores$pBD,
            bScores$pFi, bScores$pFc, bScores$cid))
cat(sprintf("  AWARE: clA=%.3f clB=%.3f AC=%.3f BD=%.3f falseInner=%.3f falseAB=%.3f CID=%.3f\n",
            aScores$pA, aScores$pB, aScores$pAC, aScores$pBD,
            aScores$pFi, aScores$pFc, aScores$cid))
cat("\nTarget: blind P(falseInner) or P(falseAB) > 0; aware recovers (A,C),(B,D).\n")
cat("\nDone. Results in", outRoot, "\n")
