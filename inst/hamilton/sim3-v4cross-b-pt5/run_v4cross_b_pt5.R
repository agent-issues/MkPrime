# run_v4cross_b_pt5.R — Sim 3 v4-cross balanced regime, PT 5-rep array sweep.
#
# Clone of run_v4cross_b_pt.R parameterised by SLURM_ARRAY_TASK_ID.
# Goal: determine whether the clade-C collapse in v4cross-b-pt rep 1
# (P=0.219 while clB/clA recovered) is stochastic mode-flipping or systematic.
#
# Each array task uses seed = seedBase + repId, output to per-rep subdir.

suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKP_LIB",
    file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/lib"))
  .libPaths(c(libPath, .libPaths()))
  library("MkPrime")
  library("TreeTools")
  library("TreeDist")
})

repId    <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
nIter    <- 100000L
nChains  <- 4L
seedBase <- 20260601L
PARAM_SET <- "v4b"  # reuse v4b balanced regime

outRoot <- file.path("/nobackup", Sys.getenv("USER"),
                     "mkp-sim3-v4cross-b-pt5", "results",
                     sprintf("rep%02d", repId))
dir.create(outRoot, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== Sim3 v4cross-b-pt5 (nChains=%d): rep%d seed=%d nIter=%d ===\n",
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
  cat(sprintf("\n--- %s chain (nChains=%d, rep%02d) ---\n", label, nChains, repId))
  logFile <- file.path(outRoot, sprintf("%s-rep%02d.log", label, repId))
  ckpFile <- file.path(outRoot, sprintf("%s-rep%02d.ckp", label, repId))
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
  saveRDS(res, file.path(outRoot, sprintf("%s-rep%02d-result.rds", label, repId)))
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

cat(sprintf("\n=== v4cross-b-pt5 summary (rep %02d, seed %d, nChains=%d) ===\n",
            repId, seedBase + repId, nChains))
cat(sprintf("  BLIND: clA=%.3f clB=%.3f clC=%.3f clD=%.3f AC=%.3f BD=%.3f falseInner=%.3f falseAB=%.3f CID=%.3f\n",
            bScores$pA, bScores$pB, bScores$pC, bScores$pD,
            bScores$pAC, bScores$pBD,
            bScores$pFi, bScores$pFc, bScores$cid))
cat(sprintf("  AWARE: clA=%.3f clB=%.3f clC=%.3f clD=%.3f AC=%.3f BD=%.3f falseInner=%.3f falseAB=%.3f CID=%.3f\n",
            aScores$pA, aScores$pB, aScores$pC, aScores$pD,
            aScores$pAC, aScores$pBD,
            aScores$pFi, aScores$pFc, aScores$cid))
cat("\nKey diagnostic: did AWARE clade C recover (was P=0.219 in single rep)?\n")
cat("\nDone. Results in", outRoot, "\n")
