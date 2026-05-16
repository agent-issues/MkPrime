# run_v4b.R — Sim 3 v4 (balanced "sweet spot" regime) blind+aware single-rep test.
#
# v4b: moderate eco stem (stemBrEco=0.08), moderate clade stems
# (stemBrClade=0.15).  Convergent signal is real but ancestry signal is
# still substantial.  Target story: blind chain pulled toward the spurious
# {A1,A2,C1,C2} bipartition; aware chain recovers true clade A and C
# monophyly and their AC sister relationship.
#
# Spurious bipartition target: {A1,A2,C1,C2} — breaks clade A and C monophyly.
# True bipartitions to recover:
#   trueClade_A   : {A1,A2,A3,A4}
#   trueClade_C   : {C1,C2,C3,C4}
#   trueSister_AC : {A1..A4,C1..C4}
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
PARAM_SET <- "v4b"

outRoot <- file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-v4b", "results")
dir.create(outRoot, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== Sim3 v4b balanced test: rep%d seed=%d nIter=%d ===\n",
            repId, seedBase + repId, nIter))
cat("Output:", outRoot, "\n")

mkpRoot <- Sys.getenv("MKP_REPO_ROOT",
  file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/MkPrime"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3-simulate.R"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3v4-helpers.R"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3v4-params.R"))

params <- SIM3V4_PARAMS[[PARAM_SET]]
cat(sprintf("Params: tipBr=%.3f stemBrEco=%.3f stemBrClade=%.3f rootBr=%.3f\n",
            params$tipBr, params$stemBrEco, params$stemBrClade, params$rootBr))
cat(sprintf("        nNeo=%d nTrans=%d phi=%g pi0=%.2f theta=%.2f\n",
            params$nNeo, params$nTrans, params$phi, params$pi0, params$theta))

tree   <- .BuildConvergentTreeV4(params$tipBr, params$stemBrEco,
                                  params$stemBrClade, params$rootBr)
eco    <- .ConvergentEcologyV4(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)
biparts <- .SimBipartitionsV4()

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
  pC   <- mean(hasBipart(tr, biparts$trueClade_C))
  pAC  <- mean(hasBipart(tr, biparts$trueSister_AC))
  pFal <- mean(hasBipart(tr, biparts$falseSister))
  cid  <- mean(as.numeric(
    TreeDist::ClusteringInfoDistance(tr, tree, normalize = TRUE)))
  cat(sprintf(
    "[%s] P(clA)=%.3f P(clC)=%.3f P(AC)=%.3f P(eco12)=%.3f CID=%.3f\n",
    label, pA, pC, pAC, pFal, cid))
  invisible(list(pA=pA, pC=pC, pAC=pAC, pFal=pFal, cid=cid))
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

cat(sprintf("\n=== v4b summary (rep %d, seed %d) ===\n", repId, seedBase + repId))
cat(sprintf("  BLIND: P(clA)=%.3f P(clC)=%.3f P(AC)=%.3f P(eco12)=%.3f CID=%.3f\n",
            bScores$pA, bScores$pC, bScores$pAC, bScores$pFal, bScores$cid))
cat(sprintf("  AWARE: P(clA)=%.3f P(clC)=%.3f P(AC)=%.3f P(eco12)=%.3f CID=%.3f\n",
            aScores$pA, aScores$pC, aScores$pAC, aScores$pFal, aScores$cid))
cat("\nTarget: blind P(eco12) high (falseSister supported); aware P(AC)>blind P(AC).\n")
cat("\nDone. Results in", outRoot, "\n")
