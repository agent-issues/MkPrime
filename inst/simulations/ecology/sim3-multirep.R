# sim3-multirep.R -------------------------------------------------------------
# Multi-replicate Sim 3 driver for paper-quality statistics.  Each
# replicate draws fresh simulated data at the Goldilocks config and runs
# both the blind and aware MCMC chains.  Per-replicate bipartition
# probabilities and CID summaries are saved; a final aggregated table is
# printed.
#
# Cost: ~30-60 min per replicate (1 chain x 200k iter, blind + aware).
# Use nReps = 4 for a local sanity check, nReps = 20+ on HPC for the
# headline figure.
#
# Args (commandArgs trailing):
#   1. nReps  (default 4)
#   2. seedBase (default 20260601)
#   3. nIter (default 200000)
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
  library("TreeDist")
})
source("inst/simulations/ecology/sim3-helpers.R")
source("inst/simulations/ecology/sim3-simulate.R")

args <- commandArgs(trailingOnly = TRUE)
nReps    <- if (length(args) >= 1) as.integer(args[1]) else 4L
seedBase <- if (length(args) >= 2) as.integer(args[2]) else 20260601L
nIter    <- if (length(args) >= 3) as.integer(args[3]) else 200000L
cat(sprintf("Sim 3 multi-rep: nReps=%d  seedBase=%d  nIter=%d\n",
            nReps, seedBase, nIter))

# Goldilocks config (from sim3-pilot.R).
nEco <- 60L; nBase <- 180L; phi <- 4; stemBr <- 0.10; rootBr <- 0.15
tree <- .BuildConvergentTree(stemBranch = stemBr, rootBranch = rootBr)
eco  <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)

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

zFull <- matrix(0L, nrow = nEco + nBase, ncol = 2)
zFull[seq_len(nEco), 2] <- 1L
cFull <- c(rep("neomorphic", nEco), rep("transformational", nBase))

results <- vector("list", nReps)
for (rep in seq_len(nReps)) {
  set.seed(seedBase + rep)
  cat(sprintf("\n=== Replicate %d / %d (seed %d) ===\n",
              rep, nReps, seedBase + rep))
  datSim <- .SimulateMkPrimeEcology(tree, edgeEc, zFull, phi = phi,
                                    type = cFull, baseRate = 0.5,
                                    rateLoss = 1)
  pdSim <- MatrixToPhyDat(datSim)
  mkdBlind <- MkPrimeData(pdSim, neomorphic = seq_len(nEco))
  ecoTipVec <- eco[mkdBlind$taxon_names]
  mkdAware <- MkPrimeData(pdSim, neomorphic = seq_len(nEco),
                          ecology = ecoTipVec)

  set.seed(1)
  randTree  <- Preorder(ape::rtree(mkdBlind$nTip,
                                   tip.label = mkdBlind$taxon_names))
  psTrees   <- TreeSearch::MaximizeParsimony(pdSim, tree = randTree,
                                             verbosity = 0)
  startTree <- Preorder(psTrees[[1]])
  startTree$edge.length <- rep(0.1, nrow(startTree$edge))

  modelBlind <- MkPrimeModel(ecologyAware = FALSE,
                             kPrimePrior = "geometric",
                             coding = "variable", expSteps = 10)
  modelAware <- MkPrimeModel(ecologyAware = TRUE,
                             magnitudeMode = "global",
                             kPrimePrior = "geometric",
                             coding = "variable", expSteps = 10)
  mcmc <- MkPrimeMCMC(nIter = nIter, nChains = 1L, nRuns = 1L,
                      thin = max(1L, nIter %/% 500L),
                      treeThin = max(1L, nIter %/% 500L),
                      minWarmup = 5000L, maxWarmup = 40000L,
                      logFile = NULL, checkpointFile = NULL)

  tagBlind <- sprintf("sim3mr-rep%02d-blind", rep)
  tagAware <- sprintf("sim3mr-rep%02d-aware", rep)
  for (tag in c(tagBlind, tagAware)) {
    for (ext in c(".log", ".ckp")) {
      f <- paste0(tag, ext)
      if (file.exists(f)) file.remove(f)
    }
  }
  m <- mcmc; m$logFile <- paste0(tagBlind, ".log")
  resBlind <- RunMkPrime(mkdBlind, tree = startTree,
                         model = modelBlind, mcmc = m)
  m <- mcmc; m$logFile <- paste0(tagAware, ".log")
  resAware <- RunMkPrime(mkdAware, tree = startTree,
                         model = modelAware, mcmc = m)

  trB <- discard(resBlind$trees)
  trA <- discard(resAware$trees)
  class(trB) <- class(trA) <- "multiPhylo"
  scoreOne <- function(trees) {
    list(
      pTrue  = mean(hasBipart(trees, trueSplit)),
      pWrong = mean(hasBipart(trees, wrongSplit)),
      cidTrue  = mean(as.numeric(TreeDist::ClusteringInfoDist(
        trees, tree, normalize = TRUE))),
      cidWrong = mean(as.numeric(TreeDist::ClusteringInfoDist(
        trees,
        Preorder(ape::read.tree(text = paste0(
          "(((A1,A2),(A3,A4)),((B1,B2),(B3,B4)),",
          "(((C1,C2),(C3,C4)),((D1,D2),(D3,D4))));"))),
        normalize = TRUE)))
    )
  }
  sB <- scoreOne(trB); sA <- scoreOne(trA)
  results[[rep]] <- list(blind = sB, aware = sA)
  cat(sprintf("  blind  P(true)=%.3f  P(wrong)=%.3f  CIDtrue=%.3f  CIDwrong=%.3f\n",
              sB$pTrue, sB$pWrong, sB$cidTrue, sB$cidWrong))
  cat(sprintf("  aware  P(true)=%.3f  P(wrong)=%.3f  CIDtrue=%.3f  CIDwrong=%.3f\n",
              sA$pTrue, sA$pWrong, sA$cidTrue, sA$cidWrong))
}

# Aggregate.
agg <- do.call(rbind, lapply(seq_along(results), function(i) {
  r <- results[[i]]
  data.frame(rep = i,
             pTrue_blind  = r$blind$pTrue,  pTrue_aware  = r$aware$pTrue,
             pWrong_blind = r$blind$pWrong, pWrong_aware = r$aware$pWrong,
             cidTrue_blind  = r$blind$cidTrue,  cidTrue_aware  = r$aware$cidTrue,
             cidWrong_blind = r$blind$cidWrong, cidWrong_aware = r$aware$cidWrong)
}))
cat("\n== Per-replicate summary ==\n")
print(agg, row.names = FALSE)
cat("\n== Across replicates (mean +/- sd) ==\n")
for (col in setdiff(colnames(agg), "rep")) {
  x <- agg[[col]]
  cat(sprintf("  %-18s : %.3f +/- %.3f\n", col, mean(x), stats::sd(x)))
}

saveRDS(list(tree = tree, eco = eco, edgeEco = edgeEc,
             nReps = nReps, seedBase = seedBase, nIter = nIter,
             results = results, agg = agg),
        sprintf("inst/simulations/ecology/sim3-multirep-n%d.rds", nReps))
cat("\nDone. Multi-rep result saved.\n")
