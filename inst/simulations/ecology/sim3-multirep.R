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

# Sim 3 v3 config: 2x chars + 3x stem (selected by parsimony-grid.R).
# Stronger phylogenetic signal than v2 Goldilocks; parsimony shows
# nTrueBase=7/8 and nFooled=8/8 (see dev/pilots/2026-05-14-sim3-v3-redesign).
nEco <- 120L; nBase <- 360L; phi <- 4; stemBr <- 0.30; rootBr <- 0.15
tipBr <- 0.5
tree <- .BuildConvergentTree(tipBranch = tipBr,
                              stemBranch = stemBr, rootBranch = rootBr)
eco  <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)

# Use root-invariant HasBipartSplits() — see sim3-scoring.R.
# The agg/results structure saved to RDS bakes pTrue/pWrong values in,
# so any existing sim3-multirep-n*.rds produced before this patch
# carries root-dependent (biased) per-rep support values. Regenerate
# from raw chains to refresh them.
source("inst/simulations/ecology/sim3-scoring.R")
trueSplit  <- c(paste0("A", 1:4), paste0("C", 1:4))
wrongSplit <- c(paste0("A", 1:4), paste0("B", 1:4))
hasBipart <- HasBipartSplits
discard <- function(x) x[seq.int(ceiling(length(x) / 4) + 1L, length(x))]

zFull <- matrix(0L, nrow = nEco + nBase, ncol = 2)
zFull[seq_len(nEco), 2] <- 1L
cFull <- c(rep("neomorphic", nEco), rep("transformational", nBase))

results <- vector("list", nReps)
for (rep in seq_len(nReps)) {
  set.seed(seedBase + rep)
  cat(sprintf("\n=== Replicate %d / %d (seed %d) ===\n",
              rep, nReps, seedBase + rep))
  # Aligned-units simulation: baseRate=1.0 matches model JC kernel;
  # normalize=TRUE + explicit pi0/theta makes simulator gammaE match
  # the model's gammaE at truth. See dev/pilots/2026-05-14-aligned-units-resim.
  datSim <- .SimulateMkPrimeEcology(tree, edgeEc, zFull, phi = phi,
                                    type = cFull, baseRate = 1.0,
                                    normalize = TRUE,
                                    pi0 = 0.75, theta = 1.0,
                                    refEcology = 0L,
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

  # Apply post-hoc relabel to canonical phi >= 1 convention. The likelihood
  # (and topology / logL gap) is invariant under (phi, theta, z) reflection,
  # but RNG-dependent which mode each rep finds. Relabelling makes
  # cross-rep phi/theta summaries interpretable. Robust to zero-sample reps
  # (e.g. warmup ate the whole run on a short config).
  ecoSummary <- tryCatch({
    resAware <- RelabelEcology(resAware)
    samp <- resAware$samples
    if (NROW(samp) == 0L && !is.null(resAware$logFile)) {
      samp <- ReadMkLog(resAware$logFile)
    }
    if (NROW(samp) == 0L) {
      stop("zero samples; rep produced no usable posterior")
    }
    samp <- samp[seq.int(ceiling(nrow(samp) / 4) + 1L, nrow(samp)), ,
                 drop = FALSE]
    phiCols   <- grep("^phi(_|$)",   colnames(samp), value = TRUE)
    thetaCols <- grep("^theta_",     colnames(samp), value = TRUE)
    list(
      phi_median   = vapply(phiCols,   function(c) median(samp[, c]),
                            numeric(1)),
      theta_median = vapply(thetaCols, function(c) median(samp[, c]),
                            numeric(1)),
      pi0_median   = median(samp[, "pi0"]),
      tl_median    = median(samp[, "tree_length"]),
      nRetained    = nrow(samp)
    )
  }, error = function(e) {
    message(sprintf("  [eco-summary skipped: %s]", conditionMessage(e)))
    list(phi_median = NA_real_, theta_median = NA_real_,
         pi0_median = NA_real_, tl_median = NA_real_, nRetained = 0L)
  })

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
  results[[rep]] <- list(blind = sB, aware = sA, eco = ecoSummary)
  cat(sprintf("  blind  P(true)=%.3f  P(wrong)=%.3f  CIDtrue=%.3f  CIDwrong=%.3f\n",
              sB$pTrue, sB$pWrong, sB$cidTrue, sB$cidWrong))
  cat(sprintf("  aware  P(true)=%.3f  P(wrong)=%.3f  CIDtrue=%.3f  CIDwrong=%.3f\n",
              sA$pTrue, sA$pWrong, sA$cidTrue, sA$cidWrong))
  cat(sprintf("  aware  phi_med=%.3f  theta_med=%.3f  pi0_med=%.3f  tl_med=%.3f\n",
              mean(ecoSummary$phi_median),
              mean(ecoSummary$theta_median),
              ecoSummary$pi0_median,
              ecoSummary$tl_median))
}

# Aggregate.
agg <- do.call(rbind, lapply(seq_along(results), function(i) {
  r <- results[[i]]
  data.frame(rep = i,
             pTrue_blind  = r$blind$pTrue,  pTrue_aware  = r$aware$pTrue,
             pWrong_blind = r$blind$pWrong, pWrong_aware = r$aware$pWrong,
             cidTrue_blind  = r$blind$cidTrue,  cidTrue_aware  = r$aware$cidTrue,
             cidWrong_blind = r$blind$cidWrong, cidWrong_aware = r$aware$cidWrong,
             phi_aware_med   = mean(r$eco$phi_median),
             theta_aware_med = mean(r$eco$theta_median),
             pi0_aware_med   = r$eco$pi0_median,
             tl_aware_med    = r$eco$tl_median)
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
