# run_rep.R — one Sim 3 multirep replicate (blind + aware) on Hamilton
#
# v4: nChains = 8 (parallel tempering, geometric heat ladder).
# Reuses the mkp-sim3-multirep-v3 library; outputs to
# /nobackup/$USER/mkp-sim3-multirep-v4/results/rep%02d/.
#
# Usage on Hamilton:
#   Rscript run_rep.R <rep_id> <n_iter> <seed_base>

suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKP_LIB",
    file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/lib"))
  .libPaths(c(libPath, .libPaths()))
  library("MkPrime")
  library("TreeTools")
  library("TreeDist")
})

args <- commandArgs(trailingOnly = TRUE)
repId    <- if (length(args) >= 1) as.integer(args[1]) else stop("rep_id required")
nIter    <- if (length(args) >= 2) as.integer(args[2]) else 100000L
seedBase <- if (length(args) >= 3) as.integer(args[3]) else 20260601L

outRoot <- file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v4",
                     "results", sprintf("rep%02d", repId))
dir.create(outRoot, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== Rep %d / seed %d / nIter %d / nChains 8 ===\n",
            repId, seedBase + repId, nIter))
cat("Output:", outRoot, "\n")

mkpRoot <- Sys.getenv("MKP_REPO_ROOT",
  file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/MkPrime"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3-helpers.R"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3-simulate.R"))

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
randTree <- Preorder(ape::rtree(mkdBlind$nTip,
                                tip.label = mkdBlind$taxon_names))
psTrees <- TreeSearch::MaximizeParsimony(pdSim, tree = randTree,
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

mcmcCommon <- function(logFile, ckpFile, nChains = 8L) {
  MkPrimeMCMC(nIter = nIter, nChains = nChains, heat = 0.2, nRuns = 1L,
              thin = max(1L, nIter %/% 500L),
              treeThin = max(1L, nIter %/% 500L),
              minWarmup = 5000L, maxWarmup = 40000L,
              logFile = logFile, checkpointFile = ckpFile)
}

run_one <- function(tag, mkd, model, nChains = 8L) {
  logFile <- file.path(outRoot, paste0(tag, "-chain.log"))
  ckpFile <- file.path(outRoot, paste0(tag, "-chain.ckp"))
  mcmc <- mcmcCommon(logFile, ckpFile, nChains = nChains)

  cat(sprintf("\n--- %s chain (nChains=%d) ---\n", tag, nChains))
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

scoreOne <- function(trees) {
  list(
    pTrue  = mean(hasBipart(trees, trueSplit)),
    pWrong = mean(hasBipart(trees, wrongSplit)),
    cidTrue = mean(as.numeric(TreeDist::ClusteringInfoDistance(
      trees, tree, normalize = TRUE)))
  )
}

trB <- discard(resBlind$trees); class(trB) <- "multiPhylo"
trA <- discard(resAware$trees); class(trA) <- "multiPhylo"
sB <- scoreOne(trB); sA <- scoreOne(trA)

summary <- list(
  repId = repId, seed = seedBase + repId, nIter = nIter,
  config = list(nEco = nEco, nBase = nBase, phi = phi,
                stemBr = stemBr, rootBr = rootBr, tipBr = tipBr),
  truthTL = truthTL,
  blind = sB, aware = sA
)
saveRDS(summary, file.path(outRoot, "summary.rds"))
cat("\n=== Per-rep scores ===\n")
cat(sprintf("blind  P(true)=%.3f  P(wrong)=%.3f  CIDtrue=%.3f\n",
            sB$pTrue, sB$pWrong, sB$cidTrue))
cat(sprintf("aware  P(true)=%.3f  P(wrong)=%.3f  CIDtrue=%.3f\n",
            sA$pTrue, sA$pWrong, sA$cidTrue))
if (!is.null(resAware$swap_rates)) {
  cat(sprintf("aware  betas: %s\n",
              paste(round(resAware$betas, 3), collapse = " -> ")))
  cat(sprintf("aware  swap rates: %s\n",
              paste(round(resAware$swap_rates, 3), collapse = " / ")))
}
cat("\nDone. Saved to", outRoot, "\n")
