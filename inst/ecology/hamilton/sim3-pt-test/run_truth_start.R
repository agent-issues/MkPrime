# run_truth_start.R — diagnostic: start aware chain at truth tree (rep 05 seed).
#
# Answers the key strategic question: is the AC bipartition actually an
# attractive posterior mode, or just high-logPost at MAP nuisance?
#
# If the chain STAYS at AC (P(true) ≈ 1, CID low) → AC is genuinely the mode
#   and the mixing problem is pure topology-find. Eco-aware proposals will help.
# If the chain DRIFTS OFF AC → AC is a conditional ridge, not the integrated
#   posterior mode. The +106-nat gap was a counterfactual artefact.
#
# Setup: identical to run_pt_test.R (rep 05 seed 20260606, v3 config, aware model)
# but startTree = truthTree with truth branch lengths; nChains = 1.

suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKP_LIB",
    file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/lib"))
  .libPaths(c(libPath, .libPaths()))
  library("MkPrime")
  library("TreeTools")
  library("TreeDist")
})

repId    <- 5L
nIter    <- 100000L
nChains  <- 1L
seedBase <- 20260601L

outRoot <- file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-truth-start",
                     "results")
dir.create(outRoot, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== Truth-start diagnostic: rep%d seed=%d nIter=%d ===\n",
            repId, seedBase + repId, nIter))
cat("Output:", outRoot, "\n")

mkpRoot <- Sys.getenv("MKP_REPO_ROOT",
  file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/MkPrime"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-helpers.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-simulate.R"))

# V3 config (matches sim3-multirep-v3 and run_pt_test.R)
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

# START AT TRUTH TREE (with truth branch lengths)
startTree <- Preorder(tree)

modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)

logFile <- file.path(outRoot, "truth-start-aware.log")
ckpFile <- file.path(outRoot, "truth-start-aware.ckp")
mcmc <- MkPrimeMCMC(nIter = nIter, nChains = nChains,
                    nRuns = 1L,
                    thin = max(1L, nIter %/% 500L),
                    treeThin = max(1L, nIter %/% 500L),
                    minWarmup = 5000L, maxWarmup = 40000L,
                    logFile = logFile, checkpointFile = ckpFile)

cat("\n--- AWARE chain (truth start, nChains=1) ---\n")
t0 <- Sys.time()
if (file.exists(ckpFile)) {
  cat("Resuming from checkpoint:", ckpFile, "\n")
  res <- ResumeMkPrime(mkdAware, mcmc = mcmc, model = modelAware)
} else {
  res <- RunMkPrime(mkdAware, tree = startTree,
                    model = modelAware, mcmc = mcmc)
}
cat(sprintf("Elapsed: %s\n", format(Sys.time() - t0)))
res <- RelabelEcology(res)
saveRDS(res, file.path(outRoot, "truth-start-result.rds"))

# Score
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

# TL trace summary
if (!is.null(res$logFile)) {
  samp <- tryCatch(ReadMkLog(res$logFile), error = function(e) NULL)
} else {
  samp <- tryCatch(as.data.frame(res$samples), error = function(e) NULL)
}
if (!is.null(samp) && "tree_length" %in% names(samp)) {
  tl <- samp$tree_length
  cat(sprintf("TL: min=%.2f  med=%.2f  max=%.2f  last100mean=%.2f  truth=%.2f\n",
              min(tl), median(tl), max(tl), mean(tail(tl, 100)), truthTL))
}

cat("\n=== Truth-start scores ===\n")
cat(sprintf("P(true AC)=%.3f  P(wrong AB)=%.3f  CID-to-truth=%.3f\n",
            pTrue, pWrong, cidTrue))
cat("Interpretation:\n")
if (pTrue > 0.5) {
  cat("  Chain STAYED near AC -> AC is genuinely attractive. Fix: eco-aware proposals.\n")
} else if (pTrue > 0.05) {
  cat("  Chain PARTIALLY retained AC -> AC has weight but weak attraction.\n")
} else {
  cat("  Chain DRIFTED AWAY from AC -> AC may be a conditional ridge, not integrated mode.\n")
  cat("  The +106-nat MAP-nuisance advantage may not reflect full posterior.\n")
}
cat(sprintf("Compare PT cold chain (parsimony start): P(true)=0.000 P(wrong)=0.186 CID=0.445\n"))
cat(sprintf("Compare 1-chain multirep rep 05:          P(true)=0.000 P(wrong)=0.462 CID=0.473\n"))

cat("\nDone. Result + log + ckp in", outRoot, "\n")
