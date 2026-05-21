# run_phi2_fixedpriors.R — phi=2 rerun with tightened priors (R5-1/R5-2/R5-3 fix).
#
# Red-team audit identified a pi0 feedback loop caused by weak priors:
#   R5-1: sigmaPhi=1.0 was too concentrated at null — raised to 1.5.
#   R5-2: thetaAlpha=thetaBeta=1 (uniform) allows drift to 0.5 saddle — raised
#          to Beta(2,2) to softly centre theta and prevent saddle trapping.
#   R5-3: pi0 prior ESS=100 (Beta(75,25)) too weak vs ~480 z-cell observations
#          — raised to Beta(360,120) (ESS=480) so the sparsity prior can resist
#          random z reorganisation driving pi0 << 0.5.
#
# These fixes are passed explicitly as MkPrimeModel() hyperparameters and do
# NOT require reinstalling the package on Hamilton.
#
# Diagnostic question: do tightened priors bring pi0 closer to truth (0.75)
# and does aware P(true AC) improve vs blind?
#
# Both blind and aware models receive the shared prior fixes. Only modelAware
# uses the ecology-specific parameters (rho0Alpha/Beta, thetaAlpha/Beta,
# sigmaPhi); they are harmless no-ops in modelBlind (ecologyAware=FALSE).

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
                     "results-fixedpriors")
dir.create(outRoot, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== phi=%g fixed-priors test: rep%d seed=%d nIter=%d ===\n",
            phiSim, repId, seedBase + repId, nIter))
cat("Output:", outRoot, "\n")

# Prior hyperparameters — R5-1/R5-2/R5-3 fixes
SIGMA_PHI    <- 1.5   # was 1.0; LogNormal sdlog loosened to reduce null pull
THETA_ALPHA  <- 2     # was 1; Beta(2,2) softly centres theta, avoids 0.5 saddle
THETA_BETA   <- 2     # was 1
RHO0_ALPHA   <- 360   # was 75; ESS=480 vs old ESS=100 (pi0 sparsity tightened)
RHO0_BETA    <- 120   # was 25; Beta(360,120) mode=0.75, ESS=480

cat(sprintf(
  "Fixed priors: sigmaPhi=%.1f  thetaAlpha=%d  thetaBeta=%d  rho0Alpha=%d  rho0Beta=%d\n",
  SIGMA_PHI, THETA_ALPHA, THETA_BETA, RHO0_ALPHA, RHO0_BETA))

mkpRoot <- Sys.getenv("MKP_REPO_ROOT",
  file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/MkPrime"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-helpers.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-simulate.R"))

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
startTree <- Preorder(TreeSearch::AdditionTree(pdSim))
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
  cat(sprintf("\n--- %s chain (phi_sim=%g, nChains=1, fixed-priors) ---\n",
              label, phiSim))
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

  # --- Ecology parameter summaries (pi0 and phi) for prior-fix verification ---
  if (isTRUE(model$ecologyAware) && !is.null(res$log)) {
    logDf <- res$log
    # pi0
    if ("pi0" %in% names(logDf)) {
      pi0samp <- discard(logDf$pi0)
      cat(sprintf("[%s] pi0: median=%.3f  mean=%.3f  (truth=0.75)\n",
                  label, median(pi0samp), mean(pi0samp)))
    }
    # phi (may be scalar or first column if vector)
    phiCols <- grep("^phi", names(logDf), value = TRUE)
    if (length(phiCols) > 0L) {
      phisamp <- discard(logDf[[phiCols[1]]])
      cat(sprintf("[%s] %s: median=%.3f  mean=%.3f  (truth=%.1f)\n",
                  label, phiCols[1], median(phisamp), mean(phisamp), phiSim))
    }
  }

  invisible(list(pTrue = pTrue, pWrong = pWrong, cid = cidTrue))
}

# modelBlind: ecology-unaware baseline.
# sigmaPhi/thetaAlpha/thetaBeta/rho0Alpha/rho0Beta are ignored when
# ecologyAware=FALSE — listed here for documentation clarity only.
modelBlind <- MkPrimeModel(ecologyAware = FALSE,
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)

# modelAware: R5-1/R5-2/R5-3 fixed priors applied explicitly.
modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10,
                           sigmaPhi   = SIGMA_PHI,
                           thetaAlpha = THETA_ALPHA,
                           thetaBeta  = THETA_BETA,
                           rho0Alpha  = RHO0_ALPHA,
                           rho0Beta   = RHO0_BETA)

bScores <- runOne(mkdBlind, modelBlind, "blind-fixedpriors")
aScores <- runOne(mkdAware, modelAware, "aware-fixedpriors")

cat(sprintf("\n=== phi=%g fixed-priors summary (rep %d, seed %d) ===\n",
            phiSim, repId, seedBase + repId))
cat(sprintf("  BLIND: P(true)=%.3f  P(wrong)=%.3f  CID=%.3f\n",
            bScores$pTrue, bScores$pWrong, bScores$cid))
cat(sprintf("  AWARE: P(true)=%.3f  P(wrong)=%.3f  CID=%.3f\n",
            aScores$pTrue, aScores$pWrong, aScores$cid))
cat(sprintf("Baseline (weak priors, job 17185756): aware P(true)=? P(wrong)=?\n"))
cat("\nPosture 2 success criterion: aware P(true) substantially > blind P(true).\n")
cat("Prior-fix criterion: aware pi0 median closer to 0.75 than weak-prior run.\n")

cat("\nDone. Results in", outRoot, "\n")
