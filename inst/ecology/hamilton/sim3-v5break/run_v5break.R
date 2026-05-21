# run_v5break.R — Sim 3 v5-break (BREAK-BLIND regime) blind+aware PT array.
#
# PURPOSE
#   Across 7+ v4 / v4-cross regimes the (correctly-scored) result has been
#   that blind essentially recovers truth at 16 tips x 200-300 chars,
#   phi = 4. The ecology confound is too weak relative to clade-stem
#   ancestry signal. The methods paper needs a configuration where blind
#   genuinely fails so aware can demonstrate rescue.
#
# DESIGN (see SIM3V4_PARAMS$v5break for full rationale and arithmetic)
#   v5break flips the signal-to-noise ratio by:
#     (i)   shrinking clade-stem ancestry (stemBrClade = 0.03);
#     (ii)  shrinking the matrix (80 chars total) so per-clade synapomorphy
#           counts do not grow with character count;
#     (iii) keeping the eco stem in the UNSATURATED band
#           (stemBrEco = 0.10, phi = 6) so spurious convergence accumulates
#           without saturating;
#     (iv)  lowering pi0 to 0.45 (55 % of chars ecology-encoded).
#
#   Predicted per-character expected synapomorphies:
#     falseInner (parallel changes on both eco stems): ~ 5.4
#     trueClade_A   (changes on clade A stem)        : ~ 2.3
#     trueSister_AC (changes on (A,C) root edge)     : ~ 3.8
#   Ratio falseInner : trueClade_A ~ 2.3:1 -- blind should genuinely
#   struggle to keep clade A (and B) monophyletic.
#
# GEOMETRY
#   v4-cross: eco-1 tips = {A1, A2, B1, B2} (NON-sister clade inner cherries).
#   The spurious bipartition {A1,A2,B1,B2} directly contradicts the true
#   ((A,C),(B,D)) clade-level topology.
#
# SCORING
#   Uses root-invariant HasBipartSplits / ScoreTreesUnrooted from
#   inst/ecology/simulations/sim3-scoring.R (NOT the legacy inline
#   prop.part-based hasBipart). Same as the rescore harness applied to
#   v4cross-b-pt5.
#
# MCMC
#   PT 4-chain (single-chain aware mode-traps on hard configs).
#   nIter = 100,000; thin and treeThin = nIter / 500 so ~500 retained
#   samples per chain.

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
seedBase <- 20260619L
PARAM_SET <- "v5break"

outRoot <- file.path("/nobackup", Sys.getenv("USER"),
                     "mkp-sim3-v5break", "results",
                     sprintf("rep%02d", repId))
dir.create(outRoot, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== Sim3 v5break PT (nChains=%d): rep%02d seed=%d nIter=%d ===\n",
            nChains, repId, seedBase + repId, nIter))
cat("Output:", outRoot, "\n")

mkpRoot <- Sys.getenv("MKP_REPO_ROOT",
  file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/MkPrime"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-simulate.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3v4-helpers.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3v4cross-helpers.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3v4-params.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-scoring.R"))

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
startTree <- Preorder(TreeSearch::AdditionTree(v4dat$pdSim))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

scoreReport <- function(treeList, label) {
  df <- ScoreTreesUnrooted(treeList, biparts, refTree = tree,
                           discardBurnin = TRUE)
  # Tidy named vector for the single-line report.
  get <- function(nm) df$value_corrected[df$metric == nm]
  pA  <- get("trueClade_A"); pB <- get("trueClade_B")
  pC  <- get("trueClade_C"); pD <- get("trueClade_D")
  pAC <- get("trueSister_AC"); pBD <- get("trueSister_BD")
  pFi <- get("falseInner"); pFc <- get("falseClade_AB")
  cid <- get("CID")
  cat(sprintf(
    "[%s] clA=%.3f clB=%.3f clC=%.3f clD=%.3f | AC=%.3f BD=%.3f | falseInner=%.3f falseAB=%.3f CID=%.3f\n",
    label, pA, pB, pC, pD, pAC, pBD, pFi, pFc, cid))
  invisible(list(df = df,
                 pA = pA, pB = pB, pC = pC, pD = pD,
                 pAC = pAC, pBD = pBD,
                 pFi = pFi, pFc = pFc, cid = cid))
}

runOne <- function(mkd, model, label) {
  cat(sprintf("\n--- %s chain (nChains=%d, rep%02d) ---\n",
              label, nChains, repId))
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
  scoreReport(res$trees, label)
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

cat(sprintf("\n=== v5break summary (rep %02d, seed %d, nChains=%d) ===\n",
            repId, seedBase + repId, nChains))
cat(sprintf("  BLIND: clA=%.3f clB=%.3f clC=%.3f clD=%.3f AC=%.3f BD=%.3f falseInner=%.3f falseAB=%.3f CID=%.3f\n",
            bScores$pA, bScores$pB, bScores$pC, bScores$pD,
            bScores$pAC, bScores$pBD,
            bScores$pFi, bScores$pFc, bScores$cid))
cat(sprintf("  AWARE: clA=%.3f clB=%.3f clC=%.3f clD=%.3f AC=%.3f BD=%.3f falseInner=%.3f falseAB=%.3f CID=%.3f\n",
            aScores$pA, aScores$pB, aScores$pC, aScores$pD,
            aScores$pAC, aScores$pBD,
            aScores$pFi, aScores$pFc, aScores$cid))
cat("\nTarget: BLIND falseInner > 0.1 OR BLIND AC < 0.5; AWARE recovers truth.\n")
cat("\nDone. Results in", outRoot, "\n")
