# run_v9.R — Sim 3 v9-induce replicate (blind + aware) on Hamilton.
#
# v9-induce tests whether the dataset configuration that the MP
# pre-screen (dev/sim-design/v9-discriminate.R, commit 8fc37ef) flagged
# as a parsimony trap also fools BLIND ML MCMC, while AWARE ML rescues
# truth.
#
# Pre-screen winner (v9c_M2_realistic): at phi=6, TL=1.48,
#   * MP recovers TRUE sister (A,C) in ~38% of reps
#   * MP groups FALSE (A,B) sister in ~62% of reps
#
# Mechanism M2 "parallel":
#   * Two ecology labels: eco-1 over the full A clade, eco-2 over the
#     full B clade; C and D stay baseline (eco-0).
#   * z matrix is nChar x 3 (baseline + eco-1 + eco-2). For ~50% of
#     chars (p_eco), z[c, 2] = 1L AND z[c, 3] = 1L — both ecologies
#     elevate the gain rate on the SAME characters.  Eco-A and eco-B
#     clades therefore drift toward state 1 in parallel.
#   * All chars `neomorphic`, rateLoss = 1.0; normalize = FALSE.
#
# Uses the auto-expSteps default (commit bc24d52): expSteps is set
# from 1.05 * parsimony(start tree) at .FinalizeModel time. We DO NOT
# pass expSteps to MkPrimeModel().
#
# Root-invariant scoring via inst/ecology/simulations/sim3-scoring.R
# (HasBipartSplits, ScoreTreesUnrooted).
#
# Usage on Hamilton:
#   Rscript run_v9.R <rep_id> <n_iter> <n_chains> <seed_base>
#
# Outputs to /nobackup/$USER/mkp-sim3-v9-induce/results/rep%02d/
#   - blind-chain.log, blind-chain.ckp, blind-result.rds
#   - aware-chain.log, aware-chain.ckp, aware-result.rds
#   - summary.rds   (named list of corrected scores; AC = true, AB = false)

suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKP_LIB",
    file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-v9-induce/lib"))
  .libPaths(c(libPath, .libPaths()))
  library("MkPrime")
  library("TreeTools")
  library("TreeDist")
})

args <- commandArgs(trailingOnly = TRUE)
repId    <- if (length(args) >= 1) as.integer(args[1]) else stop("rep_id required")
nIter    <- if (length(args) >= 2) as.integer(args[2]) else 100000L
nChains  <- if (length(args) >= 3) as.integer(args[3]) else 4L
seedBase <- if (length(args) >= 4) as.integer(args[4]) else 20260520L

outRoot <- file.path("/nobackup", Sys.getenv("USER"),
                     "mkp-sim3-v9-induce",
                     "results", sprintf("rep%02d", repId))
dir.create(outRoot, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== v9-induce rep %d / seed %d / nIter %d / nChains %d ===\n",
            repId, seedBase + repId, nIter, nChains))
cat("Output:", outRoot, "\n")

# Helper + simulator scripts -- reuse repo checkout
mkpRoot <- Sys.getenv("MKP_REPO_ROOT",
  file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-v9-induce/MkPrime"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-helpers.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-simulate.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-scoring.R"))

# --- v9 configuration (matches v9c_M2_realistic in v9-discriminate.R) -----
cfg <- list(
  tipBr = 0.04, stemBrClade = 0.045, stemBrEco = 0.18, rootBr = 0.035,
  nChar = 200L, phi = 6, p_eco = 0.5,
  rateLoss = 1.0, baseRate = 1.0
)
nChar <- cfg$nChar
phi   <- cfg$phi

# --- Build the v9 tree (clones buildTreeV9 from v9-discriminate.R) --------
buildTreeV9 <- function(tipBr, stemBrClade, stemBrEco, rootBr) {
  cladeNw <- function(prefix, tb, sb) {
    sprintf(
      "((%s1:%g,%s2:%g):%g,(%s3:%g,%s4:%g):%g):%g",
      prefix, tb, prefix, tb, tb,
      prefix, tb, prefix, tb, tb, sb
    )
  }
  cladeA <- cladeNw("A", tipBr, stemBrEco)
  cladeB <- cladeNw("B", tipBr, stemBrEco)
  cladeC <- cladeNw("C", tipBr, stemBrClade)
  cladeD <- cladeNw("D", tipBr, stemBrClade)
  newick <- sprintf("((%s,%s):%g,(%s,%s):%g);",
                    cladeA, cladeC, rootBr,
                    cladeB, cladeD, rootBr)
  TreeTools::Preorder(ape::read.tree(text = newick))
}

# Two-ecology tip assignment (clones .TwoEcologyV9): eco-1 on full A
# clade, eco-2 on full B clade. C and D stay baseline 0.
.TwoEcologyV9 <- function(tree) {
  tips <- tree$tip.label
  eco  <- integer(length(tips))
  names(eco) <- tips
  eco[grepl("^A", tips)] <- 1L
  eco[grepl("^B", tips)] <- 2L
  eco
}

tree <- buildTreeV9(cfg$tipBr, cfg$stemBrClade, cfg$stemBrEco, cfg$rootBr)
eco    <- .TwoEcologyV9(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)
truthTL <- sum(tree$edge.length)
cat(sprintf("Truth TL = %.4f (expected ~1.48 for v9c_M2_realistic)\n", truthTL))
cat("Tip ecology:\n"); print(table(eco))
cat("Edge ecology distribution:\n"); print(table(edgeEc))

# --- Build the v9 z matrix: M2 parallel mechanism -------------------------
# nChar x 3 (baseline + eco-1 + eco-2). For ~p_eco of chars, both
# eco-1 and eco-2 push toward state 1 (z = 1L on both eco columns).
set.seed(seedBase + repId)
zMat <- matrix(0L, nrow = nChar, ncol = 3L)
ecoActive <- stats::runif(nChar) < cfg$p_eco
zMat[ecoActive, 2L] <- 1L  # eco-1 elevates gain on these chars
zMat[ecoActive, 3L] <- 1L  # eco-2 elevates gain on the SAME chars
cat(sprintf("Eco-active chars: %d / %d  (p_eco = %.2f)\n",
            sum(ecoActive), nChar, cfg$p_eco))

# --- Simulate (neomorphic only; no normalization) -------------------------
cFull <- rep("neomorphic", nChar)
set.seed(seedBase + repId + 1000L)
datSim <- .SimulateMkPrimeEcology(
  tree, edgeEc, zMat, phi = phi,
  type = cFull, baseRate = cfg$baseRate,
  normalize = FALSE,
  rateLoss = cfg$rateLoss
)

# Filter constant columns (neomorphic with strong drift can create
# all-0 or all-1 columns; MP can't use them and they bloat run time).
varCols <- apply(datSim, 2L, function(col) length(unique(col)) > 1L)
datSim <- datSim[, varCols, drop = FALSE]
nVar <- ncol(datSim)
cat(sprintf("Variable chars after constant filter: %d / %d\n", nVar, nChar))

pdSim    <- MatrixToPhyDat(datSim)
mkdBlind <- MkPrimeData(pdSim, neomorphic = seq_len(nVar))
ecoTipVec <- eco[mkdBlind$taxon_names]
mkdAware <- MkPrimeData(pdSim, neomorphic = seq_len(nVar),
                        ecology = ecoTipVec)

# --- Parsimony start tree -------------------------------------------------
set.seed(1)
randTree <- Preorder(ape::rtree(mkdBlind$nTip,
                                tip.label = mkdBlind$taxon_names))
psTrees <- TreeSearch::MaximizeParsimony(pdSim, tree = randTree,
                                          verbosity = 0)
if (inherits(psTrees, "phylo")) psTrees <- list(psTrees)
startTree <- Preorder(psTrees[[1]])
startTree$edge.length <- rep(0.1, nrow(startTree$edge))
cli::cli_alert_info(
  "Parsimony start tree: TL slot will be 0.1/edge; pre-MCMC parsimony \\
   used as the expSteps anchor by .FinalizeModel."
)

# --- Models: AUTO-expSteps (do NOT pass expSteps) -------------------------
modelBlind <- MkPrimeModel(ecologyAware = FALSE,
                           kPrimePrior = "geometric",
                           coding = "variable")
modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable")
cli::cli_alert_info("Models constructed with auto-expSteps (will resolve in .FinalizeModel)")

mcmcCfg <- function(logFile, ckpFile) {
  MkPrimeMCMC(nIter = nIter, nChains = nChains, nRuns = 1L,
              thin = max(1L, nIter %/% 500L),
              treeThin = max(1L, nIter %/% 500L),
              minWarmup = 5000L, maxWarmup = 40000L,
              logFile = logFile, checkpointFile = ckpFile)
}

run_one <- function(tag, mkd, model) {
  logFile <- file.path(outRoot, paste0(tag, "-chain.log"))
  ckpFile <- file.path(outRoot, paste0(tag, "-chain.ckp"))
  mcmc <- mcmcCfg(logFile, ckpFile)

  cat(sprintf("\n--- %s chain (nChains = %d) ---\n", tag, nChains))
  t0 <- Sys.time()
  if (file.exists(ckpFile)) {
    cat("Resuming from existing checkpoint:", ckpFile, "\n")
    res <- ResumeMkPrime(checkpointFile = ckpFile, data = mkd, model = model)
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

# --- Score using root-invariant ScoreTreesUnrooted ------------------------
# AC = true sister; AB = false (parsimony trap)
biparts <- list(
  AC = c(paste0("A", 1:4), paste0("C", 1:4)),
  AB = c(paste0("A", 1:4), paste0("B", 1:4)),
  cladeA = paste0("A", 1:4),
  cladeB = paste0("B", 1:4),
  cladeC = paste0("C", 1:4),
  cladeD = paste0("D", 1:4)
)

scoreBlind <- ScoreTreesUnrooted(resBlind$trees, biparts, refTree = tree)
scoreAware <- ScoreTreesUnrooted(resAware$trees, biparts, refTree = tree)

summary <- list(
  repId = repId, seed = seedBase + repId,
  nIter = nIter, nChains = nChains,
  config = cfg,
  truthTL = truthTL,
  nVar = nVar,
  blind = scoreBlind,
  aware = scoreAware
)
saveRDS(summary, file.path(outRoot, "summary.rds"))

cat("\n=== v9-induce per-rep scores (AC = true sister, AB = false trap) ===\n")
cat("-- blind --\n")
print(scoreBlind)
cat("-- aware --\n")
print(scoreAware)
cat("\nDone. Saved to", outRoot, "\n")
