# v7-discriminate.R — pre-submission discriminating check for v7-induce.
#
# Goal: find parameters where parsimony(AB-grouped) < parsimony(truth) by
# at least ~5-10 steps across reps, at truth TL <= 2.0.
#
# Three levers vs v6: reduce nChar, increase eco confound strength, drift
# TL up modestly.

suppressPackageStartupMessages({
  library("MkPrime")
  library("TreeTools")
  library("TreeSearch")
  library("ape")
})

mkpRoot <- normalizePath("inst/simulations/ecology")
source(file.path(mkpRoot, "sim3-helpers.R"))
source(file.path(mkpRoot, "sim3-simulate.R"))
source(file.path(mkpRoot, "sim3-scoring.R"))

# Reuse v6 builder: takes (tipBr, stemBrClade, stemBrEco, rootBr)
buildTreeV7 <- function(tipBr, stemBrClade, stemBrEco, rootBr) {
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

buildAB <- function(tipBr, stemBrClade, stemBrEco, rootBr) {
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
                    cladeA, cladeB, rootBr,
                    cladeC, cladeD, rootBr)
  TreeTools::Preorder(ape::read.tree(text = newick))
}

# v7 variants
variant <- Sys.getenv("V7_VARIANT", "a")
cfg <- switch(variant,
  # User's proposed starting point
  a = list(
    name        = "v7a-proposed",
    tipBr       = 0.04,
    stemBrClade = 0.045,
    stemBrEco   = 0.20,
    rootBr      = 0.04,
    nNeo        = 25L,
    nTrans      = 75L,
    phi         = 10,
    pi0         = 0.30,
    theta       = 1.0,
    baseRate    = 1.0,
    rateLoss    = 1
  ),
  # Stronger eco
  b = list(
    name        = "v7b-stronger",
    tipBr       = 0.04,
    stemBrClade = 0.045,
    stemBrEco   = 0.25,
    rootBr      = 0.035,
    nNeo        = 25L,
    nTrans      = 75L,
    phi         = 12,
    pi0         = 0.20,
    theta       = 1.0,
    baseRate    = 1.0,
    rateLoss    = 1
  ),
  # Push to TL ~ 1.8
  c = list(
    name        = "v7c-pushTL",
    tipBr       = 0.045,
    stemBrClade = 0.05,
    stemBrEco   = 0.28,
    rootBr      = 0.03,
    nNeo        = 25L,
    nTrans      = 75L,
    phi         = 12,
    pi0         = 0.20,
    theta       = 1.0,
    baseRate    = 1.0,
    rateLoss    = 1
  ),
  # Aggressive: very short root, very long eco stems
  d = list(
    name        = "v7d-aggressive",
    tipBr       = 0.04,
    stemBrClade = 0.04,
    stemBrEco   = 0.30,
    rootBr      = 0.025,
    nNeo        = 25L,
    nTrans      = 75L,
    phi         = 15,
    pi0         = 0.15,
    theta       = 1.0,
    baseRate    = 1.0,
    rateLoss    = 1
  ),
  stop("Unknown variant: ", variant)
)
cat(sprintf("Variant: %s\n", cfg$name))

tree <- buildTreeV7(cfg$tipBr, cfg$stemBrClade, cfg$stemBrEco, cfg$rootBr)
truthTL <- sum(tree$edge.length)
cat(sprintf("Truth TL = %.4f\n", truthTL))
if (truthTL > 2.0) {
  cat(sprintf("WARNING: TL %.3f exceeds 2.0 ceiling\n", truthTL))
}

abTree <- buildAB(cfg$tipBr, cfg$stemBrClade, cfg$stemBrEco, cfg$rootBr)
eco  <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)
cat("Edge ecology assignment:\n")
print(table(edgeEc))
ecoBranchTotal <- sum(tree$edge.length[edgeEc == 1L])
cat(sprintf("Total eco-1 branch length: %.3f (%.0f%% of TL)\n",
            ecoBranchTotal, 100 * ecoBranchTotal / truthTL))

# Predicted parsimony deltas
pEcoStemChange <- (1 - exp(-cfg$phi * cfg$stemBrEco)) / 2
pParallel <- pEcoStemChange^2
ecoChars <- (cfg$nNeo + cfg$nTrans) * (1 - cfg$pi0)
expectedEcoSyn <- ecoChars * pParallel
expectedRootSyn <- cfg$rootBr * (cfg$nNeo + cfg$nTrans)
cat(sprintf(
  "Predicted: rootSyn~%.2f, P(parallel)=%.3f, ecoSyn~%.2f, delta(AB-truth)~%.2f\n",
  expectedRootSyn, pParallel, expectedEcoSyn,
  expectedRootSyn - expectedEcoSyn))

# --- Run N reps of simulate + MP search ------------------------------------
nReps <- as.integer(Sys.getenv("V7_NREPS", "8"))
seedBase <- 20260720L
nNeo <- cfg$nNeo; nTrans <- cfg$nTrans
nChar <- nNeo + nTrans
phi <- cfg$phi

results <- list()
for (rep in seq_len(nReps)) {
  set.seed(seedBase + rep)
  zFull <- matrix(0L, nrow = nChar, ncol = 2L)
  ecoActive <- stats::runif(nChar) > cfg$pi0
  zFull[ecoActive, 2L] <- 1L
  cFull <- c(rep("neomorphic", nNeo), rep("transformational", nTrans))
  datSim <- .SimulateMkPrimeEcology(
    tree, edgeEc, zFull, phi = phi,
    type = cFull, baseRate = cfg$baseRate,
    normalize = TRUE,
    pi0 = cfg$pi0, theta = cfg$theta,
    refEcology = 0L, rateLoss = cfg$rateLoss
  )

  pdSim <- MatrixToPhyDat(datSim)
  mkdBlind <- MkPrimeData(pdSim, neomorphic = seq_len(nNeo))

  parsTrue <- as.integer(TreeSearch::TreeLength(tree, pdSim))
  parsAB   <- as.integer(TreeSearch::TreeLength(abTree, pdSim))

  set.seed(rep)
  randTree <- TreeTools::Preorder(
    ape::rtree(mkdBlind$nTip, tip.label = mkdBlind$taxon_names))
  mpTrees <- TreeSearch::MaximizeParsimony(pdSim, tree = randTree,
                                            verbosity = 0)
  if (inherits(mpTrees, "phylo")) mpTrees <- list(mpTrees)
  bestTree <- mpTrees[[1]]
  parsBest <- as.integer(TreeSearch::TreeLength(bestTree, pdSim))

  acTips <- c(paste0("A", 1:4), paste0("C", 1:4))
  abTips <- c(paste0("A", 1:4), paste0("B", 1:4))
  hasAC <- HasBipartSplits(list(bestTree), acTips)
  hasAB <- HasBipartSplits(list(bestTree), abTips)

  cat(sprintf(
    "rep %d: pars(true)=%d, pars(AB)=%d, pars(MP)=%d  delta(AB-true)=%d  MP has AC=%s AB=%s\n",
    rep, parsTrue, parsAB, parsBest, parsAB - parsTrue, hasAC, hasAB))

  results[[rep]] <- list(rep = rep,
                          parsTrue = parsTrue,
                          parsAB   = parsAB,
                          parsBest = parsBest,
                          hasAC = hasAC,
                          hasAB = hasAB)
}

cat("\n=== Summary ===\n")
parsTrueV <- vapply(results, `[[`, numeric(1), "parsTrue")
parsABV   <- vapply(results, `[[`, numeric(1), "parsAB")
parsBestV <- vapply(results, `[[`, numeric(1), "parsBest")
acV <- vapply(results, `[[`, logical(1), "hasAC")
abV <- vapply(results, `[[`, logical(1), "hasAB")
cat(sprintf("Truth TL: %.3f\n", truthTL))
cat(sprintf("Mean parsimony true: %.1f\n", mean(parsTrueV)))
cat(sprintf("Mean parsimony AB:   %.1f\n", mean(parsABV)))
cat(sprintf("Mean parsimony best: %.1f\n", mean(parsBestV)))
cat(sprintf("Mean (AB - true) step difference: %.2f  (NEGATIVE = AB beats truth)\n",
            mean(parsABV - parsTrueV)))
cat(sprintf("Per-rep (AB - true): %s\n",
            paste(parsABV - parsTrueV, collapse = ", ")))
cat(sprintf("Fraction MP best has AC split: %.2f\n", mean(acV)))
cat(sprintf("Fraction MP best has AB split: %.2f\n", mean(abV)))
cat(sprintf("\nPredicted auto expSteps = 1.05 * parsimony best ~ %.1f\n",
            1.05 * mean(parsBestV)))
cat(sprintf("parsimony/(TL*nChar) ratio = %.3f\n",
            mean(parsBestV) / (truthTL * nChar)))
