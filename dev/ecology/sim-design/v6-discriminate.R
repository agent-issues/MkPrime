# v6-discriminate.R — pre-submission discriminating check for v6-realistic.
#
# Tests whether the proposed v6 parameter set produces simulated data on
# which a parsimony tree search either (a) recovers the truth (uninformative)
# or (b) finds an AB-grouped tree within a small step margin of the truth
# (blind might fail; submit). Quick local R sanity check before Hamilton.
#
# Runs N replicates, simulates each, runs MaximizeParsimony on each, reports:
#   - Truth TL (deterministic)
#   - Fitch score (parsimony) of true tree vs MP tree
#   - Whether MP tree has the AC split (truth) or AB split (false eco clade)
#   - Step difference between parsimony of true vs AB-grouped trees

suppressPackageStartupMessages({
  library("MkPrime")
  library("TreeTools")
  library("TreeSearch")
  library("ape")
})

mkpRoot <- normalizePath("inst/ecology/simulations")
source(file.path(mkpRoot, "sim3-helpers.R"))
source(file.path(mkpRoot, "sim3-simulate.R"))
source(file.path(mkpRoot, "sim3-scoring.R"))

# -- candidate parameter sets -------------------------------------------------
# Two variants to compare:
#   v6a: balanced stems (single stemBranch from the v3 helper as-written)
#   v6b: longer eco stems (requires extending the helper or building newick
#        inline). We do this inline below.

buildTreeV6 <- function(tipBr, stemBrClade, stemBrEco, rootBr) {
  # Multirep-v3 topology: ((A,C),(B,D)) with 4 clades of 4 tips each.
  # Eco-1 = full clades A and B; eco-0 = full clades C and D.
  # Eco clades get stemBrEco; non-eco clades get stemBrClade.
  cladeNw <- function(prefix, tb, sb) {
    sprintf(
      "((%s1:%g,%s2:%g):%g,(%s3:%g,%s4:%g):%g):%g",
      prefix, tb, prefix, tb, tb,
      prefix, tb, prefix, tb, tb, sb
    )
  }
  cladeA <- cladeNw("A", tipBr, stemBrEco)    # eco-1 clade
  cladeB <- cladeNw("B", tipBr, stemBrEco)    # eco-1 clade
  cladeC <- cladeNw("C", tipBr, stemBrClade)  # eco-0 clade
  cladeD <- cladeNw("D", tipBr, stemBrClade)  # eco-0 clade
  newick <- sprintf("((%s,%s):%g,(%s,%s):%g);",
                    cladeA, cladeC, rootBr,
                    cladeB, cladeD, rootBr)
  TreeTools::Preorder(ape::read.tree(text = newick))
}

# v6-realistic candidate parameters --
# Default: pass via env var V6_VARIANT (a / b / c)
variant <- Sys.getenv("V6_VARIANT", "b")
cfg <- switch(variant,
  a = list(
    name        = "v6a-mild",
    tipBr       = 0.03,
    stemBrClade = 0.04,
    stemBrEco   = 0.10,
    rootBr      = 0.03,
    nNeo        = 60L,
    nTrans      = 140L,
    phi         = 6,
    pi0         = 0.50,
    theta       = 1.0,
    baseRate    = 1.0,
    rateLoss    = 1
  ),
  b = list(
    name        = "v6b-strong-eco",
    tipBr       = 0.02,
    stemBrClade = 0.025,
    stemBrEco   = 0.14,
    rootBr      = 0.02,
    nNeo        = 60L,
    nTrans      = 140L,
    phi         = 8,
    pi0         = 0.40,
    theta       = 1.0,
    baseRate    = 1.0,
    rateLoss    = 1
  ),
  c = list(
    name        = "v6c-very-strong",
    tipBr       = 0.02,
    stemBrClade = 0.02,
    stemBrEco   = 0.18,
    rootBr      = 0.015,
    nNeo        = 80L,
    nTrans      = 200L,
    phi         = 10,
    pi0         = 0.35,
    theta       = 1.0,
    baseRate    = 1.0,
    rateLoss    = 1
  ),
  # v6d: target TL ~ 1.2 like real morphology, phi=8, moderate-strong eco,
  # but with enough characters so signal accumulates over many.
  d = list(
    name        = "v6d-balanced-realistic",
    tipBr       = 0.030,
    stemBrClade = 0.035,
    stemBrEco   = 0.15,
    rootBr      = 0.035,
    nNeo        = 80L,
    nTrans      = 220L,
    phi         = 8,
    pi0         = 0.45,
    theta       = 1.0,
    baseRate    = 1.0,
    rateLoss    = 1
  ),
  # v6e: like v6d but slightly weaker eco, more eco-coded characters
  e = list(
    name        = "v6e-numerous-mild-eco",
    tipBr       = 0.025,
    stemBrClade = 0.03,
    stemBrEco   = 0.10,
    rootBr      = 0.03,
    nNeo        = 80L,
    nTrans      = 220L,
    phi         = 6,
    pi0         = 0.30,
    theta       = 1.0,
    baseRate    = 1.0,
    rateLoss    = 1
  ),
  stop("Unknown variant: ", variant)
)
cat(sprintf("Variant: %s\n", cfg$name))

tree <- buildTreeV6(cfg$tipBr, cfg$stemBrClade, cfg$stemBrEco, cfg$rootBr)
truthTL <- sum(tree$edge.length)
cat(sprintf("Truth TL = %.4f\n", truthTL))
cat(sprintf("  4 eco-1 tips x %g = %.3f\n",
            cfg$tipBr, 4 * cfg$tipBr))  # debug per-clade
cat(sprintf("  Per-clade total (eco): %g\n",
            4 * cfg$tipBr + 2 * cfg$tipBr + cfg$stemBrEco))
cat(sprintf("  Per-clade total (non-eco): %g\n",
            4 * cfg$tipBr + 2 * cfg$tipBr + cfg$stemBrClade))

# Construct AB-grouped reference tree (false eco clade)
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
  # AB grouping: ((A,B),(C,D))
  newick <- sprintf("((%s,%s):%g,(%s,%s):%g);",
                    cladeA, cladeB, rootBr,
                    cladeC, cladeD, rootBr)
  TreeTools::Preorder(ape::read.tree(text = newick))
}
abTree <- buildAB(cfg$tipBr, cfg$stemBrClade, cfg$stemBrEco, cfg$rootBr)

eco  <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)
cat("Edge ecology assignment:\n")
print(table(edgeEc))
# Sanity: eco edge total length
ecoBranchTotal <- sum(tree$edge.length[edgeEc == 1L])
cat(sprintf("Total eco-1 branch length: %.3f (%.0f%% of TL)\n",
            ecoBranchTotal, 100 * ecoBranchTotal / truthTL))

# --- Run N reps of simulate + MP search -------------------------------------
nReps <- 5L
seedBase <- 20260601L
nNeo <- cfg$nNeo; nTrans <- cfg$nTrans
nChar <- nNeo + nTrans
phi <- cfg$phi

results <- list()
for (rep in seq_len(nReps)) {
  set.seed(seedBase + rep)
  # z is nChar x kEco (kEco = 2: ref ecology 0 has z=0; ecology 1 mixed)
  zFull <- matrix(0L, nrow = nChar, ncol = 2L)
  # Eco column 2 (eco state = 1) gets z = 1 for a fraction (1 - pi0) of chars
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

  # Parsimony of true tree
  parsTrue <- as.integer(TreeSearch::TreeLength(tree, pdSim))
  # Parsimony of AB-grouped tree (false eco clade)
  parsAB   <- as.integer(TreeSearch::TreeLength(abTree, pdSim))

  # MP search to find a best-parsimony tree (random start, several restarts)
  set.seed(rep)
  randTree <- TreeTools::Preorder(
    ape::rtree(mkdBlind$nTip, tip.label = mkdBlind$taxon_names))
  mpTrees <- TreeSearch::MaximizeParsimony(pdSim, tree = randTree,
                                            verbosity = 0)
  if (inherits(mpTrees, "phylo")) mpTrees <- list(mpTrees)
  bestTree <- mpTrees[[1]]
  parsBest <- as.integer(TreeSearch::TreeLength(bestTree, pdSim))

  # Test AC vs AB splits on the best tree
  acTips <- c(paste0("A", 1:4), paste0("C", 1:4))
  abTips <- c(paste0("A", 1:4), paste0("B", 1:4))
  hasAC <- HasBipartSplits(list(bestTree), acTips)
  hasAB <- HasBipartSplits(list(bestTree), abTips)

  cat(sprintf(
    "rep %d: pars(true)=%d, pars(AB)=%d, pars(MP best)=%d  diff(AB-true)=%d  MP has AC=%s AB=%s\n",
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
cat(sprintf("Mean parsimony true: %.1f\n", mean(parsTrueV)))
cat(sprintf("Mean parsimony AB:   %.1f\n", mean(parsABV)))
cat(sprintf("Mean parsimony best: %.1f\n", mean(parsBestV)))
cat(sprintf("Mean (AB - true) step difference: %.1f\n",
            mean(parsABV - parsTrueV)))
cat(sprintf("Fraction MP best has AC split: %.2f\n", mean(acV)))
cat(sprintf("Fraction MP best has AB split: %.2f\n", mean(abV)))

# Predicted expSteps under auto-default (1.05 * parsimony)
cat(sprintf("\nPredicted auto expSteps = 1.05 * parsimony best ~ %.1f\n",
            1.05 * mean(parsBestV)))
cat(sprintf("Truth TL = %.3f; ratio parsimony/(TL*nChar) = %.3f\n",
            truthTL, mean(parsBestV) / (truthTL * nChar)))
