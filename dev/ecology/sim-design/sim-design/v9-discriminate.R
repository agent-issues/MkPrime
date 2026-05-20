# v9-discriminate.R — neomorphic-directional pre-screen.
#
# User intuition (2026-05-20):
#   "Move to a neomorphic world, set simulations with just gain/losses,
#    with a symmetric slab so ecology A is gaining some characters at a
#    higher rate, and that ecology B is losing these same characters at
#    a higher rate. That's the cleanest framing of a convergent signal."
#
# v7/v8 used type="transformational" (symmetric mult = c(1, phi, 1/phi));
# phi only scales the OVERALL rate, never the direction. The existing
# `.SimulateMkPrimeEcology()` ALREADY contains an asymmetric neomorphic
# branch — `z = 1` gives (rate01 * phi, rate10 / phi) — that pushes the
# Markov chain toward state 1; `z = 2` reverses direction. The asymmetric
# mode was never exercised in sim3.
#
# This driver probes two mechanisms by setting the z matrix accordingly:
#
#   Mechanism M1 ("divergent"): z[c, A] = 1 AND z[c, B] = 2  for the
#     same character c. Eco-A clades drift toward state 1 (gains
#     elevated, losses suppressed); eco-B clades drift toward state 0
#     (losses elevated, gains suppressed). A and B become DISSIMILAR.
#
#   Mechanism M2 ("parallel"):  z[c, A] = 1 AND z[c, B] = 1  for the
#     same character c. Both eco-A and eco-B clades drift toward state 1.
#     A and B share derived states; potential to trick MP into pulling
#     (A,B) sister.
#
# Tree: 16-tip canonical ((A,C),(B,D)) with eco-1 on full A clade and
# eco-2 on full B clade (analogous to v6 but using two ecology labels so
# the z columns can differ between A and B). C and D stay baseline.
#
# We score:
#   * pTrue  — MP recovers (AC|BD) (true sister).
#   * pAB    — MP groups (A,B) as a clade (the convergent trap).
#   * delta  — pars(true) - pars(MP best).
#
# Cap: 5 parameter rounds across the two mechanisms; STOP if any setting
# achieves pAB >= 50 %.

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


# Build the canonical convergent tree but expose tip/stem/root knobs ----
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


# Two-ecology tip assignment: eco-1 on full A clade, eco-2 on full B clade.
.TwoEcologyV9 <- function(tree) {
  tips <- tree$tip.label
  eco  <- integer(length(tips))
  names(eco) <- tips
  eco[grepl("^A", tips)] <- 1L
  eco[grepl("^B", tips)] <- 2L
  eco
}


# Variant catalogue ---------------------------------------------------------
# Each variant pairs a mechanism (M1 divergent or M2 parallel) with a
# parameter set that respects the TL <= 2 ceiling.
variants <- list(

  # ---- M1 divergent: z[c,A]=1, z[c,B]=2 ---------------------------------
  v9a_M1_realistic = list(
    name = "v9a_M1_realistic",
    mechanism = "divergent",
    tipBr = 0.04, stemBrClade = 0.045, stemBrEco = 0.18, rootBr = 0.035,
    nChar = 200L, phi = 6, p_eco = 0.5,
    rateLoss = 1.0, baseRate = 1.0
  ),

  v9b_M1_pushed = list(
    name = "v9b_M1_pushed",
    mechanism = "divergent",
    tipBr = 0.04, stemBrClade = 0.045, stemBrEco = 0.20, rootBr = 0.035,
    nChar = 200L, phi = 8, p_eco = 0.6,
    rateLoss = 1.0, baseRate = 1.0
  ),

  # ---- M2 parallel: z[c,A] = z[c,B] = 1 ---------------------------------
  v9c_M2_realistic = list(
    name = "v9c_M2_realistic",
    mechanism = "parallel",
    tipBr = 0.04, stemBrClade = 0.045, stemBrEco = 0.18, rootBr = 0.035,
    nChar = 200L, phi = 6, p_eco = 0.5,
    rateLoss = 1.0, baseRate = 1.0
  ),

  v9d_M2_pushed = list(
    name = "v9d_M2_pushed",
    mechanism = "parallel",
    tipBr = 0.04, stemBrClade = 0.045, stemBrEco = 0.20, rootBr = 0.035,
    nChar = 200L, phi = 8, p_eco = 0.6,
    rateLoss = 1.0, baseRate = 1.0
  ),

  # ---- M2 maxed: as-aggressive-as-feasible under TL <= 2 ----------------
  v9e_M2_maxout = list(
    name = "v9e_M2_maxout",
    mechanism = "parallel",
    tipBr = 0.04, stemBrClade = 0.04, stemBrEco = 0.22, rootBr = 0.03,
    nChar = 150L, phi = 10, p_eco = 0.7,
    rateLoss = 1.0, baseRate = 1.0
  )
)

requestedVariant <- Sys.getenv("V9_VARIANT", "")
nReps    <- as.integer(Sys.getenv("V9_NREPS", "8"))
seedBase <- 20260520L


# ---- One rep ---------------------------------------------------------------
# Build z matrix per mechanism. nChar x 3 (baseline + eco-A + eco-B).
makeZ <- function(nChar, p_eco, mechanism, seed) {
  set.seed(seed)
  z <- matrix(0L, nrow = nChar, ncol = 3L)
  ecoActive <- stats::runif(nChar) < p_eco
  if (mechanism == "divergent") {
    # eco-A pushes toward 1; eco-B pushes toward 0.
    z[ecoActive, 2L] <- 1L
    z[ecoActive, 3L] <- 2L
  } else if (mechanism == "parallel") {
    # Both eco-A and eco-B push toward 1.
    z[ecoActive, 2L] <- 1L
    z[ecoActive, 3L] <- 1L
  } else stop("unknown mechanism: ", mechanism)
  z
}


runOneRep <- function(cfg, rep, tree, edgeEc) {
  zMat <- makeZ(cfg$nChar, cfg$p_eco, cfg$mechanism, seedBase + rep)
  set.seed(seedBase + rep + 1000L)
  cFull <- rep("neomorphic", cfg$nChar)
  datSim <- .SimulateMkPrimeEcology(
    tree, edgeEc, zMat, phi = cfg$phi,
    type = cFull, baseRate = cfg$baseRate,
    normalize = FALSE,
    rateLoss = cfg$rateLoss
  )
  # Constant-char filter — neomorphic with strong directional push may
  # generate columns with no variation; MP can't use them. Filter them
  # so the parsimony comparison is meaningful.
  varCols <- apply(datSim, 2L, function(col) length(unique(col)) > 1L)
  datSim <- datSim[, varCols, drop = FALSE]
  nVar <- ncol(datSim)
  pdSim <- MatrixToPhyDat(datSim)

  trueSister <- c(paste0("A", 1:4), paste0("C", 1:4))
  falseAB    <- c(paste0("A", 1:4), paste0("B", 1:4))

  parsTrue <- as.integer(TreeSearch::TreeLength(tree, pdSim))
  set.seed(rep)
  randTree <- TreeTools::Preorder(
    ape::rtree(length(tree$tip.label), tip.label = tree$tip.label))
  mpTrees <- TreeSearch::MaximizeParsimony(
    pdSim, tree = randTree, concavity = Inf,
    verbosity = 0, maxReplicates = 24L, targetHits = 10L)
  if (inherits(mpTrees, "phylo")) mpTrees <- list(mpTrees)
  class(mpTrees) <- "multiPhylo"
  parsBest <- as.integer(TreeSearch::TreeLength(mpTrees[[1]], pdSim))
  list(
    rep = rep, nVar = nVar,
    parsTrue = parsTrue, parsBest = parsBest,
    hasTrueSister = any(HasBipartSplits(mpTrees, trueSister)),
    hasFalseAB    = any(HasBipartSplits(mpTrees, falseAB))
  )
}


runVariant <- function(cfg) {
  cat(sprintf("\n=== %s  (mechanism=%s) ===\n", cfg$name, cfg$mechanism))
  cat(sprintf("  phi=%g  p_eco=%.2f  rateLoss=%g  nChar=%d\n",
              cfg$phi, cfg$p_eco, cfg$rateLoss, cfg$nChar))
  tree <- buildTreeV9(cfg$tipBr, cfg$stemBrClade, cfg$stemBrEco, cfg$rootBr)
  TL <- sum(tree$edge.length)
  cat(sprintf("  TL=%.3f\n", TL))
  if (TL > 2.0) cat(sprintf("  WARNING: TL %.3f > 2.0 ceiling\n", TL))
  eco <- .TwoEcologyV9(tree)
  edgeEc <- .AssignEdgeEcology(tree, eco)
  cat("  edge-ecology table:\n"); print(table(edgeEc))

  results <- lapply(seq_len(nReps),
                    function(r) runOneRep(cfg, r, tree, edgeEc))

  hitMatrix <- do.call(rbind, lapply(results, function(r)
    c(hasTrueSister = r$hasTrueSister,
      hasFalseAB    = r$hasFalseAB)))
  deltas <- vapply(results, function(r) r$parsTrue - r$parsBest, integer(1))
  nVars  <- vapply(results, function(r) r$nVar,    integer(1))
  pTrue  <- mean(hitMatrix[, "hasTrueSister"])
  pAB    <- mean(hitMatrix[, "hasFalseAB"])

  cat(sprintf("  variable chars per rep: %s\n",
              paste(nVars, collapse = ", ")))
  cat(sprintf("  pars(true) mean: %.1f  pars(MP) mean: %.1f\n",
              mean(vapply(results, function(r) r$parsTrue, integer(1))),
              mean(vapply(results, function(r) r$parsBest, integer(1)))))
  cat(sprintf("  delta(true - MP) per rep: %s\n",
              paste(deltas, collapse = ", ")))
  cat(sprintf("  MP contains TRUE sister: %.0f%%\n", 100 * pTrue))
  cat(sprintf("  MP contains FALSE (A,B): %.0f%%\n", 100 * pAB))

  list(name = cfg$name, mechanism = cfg$mechanism, TL = TL,
       pTrue = pTrue, pAB = pAB,
       deltas = deltas, nVars = nVars, cfg = cfg)
}

summary <- list()
varNames <- if (nzchar(requestedVariant)) requestedVariant else names(variants)
for (vn in varNames) {
  if (!vn %in% names(variants)) {
    cat("Unknown variant: ", vn, "\n"); next
  }
  res <- runVariant(variants[[vn]])
  summary[[vn]] <- res
  if (res$pAB >= 0.5) {
    cat(sprintf("\n*** TRAP FOUND at %s — pAB = %.0f%% ***\n",
                res$name, 100 * res$pAB))
    break
  }
}

cat("\n\n=== Variant summary ===\n")
for (vn in names(summary)) {
  s <- summary[[vn]]
  cat(sprintf("%-26s  mech=%-9s  TL=%.2f  pTrue=%.0f%%  pAB=%.0f%%\n",
              s$name, s$mechanism, s$TL,
              100 * s$pTrue, 100 * s$pAB))
}

saveRDS(summary,
        file.path("dev", "sim-design", "v9-discriminate-results.rds"))
cat("\nWrote dev/sim-design/v9-discriminate-results.rds\n")
