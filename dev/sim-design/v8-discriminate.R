# v8-discriminate.R — multi-ecology pre-screen for the parsimony trap.
#
# Goal: do 3 independent ecologies (each driving a different 4-tip
# convergent clade) at realistic phi/pi0 cause MaximizeParsimony to
# return ANY of the false eco clades?
#
# Design context (advisor-checked, 2026-05-20)
# ---------------------------------------------
# v6 (phi=8, pi0=0.45, single binary eco, full A+B clades, TL=1.16)
#   → MP recovered truth on 3/3 reps; deltaSteps=0; no trap.
# v7 (4 variants, phi up to 15, pi0 down to 0.15, single eco, TL ≤ 1.8)
#   → best was v7d, deltaSteps -0.25 (coin flip), past the optimum.
# multirep-v3 (TL=13.5, saturated, phi=4, single 8-tip eco)
#   → 7/8 reps MP hit false clade (saturation regime, unrealistic).
#
# Hypothesis to falsify: at realistic TL (~1.5), three independent
# ecologies — each weaker per-trap than v6/v7 but ORed together — fool
# MP at least 50 % of the time.
#
# The simulator uses ONE shared phi across ecologies; edgeEcology
# assigns each edge to a single eco label. So the multi-eco design
# encodes each "convergent clade pair" as a distinct ecology state
# (0=baseline, 1=arboreal, 2=semiaquatic, 3=fossorial); each eco gets
# its own z[, e] column. With normalize=TRUE the per-eco shrinkage
# `gammaE = pi0 + (1-pi0)*(theta*phi + (1-theta)/phi)` applies; the
# realised rate factor is then phi/gammaE.
#
# Variants probe normalize on/off and increasing aggression. Cap: 5
# variants; STOP after the first that achieves the target (>= 50 % of
# reps with MP containing any false eco clade) OR after exhausting the
# list.

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

# --- Topology builder reused from v6 (eco-1 vs eco-0 stems) -----------
# For multi-eco we lengthen the stems for *all* three eco-affected clade
# pairs, but only on the side where the false convergent group sits.
buildTreeV8 <- function(tipBr, stemBrClade, stemBrEco, rootBr) {
  cladeNw <- function(prefix, tb, sb) {
    sprintf(
      "((%s1:%g,%s2:%g):%g,(%s3:%g,%s4:%g):%g):%g",
      prefix, tb, prefix, tb, tb,
      prefix, tb, prefix, tb, tb, sb
    )
  }
  cladeA <- cladeNw("A", tipBr, stemBrEco)
  cladeB <- cladeNw("B", tipBr, stemBrEco)
  cladeC <- cladeNw("C", tipBr, stemBrEco)
  cladeD <- cladeNw("D", tipBr, stemBrEco)
  newick <- sprintf("((%s,%s):%g,(%s,%s):%g);",
                    cladeA, cladeC, rootBr,
                    cladeB, cladeD, rootBr)
  TreeTools::Preorder(ape::read.tree(text = newick))
}

# --- Multi-ecology tip-state assignment -----------------------------------
# Three ecologies; each defined over a 4-tip group spanning two clades.
#   eco-1 ("arboreal")    : inner cherries of A and B = {A1,A2,B1,B2}
#   eco-2 ("semiaquatic") : outer cherries of C and D = {C3,C4,D3,D4}
#   eco-3 ("fossorial")   : outer cherry of A + inner cherry of D
#                                              = {A3,A4,D1,D2}
# Each ecology drives a different false bipartition; truth (AC) sister
# crosses none of them.
.MultiEcologyV8 <- function(tree, twoEcoFullAB = FALSE) {
  tips <- tree$tip.label
  eco  <- integer(length(tips))
  names(eco) <- tips
  if (isTRUE(twoEcoFullAB)) {
    # v6-style strong eco on full A+B clades, plus weaker eco-2 on
    # outer cherries of C, D.
    eco[grepl("^[AB]", tips)] <- 1L
    eco[c("C3","C4","D3","D4")] <- 2L
  } else {
    eco[c("A1","A2","B1","B2")] <- 1L
    eco[c("C3","C4","D3","D4")] <- 2L
    eco[c("A3","A4","D1","D2")] <- 3L
  }
  eco
}

# --- False bipartitions induced by each ecology ---------------------------
.MultiBipartsV8 <- function(twoEcoFullAB = FALSE) {
  if (isTRUE(twoEcoFullAB)) {
    return(list(
      trueSister = c(paste0("A", 1:4), paste0("C", 1:4)),
      falseEco1  = c(paste0("A", 1:4), paste0("B", 1:4)), # 8-tip A+B
      falseEco2  = c("C3","C4","D3","D4"),
      falseEco3  = c("A3","A4","D1","D2")  # unused; left for unified API
    ))
  }
  list(
    trueSister  = c(paste0("A", 1:4), paste0("C", 1:4)),
    falseEco1   = c("A1","A2","B1","B2"),
    falseEco2   = c("C3","C4","D3","D4"),
    falseEco3   = c("A3","A4","D1","D2")
  )
}

# Variant catalogue ----------------------------------------------------
variants <- list(
  v8a = list(  # User's proposal — realistic phi, normalize=TRUE
    name = "v8a_realistic_norm",
    tipBr = 0.04, stemBrClade = 0.045, stemBrEco = 0.10, rootBr = 0.04,
    nNeo = 0L, nTrans = 200L,
    phi = 4, pi0 = 0.5, theta = 1.0,
    normalize = TRUE, baseRate = 1.0
  ),
  v8b = list(  # Same but normalize OFF → mechanism ceiling at realistic phi
    name = "v8b_realistic_unnorm",
    tipBr = 0.04, stemBrClade = 0.045, stemBrEco = 0.10, rootBr = 0.04,
    nNeo = 0L, nTrans = 200L,
    phi = 4, pi0 = 0.5, theta = 1.0,
    normalize = FALSE, baseRate = 1.0
  ),
  v8c = list(  # Lengthen eco stems, drop pi0, normalize OFF
    name = "v8c_aggressive_unnorm",
    tipBr = 0.04, stemBrClade = 0.045, stemBrEco = 0.15, rootBr = 0.035,
    nNeo = 0L, nTrans = 200L,
    phi = 6, pi0 = 0.3, theta = 1.0,
    normalize = FALSE, baseRate = 1.0
  ),
  v8d = list(  # Push hard — phi=8, very long eco stems, normalize OFF
    name = "v8d_pushed_unnorm",
    tipBr = 0.04, stemBrClade = 0.04, stemBrEco = 0.20, rootBr = 0.03,
    nNeo = 0L, nTrans = 200L,
    phi = 8, pi0 = 0.25, theta = 1.0,
    normalize = FALSE, baseRate = 1.0
  ),
  v8e = list(  # Most aggressive feasible at TL <= 2: phi=10, pi0=0.2
    name = "v8e_maxout_unnorm",
    tipBr = 0.04, stemBrClade = 0.04, stemBrEco = 0.22, rootBr = 0.025,
    nNeo = 0L, nTrans = 150L,
    phi = 10, pi0 = 0.2, theta = 1.0,
    normalize = FALSE, baseRate = 1.0
  ),
  # Reuses v6's 8-tip full-clade eco (A+B), unnormalised, plus second
  # smaller eco-2 over a 4-tip subset {C3,C4,D3,D4}. Probes whether the
  # v6 setup (only fails by 9 steps when normalised) becomes a trap when
  # unnormalised + supplemented.
  v8f = list(
    name = "v8f_v6plus_unnorm",
    tipBr = 0.04, stemBrClade = 0.045, stemBrEco = 0.18, rootBr = 0.035,
    nNeo = 0L, nTrans = 200L,
    phi = 6, pi0 = 0.35, theta = 1.0,
    normalize = FALSE, baseRate = 1.0,
    twoEcoFullAB = TRUE
  )
)

requestedVariant <- Sys.getenv("V8_VARIANT", "")
nReps <- as.integer(Sys.getenv("V8_NREPS", "8"))
seedBase <- 20260820L

# Helper: simulate one rep and return MP analysis
runOneRep <- function(cfg, rep, tree, edgeEc) {
  set.seed(seedBase + rep)
  nChar <- cfg$nNeo + cfg$nTrans
  kEco  <- max(edgeEc) + 1L  # number of eco labels including baseline
  zFull <- matrix(0L, nrow = nChar, ncol = kEco)
  # For each non-baseline ecology, independently draw which chars are
  # eco-active and force z=1L (matching theta=1.0).
  for (e in 2:kEco) {
    ecoActive <- stats::runif(nChar) > cfg$pi0
    zFull[ecoActive, e] <- 1L
  }
  cFull <- c(rep("neomorphic", cfg$nNeo),
             rep("transformational", cfg$nTrans))
  datSim <- .SimulateMkPrimeEcology(
    tree, edgeEc, zFull, phi = cfg$phi,
    type = cFull, baseRate = cfg$baseRate,
    normalize = cfg$normalize,
    pi0 = cfg$pi0, theta = cfg$theta,
    refEcology = 0L, rateLoss = 1
  )
  pdSim <- MatrixToPhyDat(datSim)
  biparts <- .MultiBipartsV8(isTRUE(cfg$twoEcoFullAB))
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
    rep = rep,
    parsTrue = parsTrue, parsBest = parsBest,
    hasTrueSister = any(HasBipartSplits(mpTrees, biparts$trueSister)),
    hasFalseEco1  = any(HasBipartSplits(mpTrees, biparts$falseEco1)),
    hasFalseEco2  = any(HasBipartSplits(mpTrees, biparts$falseEco2)),
    hasFalseEco3  = any(HasBipartSplits(mpTrees, biparts$falseEco3))
  )
}

# Iterate ----------------------------------------------------------------
runVariant <- function(cfg) {
  cat(sprintf("\n=== %s ===\n", cfg$name))
  cat(sprintf("  phi=%g  pi0=%.2f  normalize=%s  nChar=%d\n",
              cfg$phi, cfg$pi0, cfg$normalize,
              cfg$nNeo + cfg$nTrans))
  if (isTRUE(cfg$normalize)) {
    gammaE <- cfg$pi0 + (1 - cfg$pi0) *
      (cfg$theta * cfg$phi + (1 - cfg$theta) / cfg$phi)
    eff <- cfg$phi / gammaE
    cat(sprintf("  gammaE=%.3f  effective rate factor=%.3f\n", gammaE, eff))
  } else {
    cat(sprintf("  unnormalised — full rate factor = phi = %g\n", cfg$phi))
  }
  tree <- buildTreeV8(cfg$tipBr, cfg$stemBrClade, cfg$stemBrEco, cfg$rootBr)
  TL <- sum(tree$edge.length)
  cat(sprintf("  TL=%.3f\n", TL))
  if (TL > 2.0) cat(sprintf("  WARNING: TL %.3f > 2.0 ceiling\n", TL))
  eco <- .MultiEcologyV8(tree, isTRUE(cfg$twoEcoFullAB))
  edgeEc <- .AssignEdgeEcology(tree, eco)
  cat("  edge-ecology table:\n")
  print(table(edgeEc))

  results <- lapply(seq_len(nReps),
                    function(r) runOneRep(cfg, r, tree, edgeEc))

  hitMatrix <- do.call(rbind, lapply(results, function(r)
    c(hasTrueSister = r$hasTrueSister,
      hasFalseEco1  = r$hasFalseEco1,
      hasFalseEco2  = r$hasFalseEco2,
      hasFalseEco3  = r$hasFalseEco3)))
  anyFalse <- (hitMatrix[, "hasFalseEco1"] |
               hitMatrix[, "hasFalseEco2"] |
               hitMatrix[, "hasFalseEco3"])
  deltas <- vapply(results, function(r) r$parsTrue - r$parsBest, integer(1))
  pTrue  <- mean(hitMatrix[, "hasTrueSister"])
  pAnyF  <- mean(anyFalse)
  pE1    <- mean(hitMatrix[, "hasFalseEco1"])
  pE2    <- mean(hitMatrix[, "hasFalseEco2"])
  pE3    <- mean(hitMatrix[, "hasFalseEco3"])

  cat(sprintf("  pars(true) mean: %.1f  pars(MP) mean: %.1f\n",
              mean(vapply(results, function(r) r$parsTrue, integer(1))),
              mean(vapply(results, function(r) r$parsBest, integer(1)))))
  cat(sprintf("  delta(true - MP) per rep: %s\n",
              paste(deltas, collapse = ", ")))
  cat(sprintf("  MP contains TRUE sister:  %.0f%%\n", 100 * pTrue))
  cat(sprintf("  MP contains ANY falseEco: %.0f%%   (e1=%.0f%% e2=%.0f%% e3=%.0f%%)\n",
              100 * pAnyF, 100 * pE1, 100 * pE2, 100 * pE3))

  list(name = cfg$name, TL = TL,
       pTrue = pTrue, pAnyF = pAnyF,
       pE1 = pE1, pE2 = pE2, pE3 = pE3,
       deltas = deltas, cfg = cfg)
}

summary <- list()
varNames <- if (nzchar(requestedVariant)) requestedVariant else names(variants)
for (vn in varNames) {
  if (!vn %in% names(variants)) {
    cat("Unknown variant: ", vn, "\n"); next
  }
  res <- runVariant(variants[[vn]])
  summary[[vn]] <- res
  if (res$pAnyF >= 0.5) {
    cat(sprintf("\n*** TRAP FOUND at %s — pAnyFalse = %.0f%% ***\n",
                res$name, 100 * res$pAnyF))
    break
  }
}

cat("\n\n=== Variant summary ===\n")
for (vn in names(summary)) {
  s <- summary[[vn]]
  cat(sprintf("%-26s  TL=%.2f  trueSister=%.0f%%  anyFalse=%.0f%%\n",
              s$name, s$TL,
              100 * s$pTrue, 100 * s$pAnyF))
}

# Persist for inspection
saveRDS(summary,
        file.path("dev", "sim-design", "v8-discriminate-results.rds"))
cat("\nWrote dev/sim-design/v8-discriminate-results.rds\n")
