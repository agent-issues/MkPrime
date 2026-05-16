# sim3v4-helpers.R ------------------------------------------------------------
#
# Helpers for the redesigned Sim 3 ("v4") convergent-ecology study.
#
# v4 differs from v3 in three important ways:
#
# 1. Within-clade topology is PECTINATE, not balanced.
#       v3:  ((A1,A2),(A3,A4))           — balanced
#       v4:  (A4,(A3,(A2,A1)))           — pectinate; A1,A2 innermost pair
#
# 2. Ecology is assigned only to the INNERMOST pair of clades A and C:
#       eco1 tips: A1, A2, C1, C2  (4 tips total)
#       eco0 tips: A3, A4, B1..B4, C3, C4, D1..D4 (12 tips)
#    A and C are TRUE sisters (true topology ((A,C),(B,D))), so this is
#    WITHIN-clade partial ecology rather than the v3 across-clade design.
#
#    The spurious bipartition driven by convergence is
#        {A1, A2, C1, C2}  vs  rest
#    which is incompatible with the true clades A and C.
#
# 3. Characters are split into TWO arms with no mixing of types:
#       neomorphic arm       — asymmetric gain/loss CTMC (M2-NT)
#       transformational arm — symmetric k-state Mk
#    Both arms are ecology-encoded (z = 1 in ecology-1 for a fraction
#    1 - pi0 of characters); arm membership is recorded via
#    MkPrimeData(neomorphic = ...).
#
# Branch-length levers (per parameter set; see sim3v4-params.R):
#   tipBr        — within-clade tip & internal pectinate branches
#   stemBrEco    — branch subtending the innermost ecology-1 pair
#                  (the (A2,A1) cherry within clade A, and (C2,C1) within C)
#                  ** THE KEY LEVER for convergent-signal strength **
#   stemBrClade  — branch from clade stem to root (clade A, B, C, D each)
#   rootBr       — internal branch separating (A,C) from (B,D)


#' Build the v4 pectinate convergent tree.
#'
#' Within each 4-tip clade, topology is `(X4,(X3,(X2,X1)))`: X1 and X2 form
#' the innermost cherry under a stem of length `stemBrEco`; X3 is the next
#' branch out, X4 the outermost. All tip branches have length `tipBr`; the
#' single pectinate-internal branches between (X2,X1) and X3, and between
#' that node and X4, also have length `tipBr`. The whole clade hangs off a
#' stem of length `stemBrClade`. The two clade pairs (A,C) and (B,D) are
#' joined by branches of length `rootBr` either side of the root.
#'
#' True topology at the clade level: `((A, C), (B, D))`.
#'
#' @param tipBr      Positive scalar; tip and pectinate-internal branch
#'   length within each clade.
#' @param stemBrEco  Positive scalar; stem branch subtending the innermost
#'   (X2,X1) cherry (X = A, C) — the ecology-1 stem.
#' @param stemBrClade Positive scalar; stem subtending each 4-tip clade.
#' @param rootBr     Positive scalar; branch length on each side of the root.
#'
#' @return A `phylo` object in Preorder, 16 labelled tips.
#' @importFrom ape read.tree
#' @importFrom TreeTools Preorder
.BuildConvergentTreeV4 <- function(tipBr, stemBrEco, stemBrClade, rootBr) {
  stopifnot(tipBr > 0, stemBrEco > 0, stemBrClade > 0, rootBr > 0)
  # Pectinate within-clade Newick.  Innermost cherry (X2,X1) sits under a
  # stem of length stemBrEco when X is A or C (the ecology-1 clades);
  # otherwise the inner stem also uses tipBr so B and D remain symmetric
  # in branch length to A and C apart from the eco stem.
  cladeNewick <- function(prefix, innerStem) {
    sprintf(
      "(%s4:%g,(%s3:%g,(%s2:%g,%s1:%g):%g):%g):%g",
      prefix, tipBr,        # X4 tip
      prefix, tipBr,        # X3 tip
      prefix, tipBr,        # X2 tip
      prefix, tipBr,        # X1 tip
      innerStem,            # stem subtending (X2,X1)
      tipBr,                # internal pectinate branch between X3 and that node
      stemBrClade           # clade stem
    )
  }
  cladeA <- cladeNewick("A", stemBrEco)
  cladeC <- cladeNewick("C", stemBrEco)
  cladeB <- cladeNewick("B", tipBr)
  cladeD <- cladeNewick("D", tipBr)
  newick <- sprintf(
    "((%s,%s):%g,(%s,%s):%g);",
    cladeA, cladeC, rootBr, cladeB, cladeD, rootBr
  )
  tree <- ape::read.tree(text = newick)
  TreeTools::Preorder(tree)
}


#' Per-tip ecology assignment for the v4 tree.
#'
#' Tips A1, A2, C1, C2 are ecology 1 (innermost pair within clades A and C);
#' all other tips are ecology 0.
#'
#' @param tree Output of [.BuildConvergentTreeV4()].
#' @return Named integer vector of length `NTip(tree)`, values in {0L, 1L}.
.ConvergentEcologyV4 <- function(tree) {
  tips <- tree$tip.label
  eco <- ifelse(tips %in% c("A1", "A2", "C1", "C2"), 1L, 0L)
  names(eco) <- tips
  eco
}


#' Bipartitions of interest for scoring v4 chains.
#'
#' Returns a list of tip-label sets:
#'   - `trueClade_A` : `{A1,A2,A3,A4}` — the correctly-resolved clade A
#'   - `trueClade_C` : `{C1,C2,C3,C4}` — the correctly-resolved clade C
#'   - `trueSister_AC`: `{A1..A4, C1..C4}` — the true clade-level sister pair
#'   - `falseSister`  : `{A1,A2,C1,C2}` — the spurious ecology-driven grouping
#'
#' @return Named list of character vectors.
.SimBipartitionsV4 <- function() {
  list(
    trueClade_A    = paste0("A", 1:4),
    trueClade_C    = paste0("C", 1:4),
    trueSister_AC  = c(paste0("A", 1:4), paste0("C", 1:4)),
    falseSister    = c("A1", "A2", "C1", "C2")
  )
}


#' Build per-rep v4 simulation data.
#'
#' Simulates a single replicate under the parameter set `params`, with the
#' neomorphic and transformational character arms generated separately and
#' concatenated. Ecology-encoding (z = 1 on ecology-1 edges) is applied to a
#' fraction `1 - pi0` of characters in EACH arm, drawn at random;
#' remaining characters are baseline (z = 0).
#'
#' @param tree   v4 tree from [.BuildConvergentTreeV4()].
#' @param eco    Tip ecology from [.ConvergentEcologyV4()].
#' @param edgeEc Edge ecology from [.AssignEdgeEcology()] (sim3-simulate.R).
#' @param params Named list with elements:
#'   `nNeo`, `nTrans`, `phi`, `pi0`, `theta`, `baseRate` (default 1),
#'   `rateLoss` (default 1), `kStates` (default 2L, used for trans arm).
#' @param simulator Function with the signature of
#'   [.SimulateMkPrimeEcology()] (sim3-simulate.R). Supplied by the caller
#'   so this file does not need to source the simulator itself.
#'
#' @return List with elements
#'   `pdSim`     : phyDat object,
#'   `mkdBlind`  : MkPrimeData with `neomorphic = 1:nNeo`,
#'   `mkdAware`  : MkPrimeData with `neomorphic` and `ecology` set,
#'   `zNeo`, `zTrans` : the z matrices actually drawn,
#'   `truthTL`   : total tree length.
.BuildV4Data <- function(tree, eco, edgeEc, params, simulator) {
  stopifnot(is.list(params))
  nNeo     <- as.integer(params$nNeo)
  nTrans   <- as.integer(params$nTrans)
  phi      <- params$phi
  pi0      <- params$pi0
  theta    <- params$theta
  baseRate <- if (is.null(params$baseRate)) 1 else params$baseRate
  rateLoss <- if (is.null(params$rateLoss)) 1 else params$rateLoss
  kStates  <- if (is.null(params$kStates)) 2L else as.integer(params$kStates)

  # Draw z per character per arm. theta gives the fraction of
  # ecology-affected characters that are "encouraged" (z = 1) vs
  # "discouraged" (z = 2); theta = 1 means all encouraged.
  drawZ <- function(nChar) {
    z <- matrix(0L, nrow = nChar, ncol = 2L)  # kEco = 2 (eco0, eco1)
    isEco <- stats::runif(nChar) > pi0
    if (any(isEco)) {
      enc <- stats::runif(sum(isEco)) < theta
      z[isEco, 2L] <- ifelse(enc, 1L, 2L)
    }
    z
  }

  zNeo   <- drawZ(nNeo)
  zTrans <- drawZ(nTrans)

  datNeo <- simulator(tree, edgeEc, zNeo, phi = phi,
                      type = "neomorphic", baseRate = baseRate,
                      rateLoss = rateLoss, kStates = 2L,
                      normalize = TRUE, pi0 = pi0, theta = theta,
                      refEcology = 0L)
  datTrans <- simulator(tree, edgeEc, zTrans, phi = phi,
                        type = "transformational", baseRate = baseRate,
                        kStates = kStates,
                        normalize = TRUE, pi0 = pi0, theta = theta,
                        refEcology = 0L)
  datSim <- cbind(datNeo, datTrans)

  pdSim <- MkPrime::MatrixToPhyDat(datSim)
  neoIdx <- seq_len(nNeo)
  mkdBlind <- MkPrime::MkPrimeData(pdSim, neomorphic = neoIdx)
  ecoTipVec <- eco[mkdBlind$taxon_names]
  mkdAware <- MkPrime::MkPrimeData(pdSim, neomorphic = neoIdx,
                                   ecology = ecoTipVec)

  list(pdSim    = pdSim,
       mkdBlind = mkdBlind,
       mkdAware = mkdAware,
       zNeo     = zNeo,
       zTrans   = zTrans,
       truthTL  = sum(tree$edge.length))
}
