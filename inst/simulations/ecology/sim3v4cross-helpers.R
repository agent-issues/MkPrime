# sim3v4cross-helpers.R --------------------------------------------------------
#
# Helpers for the "v4-cross" variant of the Sim 3 convergent-ecology study.
#
# v4-cross differs from v4 in the placement of the ecology-1 tips:
#
#   v4       : eco-1 = {A1, A2, C1, C2}    — inner cherries of SISTER clades
#   v4-cross : eco-1 = {A1, A2, B1, B2}    — inner cherries of NON-SISTER clades
#
# True topology is unchanged: ((A, C), (B, D)).  In v4 the eco signal aligns
# with the true (A,C) inter-clade sister relationship and the confound only
# operates within clades.  In v4-cross the eco signal pulls A1,A2 toward
# B1,B2 ACROSS the true (A,C)+(B,D) divide, so the false bipartition
# {A1,A2,B1,B2} directly contradicts the true clade-level topology — the
# blind chain must choose between (A,C)+(B,D) and the spurious eco grouping.
#
# Tree-building and data-building logic is otherwise identical to v4; we
# reuse `.BuildV4Data()` from sim3v4-helpers.R.  Only the inner-stem
# placement (which clades carry `stemBrEco` vs `tipBr`) and the tip ecology
# vector differ.


#' Build the v4-cross pectinate convergent tree.
#'
#' Within-clade topology, branch lengths, and clade-level topology
#' `((A, C), (B, D))` are identical to [.BuildConvergentTreeV4()]; the only
#' change is that the eco-stem branch length `stemBrEco` is applied to the
#' innermost cherry of clades A and B (not A and C).  C and D inner stems
#' use `tipBr` for length symmetry with the v4 design.
#'
#' @param tipBr,stemBrEco,stemBrClade,rootBr See [.BuildConvergentTreeV4()].
#' @return A `phylo` object in Preorder, 16 labelled tips.
#' @importFrom ape read.tree
#' @importFrom TreeTools Preorder
.BuildConvergentTreeV4Cross <- function(tipBr, stemBrEco, stemBrClade,
                                         rootBr) {
  stopifnot(tipBr > 0, stemBrEco > 0, stemBrClade > 0, rootBr > 0)
  cladeNewick <- function(prefix, innerStem) {
    sprintf(
      "(%s4:%g,(%s3:%g,(%s2:%g,%s1:%g):%g):%g):%g",
      prefix, tipBr,
      prefix, tipBr,
      prefix, tipBr,
      prefix, tipBr,
      innerStem,
      tipBr,
      stemBrClade
    )
  }
  cladeA <- cladeNewick("A", stemBrEco)   # eco-1 inner cherry
  cladeB <- cladeNewick("B", stemBrEco)   # eco-1 inner cherry
  cladeC <- cladeNewick("C", tipBr)        # no eco
  cladeD <- cladeNewick("D", tipBr)        # no eco
  # True topology at the clade level: ((A, C), (B, D)).
  newick <- sprintf(
    "((%s,%s):%g,(%s,%s):%g);",
    cladeA, cladeC, rootBr, cladeB, cladeD, rootBr
  )
  tree <- ape::read.tree(text = newick)
  TreeTools::Preorder(tree)
}


#' Per-tip ecology assignment for the v4-cross tree.
#'
#' Tips A1, A2, B1, B2 are ecology 1 (innermost pair within clades A and B);
#' all other tips are ecology 0.
#'
#' @param tree Output of [.BuildConvergentTreeV4Cross()].
#' @return Named integer vector of length `NTip(tree)`, values in {0L, 1L}.
.ConvergentEcologyV4Cross <- function(tree) {
  tips <- tree$tip.label
  eco <- ifelse(tips %in% c("A1", "A2", "B1", "B2"), 1L, 0L)
  names(eco) <- tips
  eco
}


#' Bipartitions of interest for scoring v4-cross chains.
#'
#' Returns a list of tip-label sets:
#'   - `trueClade_A` ... `trueClade_D`: each 4-tip clade
#'   - `trueSister_AC` : `{A1..A4, C1..C4}` — true clade-level sister pair
#'   - `trueSister_BD` : `{B1..B4, D1..D4}` — true clade-level sister pair
#'   - `falseInner`    : `{A1,A2,B1,B2}` — eco-driven inner-pair grouping
#'                       (requires breaking clade A and B monophyly)
#'   - `falseClade_AB` : `{A1..A4, B1..B4}` — eco-driven clade-level grouping
#'                       (preserves clade monophyly but breaks (A,C),(B,D))
#'
#' @return Named list of character vectors.
.SimBipartitionsV4Cross <- function() {
  list(
    trueClade_A    = paste0("A", 1:4),
    trueClade_B    = paste0("B", 1:4),
    trueClade_C    = paste0("C", 1:4),
    trueClade_D    = paste0("D", 1:4),
    trueSister_AC  = c(paste0("A", 1:4), paste0("C", 1:4)),
    trueSister_BD  = c(paste0("B", 1:4), paste0("D", 1:4)),
    falseInner     = c("A1", "A2", "B1", "B2"),
    falseClade_AB  = c(paste0("A", 1:4), paste0("B", 1:4))
  )
}
