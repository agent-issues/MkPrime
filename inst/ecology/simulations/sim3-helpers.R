# sim3-helpers.R --------------------------------------------------------------
#
# Helpers for the convergent-ecology validation study (Sim 3 in the
# ecology-aware MkPrime simulation plan). Builds the canonical convergent-
# clades tree and exposes utilities for downstream simulation / scoring
# scripts.
#
# The tree has four 4-tip clades A, B, C, D arranged so that the TRUE
# topology is ((A, C), (B, D)). A and B are assigned ecology state 1; C
# and D are assigned ecology state 0. An ecology-blind Mk model that
# inflates parallel state changes in A and B is expected to incorrectly
# group A with B.
#
# Branch lengths:
#   - terminal branches inside each clade: tipBranch (default 0.5)
#   - stem branches subtending each clade:  stemBranch (default 0.05)
#     short stems so the convergent signal can plausibly dominate
#   - root branch separating ((A,C),(B,D)): rootBranch (default 0.10)

#' Build the canonical convergent-clades simulation tree.
#'
#' @param tipBranch Numeric. Length of within-clade terminal branches.
#' @param stemBranch Numeric. Length of the stem branches subtending each
#'   4-tip clade.
#' @param rootBranch Numeric. Length of the internal branch separating the
#'   two cherries-of-clades.
#'
#' @return A `phylo` object in Preorder with 16 labelled tips. Tips
#'   `A1..A4` belong to clade A, etc. Branch lengths follow the rules
#'   described in the file header.
#'
#' @importFrom ape read.tree
#' @importFrom TreeTools Preorder
.BuildConvergentTree <- function(tipBranch = 0.5,
                                 stemBranch = 0.05,
                                 rootBranch = 0.10) {
  stopifnot(tipBranch > 0, stemBranch > 0, rootBranch > 0)
  # Newick built by hand to give explicit per-branch control.
  # Within-clade: ((A1:tb,A2:tb):tb,(A3:tb,A4:tb):tb):sb
  cladeNewick <- function(prefix, tb, sb) {
    sprintf(
      "((%s1:%g,%s2:%g):%g,(%s3:%g,%s4:%g):%g):%g",
      prefix, tb, prefix, tb, tb,
      prefix, tb, prefix, tb, tb, sb
    )
  }
  cladeA <- cladeNewick("A", tipBranch, stemBranch)
  cladeC <- cladeNewick("C", tipBranch, stemBranch)
  cladeB <- cladeNewick("B", tipBranch, stemBranch)
  cladeD <- cladeNewick("D", tipBranch, stemBranch)
  # True topology: ((A,C), (B,D))
  newick <- sprintf(
    "((%s,%s):%g,(%s,%s):%g);",
    cladeA, cladeC, rootBranch, cladeB, cladeD, rootBranch
  )
  tree <- ape::read.tree(text = newick)
  TreeTools::Preorder(tree)
}


#' Per-tip ecology assignments for the convergent tree.
#'
#' Tips A* and B* are ecology 1; tips C* and D* are ecology 0.
#'
#' @param tree Output of [.BuildConvergentTree()].
#' @return Integer vector of length `NTip(tree)`, named by tip label,
#'   with values in `{0L, 1L}`.
.ConvergentEcology <- function(tree) {
  tips <- tree$tip.label
  eco <- ifelse(grepl("^[AB]", tips), 1L, 0L)
  names(eco) <- tips
  eco
}


#' Bipartitions of interest for scoring.
#'
#' Returns a list of tip-label sets representing
#'   - `trueSister`: `(A,C)` — present in the truth
#'   - `falseSister`: `(A,B)` — the convergence-driven false grouping
#'
#' @return Named list of character vectors.
.SimBipartitions <- function() {
  list(
    trueSister  = c(paste0("A", 1:4), paste0("C", 1:4)),
    falseSister = c(paste0("A", 1:4), paste0("B", 1:4))
  )
}
