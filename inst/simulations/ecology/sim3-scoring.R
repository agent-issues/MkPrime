# sim3-scoring.R --------------------------------------------------------------
#
# Root-invariant bipartition scoring for Sim 3 (v4 and v4-cross) MCMC chains.
#
# The original run_v4*.R / run_v4cross*.R harnesses scored target
# bipartitions with ape::prop.part(), which returns ROOTED partitions: the
# set of tips descending from each internal node given the current root
# placement. When MCMC trees are rooted in different positions across
# samples (they are — see v4cross-b-pt5 rep01: 8 root configurations seen,
# ~30 % of samples rooted inside the BD side), prop.part fails to
# recognise a target bipartition whenever the root lies inside the target
# tip set. The split exists in the unrooted tree but appears in
# prop.part's output as its complement.
#
# Concrete symptom: P(AC) and P(BD) — the SAME unrooted bipartition on a
# binary tree — receive different scores (e.g. 0.978 vs 0.148).
#
# Fix: use TreeTools::as.Splits() to enumerate unrooted bipartitions, and
# test target presence as
#     (row == target)  OR  (row == !target)
# in tip-label-ordered logical form.
#
# This file provides:
#   HasBipartSplits()      — root-invariant per-tree presence test
#   ScoreTreesUnrooted()   — mirrors the legacy scoreTrees() in run_v4*.R
#   ScoreTreesLegacy()     — exact copy of the original buggy scorer, for
#                            side-by-side comparison.
#   RootInsideSet()        — for diagnostic reporting: per-tree, is the root
#                            edge attached to a tip in the target set?
#
# Dependencies: TreeTools (as.Splits), ape (prop.part, root), TreeDist
# (ClusteringInfoDistance). These are already on every Hamilton library
# path used by the sim3 harnesses.


#' Root-invariant per-tree bipartition presence test.
#'
#' For each tree in `treeList`, return TRUE iff the unrooted bipartition
#' separating `splitTips` from the remaining tips is present.
#'
#' Implementation: Each tree's set of (unrooted) bipartitions is obtained
#' via `TreeTools::as.Splits()` and coerced to a logical matrix (rows =
#' splits, columns = tips, in the order of `treeList[[1]]$tip.label`). A
#' target tip-set is converted to a logical vector under the same tip
#' ordering. The split is present iff some row equals the target vector
#' OR its complement.
#'
#' @param treeList multiPhylo or list of `phylo`. All trees must share the
#'   same tip-label set; tip labels need not be in the same order across
#'   trees (Splits is computed in `tipLabels = treeList[[1]]$tip.label`
#'   order).
#' @param splitTips character vector of tip labels defining one side of
#'   the target split.
#'
#' @return logical vector of length `length(treeList)`.
HasBipartSplits <- function(treeList, splitTips) {
  if (length(treeList) == 0L) return(logical(0))
  tipLab <- treeList[[1]]$tip.label
  target <- tipLab %in% splitTips
  vapply(treeList, function(tr) {
    spl <- TreeTools::as.Splits(tr, tipLabels = tipLab)
    # as.logical(spl) gives a logical matrix with rows = splits and
    # columns = tips, columns named by tip label in `tipLab` order.
    m <- as.logical(spl)
    if (is.null(dim(m))) m <- matrix(m, nrow = 1L)
    any(apply(m, 1L, function(r) all(r == target) || all(r == !target)))
  }, logical(1))
}


#' Legacy (buggy) per-tree bipartition presence test, kept for comparison.
#'
#' Exact copy of the `hasBipart` function inlined in every run_v4*.R and
#' run_v4cross*.R script. Tests the ROOTED partition produced by
#' `ape::prop.part`; sensitive to root placement.
#'
#' @inheritParams HasBipartSplits
#' @return logical vector of length `length(treeList)`.
HasBipartLegacy <- function(treeList, splitTips) {
  if (length(treeList) == 0L) return(logical(0))
  vapply(treeList, function(tr) {
    cl   <- ape::prop.part(tr)
    tips <- attr(cl, "labels")
    splitSet <- which(tips %in% splitTips)
    any(vapply(cl, function(p) setequal(p, splitSet), logical(1)))
  }, logical(1))
}


#' Discard the first 25 % of a sample as burn-in.
#'
#' Identical to the inline `discard` function in run_v4*.R: drop the
#' first `ceiling(n / 4)` samples and return the remainder.
.DiscardBurnin <- function(x) {
  n <- length(x)
  if (n <= 1L) return(x)
  keepFrom <- ceiling(n / 4) + 1L
  if (keepFrom > n) return(x[integer(0)])
  x[seq.int(keepFrom, n)]
}


#' Score a sample of trees against a set of named bipartitions.
#'
#' Mirrors the legacy `scoreTrees()` inside run_v4*.R and run_v4cross*.R,
#' but uses [HasBipartSplits()] (root-invariant) for all per-bipartition
#' calculations, and additionally returns the matching legacy values for
#' side-by-side comparison.
#'
#' @param treeList multiPhylo (typically `res$trees` from `RunMkPrime`).
#' @param biparts  Named list of character vectors (tip-label sets) to
#'   score. Names become metric labels in the returned data frame.
#' @param refTree  Reference (true) tree for `TreeDist::ClusteringInfoDistance`.
#' @param discardBurnin Logical; if TRUE (default), drop the first 25 % of
#'   `treeList` before scoring (matches the original harness).
#'
#' @return Data frame with columns
#'   `metric` (one row per element of `biparts`, plus "CID"),
#'   `value_corrected` (Splits-based),
#'   `value_legacy`    (prop.part-based; for "CID" they are equal),
#'   `nTrees`          (number of post-burnin trees scored).
ScoreTreesUnrooted <- function(treeList, biparts, refTree,
                               discardBurnin = TRUE) {
  if (isTRUE(discardBurnin)) {
    tr <- .DiscardBurnin(treeList)
  } else {
    tr <- treeList
  }
  class(tr) <- "multiPhylo"
  nTrees <- length(tr)

  rows <- lapply(names(biparts), function(nm) {
    target <- biparts[[nm]]
    pNew <- mean(HasBipartSplits(tr, target))
    pOld <- mean(HasBipartLegacy(tr, target))
    data.frame(metric          = nm,
               value_corrected = pNew,
               value_legacy    = pOld,
               nTrees          = nTrees,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)

  cid <- if (nTrees > 0L) {
    mean(as.numeric(
      TreeDist::ClusteringInfoDistance(tr, refTree, normalize = TRUE)))
  } else {
    NA_real_
  }
  out <- rbind(out,
               data.frame(metric          = "CID",
                          value_corrected = cid,
                          value_legacy    = cid,
                          nTrees          = nTrees,
                          stringsAsFactors = FALSE))
  out
}


#' Per-tree report of root-edge attachment.
#'
#' For diagnostic / reporting use. Returns a logical vector of length
#' `length(treeList)` indicating, per tree, whether the rooting edge is
#' attached to a tip that lies inside `splitTips`. Trees with no rooting
#' (unrooted) yield NA. This is the rough bound on how much the legacy
#' `prop.part` scoring can disagree with the corrected scoring for the
#' bipartition `splitTips | rest`.
#'
#' Note: ape phylo objects always carry an explicit root in the edge
#' matrix, so "is.rooted" is determined here by whether root edge points
#' to a leaf on either the target side or the complement side.
#'
#' @param treeList multiPhylo or list of phylo.
#' @param splitTips character vector defining one side of the split.
#' @return logical vector; TRUE = root is inside `splitTips`, FALSE =
#'   inside the complement.
RootInsideSet <- function(treeList, splitTips) {
  if (length(treeList) == 0L) return(logical(0))
  vapply(treeList, function(tr) {
    # The root node is (nTip + 1) in an ape phylo. Descend the tree from
    # the root and ask which side of the root the FIRST-found tip from
    # `splitTips` lies on. Equivalent: of the two clades hanging off the
    # root, which one (if any) is a subset of splitTips?
    tipLab <- tr$tip.label
    nT <- length(tipLab)
    rootNode <- nT + 1L
    rootChildren <- tr$edge[tr$edge[, 1L] == rootNode, 2L]
    if (length(rootChildren) < 2L) return(NA)
    # First child: collect descendant tips.
    descTips <- function(node) {
      if (node <= nT) return(tipLab[node])
      kids <- tr$edge[tr$edge[, 1L] == node, 2L]
      unlist(lapply(kids, descTips))
    }
    side1 <- descTips(rootChildren[1L])
    inTarget1 <- all(side1 %in% splitTips)
    inComp1   <- all(!(side1 %in% splitTips))
    if (inTarget1) TRUE
    else if (inComp1) FALSE
    else NA  # root edge cuts ACROSS the target split (also affects scoring)
  }, logical(1))
}


#' Distribution of distinct rootings (root-position bitmasks) in a tree list.
#'
#' For each tree, encode "which tip set is the descendant of root child 1"
#' as a sorted comma-separated string. Returns the table of frequencies.
#'
#' @param treeList multiPhylo or list of phylo.
#' @return Named integer vector (counts), with names being concatenated
#'   tip-label sets defining the rooting.
RootConfigCounts <- function(treeList) {
  if (length(treeList) == 0L) return(integer(0))
  keys <- vapply(treeList, function(tr) {
    tipLab <- tr$tip.label
    nT <- length(tipLab)
    rootNode <- nT + 1L
    rootChildren <- tr$edge[tr$edge[, 1L] == rootNode, 2L]
    if (length(rootChildren) < 2L) return(NA_character_)
    descTips <- function(node) {
      if (node <= nT) return(tipLab[node])
      kids <- tr$edge[tr$edge[, 1L] == node, 2L]
      unlist(lapply(kids, descTips))
    }
    side1 <- sort(descTips(rootChildren[1L]))
    side2 <- sort(descTips(rootChildren[2L]))
    # Canonicalise: use the lexicographically smaller side as the key.
    key1 <- paste(side1, collapse = ",")
    key2 <- paste(side2, collapse = ",")
    if (key1 < key2) key1 else key2
  }, character(1))
  table(keys, useNA = "ifany")
}
