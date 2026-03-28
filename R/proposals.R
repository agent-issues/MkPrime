# MCMC proposal functions for MkPrime
#
# Scalar/simplex proposals return list(value = ..., logHastings = ...).
# Tree topology proposals return list(tree = ..., rel_br_lengths = ...,
#   logHastings = ...).
# The Hastings ratio is log q(current | proposed) - log q(proposed | current).

#' Scale proposal for positive scalars
#'
#' Proposes x' = x * exp(tuning * (U - 0.5)) where U ~ Uniform(0,1).
#' Hastings ratio: log(x'/x) = tuning * (U - 0.5).
#'
#' @param x Current value (positive scalar).
#' @param tuning Scale parameter controlling proposal width.
#' @return `list(value, logHastings)`.
#' @keywords internal
ProposeScale <- function(x, tuning = 1.0) {
  u <- runif(1)
  m <- exp(tuning * (u - 0.5))
  list(value = x * m, logHastings = log(m))
}


#' BetaSimplex proposal for simplex vectors
#'
#' Picks element `index` and a random other element, proposes redistributing
#' mass between them using a Beta draw. Other elements are unchanged.
#'
#' @param x Current simplex vector (sums to 1, all positive).
#' @param index Index of the element to perturb. If NULL, picks randomly.
#' @param tuning Concentration parameter (higher = more conservative).
#'   Effective concentration is `tuning * 2`.
#' @return `list(value, logHastings)`.
#' @keywords internal
ProposeBetaSimplex <- function(x, index = NULL, tuning = 10.0) {
  n <- length(x)
  if (n < 2L) {
    return(list(value = x, logHastings = 0))
  }

  if (is.null(index)) {
    index <- sample.int(n, 1L)
  }

  # C++ uses 0-based indexing
  beta_simplex_proposal(x, index - 1L, tuning)
}


#' BoundedIntegerWalk proposal
#'
#' Proposes x' = x + delta where delta ~ Uniform(-window, ..., window).
#' Rejects (returns current value with logHastings = -Inf) if x' < lower.
#'
#' @param x Current integer value.
#' @param lower Lower bound (inclusive).
#' @param window Half-width of the proposal window.
#' @return `list(value, logHastings)`.
#' @keywords internal
ProposeBoundedIntWalk <- function(x, lower, window = 1L) {
  delta <- sample(-window:window, 1L)
  xNew <- x + delta

  if (xNew < lower) {
    return(list(value = x, logHastings = -Inf))
  }

  # Symmetric proposal: logHastings = 0
  list(value = xNew, logHastings = 0)
}


#' NNI proposal on an unrooted binary tree
#'
#' Picks a random internal edge and swaps one subtree from each side,
#' producing a nearest-neighbor interchange. The proposal is symmetric
#' (Hastings ratio = 1).
#'
#' Delegates to C++ (\code{nni_proposal} in tree_moves.cpp) for edge
#' manipulation and canonical preorder reordering.
#'
#' @param tree A `phylo` object in canonical preorder (unrooted binary).
#' @param tree_length Current total tree length.
#' @param rel_br_lengths Current relative branch lengths (simplex).
#' @return `list(tree, rel_br_lengths, logHastings)`. The returned tree
#'   is in canonical preorder.
#' @keywords internal
ProposeNni <- function(tree, tree_length, rel_br_lengths) {
  result <- nni_proposal(tree$edge, length(tree$tip.label),
                         tree_length, rel_br_lengths)

  if (!is.finite(result$logHastings)) {
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                logHastings = -Inf))
  }

  newTree <- tree
  newTree$edge <- result$edge
  newTree$edge.length <- tree_length * result$rel_br_lengths

  list(tree = newTree, rel_br_lengths = result$rel_br_lengths,
       logHastings = result$logHastings)
}


#' SPR proposal on an unrooted binary tree
#'
#' Subtree pruning and regrafting. Picks a random edge (u->v), detaches v's
#' subtree, suppresses u (merging its parent and sibling edges), then
#' reattaches u on a random backbone edge.
#'
#' Delegates to C++ (\code{spr_proposal} in proposals.cpp) for edge
#' manipulation and canonical preorder reordering.
#'
#' The Hastings ratio includes a Jacobian correction for the edge length
#' redistribution: `log(lRegraft) - log(lMerge)`.
#'
#' @param tree A `phylo` object in canonical preorder (unrooted binary).
#' @param tree_length Current total tree length.
#' @param rel_br_lengths Current relative branch lengths (simplex).
#' @return `list(tree, rel_br_lengths, logHastings)`. The returned tree
#'   is in canonical preorder.
#' @keywords internal
ProposeSpr <- function(tree, tree_length, rel_br_lengths) {
  result <- spr_proposal(tree$edge, length(tree$tip.label),
                         tree_length, rel_br_lengths)

  if (!is.finite(result$logHastings)) {
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                logHastings = -Inf))
  }

  newTree <- tree
  newTree$edge <- result$edge
  newTree$edge.length <- tree_length * result$rel_br_lengths

  list(tree = newTree, rel_br_lengths = result$rel_br_lengths,
       logHastings = result$logHastings)
}


#' TBR topology proposal (M-053)
#'
#' Tree Bisection and Reconnection. Extends SPR by additionally re-rooting
#' the pruned subtree at a random internal edge before regrafting.
#'
#' @param tree A `phylo` object in canonical preorder (unrooted binary).
#' @param tree_length Current total tree length.
#' @param rel_br_lengths Current relative branch lengths (simplex).
#' @return `list(tree, rel_br_lengths, logHastings)`. The returned tree
#'   is in canonical preorder.
#' @keywords internal
ProposeTbr <- function(tree, tree_length, rel_br_lengths) {
  result <- tbr_proposal(tree$edge, length(tree$tip.label),
                         tree_length, rel_br_lengths)

  if (!is.finite(result$logHastings)) {
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                logHastings = -Inf))
  }

  newTree <- tree
  newTree$edge <- result$edge
  newTree$edge.length <- tree_length * result$rel_br_lengths

  list(tree = newTree, rel_br_lengths = result$rel_br_lengths,
       logHastings = result$logHastings)
}
