# MCMC proposal functions for MkPrime
#
# Scalar/simplex proposals return list(value = ..., log_hastings = ...).
# Tree topology proposals return list(tree = ..., rel_br_lengths = ...,
#   log_hastings = ...).
# The Hastings ratio is log q(current | proposed) - log q(proposed | current).

#' Scale proposal for positive scalars
#'
#' Proposes x' = x * exp(tuning * (U - 0.5)) where U ~ Uniform(0,1).
#' Hastings ratio: log(x'/x) = tuning * (U - 0.5).
#'
#' @param x Current value (positive scalar).
#' @param tuning Scale parameter controlling proposal width.
#' @return `list(value, log_hastings)`.
#' @keywords internal
propose_scale <- function(x, tuning = 1.0) {
  u <- runif(1)
  m <- exp(tuning * (u - 0.5))
  list(value = x * m, log_hastings = log(m))
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
#' @return `list(value, log_hastings)`.
#' @keywords internal
propose_beta_simplex <- function(x, index = NULL, tuning = 10.0) {
  n <- length(x)
  if (n < 2L) {
    return(list(value = x, log_hastings = 0))
  }

  if (is.null(index)) {
    index <- sample.int(n, 1L)
  }
  # Pick a random other element
  other <- sample.int(n - 1L, 1L)
  if (other >= index) other <- other + 1L

  # Current values of the two elements
  old_a <- x[index]
  old_b <- x[other]
  total <- old_a + old_b

  if (total <= 0) {
    return(list(value = x, log_hastings = 0))
  }

  # Current fraction: f = old_a / total
  old_f <- old_a / total

  # Propose new fraction from Beta centered on old_f
  alpha <- old_f * tuning + 1
  beta_param <- (1 - old_f) * tuning + 1
  new_f <- rbeta(1, alpha, beta_param)

  # Compute new values
  new_a <- new_f * total
  new_b <- (1 - new_f) * total
  x_new <- x
  x_new[index] <- new_a
  x_new[other] <- new_b

  # Hastings ratio: q(old|new) / q(new|old)
  # Forward: Beta(old_f * tuning + 1, (1-old_f) * tuning + 1) at new_f
  # Reverse: Beta(new_f * tuning + 1, (1-new_f) * tuning + 1) at old_f
  log_fwd <- dbeta(new_f, alpha, beta_param, log = TRUE)
  rev_alpha <- new_f * tuning + 1
  rev_beta <- (1 - new_f) * tuning + 1
  log_rev <- dbeta(old_f, rev_alpha, rev_beta, log = TRUE)

  list(value = x_new, log_hastings = log_rev - log_fwd)
}


#' BoundedIntegerWalk proposal
#'
#' Proposes x' = x + delta where delta ~ Uniform(-window, ..., window).
#' Rejects (returns current value with log_hastings = -Inf) if x' < lower.
#'
#' @param x Current integer value.
#' @param lower Lower bound (inclusive).
#' @param window Half-width of the proposal window.
#' @return `list(value, log_hastings)`.
#' @keywords internal
propose_bounded_int_walk <- function(x, lower, window = 1L) {
  delta <- sample(-window:window, 1L)
  x_new <- x + delta

  if (x_new < lower) {
    return(list(value = x, log_hastings = -Inf))
  }

  # Symmetric proposal: log_hastings = 0
  # (Both x→x' and x'→x have same number of valid proposals)
  list(value = x_new, log_hastings = 0)
}


#' NNI proposal on an unrooted binary tree
#'
#' Picks a random internal edge and swaps one subtree from each side,
#' producing a nearest-neighbor interchange. The proposal is symmetric
#' (Hastings ratio = 1).
#'
#' For unrooted binary trees in ape's rooted representation, the root has
#' degree 3 and all other internal nodes have degree 3 (1 parent + 2
#' children). Internal edges are those where both endpoints are internal.
#'
#' @param tree A `phylo` object in postorder (unrooted binary).
#' @param tree_length Current total tree length.
#' @param rel_br_lengths Current relative branch lengths (simplex).
#' @return `list(tree, rel_br_lengths, log_hastings)`. The returned tree
#'   is in postorder.
#' @keywords internal
propose_nni <- function(tree, tree_length, rel_br_lengths) {
  nTip <- length(tree$tip.label)
  edge <- tree$edge

  # Internal edges: both endpoints are internal nodes
  internal_rows <- which(edge[, 1] > nTip & edge[, 2] > nTip)

  if (length(internal_rows) == 0L) {
    # Too few tips for NNI (n <= 3)
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                log_hastings = -Inf))
  }

  # Pick a random internal edge
  edge_idx <- if (length(internal_rows) == 1L) {
    internal_rows
  } else {
    sample(internal_rows, 1L)
  }
  u <- edge[edge_idx, 1]
  v <- edge[edge_idx, 2]

  # v's children (always exactly 2)
  v_child_rows <- which(edge[, 1] == v)
  v_children <- edge[v_child_rows, 2]

  # u's children other than v (1 if non-root, 2 if root)
  u_sib_rows <- which(edge[, 1] == u & edge[, 2] != v)
  u_siblings <- edge[u_sib_rows, 2]

  # Pick one child of v and one sibling on u's side to swap
  c_idx <- sample.int(length(v_children), 1L)
  w_idx <- sample.int(length(u_siblings), 1L)

  c_row <- v_child_rows[c_idx]
  w_row <- u_sib_rows[w_idx]

  # Swap: change parent of selected v-child to u, parent of u-sibling to v
  new_edge <- edge
  new_edge[c_row, 1] <- u
  new_edge[w_row, 1] <- v

  # Rebuild tree and reorder to postorder
  # Edge lengths must be absolute for reorder to preserve associations
  new_tree <- tree
  new_tree$edge <- new_edge
  new_tree$edge.length <- tree_length * rel_br_lengths
  new_tree <- ape::reorder.phylo(new_tree, "postorder")

  # Recompute relative branch lengths (row order may have changed)
  new_rel_br <- new_tree$edge.length / tree_length

  list(tree = new_tree, rel_br_lengths = new_rel_br, log_hastings = 0)
}


#' SPR proposal on an unrooted binary tree
#'
#' Subtree pruning and regrafting. Picks a random edge (u->v), detaches v's
#' subtree, suppresses u (merging its parent and sibling edges), then
#' reattaches u on a random backbone edge. Node u "moves" so the edge
#' matrix stays the same size; no node renumbering needed.
#'
#' Edges adjacent to the root cannot be pruned in this R prototype (the
#' root-degree change would require node renumbering). This restricts
#' 3 out of 2n-3 edges. C++ implementation (Phase 5) will handle this.
#'
#' The Hastings ratio includes a Jacobian correction for the edge length
#' redistribution: `log(l_regraft) - log(l_merge)`, where `l_regraft` is
#' the regraft edge length and `l_merge` is the sum of the two edges
#' merged during suppression.
#'
#' @param tree A `phylo` object in postorder (unrooted binary).
#' @param tree_length Current total tree length.
#' @param rel_br_lengths Current relative branch lengths (simplex).
#' @return `list(tree, rel_br_lengths, log_hastings)`. The returned tree
#'   is in postorder.
#' @keywords internal
propose_spr <- function(tree, tree_length, rel_br_lengths) {
  nTip <- length(tree$tip.label)
  nEdge <- nrow(tree$edge)
  edge <- tree$edge
  edge_length <- tree_length * rel_br_lengths
  root <- nTip + 1L

  # Eligible prune edges: parent != root
  eligible_prune <- which(edge[, 1] != root)
  if (length(eligible_prune) == 0L) {
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                log_hastings = -Inf))
  }

  prune_row <- if (length(eligible_prune) == 1L) {
    eligible_prune
  } else {
    sample(eligible_prune, 1L)
  }
  u <- edge[prune_row, 1]
  v <- edge[prune_row, 2]

  # u's parent and v's sibling (u != root, so exactly 1 parent, 1 sibling)
  parent_row <- which(edge[, 2] == u)
  p <- edge[parent_row, 1]
  sib_row <- which(edge[, 1] == u & edge[, 2] != v)
  w <- edge[sib_row, 2]

  # Find all descendants of v (for exclusion)
  desc_v <- .descendants(v, edge, nTip)

  # Exclude: subtree edges (child in v's subtree or v itself) + adjacent to u
  subtree_and_prune <- which(edge[, 2] %in% c(v, desc_v))
  adjacent_to_u <- which(edge[, 1] == u | edge[, 2] == u)
  exclude <- union(subtree_and_prune, adjacent_to_u)
  candidate_regraft <- setdiff(seq_len(nEdge), exclude)

  if (length(candidate_regraft) == 0L) {
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                log_hastings = -Inf))
  }

  regraft_row <- if (length(candidate_regraft) == 1L) {
    candidate_regraft
  } else {
    sample(candidate_regraft, 1L)
  }
  a <- edge[regraft_row, 1]
  b <- edge[regraft_row, 2]

  tau <- runif(1)

  # Lengths for Hastings ratio
  l_regraft <- edge_length[regraft_row]
  l_merge <- edge_length[parent_row] + edge_length[sib_row]

  # --- Perform the SPR ---
  new_edge <- edge
  new_el <- edge_length

  # 1. Suppress u: (p -> u) becomes (p -> w), length absorbs (u -> w)
  new_edge[parent_row, 2] <- w
  new_el[parent_row] <- l_merge

  # 2. Insert u on regraft edge: (a -> b) becomes (a -> u)
  new_edge[regraft_row, 2] <- u
  new_el[regraft_row] <- tau * l_regraft

  # 3. Reuse sib_row for (u -> b)
  new_edge[sib_row, ] <- c(u, b)
  new_el[sib_row] <- (1 - tau) * l_regraft

  # Prune edge (u -> v) unchanged

  # Rebuild tree and reorder to postorder
  new_tree <- tree
  new_tree$edge <- new_edge
  new_tree$edge.length <- new_el
  new_tree <- ape::reorder.phylo(new_tree, "postorder")

  new_rel_br <- new_tree$edge.length / tree_length

  # Hastings: Jacobian from splitting l_regraft vs merging l_merge
  # Discrete part cancels (same |eligible| and |candidates| both ways)
  log_hastings <- log(l_regraft) - log(l_merge)

  list(tree = new_tree, rel_br_lengths = new_rel_br, log_hastings = log_hastings)
}


#' Find all descendant node IDs (BFS)
#' @keywords internal
.descendants <- function(node, edge, nTip) {
  children <- edge[edge[, 1] == node, 2]
  if (length(children) == 0L) return(integer(0))
  desc <- children
  queue <- children[children > nTip]
  while (length(queue) > 0L) {
    cur <- queue[1L]
    queue <- queue[-1L]
    ch <- edge[edge[, 1] == cur, 2]
    desc <- c(desc, ch)
    queue <- c(queue, ch[ch > nTip])
  }
  desc
}
