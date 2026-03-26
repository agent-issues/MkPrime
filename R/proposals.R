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
  # Pick a random other element
  other <- sample.int(n - 1L, 1L)
  if (other >= index) other <- other + 1L

  # Current values of the two elements
  oldA <- x[index]
  oldB <- x[other]
  total <- oldA + oldB

  if (total <= 0) {
    return(list(value = x, logHastings = 0))
  }

  # Current fraction: f = oldA / total
  oldF <- oldA / total

  # Propose new fraction from Beta centered on oldF
  alpha <- oldF * tuning + 1
  betaParam <- (1 - oldF) * tuning + 1
  newF <- rbeta(1, alpha, betaParam)

  # Compute new values
  newA <- newF * total
  newB <- (1 - newF) * total
  xNew <- x
  xNew[index] <- newA
  xNew[other] <- newB

  # Hastings ratio: q(old|new) / q(new|old)
  # Forward: Beta(oldF * tuning + 1, (1-oldF) * tuning + 1) at newF
  # Reverse: Beta(newF * tuning + 1, (1-newF) * tuning + 1) at oldF
  logFwd <- dbeta(newF, alpha, betaParam, log = TRUE)
  revAlpha <- newF * tuning + 1
  revBeta <- (1 - newF) * tuning + 1
  logRev <- dbeta(oldF, revAlpha, revBeta, log = TRUE)

  list(value = xNew, logHastings = logRev - logFwd)
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
  # (Both x→x' and x'→x have same number of valid proposals)
  list(value = xNew, logHastings = 0)
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
#' @return `list(tree, rel_br_lengths, logHastings)`. The returned tree
#'   is in postorder.
#' @keywords internal
ProposeNni <- function(tree, tree_length, rel_br_lengths) {
  nTip <- length(tree$tip.label)
  edge <- tree$edge

  # Internal edges: both endpoints are internal nodes
  internalRows <- which(edge[, 1] > nTip & edge[, 2] > nTip)

  if (length(internalRows) == 0L) {
    # Too few tips for NNI (n <= 3)
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                logHastings = -Inf))
  }

  # Pick a random internal edge
  edgeIdx <- if (length(internalRows) == 1L) {
    internalRows
  } else {
    sample(internalRows, 1L)
  }
  u <- edge[edgeIdx, 1]
  v <- edge[edgeIdx, 2]

  # v's children (always exactly 2)
  vChildRows <- which(edge[, 1] == v)
  vChildren <- edge[vChildRows, 2]

  # u's children other than v (1 if non-root, 2 if root)
  uSibRows <- which(edge[, 1] == u & edge[, 2] != v)
  uSiblings <- edge[uSibRows, 2]

  # Pick one child of v and one sibling on u's side to swap
  cIdx <- sample.int(length(vChildren), 1L)
  wIdx <- sample.int(length(uSiblings), 1L)

  cRow <- vChildRows[cIdx]
  wRow <- uSibRows[wIdx]

  # Swap: change parent of selected v-child to u, parent of u-sibling to v
  newEdge <- edge
  newEdge[cRow, 1] <- u
  newEdge[wRow, 1] <- v

  # Rebuild tree and reorder to postorder
  # Edge lengths must be absolute for reorder to preserve associations
  newTree <- tree
  newTree$edge <- newEdge
  newTree$edge.length <- tree_length * rel_br_lengths
  newTree <- TreeTools::Postorder(newTree)

  # Recompute relative branch lengths (row order may have changed)
  newRelBr <- newTree$edge.length / tree_length

  list(tree = newTree, rel_br_lengths = newRelBr, logHastings = 0)
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
#' redistribution: `log(lRegraft) - log(lMerge)`, where `lRegraft` is
#' the regraft edge length and `lMerge` is the sum of the two edges
#' merged during suppression.
#'
#' @param tree A `phylo` object in postorder (unrooted binary).
#' @param tree_length Current total tree length.
#' @param rel_br_lengths Current relative branch lengths (simplex).
#' @return `list(tree, rel_br_lengths, logHastings)`. The returned tree
#'   is in postorder.
#' @keywords internal
ProposeSpr <- function(tree, tree_length, rel_br_lengths) {
  nTip <- length(tree$tip.label)
  nEdge <- nrow(tree$edge)
  edge <- tree$edge
  edgeLength <- tree_length * rel_br_lengths
  root <- nTip + 1L

  # Eligible prune edges: parent != root
  eligiblePrune <- which(edge[, 1] != root)
  if (length(eligiblePrune) == 0L) {
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                logHastings = -Inf))
  }

  pruneRow <- if (length(eligiblePrune) == 1L) {
    eligiblePrune
  } else {
    sample(eligiblePrune, 1L)
  }
  u <- edge[pruneRow, 1]
  v <- edge[pruneRow, 2]

  # u's parent and v's sibling (u != root, so exactly 1 parent, 1 sibling)
  parentRow <- which(edge[, 2] == u)
  p <- edge[parentRow, 1]
  sibRow <- which(edge[, 1] == u & edge[, 2] != v)
  w <- edge[sibRow, 2]

  # Find all descendants of v (for exclusion)
  descV <- .Descendants(v, edge, nTip)

  # Exclude: subtree edges (child in v's subtree or v itself) + adjacent to u
  subtreeAndPrune <- which(edge[, 2] %in% c(v, descV))
  adjacentToU <- which(edge[, 1] == u | edge[, 2] == u)
  exclude <- union(subtreeAndPrune, adjacentToU)
  candidateRegraft <- setdiff(seq_len(nEdge), exclude)

  if (length(candidateRegraft) == 0L) {
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                logHastings = -Inf))
  }

  regraftRow <- if (length(candidateRegraft) == 1L) {
    candidateRegraft
  } else {
    sample(candidateRegraft, 1L)
  }
  a <- edge[regraftRow, 1]
  b <- edge[regraftRow, 2]

  tau <- runif(1)

  # Lengths for Hastings ratio
  lRegraft <- edgeLength[regraftRow]
  lMerge <- edgeLength[parentRow] + edgeLength[sibRow]

  # --- Perform the SPR ---
  newEdge <- edge
  newEl <- edgeLength

  # 1. Suppress u: (p -> u) becomes (p -> w), length absorbs (u -> w)
  newEdge[parentRow, 2] <- w
  newEl[parentRow] <- lMerge

  # 2. Insert u on regraft edge: (a -> b) becomes (a -> u)
  newEdge[regraftRow, 2] <- u
  newEl[regraftRow] <- tau * lRegraft

  # 3. Reuse sibRow for (u -> b)
  newEdge[sibRow, ] <- c(u, b)
  newEl[sibRow] <- (1 - tau) * lRegraft

  # Prune edge (u -> v) unchanged

  # Rebuild tree and reorder to postorder
  newTree <- tree
  newTree$edge <- newEdge
  newTree$edge.length <- newEl
  newTree <- TreeTools::Postorder(newTree)

  newRelBr <- newTree$edge.length / tree_length

  # Hastings: Jacobian from splitting lRegraft vs merging lMerge
  # Discrete part cancels (same |eligible| and |candidates| both ways)
  logHastings <- log(lRegraft) - log(lMerge)

  list(tree = newTree, rel_br_lengths = newRelBr, logHastings = logHastings)
}


#' Find all descendant node IDs (BFS)
#' @keywords internal
.Descendants <- function(node, edge, nTip) {
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
