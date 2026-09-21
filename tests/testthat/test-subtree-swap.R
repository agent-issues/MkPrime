# Tests for swap_subtrees / get_valid_swap_partners (M-084)
#
# swap_subtrees_cpp(edge, nTip, treeLength, relBrLengths, nodeA, nodeB)
# get_valid_swap_partners_cpp(edge, nTip, pruneNode)

# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

# 5-tip preorder tree: ((t1,t2),(t3,(t4,t5)));
# Canonical preorder renumbers internal nodes in visit order:
#   root=6, first visited subtree internal=7, etc.
.tree5 <- function() {
  tr <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,(t4:0.1,t5:0.2):0.12):0.18);")
  Preorder(tr)
}

.edge5   <- function() .tree5()$edge
.relBr5  <- function() { tr <- .tree5(); tr$edge.length / sum(tr$edge.length) }
.treeL5  <- function() sum(.tree5()$edge.length)

# Get all tips that descend from `node` in a tree given by `edge`.
.tips_below <- function(edge, node, nTip) {
  if (node <= nTip) return(node)
  tips <- integer(0)
  queue <- node
  while (length(queue) > 0) {
    cur <- queue[1L]; queue <- queue[-1L]
    children <- edge[edge[, 1L] == cur, 2L]
    for (ch in children) {
      if (ch <= nTip) tips <- c(tips, ch) else queue <- c(queue, ch)
    }
  }
  sort(tips)
}

# Find the edge-row index whose child subtree contains exactly `tips_set`.
.find_subtree_row <- function(edge, tips_set, nTip) {
  for (i in seq_len(nrow(edge))) {
    nd <- edge[i, 2L]
    if (nd <= nTip) { if (length(tips_set) == 1L && nd == tips_set) return(i); next }
    if (setequal(.tips_below(edge, nd, nTip), tips_set)) return(i)
  }
  NA_integer_
}

# ---------------------------------------------------------------------------
# Test: round-trip swap returns original tree
# ---------------------------------------------------------------------------

test_that("double swap is identity", {
  edge   <- .edge5()
  relBr  <- .relBr5()
  treeL  <- .treeL5()
  nTip   <- 5L

  # Use TIP nodes (labels stable across canonical renumbering)
  res1 <- swap_subtrees_cpp(edge, nTip, treeL, relBr, nodeA = 1L, nodeB = 4L)
  expect_equal(res1$logHastings, 0.0)

  # Swap the same tips again → identity
  res2 <- swap_subtrees_cpp(res1$edge, nTip, treeL, res1$rel_br_lengths,
                             nodeA = 1L, nodeB = 4L)
  expect_equal(res2$logHastings, 0.0)

  # Canonical preorder is unique per topology → identical edge matrices
  expect_equal(res2$edge, edge)
  expect_equal(res2$rel_br_lengths, relBr, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Test: tip swap produces correct topology
# ---------------------------------------------------------------------------

test_that("swapping two tip nodes exchanges them in the topology", {
  edge  <- .edge5()
  relBr <- .relBr5()
  treeL <- .treeL5()
  nTip  <- 5L

  # Swap tip 1 (t1) and tip 3 (t3)
  res <- swap_subtrees_cpp(edge, nTip, treeL, relBr, nodeA = 1L, nodeB = 3L)
  expect_equal(res$logHastings, 0.0)

  # After swap: t3 should be grouped with t2 (was t1's sibling)
  new_p_t3 <- res$edge[res$edge[, 2L] == 3L, 1L]
  new_p_t2 <- res$edge[res$edge[, 2L] == 2L, 1L]
  expect_equal(new_p_t3, new_p_t2)

  # t1 should now be in the other clade, next to (t4,t5) subtree
  new_p_t1 <- res$edge[res$edge[, 2L] == 1L, 1L]
  tips_under_t1_parent <- .tips_below(res$edge, new_p_t1, nTip)
  expect_true(setequal(tips_under_t1_parent, c(1L, 4L, 5L)))
})

# ---------------------------------------------------------------------------
# Test: tree length preserved after swap
# ---------------------------------------------------------------------------

test_that("tree length is preserved after swap", {
  edge  <- .edge5()
  relBr <- .relBr5()
  treeL <- .treeL5()

  res <- swap_subtrees_cpp(edge, 5L, treeL, relBr, nodeA = 7L, nodeB = 9L)

  expect_equal(sum(res$rel_br_lengths), sum(relBr), tolerance = 1e-12)
  expect_equal(sum(treeL * res$rel_br_lengths), treeL, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Test: branch lengths swap with subtrees
# ---------------------------------------------------------------------------

test_that("branch lengths swap when subtrees are exchanged", {
  edge  <- .edge5()
  relBr <- .relBr5()
  treeL <- .treeL5()
  nTip  <- 5L

  # Identify subtrees by tip content (stable across renumbering)
  rowA <- .find_subtree_row(edge, c(1L, 2L), nTip)  # (t1, t2) subtree
  rowB <- .find_subtree_row(edge, c(4L, 5L), nTip)  # (t4, t5) subtree
  origA <- relBr[rowA]
  origB <- relBr[rowB]

  res <- swap_subtrees_cpp(edge, nTip, treeL, relBr, nodeA = 7L, nodeB = 9L)

  # After swap, the subtree containing {1,2} now has the branch length from {4,5}
  newRowA <- .find_subtree_row(res$edge, c(1L, 2L), nTip)
  newRowB <- .find_subtree_row(res$edge, c(4L, 5L), nTip)

  expect_equal(res$rel_br_lengths[newRowA], origB, tolerance = 1e-14)
  expect_equal(res$rel_br_lengths[newRowB], origA, tolerance = 1e-14)
})

# ---------------------------------------------------------------------------
# Test: get_valid_swap_partners excludes descendants
# ---------------------------------------------------------------------------

test_that("get_valid_swap_partners excludes descendants of pruneNode", {
  edge <- .edge5()
  nTip <- 5L

  # Internal node 8 = (t3,(t4,t5)), so 3, 4, 5, 9 are its descendants
  partners8 <- get_valid_swap_partners_cpp(edge, nTip, pruneNode = 8L)

  # Descendants of 8: {3, 9, 4, 5}; also 8 itself excluded
  for (desc in c(3L, 9L, 4L, 5L)) {
    expect_false(desc %in% partners8,
                 label = paste("Descendant", desc, "not a valid partner of 8"))
  }
  expect_false(8L %in% partners8)
})

# ---------------------------------------------------------------------------
# Test: get_valid_swap_partners excludes ancestors
# ---------------------------------------------------------------------------

test_that("get_valid_swap_partners excludes ancestors of pruneNode", {
  edge <- .edge5()
  nTip <- 5L

  # Ancestors of node 9 (=(t4,t5)): 8 and 6 (root)
  partners9 <- get_valid_swap_partners_cpp(edge, nTip, pruneNode = 9L)

  expect_false(8L %in% partners9, label = "Ancestor 8 not a valid partner of 9")
  expect_false(6L %in% partners9, label = "Root 6 not a valid partner of 9")
  expect_false(9L %in% partners9, label = "9 not its own partner")
})

# ---------------------------------------------------------------------------
# Test: get_valid_swap_partners excludes siblings
# ---------------------------------------------------------------------------

test_that("get_valid_swap_partners excludes siblings", {
  edge <- .edge5()
  nTip <- 5L

  # Node 7 = (t1,t2), sibling is 8 = (t3,(t4,t5)) (both children of root 6)
  partners7 <- get_valid_swap_partners_cpp(edge, nTip, pruneNode = 7L)
  expect_false(8L %in% partners7, label = "Sibling 8 not a valid partner of 7")
})

# ---------------------------------------------------------------------------
# Test: get_valid_swap_partners returns non-empty set for typical node
# ---------------------------------------------------------------------------

test_that("get_valid_swap_partners returns valid candidates", {
  edge <- .edge5()
  nTip <- 5L

  # For node 7 (=(t1,t2)), valid partners are non-desc, non-anc, non-sibling:
  # All nodes: 1,2,3,4,5,6,7,8,9 (but 6=root has no parent, not a child)
  # Children in edge: 7,8,1,2,3,9,4,5
  # Exclude: desc(7)={7,1,2}, anc(7)={6}, sib(7)={8}
  # Valid: {3, 4, 5, 9}
  partners7 <- get_valid_swap_partners_cpp(edge, nTip, pruneNode = 7L)
  expect_true(length(partners7) > 0)
  for (w in partners7) {
    expect_true(w %in% c(3L, 4L, 5L, 9L),
                label = paste("Unexpected partner", w, "for node 7"))
  }
  expect_setequal(partners7, c(3L, 4L, 5L, 9L))
})

# ---------------------------------------------------------------------------
# Test: root node returns empty partner list
# ---------------------------------------------------------------------------

test_that("get_valid_swap_partners for root returns empty", {
  edge <- .edge5()
  nTip <- 5L
  root_node <- nTip + 1L  # = 6

  partners_root <- get_valid_swap_partners_cpp(edge, nTip, pruneNode = root_node)
  expect_equal(length(partners_root), 0L)
})

# ---------------------------------------------------------------------------
# Test: swap of invalid nodes (root) returns R_NegInf logHastings
# ---------------------------------------------------------------------------

test_that("swap with root node gives logHastings = -Inf", {
  edge  <- .edge5()
  relBr <- .relBr5()
  treeL <- .treeL5()
  nTip  <- 5L
  root  <- nTip + 1L  # = 6

  res <- swap_subtrees_cpp(edge, nTip, treeL, relBr, nodeA = root, nodeB = 7L)
  expect_equal(res$logHastings, -Inf)
})

# ---------------------------------------------------------------------------
# Test: 4-tip tree — verify a concrete swap
# ---------------------------------------------------------------------------

test_that("swap on 4-tip tree matches hand-computed result", {
  # ((t1,t2),(t3,t4))
  tr   <- Preorder(read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);"))
  edge <- tr$edge
  relBr <- tr$edge.length / sum(tr$edge.length)
  treeL <- sum(tr$edge.length)
  nTip  <- 4L

  # Swap tip 1 (t1) with tip 3 (t3)
  res <- swap_subtrees_cpp(edge, nTip, treeL, relBr, nodeA = 1L, nodeB = 3L)
  expect_equal(res$logHastings, 0.0)

  # Result should be ((t3,t2),(t1,t4)) in terms of tip groupings:
  # t3 and t2 share a parent, t1 and t4 share a parent
  new_p1 <- res$edge[res$edge[, 2L] == 1L, 1L]
  new_p3 <- res$edge[res$edge[, 2L] == 3L, 1L]

  new_sibling_of_t1 <- setdiff(res$edge[res$edge[, 1L] == new_p1, 2L], 1L)
  expect_equal(new_sibling_of_t1, 4L)

  new_sibling_of_t3 <- setdiff(res$edge[res$edge[, 1L] == new_p3, 2L], 3L)
  expect_equal(new_sibling_of_t3, 2L)
})
