# Tests for swap_subtrees / get_valid_swap_partners (M-084)
#
# swap_subtrees_cpp(edge, nTip, treeLength, relBrLengths, nodeA, nodeB)
# get_valid_swap_partners_cpp(edge, nTip, pruneNode)

library(ape)
library(TreeTools)

# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

# 5-tip postorder tree: ((t1,t2),(t3,(t4,t5)));
# Internal nodes: 6=root(t1,t2 | t3,(t4,t5)), 7=(t1,t2), 8=(t3,(t4,t5)), 9=(t4,t5)
.tree5 <- function() {
  tr <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,(t4:0.1,t5:0.2):0.12):0.18);")
  Postorder(tr)
}

.edge5   <- function() .tree5()$edge
.relBr5  <- function() { tr <- .tree5(); tr$edge.length / sum(tr$edge.length) }
.treeL5  <- function() sum(.tree5()$edge.length)


# ---------------------------------------------------------------------------
# Test: round-trip swap returns original tree
# ---------------------------------------------------------------------------

test_that("double swap is identity", {
  edge   <- .edge5()
  relBr  <- .relBr5()
  treeL  <- .treeL5()
  nTip   <- 5L

  # Swap nodes 7 and 9 (both internal, non-nested, non-sibling)
  res1 <- swap_subtrees_cpp(edge, nTip, treeL, relBr, nodeA = 7L, nodeB = 9L)
  expect_equal(res1$logHastings, 0.0)

  # Swap them back
  res2 <- swap_subtrees_cpp(res1$edge, nTip, treeL, res1$rel_br_lengths,
                             nodeA = 7L, nodeB = 9L)
  expect_equal(res2$logHastings, 0.0)

  # After double swap, edge and relBr should match original (up to postorder)
  # Match by sorting edge rows
  orig_sorted <- edge[order(edge[, 1L], edge[, 2L]), ]
  res2_sorted <- res2$edge[order(res2$edge[, 1L], res2$edge[, 2L]), ]
  expect_equal(orig_sorted, res2_sorted)

  orig_el  <- sort(treeL * relBr)
  res2_el  <- sort(treeL * res2$rel_br_lengths)
  expect_equal(orig_el, res2_el, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# Test: tip swap produces correct topology
# ---------------------------------------------------------------------------

test_that("swapping two tip nodes exchanges them in the topology", {
  edge  <- .edge5()
  relBr <- .relBr5()
  treeL <- .treeL5()
  nTip  <- 5L

  # Swap tip 1 (t1) and tip 3 (t3): result should have t3 where t1 was
  # and t1 where t3 was — reconstruct to check topology
  res <- swap_subtrees_cpp(edge, nTip, treeL, relBr, nodeA = 1L, nodeB = 3L)
  expect_equal(res$logHastings, 0.0)

  # Rebuild phylo and check that t1 and t3 changed places
  tr_orig <- .tree5()
  tr_swap <- tr_orig
  tr_swap$edge       <- res$edge
  tr_swap$edge.length <- treeL * res$rel_br_lengths

  # Verify: node 1's parent in original should equal node 3's parent in swapped
  orig_p1 <- edge[edge[, 2L] == 1L, 1L]  # parent of t1 in original
  swap_p3 <- res$edge[res$edge[, 2L] == 3L, 1L]  # parent of t3 after swap
  expect_equal(orig_p1, swap_p3)

  orig_p3 <- edge[edge[, 2L] == 3L, 1L]  # parent of t3 in original
  swap_p1 <- res$edge[res$edge[, 2L] == 1L, 1L]  # parent of t1 after swap
  expect_equal(orig_p3, swap_p1)
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

  rowA  <- which(edge[, 2L] == 7L)
  rowB  <- which(edge[, 2L] == 9L)
  origA <- relBr[rowA]
  origB <- relBr[rowB]

  res <- swap_subtrees_cpp(edge, 5L, treeL, relBr, nodeA = 7L, nodeB = 9L)

  # Find rows for nodes 7 and 9 in swapped result (postorder may change rows)
  new_rowA <- which(res$edge[, 2L] == 7L)
  new_rowB <- which(res$edge[, 2L] == 9L)

  expect_equal(res$rel_br_lengths[new_rowA], origB, tolerance = 1e-14)
  expect_equal(res$rel_br_lengths[new_rowB], origA, tolerance = 1e-14)
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
  # nTip=4, root=5, internal: 5, 6=(t1,t2), 7=(t3,t4)
  tr   <- Postorder(read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);"))
  edge <- tr$edge
  relBr <- tr$edge.length / sum(tr$edge.length)
  treeL <- sum(tr$edge.length)
  nTip  <- 4L

  # Swap tip 1 (t1) with tip 3 (t3)
  # Before: parent(t1)=6, parent(t3)=7
  # After:  parent(t1)=7, parent(t3)=6
  res <- swap_subtrees_cpp(edge, nTip, treeL, relBr, nodeA = 1L, nodeB = 3L)
  expect_equal(res$logHastings, 0.0)

  # Verify parents have swapped
  new_p1 <- res$edge[res$edge[, 2L] == 1L, 1L]
  new_p3 <- res$edge[res$edge[, 2L] == 3L, 1L]
  orig_p1 <- edge[edge[, 2L] == 1L, 1L]
  orig_p3 <- edge[edge[, 2L] == 3L, 1L]

  expect_equal(new_p1, orig_p3)
  expect_equal(new_p3, orig_p1)

  # Result should be ((t3,t2),(t1,t4)) in terms of tip groupings:
  # t3 and t2 share a parent, t1 and t4 share a parent
  new_sibling_of_t1 <- setdiff(res$edge[res$edge[,1] == new_p1, 2L], 1L)
  expect_equal(new_sibling_of_t1, 4L)

  new_sibling_of_t3 <- setdiff(res$edge[res$edge[,1] == new_p3, 2L], 3L)
  expect_equal(new_sibling_of_t3, 2L)
})
