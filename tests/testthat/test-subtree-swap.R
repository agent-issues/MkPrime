# Tests for get_valid_swap_partners (M-084)
#
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

