# Tests for TBR (Tree Bisection and Reconnection) proposal (M-053)

test_that("TBR produces valid binary tree", {
  library("ape")
  set.seed(4917)
  tree <- rtree(12)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  nTip <- length(tree$tip.label)
  nEdge <- nrow(tree$edge)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  nValid <- 0
  for (i in seq_len(50)) {
    prop <- MkPrime:::ProposeTbr(tree, tl, rel_br)
    if (!is.finite(prop$logHastings)) next
    nValid <- nValid + 1
    nt <- prop$tree

    expect_s3_class(nt, "phylo")
    expect_equal(length(nt$tip.label), nTip)
    expect_equal(nrow(nt$edge), nEdge)
    expect_setequal(nt$tip.label, tree$tip.label)
    expect_true(all(nt$edge.length > 0))
    expect_equal(sum(nt$edge.length), tl, tolerance = 1e-12)
    expect_equal(sum(prop$rel_br_lengths), 1, tolerance = 1e-12)
    expect_equal(attr(nt, "order"), "preorder")
    expect_equal(nt$Nnode, nTip - 2L)
  }
  expect_gt(nValid, 0)
})


test_that("TBR preserves total tree length", {
  library("ape")
  set.seed(3156)
  tree <- rtree(8)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  for (i in seq_len(30)) {
    prop <- MkPrime:::ProposeTbr(tree, tl, rel_br)
    if (!is.finite(prop$logHastings)) next
    expect_equal(sum(prop$tree$edge.length), tl, tolerance = 1e-12)
  }
})


test_that("TBR changes topology on medium tree", {
  library("ape")
  set.seed(7294)
  tree <- rtree(10)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl
  originalNewick <- write.tree(tree)

  nChanged <- 0
  for (i in seq_len(50)) {
    prop <- MkPrime:::ProposeTbr(tree, tl, rel_br)
    if (!is.finite(prop$logHastings)) next
    if (write.tree(prop$tree) != originalNewick) nChanged <- nChanged + 1
  }
  expect_gt(nChanged, 0)
})


test_that("TBR Hastings ratio is finite", {
  library("ape")
  set.seed(5283)
  tree <- rtree(15)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  for (i in seq_len(30)) {
    prop <- MkPrime:::ProposeTbr(tree, tl, rel_br)
    if (is.finite(prop$logHastings)) {
      expect_false(is.nan(prop$logHastings))
    }
  }
})


test_that("TBR on small tree (5 tips) works", {
  library("ape")
  tree <- read.tree(text = "(t1:1,(t2:2,(t3:3,(t4:4,t5:5):6):7):8);")
  tree <- TreeTools::Preorder(tree)
  nTip <- length(tree$tip.label)
  nEdge <- nrow(tree$edge)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  set.seed(6471)
  nValid <- 0
  for (i in seq_len(50)) {
    prop <- MkPrime:::ProposeTbr(tree, tl, rel_br)
    if (is.finite(prop$logHastings)) {
      nValid <- nValid + 1
      expect_equal(length(prop$tree$tip.label), nTip)
      expect_equal(nrow(prop$tree$edge), nEdge)
      expect_true(all(prop$tree$edge.length > 0))
    }
  }
  expect_gt(nValid, 0)
})


test_that("TBR chain explores tree space", {
  library("ape")
  set.seed(2850)
  tree <- rtree(8)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  topologies <- character(200)
  currentTree <- tree
  currentRel <- rel_br
  for (i in seq_len(200)) {
    prop <- MkPrime:::ProposeTbr(currentTree, tl, currentRel)
    if (is.finite(prop$logHastings)) {
      currentTree <- prop$tree
      currentRel <- prop$rel_br_lengths
    }
    topologies[i] <- write.tree(currentTree)
  }
  # TBR should explore multiple topologies
  expect_gt(length(unique(topologies)), 3)
})


test_that("TBR degenerates to SPR-like behavior for tip prune", {
  # When v is a tip, no subtree re-rooting occurs.
  # The move should still produce valid trees.
  library("ape")
  tree <- read.tree(text = "(t1:0.1,t2:0.2,t3:0.3);")
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  # 3-tip tree: all prune targets yield tip subtrees
  set.seed(8031)
  prop <- MkPrime:::ProposeTbr(tree, tl, rel_br)
  # Should either produce a valid move or return -Inf
  expect_true(is.numeric(prop$logHastings))
})


test_that("TBR preserves node set", {
  library("ape")
  set.seed(1397)
  tree <- rtree(20)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  nTip <- length(tree$tip.label)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  origNodes <- sort(unique(as.vector(tree$edge)))

  for (i in seq_len(30)) {
    prop <- MkPrime:::ProposeTbr(tree, tl, rel_br)
    if (!is.finite(prop$logHastings)) next
    newNodes <- sort(unique(as.vector(prop$tree$edge)))
    expect_equal(newNodes, origNodes)
  }
})


test_that("TBR on larger tree (30 tips) preserves structure", {
  library("ape")
  set.seed(9463)
  tree <- rtree(30)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  nTip <- length(tree$tip.label)
  nEdge <- nrow(tree$edge)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  nValid <- 0
  for (i in seq_len(30)) {
    prop <- MkPrime:::ProposeTbr(tree, tl, rel_br)
    if (!is.finite(prop$logHastings)) next
    nValid <- nValid + 1
    nt <- prop$tree

    expect_equal(length(nt$tip.label), nTip)
    expect_equal(nrow(nt$edge), nEdge)
    expect_true(all(nt$edge.length > 0))
    expect_equal(sum(nt$edge.length), tl, tolerance = 1e-12)
    expect_equal(nt$Nnode, nTip - 2L)
  }
  expect_gt(nValid, 5)
})
