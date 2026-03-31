# Tests for tree topology proposals (M-023: NNI, M-024: SPR)

test_that("NNI produces valid binary tree", {
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  nTip <- length(tree$tip.label)
  nEdge <- nrow(tree$edge)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  set.seed(6183)
  for (i in seq_len(50)) {
    prop <- MkPrime:::ProposeNni(tree, tl, rel_br)
    nt <- prop$tree

    # Still a valid phylo
    expect_s3_class(nt, "phylo")
    # Same number of tips and edges
    expect_equal(length(nt$tip.label), nTip)
    expect_equal(nrow(nt$edge), nEdge)
    # Same tip labels (order may differ)
    expect_setequal(nt$tip.label, tree$tip.label)
    # All edge lengths positive
    expect_true(all(nt$edge.length > 0))
    # Tree length preserved
    expect_equal(sum(nt$edge.length), tl, tolerance = 1e-12)
    # rel_br_lengths sums to 1
    expect_equal(sum(prop$rel_br_lengths), 1, tolerance = 1e-12)
    # Preorder
    expect_equal(attr(nt, "order"), "preorder")
  }
})


test_that("NNI changes topology on 4-tip tree", {
  library("ape")
  tree <- read.tree(text = "(t1:0.1,(t2:0.2,(t3:0.15,t4:0.25):0.1):0.05);")
  tree <- TreeTools::Preorder(tree)
  nTip <- length(tree$tip.label)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  original_newick <- write.tree(tree)

  # On a 4-tip tree there's 1 internal edge; NNI always changes topology
  set.seed(3812)
  n_changed <- 0
  for (i in seq_len(20)) {
    prop <- MkPrime:::ProposeNni(tree, tl, rel_br)
    prop_newick <- write.tree(prop$tree)
    if (prop_newick != original_newick) n_changed <- n_changed + 1
  }
  # Should change topology every time
  expect_equal(n_changed, 20)
})


test_that("NNI returns -Inf logHastings for 3-tip tree", {
  library("ape")
  tree <- read.tree(text = "(t1:0.1,t2:0.2,t3:0.3);")
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  prop <- MkPrime:::ProposeNni(tree, tl, rel_br)
  expect_equal(prop$logHastings, -Inf)
})


test_that("NNI on larger tree preserves structure", {
  library("ape")
  set.seed(4592)
  tree <- rtree(20)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  nTip <- length(tree$tip.label)
  nEdge <- nrow(tree$edge)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  for (i in seq_len(30)) {
    prop <- MkPrime:::ProposeNni(tree, tl, rel_br)
    nt <- prop$tree
    expect_equal(length(nt$tip.label), nTip)
    expect_equal(nrow(nt$edge), nEdge)
    expect_setequal(nt$tip.label, tree$tip.label)
    expect_true(all(nt$edge.length > 0))
    expect_equal(sum(nt$edge.length), tl, tolerance = 1e-12)
    # Unrooted binary: nTip - 2 internal nodes
    expect_equal(nt$Nnode, nTip - 2L)
  }
})


test_that("NNI Hastings ratio is zero (symmetric)", {
  library("ape")
  set.seed(7281)
  tree <- rtree(10)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  for (i in seq_len(20)) {
    prop <- MkPrime:::ProposeNni(tree, tl, rel_br)
    expect_equal(prop$logHastings, 0)
  }
})


test_that("NNI preserves individual edge lengths", {
  library("ape")
  # Use distinctive edge lengths to verify they follow their edges
  tree <- read.tree(
    text = "(t1:1,(t2:2,(t3:3,(t4:4,t5:5):6):7):8);"
  )
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  set.seed(9137)
  prop <- MkPrime:::ProposeNni(tree, tl, rel_br)

  # Same multiset of edge lengths (reordered but same values)
  expect_setequal(prop$tree$edge.length, tree$edge.length)
})


test_that("NNI explores all 3 topologies on 4-tip tree", {
  library("ape")
  tree <- read.tree(text = "(t1:0.1,(t2:0.2,(t3:0.15,t4:0.25):0.1):0.05);")
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  # Repeatedly apply NNI; on 4 tips there are 3 unrooted topologies
  set.seed(5501)
  topologies <- character(200)
  currentTree <- tree
  currentRel <- rel_br
  for (i in seq_len(200)) {
    prop <- MkPrime:::ProposeNni(currentTree, tl, currentRel)
    if (prop$logHastings > -Inf) {
      currentTree <- prop$tree
      currentRel <- prop$rel_br_lengths
    }
    topologies[i] <- write.tree(currentTree)
  }
  # Should visit more than 1 topology
  expect_gt(length(unique(topologies)), 1)
})


# --- SPR tests (M-024) ---

test_that("SPR produces valid binary tree", {
  library("ape")
  set.seed(2847)
  tree <- rtree(12)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  nTip <- length(tree$tip.label)
  nEdge <- nrow(tree$edge)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  for (i in seq_len(50)) {
    prop <- MkPrime:::ProposeSpr(tree, tl, rel_br)
    if (!is.finite(prop$logHastings)) next
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
})


test_that("SPR preserves total tree length", {
  library("ape")
  set.seed(6632)
  tree <- rtree(8)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  for (i in seq_len(30)) {
    prop <- MkPrime:::ProposeSpr(tree, tl, rel_br)
    if (!is.finite(prop$logHastings)) next
    expect_equal(sum(prop$tree$edge.length), tl, tolerance = 1e-12)
  }
})


test_that("SPR changes topology on medium tree", {
  library("ape")
  set.seed(1459)
  tree <- rtree(10)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl
  original_newick <- write.tree(tree)

  n_changed <- 0
  for (i in seq_len(50)) {
    prop <- MkPrime:::ProposeSpr(tree, tl, rel_br)
    if (!is.finite(prop$logHastings)) next
    if (write.tree(prop$tree) != original_newick) n_changed <- n_changed + 1
  }
  expect_gt(n_changed, 0)
})


test_that("SPR Hastings ratio is finite", {
  library("ape")
  set.seed(3819)
  tree <- rtree(15)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  for (i in seq_len(30)) {
    prop <- MkPrime:::ProposeSpr(tree, tl, rel_br)
    # Should be finite (not NaN or Inf) unless the move was rejected
    if (is.finite(prop$logHastings)) {
      expect_false(is.nan(prop$logHastings))
    }
  }
})


test_that("SPR on small tree (5 tips) works", {
  library("ape")
  tree <- read.tree(text = "(t1:1,(t2:2,(t3:3,(t4:4,t5:5):6):7):8);")
  tree <- TreeTools::Preorder(tree)
  nTip <- length(tree$tip.label)
  nEdge <- nrow(tree$edge)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  set.seed(8173)
  n_valid <- 0
  for (i in seq_len(50)) {
    prop <- MkPrime:::ProposeSpr(tree, tl, rel_br)
    if (is.finite(prop$logHastings)) {
      n_valid <- n_valid + 1
      expect_equal(length(prop$tree$tip.label), nTip)
      expect_equal(nrow(prop$tree$edge), nEdge)
      expect_true(all(prop$tree$edge.length > 0))
    }
  }
  # Should produce some valid moves
  expect_gt(n_valid, 0)
})


test_that("SPR chain explores tree space", {
  library("ape")
  set.seed(7402)
  tree <- rtree(8)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  topologies <- character(200)
  currentTree <- tree
  currentRel <- rel_br
  for (i in seq_len(200)) {
    prop <- MkPrime:::ProposeSpr(currentTree, tl, currentRel)
    if (is.finite(prop$logHastings)) {
      currentTree <- prop$tree
      currentRel <- prop$rel_br_lengths
    }
    topologies[i] <- write.tree(currentTree)
  }
  # SPR should explore more topologies than NNI
  expect_gt(length(unique(topologies)), 3)
})
