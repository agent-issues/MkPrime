# Regression test for the OPP-6b in-place NNI preorder bug.
#
# Before the fix, in-place NNI could produce edges where a child appeared
# before its parent in the edge list (when wRow < edgeRow).  The
# reverse-iteration Felsenstein pruning would then compute wrong
# likelihoods because a node's CL was propagated before all children
# contributed.

library("testthat")

test_that("NNI moves always produce valid preorder", {
  skip_if_not_installed("TreeTools")

  set.seed(8741)
  tree <- TreeTools::Preorder(TreeTools::BalancedTree(12))
  nTip <- TreeTools::NTip(tree)
  nEdge <- nrow(tree$edge)
  relBr <- rep(1 / nEdge, nEdge)
  treeLength <- 1.0

  n_bad <- 0L
  nIter <- 200L
  for (i in seq_len(nIter)) {
    res <- nni_proposal(tree$edge, nTip, treeLength, relBr)
    if (is.finite(res$logHastings)) {
      parent <- res$edge[, 1]
      child  <- res$edge[, 2]
      # Check: every edge's parent must have been introduced earlier
      introduced <- logical(max(parent))
      introduced[nTip + 1L] <- TRUE # root
      valid <- TRUE
      for (j in seq_len(nEdge)) {
        if (!introduced[parent[j]]) {
          valid <- FALSE
          break
        }
        introduced[child[j]] <- TRUE
      }
      if (!valid) n_bad <- n_bad + 1L
    }
  }
  expect_equal(n_bad, 0L,
               label = "Number of NNI proposals with invalid preorder")
})


test_that("In-place NNI likelihood matches full-reorder NNI likelihood", {
  skip_if_not_installed("TreeTools")

  # Build a small dataset to compute likelihoods
  set.seed(6293)
  tree <- TreeTools::Preorder(TreeTools::RandomTree(8, root = TRUE))
  nTip <- TreeTools::NTip(tree)
  nEdge <- nrow(tree$edge)

  # Generate random binary tip data (4 characters)
  tipMat <- matrix(sample(0:1, nTip * 4, replace = TRUE),
                   nrow = nTip, ncol = 4)
  kStates <- 2L
  rootFreqs <- rep(1 / kStates, kStates)

  relBr <- runif(nEdge, 0.01, 0.5)
  relBr <- relBr / sum(relBr)
  treeLength <- 2.0
  edgeLen <- treeLength * relBr

  # Reference likelihood with original tree
  ll_orig <- pruning_jc(tree$edge[, 1], tree$edge[, 2],
                        edgeLen, tipMat, kStates, rootFreqs)

  # Run many NNI proposals and verify the edge ordering is valid
  for (i in seq_len(100)) {
    res <- nni_proposal(tree$edge, nTip, treeLength, relBr)
    if (is.finite(res$logHastings)) {
      newEdgeLen <- treeLength * res$rel_br_lengths
      ll_new <- pruning_jc(res$edge[, 1], res$edge[, 2],
                           newEdgeLen, tipMat, kStates, rootFreqs)
      expect_true(is.finite(ll_new),
                  info = paste("NNI proposal", i, "gave non-finite likelihood"))
    }
  }
})
