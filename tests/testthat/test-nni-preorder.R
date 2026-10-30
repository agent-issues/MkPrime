# Regression test for the OPP-6b in-place NNI preorder bug.
#
# Before the fix, in-place NNI could produce edges where a child appeared
# before its parent in the edge list (when wRow < edgeRow).  The
# reverse-iteration Felsenstein pruning would then compute wrong
# likelihoods because a node's CL was propagated before all children
# contributed.

test_that("NNI moves always produce valid preorder", {
  skip_if_not_installed("TreeTools")

  set.seed(8741)
  tree <- Preorder(BalancedTree(12))
  nTip <- NTip(tree)
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

test_that("NNI proposal edge order gives the reordered tree's likelihood", {
  skip_if_not_installed("TreeTools")

  set.seed(6293)
  tree <- Preorder(RandomTree(8, root = TRUE))
  nTip <- NTip(tree)
  nEdge <- nrow(tree$edge)

  tipMat <- matrix(sample(0:1, nTip * 4, replace = TRUE),
                   nrow = nTip, ncol = 4)
  kStates <- 2L
  rootFreqs <- rep(1 / kStates, kStates)

  relBr <- runif(nEdge, 0.01, 0.5)
  relBr <- relBr / sum(relBr)
  treeLength <- 2.0

  nChecked <- 0L
  for (i in seq_len(100)) {
    res <- nni_proposal(tree$edge, nTip, treeLength, relBr)
    if (!is.finite(res$logHastings)) next
    newEdgeLen <- treeLength * res$rel_br_lengths
    llNew <- pruning_jc(res$edge[, 1], res$edge[, 2],
                         newEdgeLen, tipMat, kStates, rootFreqs)
    ref <- tree
    ref$edge <- res$edge
    ref$edge.length <- newEdgeLen
    attr(ref, "order") <- NULL
    ref <- Preorder(ref)
    llRef <- pruning_jc(ref$edge[, 1], ref$edge[, 2],
                         ref$edge.length, tipMat, kStates, rootFreqs)
    expect_equal(llNew, llRef, info = paste("NNI proposal", i))
    nChecked <- nChecked + 1L
  }
  expect_gt(nChecked, 0L)
})
