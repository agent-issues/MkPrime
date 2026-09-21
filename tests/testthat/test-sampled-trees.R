# Sampled trees must be valid phylo objects under either rooting convention.
# The state keeps the input tree's rooting, so an unrooted input yields
# 2n - 3 edges and a rooted one 2n - 2; Nnode has to follow.

.ExpectConsistentNnode <- function(tree) {
  expect_equal(tree$Nnode, nrow(tree$edge) - length(tree$tip.label) + 1L)
}

.ExpectNewickRoundTrip <- function(tree) {
  back <- ape::read.tree(text = ape::write.tree(tree))
  back <- TreeTools::Preorder(TreeTools::RenumberTips(back, tree$tip.label))
  expect_equal(back$Nnode, tree$Nnode)
  expect_equal(back$edge, tree$edge)
  expect_equal(back$edge.length, tree$edge.length, tolerance = 1e-6)
}

.SampledTreeFixture <- function(rooted) {
  set.seed(4417)
  nTip <- 7L
  tree <- TreeTools::RandomTree(nTip, root = rooted)
  tree$edge.length <- rep(0.2, nrow(tree$edge))
  mat <- matrix(sample(0:2, nTip * 12L, replace = TRUE), nTip,
                dimnames = list(tree$tip.label, NULL))
  list(tree = tree, pd = TreeTools::MatrixToPhyDat(mat))
}

test_that(".StateToTree reads Nnode from the edge matrix", {
  for (rooted in c(FALSE, TRUE)) {
    fx <- .SampledTreeFixture(rooted)
    tree <- TreeTools::Preorder(fx$tree)
    mkd <- MkPrimeData(fx$pd)
    model <- MkPrime:::.FinalizeModel(MkPrimeModel(), tree, mkd)
    statePtr <- MkPrime:::.InitMcmcChain(MkPrime:::.InitState(tree, mkd, model))

    out <- MkPrime:::.StateToTree(statePtr, tree$tip.label)
    expect_equal(out$Nnode, tree$Nnode)
    .ExpectConsistentNnode(out)
  }
})

test_that("RunMkPrime returns and writes well-formed trees from an unrooted start", {
  fx <- .SampledTreeFixture(rooted = FALSE)
  expect_equal(fx$tree$Nnode, length(fx$tree$tip.label) - 2L)
  treeFile <- tempfile(fileext = ".nwk")
  on.exit(unlink(treeFile))

  set.seed(4418)
  res <- suppressWarnings(RunMkPrime(
    fx$pd, fx$tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 300L, thin = 10L,
                       minWarmup = 50L, maxWarmup = 50L, autoTune = FALSE,
                       treeFile = treeFile)
  ))

  expect_gt(length(res$trees), 0L)
  for (tr in res$trees) {
    .ExpectConsistentNnode(tr)
    .ExpectNewickRoundTrip(tr)
  }

  written <- ape::read.tree(treeFile)
  expect_length(written, length(res$trees))
  expect_equal(ape::write.tree(written), ape::write.tree(res$trees))
})
