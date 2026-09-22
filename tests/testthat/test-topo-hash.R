# TREEMOVE-002: the `topo_hash` diagnostic column fingerprints the topology.
#
# Two ways of getting that wrong meet here.  Canonical relabelling makes the
# `parent` vector a function of the tree SHAPE alone, so hashing it collides
# every topology sharing a shape; and the in-place NNI branch leaves the edge
# list in a valid but non-canonical order, so hashing the edge list splits one
# topology across several values.  The multiset of non-trivial splits is
# immune to both.

.Hash <- function(tree) {
  tree <- TreeTools::Preorder(tree)
  compute_topo_hash(tree$edge[, 1], tree$edge[, 2], TreeTools::NTip(tree))
}


test_that("topo_hash separates topologies that share a shape", {
  lab <- paste0("t", 1:5)
  a <- TreeTools::RenumberTips(
    ape::read.tree(text = "((t1,t2),(t3,t4),t5);"), lab)
  b <- TreeTools::RenumberTips(
    ape::read.tree(text = "((t1,t3),(t2,t4),t5);"), lab)

  # These two share a canonical parent vector, so a hash of it cannot tell
  # them apart.
  expect_identical(TreeTools::Preorder(a)$edge[, 1],
                   TreeTools::Preorder(b)$edge[, 1])
  expect_false(.Hash(a) == .Hash(b))
})


test_that("topo_hash is invariant to edge order and node labelling", {
  set.seed(91L)
  tree <- TreeTools::Preorder(ape::rtree(8L, rooted = FALSE))
  h <- .Hash(tree)

  # An in-place NNI keeps parents ahead of children but abandons canonical
  # preorder; the hash must survive that, and any relabelling of the
  # internal nodes.
  shuffled <- tree
  perm <- c(seq_len(8L), 8L + sample(seq_len(tree$Nnode)))
  shuffled$edge[] <- perm[tree$edge]
  expect_identical(.Hash(shuffled), h)

  reordered <- tree
  ord <- order(reordered$edge[, 1])
  reordered$edge <- reordered$edge[ord, , drop = FALSE]
  reordered$edge.length <- reordered$edge.length[ord]
  expect_identical(compute_topo_hash(reordered$edge[, 1],
                                     reordered$edge[, 2], 8L), h)
})


test_that("topo_hash ignores where a two-child root sits", {
  # A rooted input tree stays rooted in the state, and the moves carry the
  # root about the topology; each position adds a trivial or duplicate split.
  set.seed(2L)
  tree <- TreeTools::RandomTree(10L, root = FALSE)
  rootings <- list(TreeTools::RootTree(tree, 1L), TreeTools::RootTree(tree, 5L),
                   TreeTools::RootTree(tree, c(2L, 3L)))
  expect_identical(vapply(rootings, .Hash, double(1)),
                   rep(.Hash(tree), length(rootings)))
})


test_that("topo_hash counts the topologies a chain actually visits", {
  skip_slow_tests()
  set.seed(11L)
  nTip <- 8L
  tree <- TreeTools::Preorder(TreeTools::UnrootTree(
    ape::rtree(nTip, tip.label = paste0("t", seq_len(nTip)))))
  mat <- matrix(sample(0:1, nTip * 40L, replace = TRUE), nrow = nTip,
                dimnames = list(tree$tip.label, NULL))
  mcmc <- MkPrimeMCMC(nRuns = 1L, nIter = 1000L, thin = 2L, maxWarmup = 200L,
                      minWarmup = 200L, autoTune = FALSE, maxTime = 120)
  setTimeLimit(elapsed = 400, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  set.seed(3L)
  res <- suppressWarnings(RunMkPrime(MatrixToPhyDat(mat), tree, mcmc = mcmc))

  # `res$samples` is thinned more finely than `res$trees`, so judge the hash
  # on the saved trees, against Robinson-Foulds distance as the reference.
  lab <- res$trees[[1L]]$tip.label
  trees <- lapply(res$trees, function(tr) {
    tr$Nnode <- nrow(tr$edge) - TreeTools::NTip(tr) + 1L
    TreeTools::RenumberTips(tr, lab)
  })
  hashes <- vapply(trees, .Hash, numeric(1L))
  expect_gt(length(unique(hashes)), 20L)

  # No two distinct topologies share a hash.
  reps <- trees[!duplicated(hashes)]
  expect_true(all(ape::dist.topo(structure(reps, class = "multiPhylo"),
                                 method = "PH85") > 0))

  # No topology is split across hashes.
  for (h in unique(hashes)) {
    inClass <- trees[hashes == h]
    if (length(inClass) < 2L) next
    expect_identical(
      max(ape::dist.topo(structure(head(inClass, 6L), class = "multiPhylo"),
                         method = "PH85")), 0)
  }

  expect_true("topo_hash" %in% colnames(res$samples))
})
