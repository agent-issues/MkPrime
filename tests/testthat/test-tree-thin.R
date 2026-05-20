# Tests for treeThin: differential thinning of tree samples (M-134)

test_that("treeThin = NULL gives 1:1 tree-to-scalar samples (default)", {
  tree <- .mkp_test_tree()
  pd   <- .mkp_test_pd()
  set.seed(7281)
  result <- suppressWarnings(RunMkPrime(pd, tree, fixTopology = TRUE,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L,
                       maxWarmup = 200L, minWarmup = 200L,
                       autoTune = FALSE)))
  expect_equal(length(result$trees), nrow(result$samples))
})


test_that("treeThin > thin produces fewer tree samples", {
  tree <- .mkp_test_tree()
  pd   <- .mkp_test_pd()
  set.seed(4052)
  result <- suppressWarnings(RunMkPrime(pd, tree, fixTopology = TRUE,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, treeThin = 25L,
                       maxWarmup = 200L, minWarmup = 200L,
                       autoTune = FALSE)))
  nScalar <- nrow(result$samples)
  nTrees  <- length(result$trees)
  # Tree count should be 1/5 of scalar count (treeThin / thin = 5)
  expect_equal(nTrees, nScalar %/% 5L)
  expect_true(nTrees < nScalar)
  # All trees should be valid
  expect_true(all(vapply(result$trees, inherits, logical(1), "phylo")))
})


test_that("treeThin must be a multiple of thin", {
  tree <- .mkp_test_tree()
  pd   <- .mkp_test_pd()
  expect_error(
    suppressWarnings(RunMkPrime(pd, tree, fixTopology = TRUE,
      mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 200L, thin = 3L, treeThin = 7L,
                         maxWarmup = 50L, minWarmup = 50L,
                         autoTune = FALSE))),
    "multiple"
  )
})


test_that("MkPrimeMCMC rejects non-positive treeThin", {
  expect_error(MkPrimeMCMC(treeThin = 0L), "positive")
  expect_error(MkPrimeMCMC(treeThin = -1L), "positive")
})


test_that("treeFile line count matches tree sample count with treeThin", {
  tree <- .mkp_test_tree()
  pd   <- .mkp_test_pd()
  tf   <- tempfile(fileext = ".trees")
  on.exit(unlink(tf), add = TRUE)

  set.seed(6913)
  result <- suppressWarnings(RunMkPrime(pd, tree, fixTopology = TRUE,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, treeThin = 25L,
                       maxWarmup = 200L, minWarmup = 200L,
                       autoTune = FALSE, treeFile = tf)))

  lines <- readLines(tf)
  lines <- lines[nzchar(lines)]
  expect_equal(length(lines), length(result$trees))
  # Should be fewer than scalar samples
  expect_true(length(lines) < nrow(result$samples))
})


test_that("treeThin stored in result object", {
  tree <- .mkp_test_tree()
  pd   <- .mkp_test_pd()
  set.seed(3847)
  result <- suppressWarnings(RunMkPrime(pd, tree, fixTopology = TRUE,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, treeThin = 25L,
                       maxWarmup = 200L, minWarmup = 200L,
                       autoTune = FALSE)))
  expect_equal(result$treeThin, 25L)
})


# STREAM-003 regression: brColStart must account for the 2 diagnostic
# columns (swap_cold, topo_hash) that C++ writes between the last
# hyperparam column and the per-character kPrime block.  When brColStart
# was off-by-two, the reconstructed tree's first 1-2 edge lengths were
# taken from kPrime / topo_hash cells (one of them an FNV hash ~1e16) and
# the last actual br_* were dropped.  Trees still parsed via ape, but with
# garbage edge lengths.  See dev/red-team/findings.md (TREEMOVE / STREAM-003).
test_that("streamed trees have finite, non-negative, plausible edge lengths", {
  tree <- .mkp_test_tree()
  pd   <- .mkp_test_pd()
  set.seed(2026)
  result <- suppressWarnings(RunMkPrime(pd, tree, fixTopology = TRUE,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 200L, thin = 5L,
                       maxWarmup = 100L, minWarmup = 100L,
                       autoTune = FALSE)))
  expect_true(length(result$trees) > 0L)
  for (tr in result$trees) {
    expect_true(all(is.finite(tr$edge.length)))
    expect_true(all(tr$edge.length >= 0))
    # Sanity: edge lengths should sum to a reasonable tree length (well below
    # the FNV hash scale of ~1e16 that the brColStart bug would inject).
    expect_lt(sum(tr$edge.length), 1e8)
  }
})
