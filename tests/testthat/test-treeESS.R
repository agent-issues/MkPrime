# Tests for tree-topology ESS (cherry-picked from treess)
#
# Reference values computed with treess v1.0.1 (afmagee/treess, GPL-3)
# using phangorn::RF.dist on simulated MCMC chains.

# ---- Fixtures ----

# 15x15 RF distance matrix from a simulated 5-taxon NNI chain (seed 7762)
.dmat_tiny <- matrix(c(
  0, 2, 4, 2, 0, 2, 0, 2, 2, 4, 4, 4, 4, 4, 4,
  2, 0, 2, 0, 2, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 2, 0, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  2, 0, 2, 0, 2, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
  0, 2, 4, 2, 0, 2, 0, 2, 2, 4, 4, 4, 4, 4, 4,
  2, 4, 4, 4, 2, 0, 2, 0, 2, 4, 4, 2, 4, 2, 4,
  0, 2, 4, 2, 0, 2, 0, 2, 2, 4, 4, 4, 4, 4, 4,
  2, 4, 4, 4, 2, 0, 2, 0, 2, 4, 4, 2, 4, 2, 4,
  2, 4, 4, 4, 2, 2, 2, 2, 0, 2, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 2, 0, 2, 4, 2, 4, 2,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 0, 2, 0, 2, 0,
  4, 4, 4, 4, 4, 2, 4, 2, 4, 4, 2, 0, 2, 0, 2,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 0, 2, 0, 2, 0,
  4, 4, 4, 4, 4, 2, 4, 2, 4, 4, 2, 0, 2, 0, 2,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 0, 2, 0, 2, 0
), nrow = 15, byrow = TRUE)

# 30x30 RF distance matrix from a simulated 5-taxon NNI chain (seed 4209)
.dmat_30 <- structure(c(0, 2, 0, 2, 2, 4, 4, 4, 4, 4, 2, 2, 2, 4, 4, 4, 4,
4, 2, 2, 0, 2, 2, 4, 4, 2, 0, 2, 4, 4, 2, 0, 2, 4, 4, 4, 2, 2,
4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 0, 2, 4, 4, 4, 2, 4, 4,
4, 0, 2, 0, 2, 2, 4, 4, 4, 4, 4, 2, 2, 2, 4, 4, 4, 4, 4, 2, 2,
0, 2, 2, 4, 4, 2, 0, 2, 4, 4, 2, 4, 2, 0, 2, 4, 4, 4, 4, 2, 0,
2, 0, 2, 4, 4, 4, 4, 2, 0, 2, 4, 4, 4, 2, 0, 2, 2, 4, 4, 2, 4,
2, 2, 0, 2, 4, 4, 4, 4, 2, 0, 2, 4, 4, 4, 2, 2, 0, 2, 2, 4, 4,
4, 4, 2, 2, 0, 2, 4, 4, 4, 4, 4, 2, 0, 2, 4, 4, 4, 4, 2, 4, 4,
4, 4, 2, 0, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 4, 4, 2, 4, 4, 4,
2, 0, 2, 4, 4, 4, 4, 4, 4, 2, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4,
4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 2, 0, 2, 2, 4, 4, 4, 2, 0, 2, 4,
4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2,
0, 2, 4, 4, 4, 2, 2, 0, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2,
2, 4, 4, 4, 2, 4, 4, 4, 2, 2, 0, 2, 4, 2, 0, 2, 2, 4, 4, 4, 2,
4, 4, 4, 4, 2, 2, 4, 4, 4, 4, 2, 4, 2, 0, 2, 4, 4, 4, 4, 2, 0,
2, 0, 2, 4, 4, 4, 4, 2, 0, 2, 4, 4, 4, 2, 0, 2, 2, 4, 4, 2, 4,
2, 2, 0, 2, 4, 4, 4, 4, 2, 0, 2, 4, 4, 4, 2, 2, 0, 2, 2, 4, 4,
4, 4, 2, 2, 0, 2, 4, 2, 4, 2, 0, 2, 4, 4, 4, 4, 2, 0, 2, 0, 2,
4, 4, 4, 4, 2, 0, 2, 4, 4, 4, 2, 0, 2, 2, 4, 4, 4, 4, 4, 2, 4,
4, 4, 2, 2, 0, 2, 4, 2, 0, 2, 2, 4, 4, 4, 2, 4, 4, 4, 4, 2, 2,
4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 2, 0, 2, 2, 4, 4, 4, 2, 0, 2, 4,
4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2,
0, 2, 4, 4, 4, 2, 2, 0, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2,
2, 4, 4, 4, 4, 2, 2, 4, 4, 2, 4, 4, 2, 4, 4, 4, 2, 0, 2, 2, 4,
4, 4, 4, 4, 4, 4, 4, 2, 0, 2, 4, 4, 4, 4, 2, 0, 2, 4, 4, 4, 4,
2, 4, 4, 4, 4, 2, 0, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 4, 2, 4,
2, 2, 0, 2, 4, 4, 4, 4, 2, 0, 2, 4, 4, 4, 2, 2, 0, 2, 2, 4, 4,
4, 4, 2, 2, 0, 2, 4, 2, 4, 2, 0, 2, 4, 4, 4, 4, 2, 0, 2, 0, 2,
4, 4, 4, 4, 2, 0, 2, 4, 4, 4, 2, 0, 2, 2, 4, 4, 0, 2, 0, 2, 2,
4, 4, 4, 4, 4, 2, 2, 2, 4, 4, 4, 4, 4, 2, 2, 0, 2, 2, 4, 4, 2,
0, 2, 4, 4, 2, 0, 2, 4, 4, 4, 2, 2, 4, 4, 4, 4, 4, 4, 2, 4, 4,
4, 4, 4, 2, 0, 2, 4, 4, 4, 2, 4, 4, 4, 2, 2, 2, 4, 4, 4, 4, 4,
4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 0, 2, 4, 4, 2, 4, 4,
2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
4, 4, 2, 0, 2, 4, 4, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 2,
4, 2, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 2, 0, 2, 4, 4, 4, 4, 2, 4,
2, 0, 2, 4, 4, 4, 4, 2, 0, 2, 0, 2, 4, 4, 4, 4, 2, 0, 2, 4, 4,
4, 2, 0, 2, 2, 4, 4, 0, 2, 0, 2, 2, 4, 4, 4, 4, 4, 2, 2, 2, 4,
4, 4, 4, 4, 2, 2, 0, 2, 2, 4, 4, 2, 0, 2, 4, 4, 2, 4, 2, 2, 0,
2, 4, 4, 4, 4, 2, 0, 2, 4, 4, 4, 2, 2, 0, 2, 2, 4, 4, 4, 4, 2,
2, 0, 2, 4, 4, 4, 4, 4, 2, 2, 4, 4, 2, 4, 4, 2, 4, 4, 4, 2, 0,
2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2, 0, 2, 4, 4, 4, 4, 4, 4, 4, 4,
2, 4, 4, 4, 4, 4, 4, 2, 2, 4, 4, 4, 4, 4, 2, 2, 4, 4, 4, 4, 2,
0), dim = c(30L, 30L))


# ---- Tests: FrechetCorrelationESS ----

test_that(".FrechetCorrelationESS matches treess reference (tiny)", {
  # treess v1.0.1 reference: 4.548441
  expect_equal(.FrechetCorrelationESS(.dmat_tiny, 5L), 4.548441,
               tolerance = 1e-4)
})

test_that(".FrechetCorrelationESS matches treess reference (n=30)", {
  # treess v1.0.1 reference: 13.17867
  expect_equal(.FrechetCorrelationESS(.dmat_30, 5L), 13.17867,
               tolerance = 1e-4)
})

test_that(".FrechetCorrelationESS returns 1 for constant chain", {
  dmat_const <- matrix(0, 20, 20)
  expect_equal(.FrechetCorrelationESS(dmat_const, 5L), 1)
})

test_that(".FrechetCorrelationESS returns NA for too-short chains", {
  dmat_short <- matrix(0, 5, 5)
  expect_true(is.na(.FrechetCorrelationESS(dmat_short, 5L)))
})

# ---- Tests: MedianPseudoESS ----
#
# .MedianPseudoESS now uses a C++ Geyer (1992) initial-monotone-sequence
# estimator (same class as Stan/rstan) instead of coda::effectiveSize
# (AR spectral).  Reference values differ from treess v1.0.1 which uses
# coda internally; the Geyer estimator is generally more conservative.

test_that(".MedianPseudoESS returns Geyer ESS (tiny)", {
  # Geyer initial-monotone-sequence reference: 5.6271
  expect_equal(.MedianPseudoESS(.dmat_tiny), 5.6271, tolerance = 1e-3)
})

test_that(".MedianPseudoESS returns Geyer ESS (n=30)", {
  # Geyer initial-monotone-sequence reference: 16.678
  expect_equal(.MedianPseudoESS(.dmat_30), 16.678, tolerance = 1e-2)
})

test_that(".MedianPseudoESS returns NA for a constant chain (#195)", {
  dmat_const <- matrix(0, 20, 20)
  expect_true(is.na(.MedianPseudoESS(dmat_const)))
})

test_that(".MedianPseudoESS returns NA for too-short chain", {
  dmat_short <- matrix(0, 5, 5)
  expect_true(is.na(.MedianPseudoESS(dmat_short)))
})

# ---- Tests: TreeESS wrapper ----

test_that("TreeESS default returns median only (cross-distance fast path)", {
  skip_if_not_installed("TreeDist")
  skip_if_not_installed("ape")

  tr1 <- read.tree(text = "((a,b),(c,(d,e)));")
  tr2 <- read.tree(text = "(((a,b),c),(d,e));")
  tr3 <- read.tree(text = "((a,(b,c)),(d,e));")
  trees <- c(tr1, tr2, tr3, tr1, tr2, tr3, tr1, tr2, tr3, tr1)

  result <- TreeESS(trees)
  expect_named(result, c("frechetCorrelationESS", "medianPseudoESS"))
  expect_true(is.na(result[["frechetCorrelationESS"]]))
  expect_true(is.finite(result[["medianPseudoESS"]]))
  expect_true(result[["medianPseudoESS"]] > 0)
})

test_that("TreeESS with frechet = TRUE computes both methods", {
  skip_if_not_installed("TreeDist")
  skip_if_not_installed("ape")

  tr1 <- read.tree(text = "((a,b),(c,(d,e)));")
  tr2 <- read.tree(text = "(((a,b),c),(d,e));")
  tr3 <- read.tree(text = "((a,(b,c)),(d,e));")
  trees <- c(tr1, tr2, tr3, tr1, tr2, tr3, tr1, tr2, tr3, tr1)

  result <- TreeESS(trees, frechet = TRUE)
  expect_named(result, c("frechetCorrelationESS", "medianPseudoESS"))
  expect_true(all(is.finite(result)))
  expect_true(all(result > 0))
})

test_that("median pseudo-ESS agrees between fast path and full matrix", {
  skip_if_not_installed("TreeDist")
  skip_if_not_installed("ape")

  tr1 <- read.tree(text = "((a,b),(c,(d,e)));")
  tr2 <- read.tree(text = "(((a,b),c),(d,e));")
  tr3 <- read.tree(text = "((a,(b,c)),(d,e));")
  trees <- c(tr1, tr2, tr3, tr1, tr2, tr3, tr1, tr2, tr3, tr1)

  # n=10, maxRows=200 > 10, so both paths use all rows
  fast <- TreeESS(trees, frechet = FALSE)
  full <- TreeESS(trees, frechet = TRUE)
  expect_equal(fast[["medianPseudoESS"]], full[["medianPseudoESS"]])
})

test_that("TreeESS returns finite positive number on rtree(8) posterior", {
  skip_if_not_installed("TreeDist")
  skip_if_not_installed("ape")

  set.seed(42)
  trees <- lapply(seq_len(10), function(i) ape::rtree(8))
  class(trees) <- "multiPhylo"

  result <- TreeESS(trees)
  expect_true(is.finite(result[["medianPseudoESS"]]))
  expect_true(result[["medianPseudoESS"]] > 0)
})

test_that("median pseudo-ESS is NA, not n, for a single topology (#195)", {
  skip_if_not_installed("TreeDist")
  tree <- as.phylo(0, 8)
  trees <- structure(rep(list(tree), 50), class = "multiPhylo")
  expect_true(is.na(TreeESS(trees)[["medianPseudoESS"]]))
})
