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

test_that(".MedianPseudoESS matches treess reference (tiny)", {
  skip_if_not_installed("coda")
  # treess v1.0.1 reference: 5.520146
  expect_equal(.MedianPseudoESS(.dmat_tiny), 5.520146, tolerance = 1e-4)
})

test_that(".MedianPseudoESS matches treess reference (n=30)", {
  skip_if_not_installed("coda")
  # treess v1.0.1 reference: 29.31821
  expect_equal(.MedianPseudoESS(.dmat_30), 29.31821, tolerance = 1e-4)
})

# ---- Tests: .TreeESS wrapper ----

test_that(".TreeESS returns named vector with both methods", {
  skip_if_not_installed("coda")
  skip_if_not_installed("TreeDist")
  skip_if_not_installed("ape")

  library(ape)
  # Build a small set of trees by hand
  tr1 <- read.tree(text = "((a,b),(c,(d,e)));")
  tr2 <- read.tree(text = "(((a,b),c),(d,e));")
  tr3 <- read.tree(text = "((a,(b,c)),(d,e));")
  trees <- c(tr1, tr2, tr3, tr1, tr2, tr3, tr1, tr2, tr3, tr1)

  result <- .TreeESS(trees)
  expect_named(result, c("frechetCorrelationESS", "medianPseudoESS"))
  expect_true(all(is.finite(result)))
  expect_true(all(result > 0))
})
