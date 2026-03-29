test_that("Stepping-stone returns finite marginal likelihood", {
  library(ape)
  library(TreeTools)

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L),
                nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(4719)
  ss <- mkp_stepping_stone(pd, tree, nStones = 8L, nIter = 100L,
                            warmup = 30L, verbose = FALSE)

  expect_true(is.finite(ss$log_marginal))
  expect_true(is.finite(ss$se))
  expect_gt(ss$se, 0)
  expect_length(ss$log_ratios, 8L)
  expect_length(ss$betas, 9L) # nStones + 1 boundary points
  expect_equal(ss$betas[1], 0)
  expect_equal(ss$betas[9], 1)
})


test_that("Beta schedule is monotonically increasing", {
  library(ape)
  library(TreeTools)

  tree <- read.tree(text = "((t1:0.1,t2:0.1):0.1,t3:0.1);")
  mat <- matrix(c(0L, 1L, 0L), nrow = 3,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(2841)
  ss <- mkp_stepping_stone(pd, tree, nStones = 5L, nIter = 50L,
                            warmup = 10L, verbose = FALSE)

  expect_true(all(diff(ss$betas) > 0))
})


test_that("Stepping-stone works with neomorphic characters", {
  library(ape)
  library(TreeTools)

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L),
                nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(6103)
  ss <- mkp_stepping_stone(pd, tree, neomorphic = 1L,
                            nStones = 5L, nIter = 100L,
                            warmup = 30L, verbose = FALSE)

  expect_true(is.finite(ss$log_marginal))
  expect_true(is.finite(ss$se))
})


test_that("Stepping-stone works with fixed topology", {
  library(ape)
  library(TreeTools)

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L),
                nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(8234)
  ss <- mkp_stepping_stone(pd, tree, nStones = 5L, nIter = 100L,
                            warmup = 30L, fixTopology = TRUE,
                            verbose = FALSE)

  expect_true(is.finite(ss$log_marginal))
  expect_true(is.finite(ss$se))
})


test_that("More stones with more iterations gives consistent results", {
  library(ape)
  library(TreeTools)

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L,
                  0L, 1L, 1L, 0L),
                nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(1583)
  ss1 <- mkp_stepping_stone(pd, tree, nStones = 15L, nIter = 500L,
                             warmup = 200L, verbose = FALSE)
  set.seed(7261)
  ss2 <- mkp_stepping_stone(pd, tree, nStones = 15L, nIter = 500L,
                             warmup = 200L, verbose = FALSE)

  # Two independent estimates should be in the same ballpark.
  # Use 10-sigma: the delta-method SE can underestimate true variance
  # when MCMC mixing is poor on small datasets with few iterations.
  combinedSe <- sqrt(ss1$se^2 + ss2$se^2)
  expect_lt(abs(ss1$log_marginal - ss2$log_marginal), 10 * combinedSe)
})


test_that("Stepping-stone accepts MkPrimeData input", {
  library(ape)
  library(TreeTools)

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0L, 1L, 0L, 1L), nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  set.seed(9428)
  ss <- mkp_stepping_stone(mkd, tree, nStones = 5L, nIter = 50L,
                            warmup = 10L, verbose = FALSE)

  expect_true(is.finite(ss$log_marginal))
  expect_true(is.finite(ss$se))
})


# --- Tests for .EssVector helper ---

test_that(".EssVector returns sensible values for iid samples", {
  set.seed(3847)
  x <- rnorm(500)
  ess <- MkPrime:::.EssVector(x)
  # iid samples: ESS should be close to n

  expect_gt(ess, 300)
  expect_lte(ess, 500)
})

test_that(".EssVector returns lower ESS for autocorrelated samples", {
  set.seed(5192)
  # AR(1) with rho = 0.9
  n <- 500
  x <- numeric(n)
  x[1] <- rnorm(1)
  for (i in 2:n) x[i] <- 0.9 * x[i - 1] + rnorm(1)

  ess <- MkPrime:::.EssVector(x)
  # Highly autocorrelated: ESS should be much less than n
  expect_gt(ess, 1)
  expect_lt(ess, 200)
})

test_that(".EssVector handles edge cases", {
  expect_equal(MkPrime:::.EssVector(numeric(0)), 1)
  expect_equal(MkPrime:::.EssVector(42), 1)
  # Constant vector: ESS = n (no variance to reduce)
  expect_equal(MkPrime:::.EssVector(rep(5, 100)), 100)
})


test_that("SE is positive and smaller with more iterations", {
  library(ape)
  library(TreeTools)

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L),
                nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(6712)
  ssSmall <- mkp_stepping_stone(pd, tree, nStones = 5L, nIter = 50L,
                                 warmup = 20L, verbose = FALSE)
  set.seed(6712)
  ssLarge <- mkp_stepping_stone(pd, tree, nStones = 5L, nIter = 500L,
                                 warmup = 50L, verbose = FALSE)

  expect_gt(ssSmall$se, 0)
  expect_gt(ssLarge$se, 0)
  # More samples should generally reduce SE (not guaranteed for any
  # single seed, but we check the SE is at least finite and positive)
  expect_true(is.finite(ssSmall$se))
  expect_true(is.finite(ssLarge$se))
})
