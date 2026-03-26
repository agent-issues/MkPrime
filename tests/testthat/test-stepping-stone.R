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
                            warmup = 30L, fix_topology = TRUE,
                            verbose = FALSE)

  expect_true(is.finite(ss$log_marginal))
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
  ss1 <- mkp_stepping_stone(pd, tree, nStones = 15L, nIter = 300L,
                             warmup = 100L, verbose = FALSE)
  set.seed(7261)
  ss2 <- mkp_stepping_stone(pd, tree, nStones = 15L, nIter = 300L,
                             warmup = 100L, verbose = FALSE)

  # Two independent estimates should be in the same ballpark
  combined_se <- sqrt(ss1$se^2 + ss2$se^2)
  expect_lt(abs(ss1$log_marginal - ss2$log_marginal), 5 * combined_se)
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
})
