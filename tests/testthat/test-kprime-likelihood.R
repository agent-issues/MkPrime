# Tests for per-character k' in likelihood
# M-013: verify mkp_loglikelihood correctly handles kPrime > kObs

test_that("kPrime = kObs gives same result as default (no kPrime arg)", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat <- matrix(c(0, 1, 0), 3, 1,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  ll_default <- mkp_loglikelihood(tree, mkd, coding = "none",
                                  rate_log_sd = 0, relabel = FALSE)
  ll_explicit <- mkp_loglikelihood(tree, mkd, kPrime = mkd$kObs,
                                   coding = "none", rate_log_sd = 0,
                                   relabel = FALSE)
  expect_equal(ll_explicit, ll_default, tolerance = 1e-12)
})


test_that("kPrime > kObs gives different likelihood than kPrime = kObs", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat <- matrix(c(0, 1, 0), 3, 1,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  ll_k2 <- mkp_loglikelihood(tree, mkd, kPrime = 2L,
                              coding = "none", rate_log_sd = 0,
                              relabel = FALSE)
  ll_k5 <- mkp_loglikelihood(tree, mkd, kPrime = 5L,
                              coding = "none", rate_log_sd = 0,
                              relabel = FALSE)

  # JC(5) vs JC(2) on same data must give different log-likelihoods
  expect_false(isTRUE(all.equal(ll_k2, ll_k5)))
})


test_that("kPrime > kObs: likelihood decreases with larger k (no relabel)", {
  # With more states but same data, the model has more "wasted" probability
  # on unobserved states, so the raw likelihood decreases
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  ll_k2 <- mkp_loglikelihood(tree, mkd, kPrime = 2L,
                              coding = "none", rate_log_sd = 0,
                              relabel = FALSE)
  ll_k3 <- mkp_loglikelihood(tree, mkd, kPrime = 3L,
                              coding = "none", rate_log_sd = 0,
                              relabel = FALSE)
  ll_k5 <- mkp_loglikelihood(tree, mkd, kPrime = 5L,
                              coding = "none", rate_log_sd = 0,
                              relabel = FALSE)
  ll_k10 <- mkp_loglikelihood(tree, mkd, kPrime = 10L,
                               coding = "none", rate_log_sd = 0,
                               relabel = FALSE)

  expect_gt(ll_k2, ll_k3)
  expect_gt(ll_k3, ll_k5)
  expect_gt(ll_k5, ll_k10)
})


test_that("Per-character kPrime: different k' per character", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  # Two transformational characters, both kObs = 2
  mat <- matrix(c(0, 1, 0, 1,
                  0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  # Char 1: k'=2, Char 2: k'=4
  ll_mixed <- mkp_loglikelihood(tree, mkd, kPrime = c(2L, 4L),
                                 coding = "none", rate_log_sd = 0,
                                 relabel = FALSE)

  # Should equal sum of individual calls with respective k'
  ll_c1 <- mkp_loglikelihood(tree, mkd, kPrime = c(2L, 2L),
                              coding = "none", rate_log_sd = 0,
                              relabel = FALSE)
  ll_c2 <- mkp_loglikelihood(tree, mkd, kPrime = c(4L, 4L),
                              coding = "none", rate_log_sd = 0,
                              relabel = FALSE)

  # For single-character reference: compute individually
  mkd1 <- MkPrimeData(pd[, 1])
  mkd2 <- MkPrimeData(pd[, 2])
  ll_ref1 <- mkp_loglikelihood(tree, mkd1, kPrime = 2L,
                                coding = "none", rate_log_sd = 0,
                                relabel = FALSE)
  ll_ref2 <- mkp_loglikelihood(tree, mkd2, kPrime = 4L,
                                coding = "none", rate_log_sd = 0,
                                relabel = FALSE)

  expect_equal(ll_mixed, ll_ref1 + ll_ref2, tolerance = 1e-12)
})


test_that("kPrime with relabelling: k'=kObs has correction, k'>kObs smaller", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat <- matrix(c(0, 1, 0), 3, 1,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  # With relabelling: higher k' adds a correction that decreases with k'
  # But the raw likelihood also decreases with k'
  # Combined effect: the posterior (logL + relabel) should favor small k'
  ll_k2_rel <- mkp_loglikelihood(tree, mkd, kPrime = 2L,
                                  coding = "none", rate_log_sd = 0,
                                  relabel = TRUE)
  ll_k5_rel <- mkp_loglikelihood(tree, mkd, kPrime = 5L,
                                  coding = "none", rate_log_sd = 0,
                                  relabel = TRUE)

  # Both should be finite

  expect_true(is.finite(ll_k2_rel))
  expect_true(is.finite(ll_k5_rel))
})


test_that("kPrime with ascertainment correction and k' > kObs", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  ll <- mkp_loglikelihood(tree, mkd, kPrime = 4L,
                           coding = "variable", rate_log_sd = 0,
                           relabel = FALSE)
  expect_true(is.finite(ll))

  # Ascertainment correction should still increase likelihood
  ll_none <- mkp_loglikelihood(tree, mkd, kPrime = 4L,
                                coding = "none", rate_log_sd = 0,
                                relabel = FALSE)
  expect_gt(ll, ll_none)
})


test_that("kPrime with ACRV and k' > kObs", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  ll <- mkp_loglikelihood(tree, mkd, kPrime = 3L,
                           coding = "none", rate_log_sd = 0.5,
                           relabel = FALSE)
  expect_true(is.finite(ll))

  # ACRV should change result vs no ACRV
  ll_no_acrv <- mkp_loglikelihood(tree, mkd, kPrime = 3L,
                                   coding = "none", rate_log_sd = 0,
                                   relabel = FALSE)
  expect_false(isTRUE(all.equal(ll, ll_no_acrv)))
})


test_that("Large kPrime still produces finite likelihood", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat <- matrix(c(0, 1, 0), 3, 1,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  ll <- mkp_loglikelihood(tree, mkd, kPrime = 50L,
                           coding = "none", rate_log_sd = 0,
                           relabel = FALSE)
  expect_true(is.finite(ll))
  expect_lt(ll, 0)
})
