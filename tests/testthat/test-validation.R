test_that("Standard Mk likelihood matches phangorn (2-state)", {
  skip_if_not_installed("phangorn")
  library(ape)

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat <- matrix(c(0, 1, 0), 3, 1,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- phangorn::phyDat(mat, type = "USER", levels = c("0", "1"))

  fit <- phangorn::pml(tree, pd, model = "ER")
  mkd <- MkPrimeData(pd)
  ll <- mkp_loglikelihood(tree, mkd, coding = "none",
                          rate_log_sd = 0, relabel = FALSE)

  expect_equal(ll, fit$logLik, tolerance = 1e-10)
})


test_that("Standard Mk likelihood matches phangorn (3-state, multi-char)", {
  skip_if_not_installed("phangorn")
  library(ape)

  set.seed(8832)
  tree <- rtree(6, tip.label = paste0("t", 1:6))
  tree$edge.length <- runif(length(tree$edge.length), 0.05, 0.5)
  mat <- matrix(sample(0:2, 6 * 4, replace = TRUE), 6, 4,
                dimnames = list(tree$tip.label, NULL))
  pd <- phangorn::phyDat(mat, type = "USER", levels = c("0", "1", "2"))

  fit <- phangorn::pml(tree, pd, model = "ER")
  mkd <- MkPrimeData(pd)
  ll <- mkp_loglikelihood(tree, mkd, coding = "none",
                          rate_log_sd = 0, relabel = FALSE)

  expect_equal(ll, fit$logLik, tolerance = 1e-10)
})


test_that("Standard Mk likelihood matches phangorn (10 tips, 8 chars)", {
  skip_if_not_installed("phangorn")
  library(ape)

  set.seed(5501)
  tree <- rtree(10, tip.label = paste0("t", 1:10))
  tree$edge.length <- runif(length(tree$edge.length), 0.01, 0.3)
  mat <- matrix(sample(0:1, 10 * 8, replace = TRUE), 10, 8,
                dimnames = list(tree$tip.label, NULL))
  pd <- phangorn::phyDat(mat, type = "USER", levels = c("0", "1"))

  fit <- phangorn::pml(tree, pd, model = "ER")
  mkd <- MkPrimeData(pd)
  ll <- mkp_loglikelihood(tree, mkd, coding = "none",
                          rate_log_sd = 0, relabel = FALSE)

  expect_equal(ll, fit$logLik, tolerance = 1e-10)
})


test_that("Mk' relabelling increases log-likelihood when applied", {
  library(ape)

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat <- matrix(c(0, 1, 0), 3, 1,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  ll_no_relabel <- mkp_loglikelihood(tree, mkd, coding = "none",
                                     rate_log_sd = 0, relabel = FALSE)
  ll_relabel <- mkp_loglikelihood(tree, mkd, coding = "none",
                                  rate_log_sd = 0, relabel = TRUE)

  # Relabelling correction with kPrime = kObs adds k*log(k) > 0
  expect_gt(ll_relabel, ll_no_relabel)
})


test_that("Ascertainment correction increases likelihood", {
  library(ape)

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  ll_none <- mkp_loglikelihood(tree, mkd, coding = "none",
                               rate_log_sd = 0, relabel = FALSE)
  ll_var <- mkp_loglikelihood(tree, mkd, coding = "variable",
                              rate_log_sd = 0, relabel = FALSE)

  # Variable coding divides by P(variable) < 1, so logL increases
  expect_gt(ll_var, ll_none)
})


test_that("ACRV changes likelihood vs no rate variation", {
  library(ape)

  set.seed(3742)
  tree <- rtree(6, tip.label = paste0("t", 1:6))
  tree$edge.length <- runif(length(tree$edge.length), 0.05, 0.5)
  mat <- matrix(sample(0:2, 6 * 6, replace = TRUE), 6, 6,
                dimnames = list(tree$tip.label, NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  ll_no_acrv <- mkp_loglikelihood(tree, mkd, coding = "none",
                                  rate_log_sd = 0, relabel = FALSE)
  ll_acrv <- mkp_loglikelihood(tree, mkd, coding = "none",
                               rate_log_sd = 1.0, relabel = FALSE)

  expect_false(isTRUE(all.equal(ll_no_acrv, ll_acrv)))
})


test_that("mkp_loglikelihood handles mixed partition types", {
  library(ape)

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1,   # binary, neomorphic
                  0, 1, 2, 0,    # 3-state, transformational
                  0, 1, 0, 1),   # binary, known k=3
                4, 3,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = 1L, known_states = c("3" = 3L))

  ll <- mkp_loglikelihood(tree, mkd, coding = "none",
                          rate_log_sd = 0, relabel = FALSE)
  expect_true(is.finite(ll))
  expect_lt(ll, 0)
})
