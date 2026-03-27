test_that("DiscreteLognormalRates returns nCat rates with mean 1", {
  for (sd in c(0.1, 0.5, 1, 2)) {
    rates <- MkPrime:::DiscreteLognormalRates(sd, nCat = 6L)
    expect_length(rates, 6)
    expect_equal(mean(rates), 1.0, tolerance = 1e-10)
    expect_true(all(rates > 0))
  }
})


test_that("DiscreteLognormalRates with sd=0 returns all 1s", {
  rates <- MkPrime:::DiscreteLognormalRates(0, nCat = 6L)
  expect_equal(rates, rep(1.0, 6))
})


test_that("DiscreteLognormalRates are sorted ascending", {
  rates <- MkPrime:::DiscreteLognormalRates(1.0, nCat = 6L)
  expect_equal(rates, sort(rates))
})


test_that("ACRV with uniform rates (all 1) matches non-ACRV", {
  library(ape)
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  tree <- TreeTools::Preorder(tree)
  root_freqs <- rep(1 / 3, 3)
  tips <- matrix(c(0L, 1L, 2L, 0L), ncol = 1)

  ll_no_acrv <- MkPrime:::pruning_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    tips, 3L, root_freqs
  )

  ll_acrv <- MkPrime:::pruning_jc_acrv(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    tips, 3L, root_freqs, rep(1.0, 6)
  )

  expect_equal(ll_acrv, ll_no_acrv, tolerance = 1e-12)
})


test_that("ACRV with single rate category matches non-ACRV", {
  library(ape)
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  tree <- TreeTools::Preorder(tree)
  root_freqs <- c(0.5, 0.5)
  tips <- matrix(c(0L, 1L, 0L, 1L), ncol = 1)

  ll_no_acrv <- MkPrime:::pruning_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    tips, 2L, root_freqs
  )

  ll_acrv <- MkPrime:::pruning_jc_acrv(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    tips, 2L, root_freqs, 1.0
  )

  expect_equal(ll_acrv, ll_no_acrv, tolerance = 1e-12)
})


test_that("ACRV produces different likelihood than non-ACRV for variable rates", {
  library(ape)
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  tree <- TreeTools::Preorder(tree)
  root_freqs <- rep(1 / 3, 3)
  tips <- matrix(c(0L, 1L, 2L, 0L, 0L, 0L, 1L, 1L), ncol = 2)

  rates <- MkPrime:::DiscreteLognormalRates(1.0, nCat = 6L)

  ll_no_acrv <- MkPrime:::pruning_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    tips, 3L, root_freqs
  )

  ll_acrv <- MkPrime:::pruning_jc_acrv(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    tips, 3L, root_freqs, rates
  )

  # ACRV likelihood should be >= non-ACRV (Jensen's inequality on log)
  # Actually, the mixture is more flexible, so the marginal likelihood
  # can be higher. The key point is they should differ.
  expect_false(isTRUE(all.equal(ll_acrv, ll_no_acrv)))
})


test_that("MkN ACRV matches non-ACRV with uniform rates", {
  library(ape)
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  tree <- TreeTools::Preorder(tree)
  rate_loss <- 2.0
  root_freqs <- as.numeric(MkPrime:::mkn_stationary_freqs(rate_loss))
  tips <- matrix(c(0L, 1L, 0L, 1L), ncol = 1)

  ll_no_acrv <- MkPrime:::pruning_mkn(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    tips, rate_loss, root_freqs
  )

  ll_acrv <- MkPrime:::pruning_mkn_acrv(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    tips, rate_loss, root_freqs, rep(1.0, 6)
  )

  expect_equal(ll_acrv, ll_no_acrv, tolerance = 1e-12)
})
