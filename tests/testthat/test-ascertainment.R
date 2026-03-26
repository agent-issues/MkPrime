test_that("Constant site prob is in [0, 1]", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.1):0.1,t3:0.2);")
  tree <- reorder(tree, "postorder")

  for (k in c(2, 3, 5)) {
    pconst <- MkPrime:::constant_site_prob_jc(
      tree$edge[, 1], tree$edge[, 2], tree$edge.length,
      3L, k, rep(1 / k, k), 1.0
    )
    expect_true(pconst > 0 && pconst < 1,
                label = sprintf("k=%d, pconst=%g", k, pconst))
  }
})


test_that("Constant site prob → 1 with zero-length branches", {
  library(ape)
  # Near-zero branch lengths: tips are nearly identical to root
  tree <- read.tree(text = "((t1:1e-10,t2:1e-10):1e-10,t3:1e-10);")
  tree <- reorder(tree, "postorder")

  pconst <- MkPrime:::constant_site_prob_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    3L, 2L, c(0.5, 0.5), 1.0
  )

  # With near-zero branches, almost all sites are constant

  expect_gt(pconst, 0.99)
})


test_that("Constant site prob → 1/k with very long branches", {
  library(ape)
  tree <- read.tree(text = "((t1:100,t2:100):100,t3:100);")
  tree <- reorder(tree, "postorder")

  for (k in c(2, 3, 4)) {
    pconst <- MkPrime:::constant_site_prob_jc(
      tree$edge[, 1], tree$edge[, 2], tree$edge.length,
      3L, k, rep(1 / k, k), 1.0
    )
    # With very long branches, states are independent at each tip
    # P(all same) = sum_s (1/k)^nTip = k * (1/k)^3 = 1/k^2 for 3 tips
    expected <- k * (1 / k)^3
    expect_equal(pconst, expected, tolerance = 0.01,
                 label = sprintf("k=%d", k))
  }
})


test_that("Ascertainment correction makes likelihood more negative", {
  library(ape)
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  tree <- reorder(tree, "postorder")
  root_freqs <- rep(1 / 3, 3)
  tips <- matrix(c(0L, 1L, 2L, 0L), ncol = 1)

  ll_uncorrected <- MkPrime:::pruning_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    tips, 3L, root_freqs
  )

  pconst <- MkPrime:::constant_site_prob_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 3L, root_freqs, 1.0
  )

  nChar <- 1L
  ll_corrected <- ll_uncorrected - nChar * log(1 - pconst)

  # Ascertainment correction makes likelihood LESS negative:
  # L_corr = L / P(variable), and P(variable) < 1, so L_corr > L
  expect_gt(ll_corrected, ll_uncorrected)
})


test_that("MkN constant site prob with rate_loss=1 matches JC(2)", {
  library(ape)
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  tree <- reorder(tree, "postorder")
  root_freqs <- c(0.5, 0.5)
  rates <- MkPrime:::discrete_lognormal_rates(0.5, 6L)

  pconst_jc <- MkPrime:::constant_site_prob_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 2L, root_freqs, rates
  )

  pconst_mkn <- MkPrime:::constant_site_prob_mkn(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 1.0, root_freqs, rates
  )

  expect_equal(pconst_mkn, pconst_jc, tolerance = 1e-12)
})


test_that("Constant site prob works with ACRV", {
  library(ape)
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  tree <- reorder(tree, "postorder")
  root_freqs <- rep(1 / 3, 3)

  pconst_no_acrv <- MkPrime:::constant_site_prob_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 3L, root_freqs, 1.0
  )

  rates <- MkPrime:::discrete_lognormal_rates(1.0, 6L)
  pconst_acrv <- MkPrime:::constant_site_prob_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 3L, root_freqs, rates
  )

  # ACRV should change the constant site probability
  expect_false(isTRUE(all.equal(pconst_no_acrv, pconst_acrv)))
  # But it should still be in [0, 1]
  expect_true(pconst_acrv > 0 && pconst_acrv < 1)
})
