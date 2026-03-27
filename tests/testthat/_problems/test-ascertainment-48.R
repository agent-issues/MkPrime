# Extracted from test-ascertainment.R:48

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "MkPrime", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
library(ape)
tree <- read.tree(text = "((t1:100,t2:100):100,t3:100);")
tree <- TreeTools::Preorder(tree)
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
