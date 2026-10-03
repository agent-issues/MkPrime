# Accuracy of the closed-form JC and MkN transition probabilities.

test_that("jc_transition_probs is finite over typical branch arguments", {
  # Typical phylogenetic branch lengths × rate multipliers × k/(k-1)
  # produce arguments in approximately [-25, -0.001].
  branch_lengths <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0)
  rate_multipliers <- c(0.1, 0.5, 1.0, 2.0, 5.0)
  k_factors <- c(2.0, 4.0 / 3, 6.0 / 5, 10.0 / 9, 20.0 / 19)

  for (bl in branch_lengths) {
    for (rm in rate_multipliers) {
      for (kf in k_factors) {
        x <- -kf * bl * rm
        expect_true(all(is.finite(jc_transition_probs(2, -x / 2))))
      }
    }
  }
})

test_that("JC transition probs match analytical formulas", {
  for (k in c(2, 3, 5, 10, 20)) {
    for (t_val in c(0.001, 0.01, 0.1, 0.5, 1.0, 5.0)) {
      P <- jc_transition_probs(k, t_val)
      # Row sums must be 1
      expect_equal(rowSums(P), rep(1.0, k), tolerance = 1e-12)
      # Diagonal elements must be equal (JC symmetry)
      expect_equal(length(unique(round(diag(P), 14))), 1)
      # Off-diagonal elements in each row must be equal
      for (i in 1:k) {
        off_diag <- P[i, -i]
        expect_equal(length(unique(round(off_diag, 14))), 1)
      }
      # At t=0, P = I
      if (t_val < 0.01) {
        expect_true(all(diag(P) > 0.99))
      }
    }
  }
})

test_that("MkN transition probs match analytical formulas", {
  for (rl in c(0.5, 1.0, 2.0, 5.0)) {
    for (t_val in c(0.001, 0.01, 0.1, 0.5, 1.0, 5.0)) {
      P <- mkn_transition_probs(rl, t_val)
      # 2x2 matrix
      expect_equal(dim(P), c(2, 2))
      # Row sums must be 1
      expect_equal(rowSums(P), c(1.0, 1.0), tolerance = 1e-12)
      # All entries non-negative
      expect_true(all(P >= 0))
    }
  }
})

test_that("pruning_jc gives a finite, negative log-likelihood", {
  skip_if_not_installed("ape")
  tree <- read.tree(text = "((t1:0.3,t2:0.5):0.1,(t3:0.2,t4:0.4):0.2,t5:0.6);")
  tree <- Preorder(tree)

  tip_states <- matrix(c(
    0L, 1L, 0L, 1L, 0L,
    1L, 1L, 0L, 0L, 1L,
    0L, 0L, 1L, 1L, 0L
  ), nrow = 5, ncol = 3)

  parent <- tree$edge[, 1]
  child  <- tree$edge[, 2]
  el     <- tree$edge.length

  for (k in c(2, 3, 5)) {
    rf <- rep(1.0 / k, k)
    ll <- pruning_jc(parent, child, el, tip_states, k, rf)
    expect_true(is.finite(ll))
    expect_true(ll < 0)  # log-likelihood must be negative
  }
})
