test_that("JC P(t) rows sum to 1", {
  for (k in c(2, 3, 5, 10)) {
    for (t in c(0, 0.01, 0.1, 1, 10)) {
      P <- MkPrime:::jc_transition_probs(k, t)
      expect_equal(rowSums(P), rep(1.0, k), tolerance = 1e-12,
                   label = sprintf("k=%d, t=%g", k, t))
    }
  }
})


test_that("JC P(0) is identity matrix", {
  for (k in c(2, 4, 8)) {
    P <- MkPrime:::jc_transition_probs(k, 0)
    expect_equal(P, diag(k), tolerance = 1e-15)
  }
})


test_that("JC P(Inf) converges to uniform (1/k)", {
  for (k in c(2, 3, 5)) {
    P <- MkPrime:::jc_transition_probs(k, 1000)
    expect_equal(as.vector(P), rep(1 / k, k * k), tolerance = 1e-10)
  }
})


test_that("JC P(t) is symmetric", {
  P <- MkPrime:::jc_transition_probs(4, 0.5)
  expect_equal(P, t(P), tolerance = 1e-15)
})


test_that("JC P(t) entries are non-negative", {
  P <- MkPrime:::jc_transition_probs(3, 0.001)
  expect_true(all(P >= 0))
})


test_that("JC validates inputs", {
  expect_error(MkPrime:::jc_transition_probs(1, 0.5), "k must be >= 2")
  expect_error(MkPrime:::jc_transition_probs(3, -1), "t must be >= 0")
})


test_that("JC k=2 matches Mk2 (symmetric binary)", {
  t <- 0.3
  P <- MkPrime:::jc_transition_probs(2, t)
  # For k=2: P_ii = 0.5 + 0.5*exp(-2t), P_ij = 0.5 - 0.5*exp(-2t)
  expected_diag <- 0.5 + 0.5 * exp(-2 * t)
  expected_off <- 0.5 - 0.5 * exp(-2 * t)
  expect_equal(P[1, 1], expected_diag, tolerance = 1e-15)
  expect_equal(P[1, 2], expected_off, tolerance = 1e-15)
})


# MkN tests

test_that("MkN P(t) rows sum to 1", {
  for (rl in c(0.5, 1, 2, 5)) {
    for (t in c(0, 0.01, 0.1, 1, 10)) {
      P <- MkPrime:::mkn_transition_probs(rl, t)
      expect_equal(rowSums(P), rep(1.0, 2), tolerance = 1e-12,
                   label = sprintf("rl=%g, t=%g", rl, t))
    }
  }
})


test_that("MkN P(0) is identity matrix", {
  P <- MkPrime:::mkn_transition_probs(2, 0)
  expect_equal(P, diag(2), tolerance = 1e-15)
})


test_that("MkN P(Inf) converges to stationary frequencies", {
  rate_loss <- 3
  P <- MkPrime:::mkn_transition_probs(rate_loss, 1000)
  pi_vec <- MkPrime:::mkn_stationary_freqs(rate_loss)
  # Both rows should equal stationary distribution
  expect_equal(P[1, ], as.numeric(pi_vec), tolerance = 1e-10)
  expect_equal(P[2, ], as.numeric(pi_vec), tolerance = 1e-10)
})


test_that("MkN with rate_loss=1 matches JC(2)", {
  t <- 0.5
  P_mkn <- MkPrime:::mkn_transition_probs(1, t)
  P_jc <- MkPrime:::jc_transition_probs(2, t)
  expect_equal(P_mkn, P_jc, tolerance = 1e-14)
})


test_that("MkN satisfies detailed balance", {
  rate_loss <- 2.5
  t <- 0.3
  P <- MkPrime:::mkn_transition_probs(rate_loss, t)
  pi_vec <- as.numeric(MkPrime:::mkn_stationary_freqs(rate_loss))
  # Detailed balance: π_i * P_ij = π_j * P_ji
  expect_equal(pi_vec[1] * P[1, 2], pi_vec[2] * P[2, 1], tolerance = 1e-14)
})


test_that("MkN Q has stationary-weighted mean rate 1", {
  # P00(t) = pi0 + pi1 exp(-lambda t), and q01 = pi1 lambda, q10 = pi0 lambda
  t <- 0.7
  for (rl in c(0.12, 0.5, 1, 2, 5)) {
    P <- MkPrime:::mkn_transition_probs(rl, t)
    pi_vec <- as.numeric(MkPrime:::mkn_stationary_freqs(rl))
    lambda <- -log((P[1, 1] - pi_vec[1]) / pi_vec[2]) / t
    expect_equal(2 * pi_vec[1] * pi_vec[2] * lambda, 1, tolerance = 1e-12,
                 label = sprintf("rl=%g", rl))
  }
})


test_that("MkN matches the normalised F81 of the Q-heterogeneity path", {
  # A single beta bin of 0.5 gives the F81 component the MkN frequencies;
  # F81 is normalised to mean rate 1, so the two models must coincide.
  tree <- Preorder(ape::read.tree(
    text = "((t1:0.3,t2:0.1):0.2,(t3:0.4,(t4:0.2,t5:0.5):0.1):0.3);"))
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  rates <- c(0.4, 1.6)
  for (rl in c(0.2, 1, 3)) {
    expect_equal(
      MkPrime:::het_const_site_prob(parent, child, tree$edge.length, 5L, 2L,
                                    rl, 0.5, rates),
      MkPrime:::constant_site_prob_mkn(parent, child, tree$edge.length, 5L, rl,
                                       MkPrime:::mkn_stationary_freqs(rl),
                                       rates),
      tolerance = 1e-12, label = sprintf("rl=%g", rl))
  }
})


test_that("MkN stationary frequencies sum to 1", {
  for (rl in c(0.1, 0.5, 1, 2, 10)) {
    pi_vec <- MkPrime:::mkn_stationary_freqs(rl)
    expect_equal(sum(pi_vec), 1.0, tolerance = 1e-15)
  }
})


test_that("MkN validates inputs", {
  expect_error(MkPrime:::mkn_transition_probs(0, 0.5), "rate_loss must be > 0")
  expect_error(MkPrime:::mkn_transition_probs(-1, 0.5), "rate_loss must be > 0")
  expect_error(MkPrime:::mkn_transition_probs(1, -1), "t must be >= 0")
  expect_error(MkPrime:::mkn_stationary_freqs(0), "rate_loss must be > 0")
})
