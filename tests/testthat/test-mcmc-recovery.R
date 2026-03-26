# Integration validation: posterior recovery of known parameters (M-021)
#
# These tests simulate data under known parameters and check that the
# MCMC posterior covers the true values. They are necessarily slow.

test_that("MCMC recovers tree_length from binary data (fixed topology)", {
  skip_on_cran()
  library(ape)

  # True parameters — 8 tips for reasonable signal
  set.seed(3847)
  true_tree <- rtree(8, tip.label = paste0("t", 1:8))
  true_tree$edge.length <- runif(nrow(true_tree$edge), 0.05, 0.25)
  true_tl <- sum(true_tree$edge.length)

  # Simulate binary characters under JC(2)
  nChar <- 60
  nNode <- 8 + 7  # tips + internals
  sim_mat <- matrix(NA_integer_, 8, nChar,
                    dimnames = list(true_tree$tip.label, NULL))
  true_tree_po <- reorder(true_tree, "postorder")

  for (ch in seq_len(nChar)) {
    root_state <- sample(0:1, 1)
    node_states <- integer(nNode)
    node_states[9] <- root_state

    edges <- true_tree_po$edge
    for (e in rev(seq_len(nrow(edges)))) {
      parent <- edges[e, 1]
      child <- edges[e, 2]
      t <- true_tree_po$edge.length[e]
      p_same <- 0.5 + 0.5 * exp(-2 * t)
      if (runif(1) < p_same) {
        node_states[child] <- node_states[parent]
      } else {
        node_states[child] <- 1L - node_states[parent]
      }
    }
    sim_mat[, ch] <- node_states[1:8]
  }

  variable <- apply(sim_mat, 2, function(x) length(unique(x)) > 1)
  sim_mat <- sim_mat[, variable, drop = FALSE]

  pd <- TreeTools::MatrixToPhyDat(sim_mat)
  mkd <- MkPrimeData(pd)

  # Set exp_steps near the true value for a reasonable prior
  set.seed(4592)
  result <- RunMkPrime(
    mkd, true_tree,
    model = MkPrimeModel(coding = "variable", relabel = TRUE,
                          exp_steps = true_tl),
    mcmc = MkPrimeMCMC(nIter = 15000L, thin = 10L, warmup = 7500L)
  )

  tl_samples <- result$samples[, "tree_length"]
  q01 <- quantile(tl_samples, 0.005)
  q99 <- quantile(tl_samples, 0.995)

  expect_true(
    true_tl > q01 && true_tl < q99,
    label = sprintf("tree_length %.2f in 99%% CI [%.2f, %.2f]",
                    true_tl, q01, q99)
  )
})


test_that("MCMC recovers rate_log_sd = 0 (no ACRV)", {
  skip_on_cran()
  library(ape)

  # Simulate data with no rate variation
  true_tree <- read.tree(
    text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);"
  )

  set.seed(1938)
  nChar <- 40
  sim_mat <- matrix(NA_integer_, 4, nChar,
                    dimnames = list(true_tree$tip.label, NULL))
  true_tree_po <- reorder(true_tree, "postorder")

  for (ch in seq_len(nChar)) {
    root_state <- sample(0:1, 1)
    node_states <- integer(7)
    node_states[5] <- root_state
    edges <- true_tree_po$edge
    for (e in rev(seq_len(nrow(edges)))) {
      parent <- edges[e, 1]
      child <- edges[e, 2]
      t <- true_tree_po$edge.length[e]
      p_same <- 0.5 + 0.5 * exp(-2 * t)
      if (runif(1) < p_same) {
        node_states[child] <- node_states[parent]
      } else {
        node_states[child] <- 1L - node_states[parent]
      }
    }
    sim_mat[, ch] <- node_states[1:4]
  }

  variable <- apply(sim_mat, 2, function(x) length(unique(x)) > 1)
  sim_mat <- sim_mat[, variable, drop = FALSE]

  pd <- TreeTools::MatrixToPhyDat(sim_mat)

  set.seed(6614)
  result <- RunMkPrime(
    pd, true_tree,
    model = MkPrimeModel(coding = "variable", exp_steps = 10),
    mcmc = MkPrimeMCMC(nIter = 10000L, thin = 10L, warmup = 5000L)
  )

  # rate_log_sd posterior should be concentrated near 0 (small values)
  rls_samples <- result$samples[, "rate_log_sd"]
  expect_lt(median(rls_samples), 2.0,
            label = "rate_log_sd median should be moderate when true = 0")
})


test_that("MCMC: k' stays near kObs for simple binary data", {
  skip_on_cran()
  library(ape)

  # Simulate binary (k=2) characters — k' should stay at 2
  true_tree <- read.tree(
    text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);"
  )

  set.seed(2795)
  nChar <- 20
  sim_mat <- matrix(NA_integer_, 4, nChar,
                    dimnames = list(true_tree$tip.label, NULL))
  true_tree_po <- reorder(true_tree, "postorder")

  for (ch in seq_len(nChar)) {
    root_state <- sample(0:1, 1)
    node_states <- integer(7)
    node_states[5] <- root_state
    edges <- true_tree_po$edge
    for (e in rev(seq_len(nrow(edges)))) {
      parent <- edges[e, 1]
      child <- edges[e, 2]
      t <- true_tree_po$edge.length[e]
      p_same <- 0.5 + 0.5 * exp(-2 * t)
      if (runif(1) < p_same) {
        node_states[child] <- node_states[parent]
      } else {
        node_states[child] <- 1L - node_states[parent]
      }
    }
    sim_mat[, ch] <- node_states[1:4]
  }

  variable <- apply(sim_mat, 2, function(x) length(unique(x)) > 1)
  sim_mat <- sim_mat[, variable, drop = FALSE]

  pd <- TreeTools::MatrixToPhyDat(sim_mat)

  set.seed(8103)
  result <- RunMkPrime(
    pd, true_tree,
    model = MkPrimeModel(coding = "variable", exp_steps = 10),
    mcmc = MkPrimeMCMC(nIter = 8000L, thin = 10L, warmup = 4000L)
  )

  # k' samples should be overwhelmingly 2 (true value)
  kp_cols <- grep("^kPrime_", colnames(result$samples))
  if (length(kp_cols)) {
    kp_mean <- mean(result$samples[, kp_cols])
    expect_lt(kp_mean, 2.5,
              label = "Mean k' should be close to 2 for true binary data")
  }
})
