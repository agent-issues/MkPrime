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

  # Set expSteps near the true value for a reasonable prior
  set.seed(4592)
  result <- RunMkPrime(
    mkd, true_tree,
    model = MkPrimeModel(coding = "variable", relabel = TRUE,
                          expSteps = true_tl),
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
    model = MkPrimeModel(coding = "variable", expSteps = 10),
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
    model = MkPrimeModel(coding = "variable", expSteps = 10),
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


# --- M-028: Tree topology recovery ---

test_that("MCMC with topology moves recovers true tree from simulated data", {
  skip_on_cran()
  library(ape)

  # True tree with clear structure: ((t1,t2),(t3,(t4,t5)))
  true_tree <- read.tree(
    text = "((t1:0.15,t2:0.15):0.2,(t3:0.1,(t4:0.1,t5:0.1):0.15):0.1);"
  )
  true_tree <- unroot(true_tree)
  true_tree <- reorder.phylo(true_tree, "postorder")
  nTip <- 5L
  nNode <- true_tree$Nnode
  nTotal <- nTip + nNode

  # Simulate 80 binary characters under JC(2)
  set.seed(2651)
  nChar <- 80
  true_tree_po <- true_tree
  sim_mat <- matrix(NA_integer_, nTip, nChar,
                    dimnames = list(true_tree_po$tip.label, NULL))

  for (ch in seq_len(nChar)) {
    root <- nTip + 1L
    node_states <- integer(nTotal)
    node_states[root] <- sample(0:1, 1)

    edges <- true_tree_po$edge
    for (e in rev(seq_len(nrow(edges)))) {
      parent <- edges[e, 1]
      child <- edges[e, 2]
      bl <- true_tree_po$edge.length[e]
      # JC(2) transition probability
      p_change <- 0.5 * (1 - exp(-2 * bl))
      if (runif(1) < p_change) {
        node_states[child] <- 1L - node_states[parent]
      } else {
        node_states[child] <- node_states[parent]
      }
    }
    sim_mat[, ch] <- node_states[seq_len(nTip)]
  }

  # Remove invariant characters
  variable <- apply(sim_mat, 2, function(x) length(unique(x)) > 1)
  sim_mat <- sim_mat[, variable, drop = FALSE]

  pd <- TreeTools::MatrixToPhyDat(sim_mat)

  # Start from a DIFFERENT tree to test topology search
  set.seed(7734)
  start_tree <- rtree(nTip, tip.label = true_tree$tip.label)
  start_tree <- unroot(start_tree)

  result <- RunMkPrime(
    pd, start_tree,
    model = MkPrimeModel(coding = "variable"),
    mcmc = MkPrimeMCMC(nIter = 20000L, thin = 20L, warmup = 10000L)
  )

  # Check: posterior trees should be close to the true tree.
  rf_dists <- vapply(result$trees, function(tr) {
    TreeDist::RobinsonFoulds(tr, true_tree)
  }, numeric(1))

  # At least some posterior trees should match the true topology (RF = 0)
  # or be within 2 RF moves
  expect_true(
    any(rf_dists <= 2),
    label = "At least one posterior tree within RF 2 of true tree"
  )
  # Median RF should be small (for a 5-tip tree, max RF = 4)
  expect_lte(median(rf_dists), 2)
})
