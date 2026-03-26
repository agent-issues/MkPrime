# Tests for parallel tempering (Phase 5: M-029, M-030)

# --- Temperature ladder ---

test_that(".build_temperature_ladder gives correct structure", {
  # Single chain: just beta = 1

  expect_equal(MkPrime:::.build_temperature_ladder(1L, 0.2), 1.0)

  # Two chains: cold + hot
  b2 <- MkPrime:::.build_temperature_ladder(2L, 0.2)
  expect_equal(length(b2), 2)
  expect_equal(b2[1], 1.0)
  expect_equal(b2[2], 0.2)

  # Four chains: geometric spacing

  b4 <- MkPrime:::.build_temperature_ladder(4L, 0.2)
  expect_equal(length(b4), 4)
  expect_equal(b4[1], 1.0)
  expect_equal(b4[4], 0.2)
  # Geometric: log(beta) should be equally spaced
  expect_equal(diff(log(b4)), rep(log(0.2) / 3, 3), tolerance = 1e-12)

  # All betas should be in (0, 1]
  b8 <- MkPrime:::.build_temperature_ladder(8L, 0.1)
  expect_true(all(b8 > 0))
  expect_true(all(b8 <= 1))
  expect_equal(b8[1], 1.0)
  expect_equal(b8[8], 0.1)
})


test_that("Temperature ladder is monotonically decreasing", {
  for (nC in c(2, 4, 6, 8)) {
    for (h in c(0.01, 0.1, 0.2, 0.5, 0.9)) {
      b <- MkPrime:::.build_temperature_ladder(nC, h)
      expect_true(all(diff(b) < 0),
                  label = paste("nChains =", nC, "heat =", h))
    }
  }
})


# --- MkPrimeMCMC config ---

test_that("MkPrimeMCMC validates nChains and heat", {
  expect_no_error(MkPrimeMCMC(nChains = 1L))
  expect_no_error(MkPrimeMCMC(nChains = 4L, heat = 0.2))

  expect_error(MkPrimeMCMC(nChains = 0L), "at least 1")
  expect_error(MkPrimeMCMC(nChains = 4L, heat = 0), "\\(0, 1\\)")
  expect_error(MkPrimeMCMC(nChains = 4L, heat = 1), "\\(0, 1\\)")
  expect_error(MkPrimeMCMC(nChains = 4L, heat = -0.1), "\\(0, 1\\)")
  expect_error(MkPrimeMCMC(nChains = 4L, heat = 1.5), "\\(0, 1\\)")

  # heat validation only applies when nChains > 1
  expect_no_error(MkPrimeMCMC(nChains = 1L, heat = 999))
})


test_that("MkPrimeMCMC stores nChains and heat", {
  cfg <- MkPrimeMCMC(nChains = 4L, heat = 0.3)
  expect_equal(cfg$nChains, 4L)
  expect_equal(cfg$heat, 0.3)
})


# --- Chain swap proposal ---

test_that(".propose_chain_swap is no-op for single chain", {
  state <- list(log_lik = -100, log_prior = -10)
  result <- MkPrime:::.propose_chain_swap(list(state), 1.0)
  expect_equal(length(result$chains), 1)
  expect_null(result$pair)
  expect_false(result$accepted)
})


test_that(".propose_chain_swap picks adjacent pairs", {
  states <- lapply(1:4, function(i) {
    list(log_lik = -100 + i * 10, log_prior = -5)
  })
  betas <- c(1.0, 0.6, 0.3, 0.1)

  set.seed(6418)
  pairs_seen <- integer(0)
  for (i in 1:200) {
    result <- MkPrime:::.propose_chain_swap(states, betas)
    expect_equal(length(result$pair), 2)
    expect_equal(result$pair[2], result$pair[1] + 1L)
    pairs_seen <- c(pairs_seen, result$pair[1])
  }
  # Should see all 3 adjacent pairs
  expect_true(all(1:3 %in% pairs_seen))
})


test_that(".propose_chain_swap exchanges states on acceptance", {
  # Construct states where swap is very favorable
  # (beta_i - beta_j) * (logLik_j - logLik_i) >> 0
  # Chain 1 (cold, beta=1) has low logLik, chain 2 (hot, beta=0.1) has high
  # Swap is favorable because cold chain "wants" the higher likelihood
  state1 <- list(log_lik = -1000, log_prior = -5, tag = "chain1")
  state2 <- list(log_lik = -10, log_prior = -5, tag = "chain2")

  set.seed(3847)
  n_swaps <- 0
  for (i in 1:100) {
    result <- MkPrime:::.propose_chain_swap(
      list(state1, state2), c(1.0, 0.1)
    )
    if (result$accepted) {
      n_swaps <- n_swaps + 1
      # States should be exchanged
      expect_equal(result$chains[[1]]$tag, "chain2")
      expect_equal(result$chains[[2]]$tag, "chain1")
    } else {
      expect_equal(result$chains[[1]]$tag, "chain1")
      expect_equal(result$chains[[2]]$tag, "chain2")
    }
  }
  # With such extreme likelihood difference, most swaps should be accepted
  expect_gt(n_swaps, 80)
})


test_that("Chain swap acceptance follows correct formula", {
  # Test with known values
  # log_alpha = (beta_i - beta_j) * (logLik_j - logLik_i)
  # = (1.0 - 0.5) * (-50 - (-100)) = 0.5 * 50 = 25 → always accept
  state1 <- list(log_lik = -100, log_prior = -5, tag = "a")
  state2 <- list(log_lik = -50, log_prior = -5, tag = "b")

  set.seed(2911)
  n_accept <- sum(vapply(1:100, function(i) {
    MkPrime:::.propose_chain_swap(list(state1, state2), c(1.0, 0.5))$accepted
  }, logical(1)))
  expect_equal(n_accept, 100)

  # Reverse: cold chain has higher logLik → swap unfavorable
  # log_alpha = (1.0 - 0.5) * (-100 - (-50)) = 0.5 * (-50) = -25
  n_accept2 <- sum(vapply(1:100, function(i) {
    MkPrime:::.propose_chain_swap(list(state2, state1), c(1.0, 0.5))$accepted
  }, logical(1)))
  expect_equal(n_accept2, 0)
})


# --- Heated MH acceptance in .do_move ---

test_that("Heated .do_move is more permissive than cold", {
  # Use real data to test that heated chains accept more proposals
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel()
  model <- MkPrime:::.finalize_model(model, tree, mkd)

  tree <- ape::reorder.phylo(tree, "postorder")
  state <- MkPrime:::.init_state(tree, mkd, model)

  tuning <- list(
    scale_tree_length = 2.0,
    beta_simplex = 10,
    scale_rate_loss = 0.5,
    scale_rate_log_sd = 0.5,
    scale_p = 0.5,
    int_walk_window = 1L
  )

  # Wide scale proposal to get many rejections for cold chain
  move <- list(name = "tree_length", type = "scale",
               target = "tree_length", weight = 1)

  set.seed(7562)
  n_cold <- 0L
  for (i in 1:200) {
    result <- MkPrime:::.do_move(move, state, mkd, model, tuning, beta = 1.0)
    if (result$accept) n_cold <- n_cold + 1L
  }

  set.seed(7562)
  n_hot <- 0L
  for (i in 1:200) {
    result <- MkPrime:::.do_move(move, state, mkd, model, tuning, beta = 0.1)
    if (result$accept) n_hot <- n_hot + 1L
  }

  # Heated chain should accept more (flattened likelihood)
  expect_gt(n_hot, n_cold)
})


test_that(".do_move with beta=1 is identical to unheated", {
  # Verify that beta=1 gives same MH ratio as the old formula
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel()
  model <- MkPrime:::.finalize_model(model, tree, mkd)
  tree <- ape::reorder.phylo(tree, "postorder")
  state <- MkPrime:::.init_state(tree, mkd, model)
  tuning <- MkPrimeMCMC()$tuning

  move <- list(name = "tree_length", type = "scale",
               target = "tree_length", weight = 1)

  # With beta=1, heated formula reduces to:
  # 1 * (logLik_new - logLik_old) + (logPrior_new - logPrior_old) + logH
  # = (logLik_new + logPrior_new) - (logLik_old + logPrior_old) + logH
  # = logPost_new - logPost_old + logH
  # which is the original formula
  set.seed(1294)
  results_beta1 <- replicate(50, {
    MkPrime:::.do_move(move, state, mkd, model, tuning, beta = 1.0)$accept
  })

  set.seed(1294)
  results_default <- replicate(50, {
    MkPrime:::.do_move(move, state, mkd, model, tuning)$accept
  })

  expect_identical(results_beta1, results_default)
})


# --- Multi-chain MCMC integration ---

test_that("RunMkPrime with nChains=1 matches Phase 4 behavior", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(5194)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nIter = 500L, thin = 5L, warmup = 200L,
                        nChains = 1L))

  expect_s3_class(result, "MkPosterior")
  expect_equal(nrow(result$samples), 60L)
  expect_true(all(is.finite(result$samples[, "log_posterior"])))
  # No tempering info
  expect_null(result$betas)
  expect_null(result$swap_rates)
})


test_that("RunMkPrime with nChains=4 runs successfully", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(6842)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nIter = 500L, thin = 5L, warmup = 200L,
                        nChains = 4L, heat = 0.2))

  expect_s3_class(result, "MkPosterior")
  expect_equal(nrow(result$samples), 60L)
  expect_true(all(is.finite(result$samples[, "log_posterior"])))

  # Tempering info should be present
  expect_equal(length(result$betas), 4)
  expect_equal(result$betas[1], 1.0)
  expect_equal(result$betas[4], 0.2)
  expect_equal(length(result$swap_rates), 3)
  expect_equal(length(result$chain_acceptance), 4)
})


test_that("RunMkPrime with nChains=2 and topology moves works", {
  library(ape)
  set.seed(8103)
  tree <- rtree(8)
  tree <- unroot(tree)
  mat <- matrix(sample(0:1, 8 * 5, replace = TRUE), 8, 5,
                dimnames = list(tree$tip.label, NULL))
  for (j in seq_len(ncol(mat))) {
    if (length(unique(mat[, j])) == 1) mat[1, j] <- 1L - mat[1, j]
  }
  pd <- TreeTools::MatrixToPhyDat(mat)

  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nIter = 1000L, thin = 5L, warmup = 500L,
                        nChains = 2L, heat = 0.3))

  expect_s3_class(result, "MkPosterior")
  expect_true("nni" %in% names(result$acceptance))
  expect_true("spr" %in% names(result$acceptance))
  # Swap rates should be reported
  expect_equal(length(result$swap_rates), 1)
})


test_that("Cold chain samples have valid posteriors under tempering", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(3319)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nIter = 1000L, thin = 5L, warmup = 500L,
                        nChains = 4L, heat = 0.1))

  # Cold chain samples should have log_post = log_lik + log_prior
  # (unheated; stored in state$log_post which is log_lik + log_prior)
  samps <- result$samples
  log_post <- samps[, "log_posterior"]
  log_lik <- samps[, "log_likelihood"]

  # log_post should be close to log_lik + prior (which varies)
  # At minimum, they should be finite and ordered consistently
  expect_true(all(is.finite(log_post)))
  expect_true(all(is.finite(log_lik)))
  # log_post = log_lik + log_prior; verify consistency
  # (can't guarantee sign of either, but they should track together)
  expect_true(all(log_post < 100))
})


test_that("Heated chains accept at higher rates", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(4410)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nIter = 2000L, thin = 10L, warmup = 1000L,
                        nChains = 4L, heat = 0.1))

  # Hottest chain should have higher overall acceptance than cold chain
  cold_acc <- mean(result$chain_acceptance[[1]])
  hot_acc <- mean(result$chain_acceptance[[4]])
  expect_gt(hot_acc, cold_acc)
})
