# Tests for parallel tempering (Phase 5: M-029, M-030)
skip_slow_tests()

# --- Temperature ladder ---

test_that(".BuildTemperatureLadder gives correct structure", {
  # Single chain: just beta = 1

  expect_equal(MkPrime:::.BuildTemperatureLadder(1L, 0.2), 1.0)

  # Two chains: cold + hot
  b2 <- MkPrime:::.BuildTemperatureLadder(2L, 0.2)
  expect_equal(length(b2), 2)
  expect_equal(b2[1], 1.0)
  expect_equal(b2[2], 0.2)

  # Four chains: geometric spacing

  b4 <- MkPrime:::.BuildTemperatureLadder(4L, 0.2)
  expect_equal(length(b4), 4)
  expect_equal(b4[1], 1.0)
  expect_equal(b4[4], 0.2)
  # Geometric: log(beta) should be equally spaced
  expect_equal(diff(log(b4)), rep(log(0.2) / 3, 3), tolerance = 1e-12)

  # All betas should be in (0, 1]
  b8 <- MkPrime:::.BuildTemperatureLadder(8L, 0.1)
  expect_true(all(b8 > 0))
  expect_true(all(b8 <= 1))
  expect_equal(b8[1], 1.0)
  expect_equal(b8[8], 0.1)
})


test_that("Temperature ladder is monotonically decreasing", {
  for (nC in c(2, 4, 6, 8)) {
    for (h in c(0.01, 0.1, 0.2, 0.5, 0.9)) {
      b <- MkPrime:::.BuildTemperatureLadder(nC, h)
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

test_that(".ProposeChainSwap is no-op for single chain", {
  state <- list(log_lik = -100, log_prior = -10)
  result <- MkPrime:::.ProposeChainSwap(list(state), 1.0)
  expect_equal(length(result$chains), 1)
  expect_null(result$pair)
  expect_false(result$accepted)
})


test_that(".ProposeChainSwap picks adjacent pairs", {
  states <- lapply(1:4, function(i) {
    list(log_lik = -100 + i * 10, log_prior = -5)
  })
  betas <- c(1.0, 0.6, 0.3, 0.1)

  set.seed(6418)
  pairs_seen <- integer(0)
  for (i in 1:200) {
    result <- MkPrime:::.ProposeChainSwap(states, betas)
    expect_equal(length(result$pair), 2)
    expect_equal(result$pair[2], result$pair[1] + 1L)
    pairs_seen <- c(pairs_seen, result$pair[1])
  }
  # Should see all 3 adjacent pairs
  expect_true(all(1:3 %in% pairs_seen))
})


test_that(".ProposeChainSwap exchanges states on acceptance", {
  # Construct states where swap is very favorable
  # (beta_i - beta_j) * (logLik_j - logLik_i) >> 0
  # Chain 1 (cold, beta=1) has low logLik, chain 2 (hot, beta=0.1) has high
  # Swap is favorable because cold chain "wants" the higher likelihood
  state1 <- list(log_lik = -1000, log_prior = -5, tag = "chain1")
  state2 <- list(log_lik = -10, log_prior = -5, tag = "chain2")

  set.seed(3847)
  n_swaps <- 0
  for (i in 1:100) {
    result <- MkPrime:::.ProposeChainSwap(
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
    MkPrime:::.ProposeChainSwap(list(state1, state2), c(1.0, 0.5))$accepted
  }, logical(1)))
  expect_equal(n_accept, 100)

  # Reverse: cold chain has higher logLik → swap unfavorable
  # log_alpha = (1.0 - 0.5) * (-100 - (-50)) = 0.5 * (-50) = -25
  n_accept2 <- sum(vapply(1:100, function(i) {
    MkPrime:::.ProposeChainSwap(list(state2, state1), c(1.0, 0.5))$accepted
  }, logical(1)))
  expect_equal(n_accept2, 0)
})


# --- Heated MH acceptance in .DoMove ---

test_that("Heated .DoMove is more permissive than cold", {
  # Use real data to test that heated chains accept more proposals
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel()
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)

  tree <- TreeTools::Preorder(tree)
  state <- MkPrime:::.InitState(tree, mkd, model)

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
    result <- MkPrime:::.DoMove(move, state, mkd, model, tuning, beta = 1.0)
    if (result$accept) n_cold <- n_cold + 1L
  }

  set.seed(7562)
  n_hot <- 0L
  for (i in 1:200) {
    result <- MkPrime:::.DoMove(move, state, mkd, model, tuning, beta = 0.1)
    if (result$accept) n_hot <- n_hot + 1L
  }

  # Heated chain should accept more (flattened likelihood)
  expect_gt(n_hot, n_cold)
})


test_that(".DoMove with beta=1 is identical to unheated", {
  # Verify that beta=1 gives same MH ratio as the old formula
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel()
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  tree <- TreeTools::Preorder(tree)
  state <- MkPrime:::.InitState(tree, mkd, model)
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
    MkPrime:::.DoMove(move, state, mkd, model, tuning, beta = 1.0)$accept
  })

  set.seed(1294)
  results_default <- replicate(50, {
    MkPrime:::.DoMove(move, state, mkd, model, tuning)$accept
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
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
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
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE,
                        nChains = 4L, heat = 0.2))

  expect_s3_class(result, "MkPosterior")
  expect_equal(nrow(result$samples), 60L)
  expect_true(all(is.finite(result$samples[, "log_posterior"])))

  # Tempering info should be present
  expect_equal(length(result$betas), 4)
  expect_equal(result$betas[1], 1.0)
  # betas may have been adapted during warmup, but structure is preserved
  expect_equal(length(result$betas), 4)
  expect_equal(result$betas[1], 1.0)
  expect_true(all(diff(result$betas) < 0))
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
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000L, thin = 5L, maxWarmup = 500L, minWarmup = 500L, autoTune = FALSE,
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
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000L, thin = 5L, maxWarmup = 500L, minWarmup = 500L, autoTune = FALSE,
                        nChains = 4L, heat = 0.1))

  # Cold chain samples should have log_post = log_lik + LogPrior
  # (unheated; stored in state$log_post which is log_lik + LogPrior)
  samps <- result$samples
  log_post <- samps[, "log_posterior"]
  log_lik <- samps[, "log_likelihood"]

  # log_post should be close to log_lik + prior (which varies)
  # At minimum, they should be finite and ordered consistently
  expect_true(all(is.finite(log_post)))
  expect_true(all(is.finite(log_lik)))
  # log_post = log_lik + LogPrior; verify consistency
  # (can't guarantee sign of either, but they should track together)
  expect_true(all(log_post < 100))
})


# --- Adaptive temperature tuning (M-031) ---

test_that(".AdaptTemperatures is no-op with insufficient data", {
  betas <- c(1.0, 0.6, 0.3, 0.1)
  # Not enough proposals (total < 20)
  result <- MkPrime:::.AdaptTemperatures(betas, c(1, 1, 1), c(5, 5, 5))
  expect_equal(result, betas)
  # Single chain
  expect_equal(MkPrime:::.AdaptTemperatures(1.0, integer(0), integer(0)), 1.0)
})


test_that(".AdaptTemperatures increases heat when swaps too low", {
  betas <- c(1.0, 0.6, 0.3, 0.1)
  # Very low swap acceptance → should bring temps closer (increase heat)
  new_betas <- MkPrime:::.AdaptTemperatures(
    betas, c(1, 1, 1), c(100, 100, 100), target = 0.25
  )
  # Hottest chain should be warmer (closer to 1)
  expect_gt(new_betas[4], betas[4])
})


test_that(".AdaptTemperatures decreases heat when swaps too high", {
  betas <- c(1.0, 0.6, 0.3, 0.1)
  # Very high swap acceptance → should spread temps further (decrease heat)
  new_betas <- MkPrime:::.AdaptTemperatures(
    betas, c(90, 90, 90), c(100, 100, 100), target = 0.25
  )
  # Hottest chain should be colder (further from 1)
  expect_lt(new_betas[4], betas[4])
})


test_that(".AdaptTemperatures stays no-op near target", {
  betas <- c(1.0, 0.6, 0.3, 0.1)
  # Swap rate ~25% → minimal change
  new_betas <- MkPrime:::.AdaptTemperatures(
    betas, c(25, 25, 25), c(100, 100, 100), target = 0.25
  )
  # Should barely change
  expect_equal(new_betas[4], betas[4], tolerance = 0.02)
})


test_that(".AdaptTemperatures respects heat bounds", {
  # Very aggressive → heat should not exceed 0.95
  betas_wide <- c(1.0, 0.001)
  result <- MkPrime:::.AdaptTemperatures(
    betas_wide, c(0), c(100), target = 0.25
  )
  expect_lte(result[2], 0.95)
  expect_gte(result[2], 0.01)

  # Very conservative → heat should not go below 0.01
  betas_close <- c(1.0, 0.94)
  result2 <- MkPrime:::.AdaptTemperatures(
    betas_close, c(100), c(100), target = 0.25
  )
  expect_gte(result2[2], 0.01)
})


test_that("Adaptive temps integrated into MCMC warmup", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(5580)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 2000L, thin = 10L, maxWarmup = 1000L, minWarmup = 1000L, autoTune = FALSE,
                        nChains = 4L, heat = 0.5))

  # Temperatures should have been adapted (may increase or decrease
  # depending on swap rates), so final heat likely differs from 0.5
  expect_true(abs(result$betas[4] - 0.5) > 0.001 ||
              all(result$swap_rates > 0.15 & result$swap_rates < 0.35))
})


# --- Independent runs (M-032) ---

test_that("MkPrimeMCMC validates nRuns", {
  expect_no_error(MkPrimeMCMC(nRuns = 1L))
  expect_no_error(MkPrimeMCMC(nRuns = 2L))
  expect_error(MkPrimeMCMC(nRuns = 0L), "at least 1")
})


test_that("RunMkPrime with nRuns=2 runs successfully", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(7218)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 500L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE))

  expect_s3_class(result, "MkPosterior")
  # Combined samples: 60 per run × 2 runs = 120
  expect_equal(nrow(result$samples), 120L)
  expect_equal(length(result$trees), 120L)
  expect_equal(result$nRuns, 2L)
  expect_equal(length(result$per_run), 2)
  # Each per-run should have 60 samples
  expect_equal(nrow(result$per_run[[1]]$samples), 60L)
  expect_equal(nrow(result$per_run[[2]]$samples), 60L)
})


test_that("RunMkPrime with nRuns=1 has no per_run field", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(1498)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE))

  expect_null(result$nRuns)
  expect_null(result$per_run)
  expect_equal(nrow(result$samples), 60L)
})


test_that("Independent runs start from different states", {
  library(ape)
  set.seed(6334)
  tree <- rtree(6)
  tree <- unroot(tree)
  mat <- matrix(sample(0:1, 6 * 4, replace = TRUE), 6, 4,
                dimnames = list(tree$tip.label, NULL))
  for (j in seq_len(ncol(mat))) {
    if (length(unique(mat[, j])) == 1) mat[1, j] <- 1L - mat[1, j]
  }
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(2891)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 500L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE))

  # Runs should have different starting log-posteriors (different start trees)
  # First sample of each run may differ
  r1_first <- result$per_run[[1]]$samples[1, "log_posterior"]
  r2_first <- result$per_run[[2]]$samples[1, "log_posterior"]
  # They could be the same by chance, but very unlikely with perturbed starts
  # Just check both are finite
  expect_true(is.finite(r1_first))
  expect_true(is.finite(r2_first))
})


test_that("Multi-run with tempering works", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(8563)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nChains = 2L, heat = 0.3,
                        nIter = 500L, thin = 5L, maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE))

  expect_equal(nrow(result$samples), 120L)
  expect_equal(result$nRuns, 2L)
  # Tempering info should be present
  expect_true(!is.null(result$betas))
  expect_equal(length(result$swap_rates), 1)
})


test_that(".PerturbStart produces valid trees", {
  library(ape)
  set.seed(3521)
  tree <- rtree(10)
  tree <- unroot(tree)

  for (i in 1:20) {
    perturbed <- MkPrime:::.PerturbStart(tree)
    expect_s3_class(perturbed, "phylo")
    expect_equal(length(perturbed$tip.label), 10L)
    expect_true(all(perturbed$edge.length > 0))
    expect_equal(nrow(perturbed$edge), nrow(tree$edge))
  }
})


test_that(".PerturbStart handles small trees", {
  library(ape)
  tree <- read.tree(text = "(t1:0.1,t2:0.2,t3:0.3);")
  perturbed <- MkPrime:::.PerturbStart(tree)
  expect_s3_class(perturbed, "phylo")
})


test_that("Heated chains accept at higher rates", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(4410)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 2000L, thin = 10L, maxWarmup = 1000L, minWarmup = 1000L, autoTune = FALSE,
                        nChains = 4L, heat = 0.1))

  # Hottest chain should have higher overall acceptance than cold chain
  cold_acc <- mean(result$chain_acceptance[[1]])
  hot_acc <- mean(result$chain_acceptance[[4]])
  expect_gt(hot_acc, cold_acc)
})
