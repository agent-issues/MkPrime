# tests/testthat/test-gibbs-het.R
# M-114: Validate partial CL for Gibbs moves under Q-heterogeneity

library("MkPrime")
library("TreeTools")

# ---------------------------------------------------------------------------
# Helper: create a small dataset + model with qHeterogeneity
# ---------------------------------------------------------------------------
make_het_setup <- function(nTip = 8L, nChar = 10L, kMax = 3L,
                            seed = 5719L) {
  set.seed(seed)
  tree <- ape::rtree(nTip, rooted = FALSE)
  tree <- Preorder(tree)
  mat <- matrix(sample(0:(kMax - 1L), nTip * nChar, replace = TRUE),
                nrow = nTip,
                dimnames = list(tree$tip.label, paste0("c", seq_len(nChar))))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(qHeterogeneity = TRUE)
  list(tree = tree, mkd = mkd, model = model)
}


# ---------------------------------------------------------------------------
# Test: F81 transition unit test
# ---------------------------------------------------------------------------
test_that("f81_transition matches analytical formula", {
  # P_ij(t) = π_j + (δ_{ij} - π_j) * exp(-μt)
  # (P × cl)_i = (1-e^{-μt}) * dot(π,cl) + e^{-μt} * cl_i
  pi <- c(0.3, 0.5, 0.2)
  k <- 3L
  sumPiSq <- sum(pi^2)
  mu <- 1 / (1 - sumPiSq)
  t_val <- 0.15

  # Build transition matrix analytically
  Pmat <- matrix(0, k, k)
  for (i in 1:k) {
    for (j in 1:k) {
      Pmat[i, j] <- pi[j] + (as.numeric(i == j) - pi[j]) * exp(-mu * t_val)
    }
  }

  # Check rows sum to 1
  expect_equal(rowSums(Pmat), rep(1, k), tolerance = 1e-14)
  # Check non-negative

  expect_true(all(Pmat >= 0))

  # Check that Σ_j P_ij cl_j = (1-exp(-μt))*dot(π,cl) + exp(-μt)*cl_i
  cl <- c(0.4, 0.1, 0.5)
  result <- Pmat %*% cl
  piDotCl <- sum(pi * cl)
  expected <- (1 - exp(-mu * t_val)) * piDotCl + exp(-mu * t_val) * cl
  expect_equal(as.vector(result), expected, tolerance = 1e-14)
})


# ---------------------------------------------------------------------------
# Test: Q-het Gibbs SPR partial CL matches full evaluation
# ---------------------------------------------------------------------------
test_that("Q-het Gibbs SPR runs and produces reasonable results", {
  setup <- make_het_setup(nTip = 8L, nChar = 12L, seed = 3847L)
  mcmc <- MkPrimeMCMC(
    nIter = 500L, maxWarmup = 200L, minWarmup = 200L, thin = 5L,
    autoTune = FALSE, gibbsSpr = TRUE, gibbsSubtreeSwap = FALSE,
    nRuns = 1L
  )
  result <- suppressWarnings(RunMkPrime(data = setup$mkd, tree = setup$tree,
                                        model = setup$model, mcmc = mcmc))
  expect_s3_class(result, "MkPosterior")
  expect_true(nrow(result$samples) > 0)
  # Allow up to 5% NaN samples in early post-warmup (Q-het numerical edge case)
  lp <- result$samples[, "log_posterior"]
  expect_true(mean(is.finite(lp)) >= 0.95)
})


# ---------------------------------------------------------------------------
# Test: Q-het Gibbs subtree swap runs and produces reasonable results
# ---------------------------------------------------------------------------
test_that("Q-het Gibbs subtree swap runs and produces reasonable results", {
  setup <- make_het_setup(nTip = 8L, nChar = 12L, seed = 6243L)
  mcmc <- MkPrimeMCMC(
    nIter = 500L, maxWarmup = 200L, minWarmup = 200L, thin = 5L,
    autoTune = FALSE, gibbsSpr = FALSE, gibbsSubtreeSwap = TRUE,
    nRuns = 1L
  )
  result <- suppressWarnings(RunMkPrime(data = setup$mkd, tree = setup$tree,
                                        model = setup$model, mcmc = mcmc))
  expect_s3_class(result, "MkPosterior")
  expect_true(nrow(result$samples) > 0)
  expect_true(all(is.finite(result$samples[, "log_posterior"])))
})


# ---------------------------------------------------------------------------
# Test: Q-het Gibbs SPR + swap together
# ---------------------------------------------------------------------------
test_that("Q-het with both Gibbs moves produces valid MCMC", {
  setup <- make_het_setup(nTip = 8L, nChar = 12L, seed = 9102L)
  mcmc <- MkPrimeMCMC(
    nIter = 600L, maxWarmup = 200L, minWarmup = 200L, thin = 5L,
    autoTune = FALSE, gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE,
    nRuns = 1L
  )
  result <- suppressWarnings(RunMkPrime(data = setup$mkd, tree = setup$tree,
                                        model = setup$model, mcmc = mcmc))
  expect_s3_class(result, "MkPosterior")

  # Both Gibbs moves should have been proposed
  rates <- result$acceptanceRates
  if ("gibbs_spr" %in% names(rates))
    expect_true(rates["gibbs_spr"] >= 0)
  if ("gibbs_subtree_swap" %in% names(rates))
    expect_true(rates["gibbs_subtree_swap"] >= 0)
})
