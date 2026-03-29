test_that("MkPrimeMCMC() defaults to thin = 'auto'", {
  mcmc <- MkPrimeMCMC()
  expect_identical(mcmc$thin, "auto")
})

test_that("MkPrimeMCMC() accepts explicit integer thin", {
  mcmc <- MkPrimeMCMC(thin = 5L)
  expect_identical(mcmc$thin, 5L)

  mcmc2 <- MkPrimeMCMC(thin = 20)
  expect_identical(mcmc2$thin, 20L)
})

test_that("MkPrimeMCMC() rejects invalid thin values", {
  expect_error(MkPrimeMCMC(thin = 0L), "positive integer")
  expect_error(MkPrimeMCMC(thin = -1L), "positive integer")
  expect_error(suppressWarnings(MkPrimeMCMC(thin = "bogus")), "positive integer")
})

test_that("auto thin resolves to number of active moves", {
  skip_if_not_installed("TreeTools")
  library(TreeTools)
  set.seed(4271)
  nTip <- 8L
  nChar <- 10L
  mat <- matrix(sample(0:2, nTip * nChar, replace = TRUE), nTip, nChar,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  tree <- TreeTools::NJTree(pd, edgeLengths = TRUE)

  # Default config: gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE, tbr = TRUE
  # Expected moves: tree_length, branch_lengths, nni, spr,
  #   gibbs_spr, gibbs_subtree_swap, tbr, kPrime, p, rate_log_sd = 10
  mcmc_auto <- suppressWarnings(MkPrimeMCMC(
    nIter = 1500L, maxWarmup = 200L, minWarmup = 200L,
    autoTune = FALSE, nRuns = 1L
  ))
  nMoves <- 10L

  # Explicit thin matching the expected auto resolution
  mcmc_explicit <- suppressWarnings(MkPrimeMCMC(
    nIter = 1500L, maxWarmup = 200L, minWarmup = 200L,
    autoTune = FALSE, nRuns = 1L, thin = nMoves
  ))

  res_auto     <- RunMkPrime(pd, tree, mcmc = mcmc_auto)
  res_explicit <- RunMkPrime(pd, tree, mcmc = mcmc_explicit)

  expect_equal(nrow(res_auto$samples), nrow(res_explicit$samples))
  expect_gt(nrow(res_auto$samples), 0L)
})

test_that("auto thin changes with move configuration", {
  skip_if_not_installed("TreeTools")
  library(TreeTools)
  set.seed(6183)
  nTip <- 8L
  nChar <- 10L
  mat <- matrix(sample(0:2, nTip * nChar, replace = TRUE), nTip, nChar,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  tree <- TreeTools::NJTree(pd, edgeLengths = TRUE)

  # Minimal moves: no Gibbs, no TBR
  # Expected: tree_length, branch_lengths, nni, spr,
  #   kPrime, p, rate_log_sd = 7
  mcmc_min <- suppressWarnings(MkPrimeMCMC(
    nIter = 1500L, maxWarmup = 200L, minWarmup = 200L,
    autoTune = FALSE, nRuns = 1L,
    gibbsSpr = FALSE, gibbsSubtreeSwap = FALSE, tbr = FALSE
  ))

  # Full moves: gibbs + TBR
  # Expected: tree_length, branch_lengths, nni, spr,
  #   gibbs_spr, gibbs_subtree_swap, tbr, kPrime, p, rate_log_sd = 10
  mcmc_max <- suppressWarnings(MkPrimeMCMC(
    nIter = 1500L, maxWarmup = 200L, minWarmup = 200L,
    autoTune = FALSE, nRuns = 1L,
    gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE, tbr = TRUE
  ))

  res_min <- RunMkPrime(pd, tree, mcmc = mcmc_min)
  res_max <- RunMkPrime(pd, tree, mcmc = mcmc_max)

  # Fewer moves = thinner thinning = more samples per iteration
  expect_gt(nrow(res_min$samples), nrow(res_max$samples))
  expect_gt(nrow(res_min$samples), 0L)
  expect_gt(nrow(res_max$samples), 0L)
})
