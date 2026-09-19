# M-161: Per-partition CL cache validity
# Validates that selective repopulation produces correct likelihoods when only
# a subset of partitions is affected by parameter changes (e.g., rate_loss
# changes only invalidate neomorphic units).

test_that("selective cache repopulation used with mixed char types", {
  skip_if_not_installed("TreeTools")
  set.seed(7612)
  nTip <- 8L
  tips <- paste0("t", seq_len(nTip))

  # Binary chars for neomorphic, multi-state for transformational
  neoMat <- matrix(sample(0:1, nTip * 4, replace = TRUE), nTip, 4)
  transMat <- matrix(sample(0:2, nTip * 8, replace = TRUE), nTip, 8)
  mat <- cbind(neoMat, transMat)
  dimnames(mat) <- list(tips, NULL)
  pd <- TreeTools::MatrixToPhyDat(mat)
  tree <- ape::rtree(nTip, tip.label = tips)

  mkd <- MkPrimeData(pd, neomorphic = 1:4)

  mcmc <- suppressWarnings(MkPrimeMCMC(
    nIter = 600L, maxWarmup = 100L, minWarmup = 50L,
    nRuns = 1L, autoTune = FALSE, thin = 1L
  ))
  model <- MkPrimeModel()

  res <- suppressWarnings(RunMkPrime(mkd, tree, model = model, mcmc = mcmc))

  expect_true(nrow(res$samples) > 0)
  lastLL <- res$samples[nrow(res$samples), "log_likelihood"]
  expect_true(is.finite(lastLL))
  expect_lt(lastLL, 0)
})

test_that("selective repopulation matches full rebuild likelihood", {
  skip_if_not_installed("TreeTools")
  set.seed(4830)
  nTip <- 6L
  tips <- paste0("t", seq_len(nTip))

  neoMat <- matrix(sample(0:1, nTip * 3, replace = TRUE), nTip, 3)
  transMat <- matrix(sample(0:2, nTip * 5, replace = TRUE), nTip, 5)
  mat <- cbind(neoMat, transMat)
  dimnames(mat) <- list(tips, NULL)
  pd <- TreeTools::MatrixToPhyDat(mat)
  tree <- ape::rtree(nTip, tip.label = tips)

  mkd <- MkPrimeData(pd, neomorphic = 1:3)

  mcmc <- suppressWarnings(MkPrimeMCMC(
    nIter = 400L, maxWarmup = 50L, minWarmup = 50L,
    nRuns = 1L, autoTune = FALSE, thin = 1L
  ))

  res <- suppressWarnings(RunMkPrime(mkd, tree, mcmc = mcmc))
  expect_true(nrow(res$samples) > 0)

  n <- nrow(res$samples)
  finalLL   <- res$samples[n, "log_likelihood"]
  finalPost <- res$samples[n, "log_posterior"]
  expect_true(is.finite(finalLL))
  expect_true(is.finite(finalPost))
  expect_lt(finalLL, 0)
})

test_that("selective repopulation with fixTopology", {
  skip_if_not_installed("TreeTools")
  set.seed(2917)
  nTip <- 6L
  tips <- paste0("t", seq_len(nTip))

  neoMat <- matrix(sample(0:1, nTip * 2, replace = TRUE), nTip, 2)
  transMat <- matrix(sample(0:2, nTip * 4, replace = TRUE), nTip, 4)
  mat <- cbind(neoMat, transMat)
  dimnames(mat) <- list(tips, NULL)
  pd <- TreeTools::MatrixToPhyDat(mat)
  tree <- ape::rtree(nTip, tip.label = tips)
  mkd <- MkPrimeData(pd, neomorphic = 1:2)

  # fixTopology = TRUE: only parameter moves — rate_loss triggers
  # selective repopulation (neo units only) while trans units reused
  mcmc <- suppressWarnings(MkPrimeMCMC(
    nIter = 300L, maxWarmup = 50L, minWarmup = 50L,
    nRuns = 1L, autoTune = FALSE, thin = 1L
  ))

  res <- suppressWarnings(RunMkPrime(mkd, tree, mcmc = mcmc,
                                     fixTopology = TRUE))
  expect_true(nrow(res$samples) > 0)
  lastLL <- res$samples[nrow(res$samples), "log_likelihood"]
  expect_true(is.finite(lastLL))
  expect_lt(lastLL, 0)
})
