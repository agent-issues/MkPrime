# Tests for adaptive tree ESS in convergence monitoring and tuning.
# These tests mock MCMC state directly and do not require full MCMC runs.

test_that("minTreeEss parameter accepted by MkPrimeMCMC", {
  mcmc <- suppressWarnings(MkPrimeMCMC(nIter = 100L, minTreeEss = 200))
  expect_equal(mcmc$minTreeEss, 200)

  mcmc2 <- suppressWarnings(MkPrimeMCMC(nIter = 100L))
  expect_null(mcmc2$minTreeEss)
})


test_that(".ComputeTreeEssInLoop returns NA when no trees saved", {
  run <- list(tree_saved_idx = 0L, tree_samples = list())
  result <- MkPrime:::.ComputeTreeEssInLoop(list(run), 500L, FALSE)
  expect_true(is.na(result))
})


test_that(".ComputeTreeEssInLoop returns NA when too few trees", {
  trees <- replicate(3, ape::rtree(5, rooted = FALSE), simplify = FALSE)
  run <- list(tree_saved_idx = 3L, tree_samples = trees)
  result <- MkPrime:::.ComputeTreeEssInLoop(list(run), 500L, FALSE)
  expect_true(is.na(result))
})


test_that(".ComputeTreeEssInLoop returns finite ESS with sufficient trees", {
  skip_if_not_installed("TreeDist")
  set.seed(6718)
  trees <- replicate(50, ape::rtree(10, rooted = FALSE), simplify = FALSE)
  run <- list(tree_saved_idx = 50L, tree_samples = trees)
  result <- MkPrime:::.ComputeTreeEssInLoop(list(run), 500L, FALSE)
  expect_true(is.finite(result))
  expect_true(result > 0)
})


test_that(".ComputeTreeEssInLoop subsamples to maxPerRun", {
  skip_if_not_installed("TreeDist")
  set.seed(4207)
  nTrees <- 200L
  trees <- replicate(nTrees, ape::rtree(8, rooted = FALSE), simplify = FALSE)
  run <- list(tree_saved_idx = nTrees, tree_samples = trees)

  # Fine (1000) vs coarse (50) — both should return finite
  fine   <- MkPrime:::.ComputeTreeEssInLoop(list(run), 1000L, FALSE)
  coarse <- MkPrime:::.ComputeTreeEssInLoop(list(run), 50L, FALSE)
  expect_true(is.finite(fine))
  expect_true(is.finite(coarse))
})


test_that(".ComputeTreeEssInLoop handles streaming mode with NULL slots", {
  skip_if_not_installed("TreeDist")
  set.seed(3384)
  # Pre-allocated list with some NULL slots at the end
  trees <- replicate(30, ape::rtree(8, rooted = FALSE), simplify = FALSE)
  trees <- c(trees, vector("list", 20))

  run <- list(tree_samples = trees)  # streaming: no saved_idx
  result <- MkPrime:::.ComputeTreeEssInLoop(list(run), 500L, TRUE)
  expect_true(is.finite(result))
})


test_that(".CheckConvergence adaptive tier: skip when scalars far off", {
  set.seed(9201)
  nSamp <- 50L
  trees <- replicate(nSamp, ape::rtree(10, rooted = FALSE), simplify = FALSE)

  # Highly autocorrelated scalar -> low ESS, far below threshold
  run <- list(
    saved_idx = nSamp,
    tree_saved_idx = nSamp,
    samples = matrix(
      cumsum(rnorm(nSamp)),
      nSamp, 1,
      dimnames = list(NULL, "log_posterior")
    ),
    tree_samples = trees
  )

  mcmc <- list(
    minEss = 200, maxRhat = NULL, minTreeEss = 100,
    warmup = 0L, maxWarmup = 0L
  )

  result <- MkPrime:::.CheckConvergence(list(run), "log_posterior", mcmc)
  expect_equal(result$treeEssPrecision, "skip")
  expect_true(is.na(result$treeEss))
  expect_false(result$converged)
})


test_that(".CheckConvergence adaptive tier: fine when scalars converged", {
  skip_if_not_installed("TreeDist")
  set.seed(3812)
  nSamp <- 100L
  trees <- replicate(nSamp, ape::rtree(10, rooted = FALSE), simplify = FALSE)

  # Independent samples -> high ESS, exceeds minEss
  run <- list(
    saved_idx = nSamp,
    tree_saved_idx = nSamp,
    samples = matrix(
      rnorm(nSamp),
      nSamp, 1,
      dimnames = list(NULL, "log_posterior")
    ),
    tree_samples = trees
  )

  mcmc <- list(
    minEss = 50, maxRhat = NULL, minTreeEss = 100,
    warmup = 0L, maxWarmup = 0L
  )

  result <- MkPrime:::.CheckConvergence(list(run), "log_posterior", mcmc)
  # Scalar ESS > minEss -> should compute at fine (or upgrade from coarse)
  expect_true(result$treeEssPrecision %in% c("fine", "coarse"))
  expect_true(is.finite(result$treeEss))
})


test_that(".CheckConvergence includes treeEss in convergence decision", {
  skip_if_not_installed("TreeDist")
  set.seed(5591)
  nSamp <- 100L
  trees <- replicate(nSamp, ape::rtree(10, rooted = FALSE), simplify = FALSE)

  run <- list(
    saved_idx = nSamp,
    tree_saved_idx = nSamp,
    samples = matrix(
      rnorm(nSamp),
      nSamp, 1,
      dimnames = list(NULL, "log_posterior")
    ),
    tree_samples = trees
  )

  # minTreeEss very high -> should prevent convergence
  mcmc_high <- list(minEss = 10, maxRhat = NULL, minTreeEss = 10000,
                    warmup = 0L, maxWarmup = 0L)
  result_high <- MkPrime:::.CheckConvergence(list(run), "log_posterior",
                                              mcmc_high)
  expect_false(result_high$converged)

  # minTreeEss very low -> should allow convergence
  mcmc_low <- list(minEss = 10, maxRhat = NULL, minTreeEss = 1,
                   warmup = 0L, maxWarmup = 0L)
  result_low <- MkPrime:::.CheckConvergence(list(run), "log_posterior",
                                             mcmc_low)
  expect_true(result_low$converged)
})


test_that(".CheckConvergence returns treeEss fields when minTreeEss NULL", {
  nSamp <- 30L
  run <- list(
    saved_idx = nSamp,
    tree_saved_idx = 0L,
    samples = matrix(
      rnorm(nSamp),
      nSamp, 1,
      dimnames = list(NULL, "log_posterior")
    ),
    tree_samples = list()
  )

  mcmc <- list(minEss = 10, maxRhat = NULL, minTreeEss = NULL,
               warmup = 0L, maxWarmup = 0L)
  result <- MkPrime:::.CheckConvergence(list(run), "log_posterior", mcmc)
  expect_true(is.na(result$treeEss))
  expect_equal(result$treeEssPrecision, "skip")
})


test_that(".CheckConvergenceFromLogs returns treeEss = NA with skip", {
  # Log-based convergence can't compute tree ESS (no trees on disk)
  skip_if_not_installed("coda")
  tmpFile <- tempfile(fileext = ".log")
  on.exit(unlink(tmpFile), add = TRUE)

  nSamp <- 30L
  mat <- matrix(rnorm(nSamp), nSamp, 1,
                dimnames = list(NULL, "log_posterior"))
  write.table(mat, tmpFile, sep = "\t", row.names = FALSE, quote = FALSE)

  mcmc <- list(minEss = 10, maxRhat = NULL, minTreeEss = 100,
               warmup = 0L, maxWarmup = 0L)
  result <- MkPrime:::.CheckConvergenceFromLogs(tmpFile, "log_posterior", mcmc)
  if (!is.null(result)) {
    expect_true(is.na(result$treeEss))
    expect_equal(result$treeEssPrecision, "skip")
  }
})


test_that(".MinEssPerSec includes tree ESS when provided", {
  skip_if_not_installed("TreeDist")
  set.seed(2847)

  nSamp <- 50L
  mat <- matrix(rnorm(nSamp * 2), nSamp, 2,
                dimnames = list(NULL, c("log_posterior", "tree_length")))

  # Without trees
  ess_scalar <- MkPrime:::.MinEssPerSec(mat, 1.0)

  # With identical trees (ESS ~ 1) — should pull min-ESS down
  identicalTrees <- replicate(nSamp, ape::rtree(8, rooted = FALSE),
                              simplify = FALSE)
  # Make them all the same tree to get low tree ESS
  singleTree <- identicalTrees[[1]]
  sameTrees <- replicate(nSamp, singleTree, simplify = FALSE)

  ess_with_trees <- MkPrime:::.MinEssPerSec(mat, 1.0,
                                             tuningTrees = sameTrees)
  # Tree ESS for identical trees is ~1, so it should pull down the min
  expect_true(ess_with_trees <= ess_scalar)
})


test_that(".TickerSummaryStr includes tree ESS when present", {
  diagCheck <- list(
    minEss = 150,
    maxRhat = NA_real_,
    treeEss = 75,
    treeEssPrecision = "coarse"
  )
  s <- MkPrime:::.TickerSummaryStr(diagCheck)
  # Should contain "treeESS" somewhere in the output
  expect_match(cli::ansi_strip(s), "treeESS")
})


test_that(".TickerSummaryStr omits tree ESS when NA", {
  diagCheck <- list(
    minEss = 150,
    maxRhat = NA_real_,
    treeEss = NA_real_,
    treeEssPrecision = "skip"
  )
  s <- MkPrime:::.TickerSummaryStr(diagCheck)
  expect_false(grepl("treeESS", cli::ansi_strip(s)))
})
