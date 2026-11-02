test_that(".FormatEssPerSec never shows a negative or missing rate", {
  expect_identical(MkPrime:::.FormatEssPerSec(-3), "\u2014")
  expect_identical(MkPrime:::.FormatEssPerSec(0), "\u2014")
  expect_identical(MkPrime:::.FormatEssPerSec(NA_real_), "\u2014")
  expect_identical(MkPrime:::.FormatEssPerSec(Inf), "\u2014")
  expect_identical(MkPrime:::.FormatEssPerSec(c(1, 2)), "\u2014")
  expect_identical(MkPrime:::.FormatEssPerSec(1.234), "1.23")
})

test_that(".RecentRunTrees pools the latest stored trees, skipping empty slots", {
  trs <- lapply(1:20, function(i) TreeTools::BalancedTree(6))
  run <- list(tree_samples = c(trs, vector("list", 5L)), tree_saved_idx = 20L)
  expect_length(MkPrime:::.RecentRunTrees(list(run)), 12L)
  expect_length(MkPrime:::.RecentRunTrees(list(run, run)), 12L)
  expect_null(MkPrime:::.RecentRunTrees(list(list(tree_saved_idx = 0L))))
  expect_null(MkPrime:::.RecentRunTrees(list(list(tree_samples = NULL))))
})

test_that(".BuildProgressInfo carries recentTrees and treeEss for the plot", {
  trs <- lapply(1:5, function(i) TreeTools::BalancedTree(6))
  run <- list(saved_idx = 0L, tree_samples = trs, tree_saved_idx = 5L,
              chainStates = list(NULL))
  testthat::local_mocked_bindings(
    get_mcmc_state = function(...) list(logLik = 0, logPrior = 0,
                                        treeLength = 1, rateLogSd = 1, p = 0.5),
    .package = "MkPrime"
  )
  mcmc <- list(nIter = 10L, warmup = 1L, nChains = 1L, logFile = NULL)
  info <- MkPrime:::.BuildProgressInfo(list(run), 5L, mcmc, proc.time()["elapsed"],
                                       0.2, "p", treeEss = 42)
  expect_length(info$recentTrees, 5L)
  expect_equal(info$treeEss, 42)
  info0 <- MkPrime:::.BuildProgressInfo(list(run), 5L, mcmc,
                                        proc.time()["elapsed"], 0.2, "p")
  expect_true(is.na(info0$treeEss))
})
