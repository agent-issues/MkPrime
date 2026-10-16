test_that(".FormatIterLabel uses k / M suffixes", {
  expect_identical(
    MkPrime:::.FormatIterLabel(c(500, 1000, 100000, 1100000, 2e6)),
    c("500", "1k", "100k", "1.1M", "2M")
  )
  expect_identical(MkPrime:::.FormatIterLabel(NA_real_), NA_character_)
})

.FakeTrees <- function(n = 15L) {
  set.seed(1)
  structure(lapply(seq_len(n), function(i) ape::rtree(8, br = NULL)),
            class = "multiPhylo")
}

test_that(".RecentConsensus returns a tree from the last 12", {
  skip_if_not_installed("TreeTools")
  trees <- lapply(seq_len(15), function(i) TreeTools::BalancedTree(8))
  cons <- MkPrime:::.RecentConsensus(trees, nMax = 12L)
  expect_s3_class(cons, "phylo")
  expect_equal(TreeTools::NTip(cons), 8L)
  expect_null(MkPrime:::.RecentConsensus(list()))
})

test_that("MkpTracePlot draws topology panel and tree ESS", {
  skip_if_not_installed("TreeTools")
  MkPrime:::.ResetEssHistory()
  MkPrime:::.ResetWarmupHistory()
  env <- MkPrime:::.tracePlotEnv
  env$lastIter <- 0L
  set.seed(2)
  samp <- cbind(log_posterior = rnorm(200), tree_length = rnorm(200, 5))
  info <- list(
    iter = 1000L, nIter = 5000L, warmup = 100L, inWarmup = FALSE,
    phase = "Sample", nRuns = 1L, nChains = 1L,
    runSamples = list(samp), currentState = NULL,
    recentAcceptance = 0.3, elapsed = 12, paramNames = colnames(samp),
    recentTrees = lapply(seq_len(15), function(i) TreeTools::BalancedTree(8)),
    treeEss = 42
  )
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(MkpTracePlot(info))
  expect_equal(
    unname(MkPrime:::.tracePlotEnv$essHistory[[1]]["tree_ess"]), 42
  )
  info$recentTrees <- NULL
  info$treeEss <- NULL
  info$iter <- 2000L
  expect_no_error(MkpTracePlot(info))
})
