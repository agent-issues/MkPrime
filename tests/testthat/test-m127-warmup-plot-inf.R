# Regression test for M-127: .PlotWarmupPanel() crash when nIter = Inf
#
# .PlotWarmupPanel() passed Inf to xlim, causing
# "Error in plot.window(...) : need finite 'xlim' values"

test_that(".PlotWarmupPanel handles nIter = Inf", {
  # Seed the warmup history so .PlotWarmupPanel has data to plot
  MkPrime:::.ResetEssHistory()
  env <- MkPrime:::.tracePlotEnv
  env$warmupIter    <- c(100L, 200L, 300L)
  env$warmupLogPost <- matrix(c(-500, -490, -480, -510, -495, -485),
                              ncol = 2)

  pdf(NULL)  # null device — no file written
  on.exit(dev.off(), add = TRUE)

  # This used to error: "need finite 'xlim' values"
  expect_no_error(
    MkPrime:::.PlotWarmupPanel(nRuns = 2L,
                               colors = c("blue", "red"),
                               iter = 300L,
                               nIter = Inf)
  )
})

test_that(".PlotWarmupPanel still works with finite nIter", {
  MkPrime:::.ResetEssHistory()
  env <- MkPrime:::.tracePlotEnv
  env$warmupIter    <- c(100L, 200L)
  env$warmupLogPost <- matrix(c(-500, -490), ncol = 1)

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_no_error(
    MkPrime:::.PlotWarmupPanel(nRuns = 1L,
                               colors = "blue",
                               iter = 200L,
                               nIter = 10000L)
  )
})
