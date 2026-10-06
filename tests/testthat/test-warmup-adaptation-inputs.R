# What the run loop feeds its adaptation controllers (#79, #83).

.AdaptInputsRun <- function(mcmc) {
  tree <- read.tree(text = paste0(
    "(((t1:0.05,t2:0.07):0.04,(t3:0.06,t4:0.05):0.03):0.05,",
    "((t5:0.04,t6:0.06):0.05,(t7:0.05,t8:0.04):0.06):0.04);"
  ))
  mat <- matrix(c(0, 0, 1, 1, 0, 1, 0, 1,
                  0, 1, 1, 0, 1, 0, 0, 1,
                  0, 1, 2, 2, 1, 1, 1, 0,
                  1, 1, 0, 1, 0, 0, 1, 0),
                nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  allow_warning(
    RunMkPrime(MatrixToPhyDat(mat), tree, neomorphic = 1L, mcmc = mcmc),
    "without stabilisation|minimum warmup|tuning"
  )
}

test_that("warmup step sizes and slice widths read one batch's counts (#79)", {
  skip_on_cran()
  seen <- new.env()
  seen$tuning <- seen$slice <- numeric(0)
  adaptTuning <- MkPrime:::.AdaptTuning
  adaptSlice <- MkPrime:::.AdaptSliceWidths
  local_mocked_bindings(
    .AdaptTuning = function(tuning, acceptCount, proposeCount, moves) {
      seen$tuning <- c(seen$tuning, sum(proposeCount))
      adaptTuning(tuning, acceptCount, proposeCount, moves)
    },
    .AdaptSliceWidths = function(tuning, proposeCount, sliceExpCount, moves) {
      seen$slice <- c(seen$slice, sum(proposeCount))
      adaptSlice(tuning, proposeCount, sliceExpCount, moves)
    }
  )
  set.seed(79)
  .AdaptInputsRun(MkPrimeMCMC(nIter = 1600L, thin = 10L, minWarmup = 1500L,
                              maxWarmup = 1500L, autoTune = FALSE, nRuns = 1L,
                              nChains = 1L, nCore = 1L, maxTime = 60))
  # Three 500-iteration warmup batches, one proposal per iteration.
  expect_equal(seen$tuning, c(500, 500, 500))
  expect_equal(seen$slice, c(500, 500, 500))
})

test_that("joint-proposal correlations freeze before the bandit (#83, #301)", {
  skip_on_cran()
  seen <- new.env()
  seen$calls <- character(0)
  seen$rows <- integer(0)
  estimate <- MkPrime:::.EstimateJointRhos
  minEssRate <- MkPrime:::.MinEssRate
  local_mocked_bindings(
    .EstimateJointRhos = function(samples, hasNeo, prior = NULL) {
      seen$calls <- c(seen$calls, "rho")
      seen$rows <- c(seen$rows, nrow(samples))
      estimate(samples, hasNeo, prior)
    },
    .MinEssRate = function(...) {
      seen$calls <- c(seen$calls, "bandit")
      minEssRate(...)
    }
  )
  # thin = 1 fills a 500-sample Tuning window in a single batch.
  set.seed(83)
  .AdaptInputsRun(MkPrimeMCMC(nIter = 4000L, thin = 1L, minWarmup = 1000L,
                              maxWarmup = 1000L, autoTune = TRUE,
                              tuningBudget = 1000L, tuningRounds = 1L,
                              nRuns = 1L, nChains = 1L, nCore = 1L,
                              maxTime = 60))
  # Warmup estimates once per batch, from one state snapshot per batch; Tuning
  # estimates once, from its first window, before any window is scored.
  expect_equal(seen$rows, c(1L, 2L, 500L))
  expect_equal(seen$calls[1:4], c("rho", "rho", "rho", "bandit"))
  expect_false("rho" %in% seen$calls[-(1:3)])
})


test_that("Tuning ends within half the iterations left, so Sample runs (#402)", {
  skip_on_cran()
  skip_under_memcheck()
  set.seed(402)
  res <- .AdaptInputsRun(MkPrimeMCMC(nIter = 5000L, nRuns = 1L, nCore = 1L,
                                     maxTime = 120))
  # Warmup ends by maxWarmup = 2500, so Tuning stops within a 500-iteration
  # batch of 1250 iterations, leaving Sample at least 750.
  expect_gte(nrow(res$samples) * res$mcmc$thin, 750L - res$mcmc$thin)
})
