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

test_that("warmup controllers read the counts since their last update (#79, #304)", {
  skip_on_cran()
  seen <- new.env()
  seen$tuning <- seen$slice <- list()
  adaptTuning <- MkPrime:::.AdaptTuning
  adaptSlice <- MkPrime:::.AdaptSliceWidths
  local_mocked_bindings(
    .AdaptTuning = function(tuning, acceptCount, proposeCount, moves) {
      seen$tuning[[length(seen$tuning) + 1L]] <- proposeCount
      seen$moves <- moves
      adaptTuning(tuning, acceptCount, proposeCount, moves)
    },
    .AdaptSliceWidths = function(tuning, proposeCount, sliceExpCount, moves) {
      seen$slice[[length(seen$slice) + 1L]] <- proposeCount
      adaptSlice(tuning, proposeCount, sliceExpCount, moves)
    }
  )
  set.seed(79)
  .AdaptInputsRun(MkPrimeMCMC(nIter = 3100L, thin = 10L, minWarmup = 3000L,
                              maxWarmup = 3000L, autoTune = FALSE, nRuns = 1L,
                              nChains = 1L, nCore = 1L, maxTime = 60))
  # Six 500-iteration warmup batches, one proposal per iteration.
  expect_length(seen$tuning, 6L)
  expect_identical(seen$slice, seen$tuning)
  expect_equal(sum(seen$tuning[[1]]), 500)

  # A move whose window was spent starts afresh: its count is one batch's,
  # so a cumulative count would show here as more than 500 in all.
  spent <- .TuningWindowSpent(seen$tuning[[1]], seen$moves)
  expect_true(any(spent))
  expect_true(any(!spent))
  # A move below the gate carries its count into the next batch rather than
  # discarding it, so a move proposed a few times a batch still adapts.
  for (k in 2:6) {
    prev <- seen$tuning[[k - 1L]]
    carried <- !.TuningWindowSpent(prev, seen$moves)
    expect_true(all(seen$tuning[[k]][carried] >= prev[carried]))
    expect_equal(sum(seen$tuning[[k]]) - sum(prev[carried]), 500)
  }
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


test_that("a run whose Tuning cannot end stops at once, and says so (#402)", {
  skip_on_cran()
  skip_under_memcheck()
  seen <- new.env()
  seen$warnings <- character(0)
  set.seed(402)
  res <- withCallingHandlers(
    .AdaptInputsRun(MkPrimeMCMC(nIter = 5000L, nRuns = 1L, nCore = 1L,
                                maxTime = 120)),
    warning = function(w) {
      seen$warnings <- c(seen$warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # Warmup ends at maxWarmup = 2500; Tuning's first round cannot end in the
  # 2500 left, so the run stops there rather than tune to nIter.
  expect_equal(nrow(res$samples), 0L)
  expect_identical(res$stop_reason, "too_short")
  expect_equal(res$actual_iter, 2500L)
  expect_match(seen$warnings,
               "Tuning needs \\d+ more iterations .* only 2500 remain",
               all = FALSE)
})


test_that("a default run too short to tune stops, rather than skip Tuning (#402)", {
  skip_on_cran()
  set.seed(4021)
  expect_warning(
    res <- .AdaptInputsRun(MkPrimeMCMC(nIter = 1000L, nRuns = 1L, nCore = 1L,
                                       maxTime = 60)),
    "Tuning needs \\d+ more iterations"
  )
  expect_identical(res$stop_reason, "too_short")
  expect_equal(nrow(res$samples), 0L)
})
