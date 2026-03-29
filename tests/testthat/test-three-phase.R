# Tests for three-phase MCMC architecture: Warmup / Tuning / Sample
#
# Covers: .CheckStabilisation(), .MinEssPerSec(), .PerturbMoveWeights(),
# MkPrimeMCMC() three-phase parameters, and phase state machine.


# ==========================================================================
# .CheckStabilisation()
# ==========================================================================

test_that("CheckStabilisation returns FALSE with insufficient data", {
  result <- MkPrime:::.CheckStabilisation(
    logPostHistory = rnorm(5),
    nStableConsecutive = 0L,
    windowSize = 10L
  )
  expect_false(result$stable)
  expect_equal(result$nStableConsecutive, 0L)
})

test_that("CheckStabilisation detects stationary series", {
  set.seed(7461)
  # Flat series should quickly stabilise
  logPost <- rnorm(100, mean = -500, sd = 1)
  nStable <- 0L
  for (i in seq(20, 100, by = 1)) {
    result <- MkPrime:::.CheckStabilisation(
      logPost[1:i], nStable, windowSize = 10L,
      zThreshold = 1.5, nStableRequired = 3L
    )
    nStable <- result$nStableConsecutive
    if (result$stable) break
  }
  expect_true(result$stable)
})

test_that("CheckStabilisation rejects trending series", {
  # Monotonically increasing — never stable
  logPost <- seq(-1000, -500, length.out = 100)
  nStable <- 0L
  stable <- FALSE
  for (i in seq(20, 100, by = 1)) {
    result <- MkPrime:::.CheckStabilisation(
      logPost[1:i], nStable, windowSize = 10L,
      zThreshold = 1.5, nStableRequired = 3L
    )
    nStable <- result$nStableConsecutive
    if (result$stable) {
      stable <- TRUE
      break
    }
  }
  expect_false(stable)
})

test_that("CheckStabilisation resets counter across a level shift", {
  set.seed(3928)
  # Build up consecutive stable count, then insert a level shift.
  # The z-score on a window spanning the shift should reset the counter.
  flat1 <- rnorm(30, -500, 0.5)
  shift <- rnorm(10, -400, 0.5)
  flat2 <- rnorm(30, -400, 0.5)
  logPost <- c(flat1, shift, flat2)

  # Run up to the shift boundary with 2 consecutive stable checks
  nStable <- 2L
  # Now evaluate a window that spans the shift (index ~30-40)
  result <- MkPrime:::.CheckStabilisation(
    logPost[1:40], nStable, windowSize = 10L,
    zThreshold = 1.5, nStableRequired = 3L
  )
  # The level shift should reset the counter

  expect_equal(result$nStableConsecutive, 0L)
  expect_false(result$stable)
})

test_that("CheckStabilisation handles zero-variance windows", {
  # All identical values
  logPost <- rep(-500, 30)
  result <- MkPrime:::.CheckStabilisation(
    logPost, nStableConsecutive = 2L,
    windowSize = 10L, nStableRequired = 3L
  )
  expect_true(result$stable)
  expect_equal(result$nStableConsecutive, 3L)
})


# ==========================================================================
# .MinEssPerSec()
# ==========================================================================

test_that("MinEssPerSec returns NA with too few samples", {
  mat <- matrix(rnorm(18), nrow = 3, ncol = 6)
  colnames(mat) <- paste0("p", 1:6)
  result <- MkPrime:::.MinEssPerSec(mat, 1.0)
  expect_true(is.na(result))
})

test_that("MinEssPerSec returns finite value for valid input", {
  skip_if_not_installed("coda")
  set.seed(5192)
  n <- 200
  mat <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  colnames(mat) <- c("log_posterior", "tree_length", "rate_log_sd")
  result <- MkPrime:::.MinEssPerSec(mat, wallTimeSec = 2.0)
  expect_true(is.finite(result))
  expect_gt(result, 0)
})

test_that("MinEssPerSec excludes kPrime columns", {
  skip_if_not_installed("coda")
  set.seed(8173)
  n <- 200
  mat <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  colnames(mat) <- c("log_posterior", "tree_length", "kPrime_1", "kPrime_2")
  result <- MkPrime:::.MinEssPerSec(mat, wallTimeSec = 2.0)
  expect_true(is.finite(result))
})

test_that("MinEssPerSec returns NA for zero wall time", {
  mat <- matrix(rnorm(200), nrow = 20, ncol = 10)
  colnames(mat) <- paste0("p", 1:10)
  expect_true(is.na(MkPrime:::.MinEssPerSec(mat, 0)))
})


# ==========================================================================
# .PerturbMoveWeights()
# ==========================================================================

test_that("PerturbMoveWeights produces valid weight vectors", {
  set.seed(2847)
  w <- c(nni = 0.25, spr = 0.25, tree_length = 0.25, branch_lengths = 0.25)
  candidates <- MkPrime:::.PerturbMoveWeights(
    w, pinnedWeights = NULL, moveNames = names(w),
    nPerturbations = 5L
  )
  expect_length(candidates, 5L)
  for (cand in candidates) {
    expect_equal(sum(cand), 1.0, tolerance = 1e-10)
    expect_true(all(cand >= 0))
    expect_equal(names(cand), names(w))
  }
})

test_that("PerturbMoveWeights respects pinned weights", {
  set.seed(4619)
  w <- c(nni = 0.3, spr = 0.2, tree_length = 0.3, branch_lengths = 0.2)
  pinned <- c(nni = 0.3)
  candidates <- MkPrime:::.PerturbMoveWeights(
    w, pinnedWeights = pinned, moveNames = names(w),
    nPerturbations = 5L
  )
  for (cand in candidates) {
    expect_equal(unname(cand["nni"]), 0.3, tolerance = 1e-10)
    expect_equal(sum(cand), 1.0, tolerance = 1e-10)
  }
})

test_that("PerturbMoveWeights returns empty list with < 2 free moves", {
  w <- c(nni = 0.5, spr = 0.5)
  pinned <- c(nni = 0.5)
  candidates <- MkPrime:::.PerturbMoveWeights(
    w, pinnedWeights = pinned, moveNames = names(w),
    nPerturbations = 3L
  )
  expect_length(candidates, 0L)
})


# ==========================================================================
# MkPrimeMCMC() three-phase parameters
# ==========================================================================

test_that("MkPrimeMCMC stores three-phase defaults", {
  mcmc <- MkPrimeMCMC(nIter = 10000L)
  expect_equal(mcmc$minWarmup, 2000L)
  expect_equal(mcmc$maxWarmup, 5000L)  # nIter / 2
  expect_true(mcmc$autoTune)
  expect_equal(mcmc$tuningBudget, 10000L)
  expect_equal(mcmc$tuningRounds, 5L)
})

test_that("MkPrimeMCMC handles Inf nIter defaults", {
  mcmc <- MkPrimeMCMC(nIter = Inf)
  expect_equal(mcmc$maxWarmup, 50000L)
  expect_equal(mcmc$warmup, 50000L)
})

test_that("MkPrimeMCMC deprecated warmup maps to maxWarmup", {
  expect_warning(
    mcmc <- MkPrimeMCMC(nIter = 10000L, warmup = 3000L),
    "deprecated"
  )
  expect_equal(mcmc$maxWarmup, 3000L)
  expect_equal(mcmc$warmup, 3000L)
})

test_that("MkPrimeMCMC warns when minWarmup > maxWarmup", {
  expect_warning(
    mcmc <- MkPrimeMCMC(nIter = 10000L, minWarmup = 8000L, maxWarmup = 3000L),
    "exceeds"
  )
  expect_equal(mcmc$minWarmup, 3000L)
})

test_that("MkPrimeMCMC rejects invalid tuningBudget", {
  expect_error(
    MkPrimeMCMC(nIter = 10000L, tuningBudget = 0L),
    "tuningBudget"
  )
})

test_that("MkPrimeMCMC rejects invalid tuningRounds", {
  expect_error(
    MkPrimeMCMC(nIter = 10000L, tuningRounds = 0L),
    "tuningRounds"
  )
})

test_that("MkPrimeMCMC autoTune = FALSE skips validation", {
  # Should not error even with zero budget when autoTune is FALSE
  mcmc <- MkPrimeMCMC(nIter = 10000L, autoTune = FALSE, tuningBudget = 0L)
  expect_false(mcmc$autoTune)
})
