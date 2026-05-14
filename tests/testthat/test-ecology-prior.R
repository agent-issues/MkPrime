# Tests for the asymmetric-slab ecology prior (v2).
#
# Pins numerical equality between cpp_log_prior (C++) and LogPrior (R) under a
# non-trivial state with mixed none/encouraged/discouraged z-cells and
# distinct theta_e per non-reference ecology.

library("TreeTools")


# --- Helper: ecology-aware fixture with kEco = 3 -----------------------------
# 6 tips, 5 transformational characters, 1 ecology character with 3 states.
# Reference ecology will be selected by edge-mass inside RunMkPrime; we drive
# the prior directly via .InitMcmcData/.InitMcmcChain so the reference choice
# is whatever the orchestrator computes from the tree + ecology vector.

.MakeAsymPriorFixture <- function() {
  set.seed(2026L)
  tips <- paste0("t", 1:6)
  # 5 transformational characters + 1 ecology character (3 states).
  mat <- matrix(c(
    0, 1, 2, 0, 1, 2,
    1, 0, 1, 2, 2, 0,
    2, 1, 0, 1, 0, 2,
    0, 1, 0, 1, 0, 1,
    1, 0, 1, 0, 1, 0,
    0, 0, 1, 1, 2, 2
  ), nrow = 6L, ncol = 6L, byrow = FALSE,
     dimnames = list(tips, NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, ecology = 6L)
  tree <- Preorder(ape::rtree(6L, tip.label = tips))
  model <- MkPrimeModel(
    ecologyAware = TRUE,
    kPrimePrior = "geometric",
    expSteps = 10L,
    treeLengthRate = 0.5,
    thetaAlpha = 2.5,
    thetaBeta = 1.7
  )
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  list(mkd = mkd, tree = tree, model = model)
}


# --- C++ vs R agreement (default initial state) -----------------------------

test_that("cpp_log_prior matches LogPrior at the .InitState seed", {
  f <- .MakeAsymPriorFixture()
  state <- MkPrime:::.InitState(f$tree, f$mkd, f$model)

  rLogPrior <- LogPrior(state, f$model, f$mkd)
  expect_true(is.finite(rLogPrior))

  mcmcData <- MkPrime:::.InitMcmcData(f$mkd, f$model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  cppLogPrior <- eval_log_prior_cpp(mcmcData, statePtr)

  expect_equal(cppLogPrior, rLogPrior, tolerance = 1e-10)
})


# --- C++ vs R agreement (perturbed state with mixed z) -----------------------

test_that("cpp_log_prior matches LogPrior under asymmetric z and theta", {
  f <- .MakeAsymPriorFixture()
  state <- MkPrime:::.InitState(f$tree, f$mkd, f$model)

  # Inject a mixed z-matrix and distinct theta per non-reference ecology.
  # state$z has dimensions nChar x (kEco - 1) = 6 x 2 (since kEco = 3).
  nCharZ <- nrow(state$z)
  nThetaZ <- ncol(state$z)
  expect_equal(nThetaZ, f$mkd$kEcology - 1L)
  expect_gt(nThetaZ, 0L)

  # Hand-built z patterns with imbalanced enc/disc/none counts per column.
  # Column 1: 2 none, 2 enc, 1 disc. Column 2: 1 none, 1 enc, 3 disc.
  zNew <- matrix(0L, nrow = nCharZ, ncol = nThetaZ)
  zNew[, 1L] <- c(0L, 0L, 1L, 1L, 2L)
  zNew[, 2L] <- c(0L, 1L, 2L, 2L, 2L)
  state$z <- zNew

  # Pick distinct, non-symmetric theta values so the asymmetric-slab kernel
  # really exercises log(theta) and log(1 - theta).
  state$theta <- c(0.30, 0.78)
  state$pi0 <- 0.42

  rLogPrior <- LogPrior(state, f$model, f$mkd)
  expect_true(is.finite(rLogPrior))

  mcmcData <- MkPrime:::.InitMcmcData(f$mkd, f$model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  cppLogPrior <- eval_log_prior_cpp(mcmcData, statePtr)

  expect_equal(cppLogPrior, rLogPrior, tolerance = 1e-10)
})


# --- Slab really is asymmetric (varies with theta) ---------------------------

test_that("LogPrior depends on theta when slab cells are present", {
  f <- .MakeAsymPriorFixture()
  state <- MkPrime:::.InitState(f$tree, f$mkd, f$model)

  # Force one enc and one disc cell in each non-reference column so changing
  # theta moves the log-prior.
  nCharZ <- nrow(state$z)
  nThetaZ <- ncol(state$z)
  zNew <- matrix(0L, nrow = nCharZ, ncol = nThetaZ)
  for (j in seq_len(nThetaZ)) {
    zNew[1L, j] <- 1L  # enc
    zNew[2L, j] <- 2L  # disc
  }
  state$z <- zNew

  state$theta <- rep(0.5, nThetaZ)
  lpSym <- LogPrior(state, f$model, f$mkd)

  state$theta <- rep(0.9, nThetaZ)
  lpAsym <- LogPrior(state, f$model, f$mkd)

  # Both finite, and they must differ: symmetric (0.5) vs asymmetric (0.9)
  # changes the slab log-density by (nEnc - nDisc) * (log theta - log(1-theta))
  # plus the Beta(thetaAlpha, thetaBeta) prior change.
  expect_true(is.finite(lpSym))
  expect_true(is.finite(lpAsym))
  expect_false(isTRUE(all.equal(lpSym, lpAsym)))
})
