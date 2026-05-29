# EBE Phase 2b: the z prior/Gibbs are restricted to NEOMORPHIC characters.
#
# Under EBE, ecology tilts only neomorphic equilibria; transformational/known
# characters get no ecology effect in the likelihood (R8).  Their z rows must
# therefore NOT contribute to the shared pi0/theta (sparsity/balance) prior —
# otherwise each transformational character injects a phantom z observation that
# biases the pi0 posterior and inflates the "fraction flagged" readout.
#
# This is a RED-BEFORE-GREEN discriminator (advisor 2026-05-29): it MUST fail on
# pre-Phase-2b HEAD (transformational z rows counted) and pass after the
# neomorphic-mask edits in LogPrior (R), cpp_log_prior (C++), and
# gibbs_z_sweep_impl (C++).  A neo-only matrix would pass both ways and prove
# nothing, so the fixture is deliberately MIXED (neomorphic + transformational).

library("TreeTools")


# --- Mixed fixture: 2 neomorphic + 3 transformational + 1 ecology (kEco = 3) --
.MakeMixedNeoTransFixture <- function() {
  set.seed(2026L)
  tips <- paste0("t", 1:6)
  # chars 1-2 neomorphic binary (0/1, variable); chars 3-5 transformational;
  # char 6 = ecology (3 states, each twice).
  mat <- matrix(c(
    0, 0, 1, 1, 0, 1,   # 1 neomorphic
    1, 0, 0, 1, 1, 0,   # 2 neomorphic
    0, 1, 2, 0, 1, 2,   # 3 transformational
    2, 1, 0, 1, 0, 2,   # 4 transformational
    0, 1, 0, 1, 0, 1,   # 5 transformational (binary, still trans)
    0, 0, 1, 1, 2, 2    # 6 ecology
  ), nrow = 6L, ncol = 6L, byrow = FALSE, dimnames = list(tips, NULL))
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = c(1L, 2L), ecology = 6L)
  tree  <- Preorder(ape::rtree(6L, tip.label = tips))
  model <- MkPrimeModel(
    ecologyAware = TRUE, kPrimePrior = "geometric", expSteps = 10L,
    treeLengthRate = 0.5, thetaAlpha = 2.5, thetaBeta = 1.7
  )
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  list(mkd = mkd, tree = tree, model = model)
}

# Evaluate both the R and C++ log-priors for a given R state.
.bothLogPriors <- function(state, f) {
  r   <- LogPrior(state, f$model, f$mkd)
  dat <- MkPrime:::.InitMcmcData(f$mkd, f$model)
  ptr <- MkPrime:::.InitMcmcChain(state)
  list(r = r, cpp = eval_log_prior_cpp(dat, ptr))
}


test_that("fixture really is mixed (has both neomorphic and transformational)", {
  f <- .MakeMixedNeoTransFixture()
  expect_gt(sum(f$mkd$type == "neomorphic"), 0L)
  expect_gt(sum(f$mkd$type == "transformational"), 0L)
  expect_equal(f$mkd$kEcology, 3L)
})


test_that("log-prior is INVARIANT to transformational z rows (R and C++)", {
  f <- .MakeMixedNeoTransFixture()
  state <- MkPrime:::.InitState(f$tree, f$mkd, f$model)

  neoRows   <- which(f$mkd$type == "neomorphic")
  transRows <- which(f$mkd$type == "transformational")
  expect_gt(length(neoRows), 0L)
  expect_gt(length(transRows), 0L)

  nThetaZ <- ncol(state$z)
  state$theta <- if (nThetaZ == 2L) c(0.30, 0.78) else rep(0.5, nThetaZ)
  state$pi0   <- 0.42

  # Baseline: a fixed nontrivial pattern on the NEOMORPHIC rows; trans rows = 0.
  zBase <- matrix(0L, nrow = nrow(state$z), ncol = nThetaZ)
  zBase[neoRows[1L], ] <- 1L                       # encouraged
  if (length(neoRows) > 1L) zBase[neoRows[2L], ] <- 2L  # discouraged

  stateA <- state; stateA$z <- zBase
  # Perturbed: identical NEOMORPHIC rows, but transformational rows set != 0.
  zPerturb <- zBase
  zPerturb[transRows, 1L] <- rep_len(c(1L, 2L), length(transRows))
  if (nThetaZ > 1L) zPerturb[transRows, 2L] <- rep_len(c(2L, 1L), length(transRows))
  stateB <- state; stateB$z <- zPerturb

  lpA <- .bothLogPriors(stateA, f)
  lpB <- .bothLogPriors(stateB, f)

  expect_true(is.finite(lpA$r) && is.finite(lpB$r))
  # The substance: transformational z rows must not move the prior.
  expect_equal(lpB$r,   lpA$r,   tolerance = 1e-10)   # R   LogPrior
  expect_equal(lpB$cpp, lpA$cpp, tolerance = 1e-10)   # C++ cpp_log_prior
  # And R <-> C++ stay in sync on both states (sync discipline).
  expect_equal(lpA$cpp, lpA$r, tolerance = 1e-10)
  expect_equal(lpB$cpp, lpB$r, tolerance = 1e-10)
})


test_that("log-prior DOES depend on neomorphic z rows (positive control)", {
  # Guards against a degenerate pass where z is ignored entirely: the z prior
  # must still be active for neomorphic characters.
  f <- .MakeMixedNeoTransFixture()
  state <- MkPrime:::.InitState(f$tree, f$mkd, f$model)
  neoRows <- which(f$mkd$type == "neomorphic")
  nThetaZ <- ncol(state$z)
  state$theta <- if (nThetaZ == 2L) c(0.30, 0.78) else rep(0.5, nThetaZ)
  state$pi0   <- 0.42

  zNone <- matrix(0L, nrow = nrow(state$z), ncol = nThetaZ)
  stateNone <- state; stateNone$z <- zNone

  zNeo <- zNone
  zNeo[neoRows[1L], ] <- 1L
  if (length(neoRows) > 1L) zNeo[neoRows[2L], ] <- 2L
  stateNeo <- state; stateNeo$z <- zNeo

  lpNone <- .bothLogPriors(stateNone, f)
  lpNeo  <- .bothLogPriors(stateNeo,  f)

  expect_false(isTRUE(all.equal(lpNeo$r, lpNone$r)))         # R: neo z matters
  expect_equal(lpNeo$cpp, lpNeo$r, tolerance = 1e-10)        # and stays in sync
  expect_equal(lpNone$cpp, lpNone$r, tolerance = 1e-10)
})
