library("TreeTools")

# Helper: build a small phyDat with an ecology column at position 4.
.MakeEcoData <- function() {
  mat <- matrix(c(
      0, 1, 0, 1, 0,   # char 1
      0, 0, 1, 1, 1,   # char 2
      0, 1, 1, 2, 2,   # char 3
      0, 0, 1, 1, 2    # char 4 (ecology, 3 states)
    ),
    nrow = 5, ncol = 4,
    dimnames = list(paste0("t", 1:5), NULL)
  )
  MatrixToPhyDat(mat)
}


test_that("MkPrimeData ignores ecology when not supplied", {
  pd <- .MakeEcoData()
  mkd <- MkPrimeData(pd)
  expect_null(mkd$ecology)
  expect_null(mkd$kEcology)
  expect_equal(mkd$nChar, 4L)
})


test_that("MkPrimeData extracts and drops an ecology column index", {
  pd <- .MakeEcoData()
  mkd <- MkPrimeData(pd, ecology = 4)
  expect_equal(mkd$nChar, 3L)
  expect_equal(mkd$ecology, c(0L, 0L, 1L, 1L, 2L))
  expect_equal(mkd$kEcology, 3L)
})


test_that("MkPrimeData accepts a taxon-aligned ecology vector", {
  pd <- .MakeEcoData()
  eco <- c(t1 = 0L, t2 = 0L, t3 = 1L, t4 = 1L, t5 = 2L)
  mkd <- MkPrimeData(pd, ecology = eco)
  expect_equal(mkd$nChar, 4L)
  expect_equal(mkd$ecology, c(0L, 0L, 1L, 1L, 2L))
  expect_equal(mkd$kEcology, 3L)
})


test_that("MkPrimeData reorders ecology by taxon names", {
  pd <- .MakeEcoData()
  eco <- c(t5 = 2L, t1 = 0L, t3 = 1L, t2 = 0L, t4 = 1L)
  mkd <- MkPrimeData(pd, ecology = eco)
  expect_equal(mkd$ecology, c(0L, 0L, 1L, 1L, 2L))
})


test_that("MkPrimeData remaps ecology states to contiguous 0-indexed", {
  pd <- .MakeEcoData()
  # Codes 1, 3, 5 should remap to 0, 1, 2
  eco <- c(t1 = 1L, t2 = 1L, t3 = 3L, t4 = 3L, t5 = 5L)
  mkd <- MkPrimeData(pd, ecology = eco)
  expect_equal(mkd$ecology, c(0L, 0L, 1L, 1L, 2L))
  expect_equal(mkd$kEcology, 3L)
})


test_that("MkPrimeData errors on invariant ecology", {
  pd <- .MakeEcoData()
  expect_error(
    MkPrimeData(pd, ecology = rep(0L, 5)),
    regexp = "must vary"
  )
})


test_that("MkPrimeData errors on out-of-range ecology column index", {
  pd <- .MakeEcoData()
  expect_error(
    MkPrimeData(pd, ecology = 99),
    regexp = "must be between 1 and"
  )
})


test_that("MkPrimeData errors on ecology vector length mismatch", {
  pd <- .MakeEcoData()
  expect_error(
    MkPrimeData(pd, ecology = c(0L, 1L)),
    regexp = "length"
  )
})


test_that("MkPrimeData remaps neomorphic indices around the ecology column", {
  pd <- .MakeEcoData()
  mkd <- MkPrimeData(pd, neomorphic = c(1L, 2L), ecology = 4)
  # Columns 1 and 2 keep their indices; column 4 was removed; column 3 remains.
  expect_equal(mkd$nChar, 3L)
  expect_equal(mkd$type[1:2], c("neomorphic", "neomorphic"))
  expect_equal(mkd$type[3], "transformational")
})


test_that("MkPrimeData rejects neomorphic / ecology column overlap", {
  pd <- .MakeEcoData()
  expect_error(
    MkPrimeData(pd, neomorphic = 4L, ecology = 4),
    regexp = "neomorphic"
  )
})


test_that("MkPrimeData honours user-supplied kEcology", {
  pd <- .MakeEcoData()
  eco <- c(t1 = 0L, t2 = 0L, t3 = 1L, t4 = 1L, t5 = 1L)
  mkd <- MkPrimeData(pd, ecology = eco, kEcology = 4L)
  expect_equal(mkd$kEcology, 4L)
})


test_that("MkPrimeData rejects kEcology below observed", {
  pd <- .MakeEcoData()
  eco <- c(t1 = 0L, t2 = 0L, t3 = 1L, t4 = 1L, t5 = 2L)
  expect_error(
    MkPrimeData(pd, ecology = eco, kEcology = 2L),
    regexp = "at least"
  )
})


test_that("MkPrimeModel accepts ecology-aware options with sane defaults", {
  m <- MkPrimeModel(ecologyAware = TRUE)
  expect_true(m$ecologyAware)
  expect_identical(m$magnitudeMode, "global")
  expect_equal(m$rho0Alpha, 75)
  expect_equal(m$rho0Beta, 25)
  expect_equal(m$sigmaPhi, 1)
  expect_equal(m$thetaAlpha, 1)
  expect_equal(m$thetaBeta, 1)
  expect_equal(m$gibbsZEvery, 50L)
})


test_that("MkPrimeModel(ecologyAware = FALSE) leaves defaults intact", {
  m <- MkPrimeModel()
  expect_false(m$ecologyAware)
  expect_identical(m$magnitudeMode, "global")
})


test_that("MkPrimeModel rejects invalid ecology hyperparameters", {
  expect_error(MkPrimeModel(ecologyAware = TRUE, rho0Alpha = 0))
  expect_error(MkPrimeModel(ecologyAware = TRUE, rho0Beta = -1))
  expect_error(MkPrimeModel(ecologyAware = TRUE, sigmaPhi = 0))
  expect_error(MkPrimeModel(ecologyAware = TRUE, gibbsZEvery = 0))
  expect_error(
    MkPrimeModel(ecologyAware = TRUE, magnitudeMode = "bogus"),
    regexp = "should be one of"
  )
})


test_that(".FinalizeModel errors when ecologyAware but no ecology in mkd", {
  pd <- .MakeEcoData()
  mkd <- MkPrimeData(pd)  # no ecology
  # expSteps set so .FinalizeModel doesn't call .FitchScore on a synthetic tree
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10)
  fakeTree <- ape::rtree(5, tip.label = paste0("t", 1:5))
  expect_error(
    MkPrime:::.FinalizeModel(model, fakeTree, mkd),
    regexp = "requires ecology data"
  )
})


test_that("LogPrior boundary checks: phi must be positive", {
  pd <- .MakeEcoData()
  mkd <- MkPrimeData(pd, ecology = 4)
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "logseries")
  state <- list(
    tree_length = 1,
    rate_log_sd = 0.1,
    rel_br_lengths = rep(1 / 7, 7),
    rate_loss = 1,
    kPrime = mkd$kObs,
    phi = -0.1,
    pi0 = 0.7,
    z = matrix(0L, nrow = mkd$nChar, ncol = mkd$kEcology)
  )
  expect_equal(MkPrime:::LogPrior(state, model, mkd), -Inf)
})


test_that("LogPrior boundary checks: pi0 must be in (0, 1)", {
  pd <- .MakeEcoData()
  mkd <- MkPrimeData(pd, ecology = 4)
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "logseries")
  state <- list(
    tree_length = 1,
    rate_log_sd = 0.1,
    rel_br_lengths = rep(1 / 7, 7),
    rate_loss = 1,
    kPrime = mkd$kObs,
    phi = 1,
    pi0 = 1.5,
    z = matrix(0L, nrow = mkd$nChar, ncol = mkd$kEcology)
  )
  expect_equal(MkPrime:::LogPrior(state, model, mkd), -Inf)
})


test_that("LogPrior returns finite for well-formed ecology state", {
  pd <- .MakeEcoData()
  mkd <- MkPrimeData(pd, ecology = 4)
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "logseries")
  # v2: state requires `theta` (length kEcology - 1) and z is nChar x (kEco-1).
  state <- list(
    tree_length = 1,
    rate_log_sd = 0.1,
    rel_br_lengths = rep(1 / 7, 7),
    rate_loss = 1,
    kPrime = mkd$kObs,
    phi = 1,
    pi0 = 0.7,
    theta = rep(0.5, mkd$kEcology - 1L),
    z = matrix(0L, nrow = mkd$nChar, ncol = mkd$kEcology - 1L)
  )
  expect_true(is.finite(MkPrime:::LogPrior(state, model, mkd)))
})
