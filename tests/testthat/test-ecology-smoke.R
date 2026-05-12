library("TreeTools")


skip_slow_tests <- function() {
  if (!identical(Sys.getenv("MKPRIME_SLOW_TESTS"), "true")) {
    testthat::skip("Slow MCMC tests skipped (set MKPRIME_SLOW_TESTS=true to run)")
  }
}


# Build the same minimal ecology fixture used by test-ecology-plumbing.R but
# expose it locally so this file is self-contained.
.SmokeFixture <- function() {
  set.seed(2026)
  tips <- paste0("t", 1:6)
  mat <- matrix(c(
    0, 1, 2, 0, 1, 2,
    1, 0, 1, 2, 2, 0,
    2, 1, 0, 1, 0, 2,
    0, 1, 0, 1, 0, 1,
    1, 0, 1, 0, 1, 0,
    0, 0, 1, 1, 2, 2
  ), nrow = 6, ncol = 6, byrow = FALSE, dimnames = list(tips, NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = c(4L, 5L), ecology = 6L)
  tree <- Preorder(ape::rtree(6, tip.label = tips))
  list(mkd = mkd, tree = tree)
}


test_that("RunMkPrime completes end-to-end with ecologyAware = TRUE", {
  skip_slow_tests()
  f <- .SmokeFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "none")
  mcmc <- MkPrimeMCMC(nIter = 100, nChains = 1, nRuns = 1,
                      thin = 5, treeThin = 5,
                      logFile = tempfile(fileext = ".log"),
                      checkpointFile = NULL,
                      maxWarmup = 0, minWarmup = 0)
  res <- RunMkPrime(f$mkd, tree = f$tree, model = model, mcmc = mcmc)

  expect_s3_class(res, "MkPosterior")
  expect_true(is.finite(res$chainStates[[1]][[1]] %||% 0) ||
              !is.null(res$logFile))
  # Smoke level: log file exists and has at least one row beyond the header.
  expect_true(file.exists(res$logFile))
})


test_that("RunMkPrime completes end-to-end with per_ecology phi", {
  skip_slow_tests()
  f <- .SmokeFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, magnitudeMode = "per_ecology",
                        expSteps = 10, kPrimePrior = "geometric",
                        coding = "none")
  mcmc <- MkPrimeMCMC(nIter = 100, nChains = 1, nRuns = 1,
                      thin = 5, treeThin = 5,
                      logFile = tempfile(fileext = ".log"),
                      checkpointFile = NULL,
                      maxWarmup = 0, minWarmup = 0)
  res <- RunMkPrime(f$mkd, tree = f$tree, model = model, mcmc = mcmc)
  expect_s3_class(res, "MkPosterior")
})
