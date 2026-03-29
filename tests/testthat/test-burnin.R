# Tests for burnin selection (M-066)
skip_slow_tests()

# Helper: create a small test posterior
.MakeTestPosterior <- function(nRuns = 1L, nIter = 200L, maxWarmup = 50L,
                               seed = 7042) {
  set.seed(seed)
  tree <- ape::read.tree(
    text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2,t5:0.25);"
  )
  mat <- matrix(c(0, 1, 0, 1, 2,
                  0, 0, 1, 1, 0,
                  1, 1, 0, 0, 1,
                  0, 1, 2, 0, 1), 5, 4,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  suppressMessages(RunMkPrime(
    pd, tree,
    mcmc = MkPrimeMCMC(nIter = nIter, maxWarmup = maxWarmup,
                       minWarmup = maxWarmup, thin = 1L, nChains = 1L,
                       nRuns = nRuns, autoTune = FALSE),
    fixTopology = TRUE
  ))
}


test_that("SetBurnin stores burnin count", {
  result <- .MakeTestPosterior(nRuns = 1L, seed = 7042)

  expect_null(result$burnin)

  r2 <- SetBurnin(result, 10)
  expect_equal(r2$burnin, 10L)

  # Fraction
  nSamp <- nrow(result$samples)
  r3 <- SetBurnin(result, 0.2)
  expect_equal(r3$burnin, floor(0.2 * nSamp))
})


test_that("SetBurnin rejects invalid values", {
  result <- .MakeTestPosterior(nRuns = 1L, seed = 2918)

  expect_error(SetBurnin(result, nrow(result$samples)),
               "must be less than")
  expect_error(SetBurnin("not_posterior", 10))
})


test_that(".PostBurninData filters correctly (single run)", {
  result <- .MakeTestPosterior(nRuns = 1L, seed = 5531)

  nOrig <- nrow(result$samples)

  pb0 <- MkPrime:::.PostBurninData(result)
  expect_equal(nrow(pb0$samples), nOrig)
  expect_equal(length(pb0$trees), nOrig)

  r2 <- SetBurnin(result, 20)
  pb20 <- MkPrime:::.PostBurninData(r2)
  expect_equal(nrow(pb20$samples), nOrig - 20)
  expect_equal(length(pb20$trees), nOrig - 20)
})


test_that(".PostBurninData filters correctly (multi-run)", {
  result <- .MakeTestPosterior(nRuns = 2L, seed = 8164)

  nPerRun <- nrow(result$per_run[[1]]$samples)
  nTotal <- nrow(result$samples)
  expect_equal(nTotal, 2 * nPerRun)

  r2 <- SetBurnin(result, 10)
  pb <- MkPrime:::.PostBurninData(r2)
  expect_equal(nrow(pb$samples), 2 * (nPerRun - 10))
  expect_equal(length(pb$trees), 2 * (nPerRun - 10))
  expect_equal(nrow(pb$per_run[[1]]$samples), nPerRun - 10)
  expect_equal(nrow(pb$per_run[[2]]$samples), nPerRun - 10)
})


test_that("summary respects burnin", {
  result <- .MakeTestPosterior(nRuns = 1L, seed = 4503)

  s_full <- summary(result)
  r2 <- SetBurnin(result, 50)
  s_burn <- summary(r2)

  expect_equal(nrow(s_full), nrow(s_burn))
  expect_true(all(c("parameter", "mean", "median") %in% names(s_full)))
  expect_true(all(c("parameter", "mean", "median") %in% names(s_burn)))
})


test_that("AutoBurnin runs without error (multi-run)", {
  skip_if_not_installed("coda")
  result <- .MakeTestPosterior(nRuns = 2L, nIter = 400L, maxWarmup = 100L,
                               seed = 6891)

  r_auto <- suppressMessages(AutoBurnin(result))
  expect_true(!is.null(r_auto$burnin))
  expect_true(r_auto$burnin >= 0L)
  expect_true(!is.null(r_auto$auto_burnin_results))

  pb <- MkPrime:::.PostBurninData(r_auto)
  expect_true(nrow(pb$samples) <= nrow(result$samples))
})


test_that("AutoBurnin runs for single-run", {
  skip_if_not_installed("coda")
  result <- .MakeTestPosterior(nRuns = 1L, nIter = 300L, maxWarmup = 80L,
                               seed = 3047)

  r_auto <- suppressMessages(AutoBurnin(result))
  expect_true(!is.null(r_auto$burnin))
  expect_true(r_auto$burnin >= 0L)
})


test_that("ConvergenceDiagnostics respects burnin", {
  skip_if_not_installed("coda")
  result <- .MakeTestPosterior(nRuns = 2L, nIter = 300L, maxWarmup = 80L,
                               seed = 1756)

  d_full <- ConvergenceDiagnostics(result)
  r2 <- SetBurnin(result, 30)
  d_burn <- ConvergenceDiagnostics(r2)

  expect_true(d_burn$nSamples < d_full$nSamples)
  expect_equal(d_burn$burnin, 30L)
})
