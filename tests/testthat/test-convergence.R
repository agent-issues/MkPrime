# Tests for convergence monitoring (Phase 5: M-033)

test_that("ConvergenceDiagnostics works with single run", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(5194)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, warmup = 200L))

  diag <- ConvergenceDiagnostics(result)

  expect_true(is.list(diag))
  expect_true("ess" %in% names(diag))
  expect_true("minEss" %in% names(diag))
  expect_true(is.numeric(diag$ess))
  expect_true(all(diag$ess > 0, na.rm = TRUE))
  expect_equal(diag$minEss, min(diag$ess))

  # No PSRF for single run
  expect_null(diag$psrf)
  expect_true(is.na(diag$maxPsrf))
})


test_that("ConvergenceDiagnostics works with multiple runs", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(2204)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 1000L, thin = 5L, warmup = 500L))

  diag <- ConvergenceDiagnostics(result)

  expect_true(is.numeric(diag$ess))
  expect_true(all(diag$ess > 0, na.rm = TRUE))

  # PSRF should be present
  expect_true(!is.null(diag$psrf))
  expect_true(is.numeric(diag$psrf))
  expect_true(all(diag$psrf > 0, na.rm = TRUE))
  expect_true(is.numeric(diag$maxPsrf))
  expect_false(is.na(diag$maxPsrf))
})


test_that("ESS is reasonable for short chains", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(4487)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 2000L, thin = 5L,
                        warmup = 1000L))

  diag <- ConvergenceDiagnostics(result)

  # 200 post-warmup samples; ESS should be between 1 and 200
  expect_true(all(diag$ess <= 200, na.rm = TRUE))
  expect_true(all(diag$ess > 0, na.rm = TRUE))
})


test_that("PSRF is near 1 for converged chains", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(9108)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 2000L, thin = 5L,
                        warmup = 1000L))

  diag <- ConvergenceDiagnostics(result)

  # For a simple problem with moderate iterations, PSRF should be < 2
  # (convergence gets better with more iterations; this is a basic check)
  expect_true(all(diag$psrf < 2.0, na.rm = TRUE))
})


test_that(".KeyParamCols excludes branch lengths", {
  nms <- c("log_posterior", "log_likelihood", "tree_length",
           "rate_loss", "rate_log_sd", "p", "kPrime_1",
           "br_1", "br_2", "br_3")
  mat <- matrix(0, 2, length(nms), dimnames = list(NULL, nms))

  key <- MkPrime:::.KeyParamCols(mat)
  expect_true(all(nms[key] %in% c("log_posterior", "log_likelihood",
                                    "tree_length", "rate_loss",
                                    "rate_log_sd", "p", "kPrime_1")))
  expect_false(any(grepl("br_", nms[key])))
})


test_that("ConvergenceDiagnostics rejects non-MkPosterior", {
  expect_error(ConvergenceDiagnostics(list()), "MkPosterior")
})
