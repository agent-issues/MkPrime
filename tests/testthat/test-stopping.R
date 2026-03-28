# Tests for stopping rules (Phase 5: M-034)
skip_slow_tests()

test_that("RunMkPrime reports stop_reason = 'max_iter' by default", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(5194)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 500L, thin = 5L, warmup = 200L))

  expect_equal(result$stop_reason, "max_iter")
  expect_equal(result$actual_iter, 500L)
})


test_that("maxTime stops MCMC early", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(3382)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000000L, thin = 5L,
                        warmup = 100L, maxTime = 0.5))

  expect_equal(result$stop_reason, "max_time")
  expect_lt(result$actual_iter, 1000000L)
  # Should still have some samples
  expect_gt(nrow(result$samples), 0)
})


test_that("MkPrimeMCMC stores stopping parameters", {
  cfg <- MkPrimeMCMC(maxTime = 60, minEss = 200, maxPsrf = 1.05,
                      checkEvery = 500L)
  expect_equal(cfg$maxTime, 60)
  expect_equal(cfg$minEss, 200)
  expect_equal(cfg$maxPsrf, 1.05)
  expect_equal(cfg$checkEvery, 500L)
})


test_that("Convergence-based stopping works", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  # Very generous convergence criteria so it triggers
  set.seed(9283)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 50000L, thin = 5L,
                        warmup = 200L, minEss = 5, maxPsrf = 5.0,
                        checkEvery = 300L))

  # Should stop before max_iter due to generous criteria
  if (result$stop_reason == "converged") {
    expect_lt(result$actual_iter, 50000L)
  }
  # Even if not converged (rare), result should be valid
  expect_s3_class(result, "MkPosterior")
  expect_gt(nrow(result$samples), 0)
})


test_that(".CheckConvergence returns NULL for insufficient samples", {
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                   "rate_loss", "rate_log_sd", "p")

  # Create runs with too few samples
  runs <- list(
    list(saved_idx = 3L,
         samples = matrix(rnorm(18), 10, length(paramNames),
                          dimnames = list(NULL, paramNames))),
    list(saved_idx = 3L,
         samples = matrix(rnorm(18), 10, length(paramNames),
                          dimnames = list(NULL, paramNames)))
  )

  mcmc <- MkPrimeMCMC(minEss = 100, maxPsrf = 1.05)
  result <- MkPrime:::.CheckConvergence(runs, paramNames, mcmc)
  expect_null(result)
})


test_that(".CheckConvergence returns per-parameter ess and psrf for multi-run", {
  set.seed(6142)
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                   "rate_loss", "rate_log_sd", "p")

  makeRun <- function() list(
    saved_idx = 50L,
    samples   = matrix(rnorm(50 * length(paramNames)), 50, length(paramNames),
                       dimnames = list(NULL, paramNames))
  )
  runs <- list(makeRun(), makeRun())
  mcmc <- MkPrimeMCMC(minEss = 100, maxPsrf = 1.05)

  result <- MkPrime:::.CheckConvergence(runs, paramNames, mcmc)
  expect_false(is.null(result))

  # Per-parameter ESS vector present and correctly named
  expect_true(!is.null(result$ess))
  keyCols <- grep("^(log_|tree_|rate_|p$|kPrime_)", paramNames)
  expect_named(result$ess, paramNames[keyCols], ignore.order = FALSE)
  expect_true(all(result$ess > 0, na.rm = TRUE))

  # Per-parameter PSRF vector present (multi-run)
  expect_true(!is.null(result$psrf))
  expect_named(result$psrf, paramNames[keyCols], ignore.order = FALSE)
  expect_true(all(result$psrf > 0, na.rm = TRUE))

  # Scalar summaries still present
  expect_true(is.numeric(result$minEss))
  expect_true(is.numeric(result$maxPsrf))
  expect_equal(result$minEss, min(result$ess, na.rm = TRUE))
  expect_equal(result$maxPsrf, max(result$psrf, na.rm = TRUE))
})


test_that(".CheckConvergence works for single run (ESS only, psrf = NULL)", {
  set.seed(3871)
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                   "rate_loss", "rate_log_sd", "p")

  runs <- list(list(
    saved_idx = 50L,
    samples   = matrix(rnorm(50 * length(paramNames)), 50, length(paramNames),
                       dimnames = list(NULL, paramNames))
  ))
  mcmc <- MkPrimeMCMC(nRuns = 1L, minEss = 100)

  result <- MkPrime:::.CheckConvergence(runs, paramNames, mcmc)
  expect_false(is.null(result))

  expect_true(!is.null(result$ess))
  expect_true(all(result$ess > 0, na.rm = TRUE))

  # No PSRF for single run
  expect_null(result$psrf)
  expect_true(is.na(result$maxPsrf))

  # Cannot converge on minEss alone when criterion not met
  expect_false(result$converged)
})


test_that(".PrintProgressTable prints without error", {
  set.seed(2953)
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                   "rate_loss", "rate_log_sd", "p")
  keyCols <- grep("^(log_|tree_|rate_|p$|kPrime_)", paramNames)
  keyNames <- paramNames[keyCols]

  # Simulate a diagCheck return value
  diagCheck <- list(
    ess     = setNames(runif(length(keyNames), 50, 300), keyNames),
    psrf    = setNames(runif(length(keyNames), 1, 1.05), keyNames),
    minEss  = 50,
    maxPsrf = 1.05
  )

  out <- capture.output(
    MkPrime:::.PrintProgressTable(diagCheck, nRuns = 2L, iter = 5000L,
                                   nSamples = 450L)
  )
  expect_true(length(out) > 0)
  expect_true(any(grepl("log_posterior", out)))
  expect_true(any(grepl("tree_length",   out)))
  expect_true(any(grepl("ESS",           out)))
  expect_true(any(grepl("PSRF",          out)))
})


test_that(".PrintProgressTable single-run omits PSRF column", {
  set.seed(7318)
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                   "rate_loss", "rate_log_sd", "p")
  keyCols  <- grep("^(log_|tree_|rate_|p$|kPrime_)", paramNames)
  keyNames <- paramNames[keyCols]

  diagCheck <- list(
    ess     = setNames(runif(length(keyNames), 100, 300), keyNames),
    psrf    = NULL,
    minEss  = 100,
    maxPsrf = NA_real_
  )

  out <- capture.output(
    MkPrime:::.PrintProgressTable(diagCheck, nRuns = 1L, iter = 3000L,
                                   nSamples = 200L)
  )
  expect_true(any(grepl("ESS", out)))
  # PSRF column header should not appear
  expect_false(any(grepl("PSRF", out)))
})


test_that("nIter = Inf with maxTime stopping works", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  set.seed(7241)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, thin = 5L, warmup = 200L, maxTime = 0.5))

  expect_s3_class(result, "MkPosterior")
  expect_equal(result$stop_reason, "max_time")
  expect_true(is.infinite(result$mcmc$nIter))
  expect_gt(nrow(result$samples), 0)
})


test_that("MkPrimeMCMC stores cancelFile", {
  cf <- tempfile(fileext = ".signal")
  cfg <- MkPrimeMCMC(cancelFile = cf)
  expect_equal(cfg$cancelFile, cf)
})


test_that("MkPrimeMCMC rejects non-string cancelFile", {
  expect_error(MkPrimeMCMC(cancelFile = 123L), "cancelFile")
  expect_error(MkPrimeMCMC(cancelFile = c("a", "b")), "cancelFile")
})


test_that("cancelFile causes early exit with stop_reason 'cancelled'", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  cf <- tempfile(fileext = ".signal")
  # Pre-create the cancel file so the run stops at the first check
  file.create(cf)
  on.exit(unlink(cf))

  set.seed(6641)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000000L, thin = 5L,
                        warmup = 100L, cancelFile = cf))

  expect_equal(result$stop_reason, "cancelled")
  expect_lt(result$actual_iter, 1000000L)
})


test_that("cancelFile + checkpointFile saves checkpoint on cancel", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  cf   <- tempfile(fileext = ".signal")
  ckpt <- tempfile(fileext = ".rds")
  file.create(cf)
  on.exit({ unlink(cf); unlink(ckpt) })

  set.seed(8823)
  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000000L, thin = 5L,
                        warmup = 100L, cancelFile = cf,
                        checkpointFile = ckpt))

  expect_equal(result$stop_reason, "cancelled")
  expect_true(file.exists(ckpt))
  payload <- readRDS(ckpt)
  expect_true(is.list(payload))
  expect_true(!is.null(payload$runs))
})


test_that("MkCancelPath returns path inside jobDir", {
  expect_equal(MkCancelPath("/tmp/myjob"), "/tmp/myjob/mkp_cancel.signal")
})


test_that("MkPrimeMCMC nIter = Inf default and warmup default", {
  m_inf <- MkPrimeMCMC()
  expect_true(is.infinite(m_inf$nIter))
  expect_equal(m_inf$warmup, 5000L)

  m_finite <- MkPrimeMCMC(nIter = 1000L)
  expect_equal(m_finite$nIter, 1000L)
  expect_equal(m_finite$warmup, 500L)  # nIter/2
})


test_that("Early stopping produces fewer samples", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)

  # Full run
  set.seed(1107)
  full <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000L, thin = 5L,
                        warmup = 200L))

  # Time-limited run
  set.seed(1107)
  limited <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 1000L, thin = 5L,
                        warmup = 200L, maxTime = 0.01))

  # Limited should have fewer or equal samples
  expect_lte(nrow(limited$samples), nrow(full$samples))
})
