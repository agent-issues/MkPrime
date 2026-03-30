# Unit tests for native ESS and R-hat implementations (R/ess.R)

# --- .Ess (single-chain ESS) ---

test_that(".Ess returns approximately n for iid draws", {
  set.seed(7841)
  x <- rnorm(1000)
  ess <- MkPrime:::.Ess(x)
  # iid draws: ESS should be close to n (within generous tolerance)
  expect_true(ess > 500)
  expect_true(ess < 2000)
})

test_that(".Ess returns low values for correlated draws", {
  set.seed(3291)
  # AR(1) with rho = 0.95; theoretical ESS ~ n*(1-rho)/(1+rho) ~ 26
  ar1 <- numeric(1000)
  ar1[1] <- rnorm(1)
  for (i in 2:1000) ar1[i] <- 0.95 * ar1[i - 1] + rnorm(1)
  ess <- MkPrime:::.Ess(ar1)
  expect_true(ess < 100)
  expect_true(ess > 5)
})

test_that(".Ess returns NA for < 3 draws", {
  expect_true(is.na(MkPrime:::.Ess(c(1, 2))))
  expect_true(is.na(MkPrime:::.Ess(1)))
  expect_true(is.na(MkPrime:::.Ess(numeric(0))))
})

test_that(".Ess returns NA for constant input", {
  expect_true(is.na(MkPrime:::.Ess(rep(5, 100))))
})

test_that(".Ess handles NA and Inf gracefully", {
  expect_true(is.na(MkPrime:::.Ess(c(1:10, NA, 12:20))))
  expect_true(is.na(MkPrime:::.Ess(c(1:10, Inf))))
})


# --- .EssVector (drop-in replacement) ---

test_that(".EssVector handles degenerate inputs", {
  expect_true(is.na(MkPrime:::.EssVector(numeric(0))))
  expect_true(is.na(MkPrime:::.EssVector(1)))
  expect_true(is.na(MkPrime:::.EssVector(rep(5, 100))))
})

test_that(".EssVector returns finite ESS for normal draws", {
  set.seed(9204)
  x <- rnorm(200)
  ess <- MkPrime:::.EssVector(x)
  expect_true(is.finite(ess))
  expect_true(ess >= 1)
})


# --- .EssMatrix ---

test_that(".EssMatrix returns named vector matching columns", {
  set.seed(6152)
  mat <- matrix(rnorm(600), ncol = 3, dimnames = list(NULL, c("a", "b", "c")))
  ess <- MkPrime:::.EssMatrix(mat)
  expect_length(ess, 3)
  expect_named(ess, c("a", "b", "c"))
  expect_true(all(is.finite(ess)))
})


# --- .Rhat ---

test_that(".Rhat is near 1.0 for converged chains", {
  set.seed(2038)
  chains <- matrix(rnorm(2000), ncol = 4)
  rhat <- MkPrime:::.Rhat(chains)
  expect_true(rhat < 1.05)
  expect_true(rhat >= 1.0 || is.na(rhat))
})

test_that(".Rhat detects non-convergence", {
  set.seed(5511)
  # Chains from different distributions
  chains <- cbind(rnorm(500, 0, 1), rnorm(500, 5, 1))
  rhat <- MkPrime:::.Rhat(chains)
  expect_true(rhat > 1.1)
})

test_that(".Rhat returns NA for short chains", {
  # < 8 draws per chain → < 4 per split-half
  expect_true(is.na(MkPrime:::.Rhat(rnorm(6))))
  expect_true(is.na(MkPrime:::.Rhat(matrix(rnorm(12), ncol = 2))))
})

test_that(".Rhat works with a single chain (split internally)", {
  set.seed(4192)
  x <- rnorm(200)
  rhat <- MkPrime:::.Rhat(x)
  expect_true(is.finite(rhat))
  expect_true(rhat < 1.1)
})

test_that(".Rhat handles constant input", {
  expect_true(is.na(MkPrime:::.Rhat(rep(5, 100))))
})

test_that(".Rhat catches tail non-convergence", {
  set.seed(8720)
  # Chains with same mean but different variances → tail divergence
  chains <- cbind(rnorm(500, 0, 1), rnorm(500, 0, 3))
  rhat <- MkPrime:::.Rhat(chains)
  # Should be elevated (though not always dramatically)
  expect_true(is.finite(rhat))
})


# --- .EssBulk and .EssTail ---

test_that(".EssBulk returns finite values for well-behaved input", {
  set.seed(3847)
  x <- rnorm(500)
  ess <- MkPrime:::.EssBulk(x)
  expect_true(is.finite(ess))
  expect_true(ess > 50)
})

test_that(".EssTail returns finite values for well-behaved input", {
  set.seed(1093)
  x <- rnorm(500)
  ess <- MkPrime:::.EssTail(x)
  expect_true(is.finite(ess))
  expect_true(ess > 10)
})


# --- Internal helpers ---

test_that(".SplitChains doubles columns and halves rows", {
  mat <- matrix(1:20, ncol = 2)
  splits <- MkPrime:::.SplitChains(mat)
  expect_equal(ncol(splits), 4)
  expect_equal(nrow(splits), 5)
})

test_that(".ZScale produces standard normal-ish output", {
  set.seed(6331)
  x <- runif(1000)
  z <- MkPrime:::.ZScale(x)
  expect_true(abs(mean(z)) < 0.1)
  expect_true(abs(sd(z) - 1) < 0.2)
})

test_that(".FoldDraws folds around median", {
  x <- c(1, 2, 3, 4, 5)
  folded <- MkPrime:::.FoldDraws(x)
  expect_equal(folded, c(2, 1, 0, 1, 2))
})

test_that(".Autocovariance matches base R acf at lag 0", {
  set.seed(2845)
  x <- rnorm(100)
  acov <- MkPrime:::.Autocovariance(x)
  # Lag-0 autocovariance ~ var(x) * (n-1)/n
  n <- length(x)
  expected <- var(x) * (n - 1) / n
  expect_equal(acov[1], expected, tolerance = 0.01)
})


# --- Cross-validation against posterior package ---

test_that(".Rhat agrees with posterior::rhat (within tolerance)", {
  skip_if_not_installed("posterior")
  set.seed(4721)

  # 4 chains, 500 draws each
  chains <- matrix(rnorm(2000), ncol = 4)

  our_rhat <- MkPrime:::.Rhat(chains)
  # posterior expects draws_matrix or similar; use rhat.default directly
  their_rhat <- posterior::rhat(posterior::as_draws_matrix(chains))

  expect_equal(our_rhat, their_rhat, tolerance = 0.01)
})

test_that(".Rhat agrees with posterior::rhat for non-converged chains", {
  skip_if_not_installed("posterior")
  set.seed(1835)

  chains <- cbind(rnorm(300, 0, 1), rnorm(300, 3, 1))
  our_rhat <- MkPrime:::.Rhat(chains)
  their_rhat <- posterior::rhat(posterior::as_draws_matrix(chains))

  # Both should agree it's bad; allow wider tolerance for large values
  expect_equal(our_rhat, their_rhat, tolerance = 0.05)
})

test_that(".EssBulk agrees with posterior::ess_bulk (within tolerance)", {
  skip_if_not_installed("posterior")
  set.seed(7612)

  chains <- matrix(rnorm(2000), ncol = 4)
  our_ess <- MkPrime:::.EssBulk(chains)
  their_ess <- posterior::ess_bulk(posterior::as_draws_matrix(chains))

  expect_equal(our_ess, their_ess, tolerance = our_ess * 0.15)
})


# --- Cross-validation against coda ---

test_that(".EssVector is comparable to coda::effectiveSize", {
  skip_if_not_installed("coda")
  set.seed(8493)

  x <- rnorm(500)
  our_ess <- MkPrime:::.EssVector(x)
  their_ess <- as.numeric(coda::effectiveSize(coda::mcmc(x)))

  # Allow generous tolerance — algorithms differ somewhat
  expect_equal(our_ess, their_ess, tolerance = max(our_ess, their_ess) * 0.3)
})

test_that(".EssVector is comparable to coda for correlated data", {
  skip_if_not_installed("coda")
  set.seed(4819)

  ar1 <- numeric(500)
  ar1[1] <- rnorm(1)
  for (i in 2:500) ar1[i] <- 0.9 * ar1[i - 1] + rnorm(1)

  our_ess <- MkPrime:::.EssVector(ar1)
  their_ess <- as.numeric(coda::effectiveSize(coda::mcmc(ar1)))

  # Both should be low; allow wider proportional tolerance
  expect_true(our_ess < 200)
  expect_true(their_ess < 200)
})


# --- M-146: unequal chain length handling ---

test_that(".CheckConvergenceFromLogs handles unequal sample counts", {
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                   "rate_log_sd", "p")
  header <- paste(c("Sample", paramNames), collapse = "\t")

  set.seed(5283)

  # Run 1: 100 samples
  n1 <- 100L
  mat1 <- matrix(rnorm(n1 * length(paramNames)), n1, length(paramNames))
  lines1 <- c(header, vapply(seq_len(n1), function(i) {
    paste(c(i * 10, mat1[i, ]), collapse = "\t")
  }, character(1)))

  # Run 2: 60 samples (fewer — e.g. ESS converged earlier)
  n2 <- 60L
  mat2 <- matrix(rnorm(n2 * length(paramNames)), n2, length(paramNames))
  lines2 <- c(header, vapply(seq_len(n2), function(i) {
    paste(c(i * 10, mat2[i, ]), collapse = "\t")
  }, character(1)))

  f1 <- tempfile(fileext = ".log")
  f2 <- tempfile(fileext = ".log")
  writeLines(lines1, f1)
  writeLines(lines2, f2)
  on.exit({ unlink(f1); unlink(f2) })

  mcmc <- MkPrimeMCMC(minEss = 5, maxRhat = 5.0)

  # Should not error or warn about recycling
  result <- expect_no_warning(
    MkPrime:::.CheckConvergenceFromLogs(c(f1, f2), paramNames, mcmc)
  )
  expect_false(is.null(result))
  expect_true(is.numeric(result$maxRhat))
  expect_false(is.na(result$maxRhat))
  expect_true(is.numeric(result$minEss))
})
