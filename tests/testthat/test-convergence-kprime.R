# Tests for M-098: kPrime excluded from convergence criteria

test_that("kPrime ESS does not affect minEss in convergence filter", {
  # Simulate a diagCheck with low kPrime ESS and high scalar ESS
  fakeEss <- c(log_posterior = 250, log_likelihood = 250, tree_length = 200,
               rate_log_sd = 300, p = 180,
               kPrime_1 = 5, kPrime_2 = 3, kPrime_3 = 8)

  # The convergence filter used in .CheckConvergence (M-098)
  isConv <- !grepl("^kPrime_", names(fakeEss)) & names(fakeEss) != "log_likelihood"
  minEss <- min(fakeEss[isConv], na.rm = TRUE)

  # Without M-098, minEss would be 3 (kPrime_2). With M-098, it's 180 (p).
  expect_equal(minEss, 180)
})


test_that("ConvergenceDiagnostics minEss ignores kPrime", {

  # Build a mock MkPosterior with fake samples (avoids running MCMC)
  set.seed(4821)
  n <- 200L
  samples <- cbind(
    log_posterior  = rnorm(n, -120, 2),
    log_likelihood = rnorm(n, -115, 2),
    tree_length    = rlnorm(n, 1, 0.3),
    rate_log_sd    = rlnorm(n, 0, 0.2),
    p              = rbeta(n, 5, 5),
    kPrime_1       = sample(2:4, n, replace = TRUE),
    kPrime_2       = rep(2L, n)
  )

  posterior <- structure(
    list(
      samples = samples,
      trees = NULL,
      acceptance = list(accepted = 50L, proposed = 200L),
      model = NULL, data = NULL, mcmc = NULL,
      warmup = 0L, tuning = NULL,
      nRuns = 1L, burnin = 0L
    ),
    class = "MkPosterior"
  )

  diag <- ConvergenceDiagnostics(posterior)

  # kPrime should still appear in the ess vector
  kpNames <- grep("^kPrime_", names(diag$ess), value = TRUE)
  expect_length(kpNames, 2L)

  # minEss should equal the min of non-kPrime, non-log_likelihood params only
  scalarEss <- diag$ess[!grepl("^kPrime_", names(diag$ess)) &
                         names(diag$ess) != "log_likelihood"]
  expect_equal(diag$minEss, min(scalarEss, na.rm = TRUE))

  # kPrime_2 has zero variance → ESS = NA. This must NOT make minEss = NA.
  expect_true(is.finite(diag$minEss))
})


test_that("ConvergenceDiagnostics maxRhat ignores kPrime", {

  # Two-run mock with kPrime that has high R-hat
  set.seed(6193)
  n <- 100L
  make_run <- function(offset) {
    cbind(
      log_posterior  = rnorm(n, -120 + offset, 2),
      log_likelihood = rnorm(n, -115 + offset, 2),
      tree_length    = rlnorm(n, 1 + offset * 0.01, 0.3),
      rate_log_sd    = rlnorm(n, 0, 0.2),
      p              = rbeta(n, 5, 5),
      # kPrime deliberately different across runs → high R-hat
      kPrime_1       = sample(2:4, n, replace = TRUE, prob = c(0.1, 0.1, 0.8) +
                                offset * c(0.3, -0.05, -0.25))
    )
  }

  run1 <- make_run(0)
  run2 <- make_run(1)

  posterior <- structure(
    list(
      samples = rbind(run1, run2),
      trees = NULL,
      acceptance = list(accepted = 50L, proposed = 200L),
      model = NULL, data = NULL, mcmc = NULL,
      warmup = 0L, tuning = NULL,
      nRuns = 2L, burnin = 0L,
      per_run = list(
        list(samples = run1, trees = NULL, acceptance = list()),
        list(samples = run2, trees = NULL, acceptance = list())
      )
    ),
    class = "MkPosterior"
  )

  diag <- ConvergenceDiagnostics(posterior)

  # maxRhat should reflect only scalar parameters, not kPrime
  if (!is.null(diag$rhat)) {
    scalarRhat <- diag$rhat[!grepl("^kPrime_", names(diag$rhat)) &
                             names(diag$rhat) != "log_likelihood"]
    expect_equal(diag$maxRhat, max(scalarRhat, na.rm = TRUE))
  }
})


test_that(".PrintProgressTable still displays kPrime row", {
  ess <- c(log_posterior = 250, tree_length = 180, rate_log_sd = 300,
           p = 200, kPrime_1 = 10, kPrime_2 = 15, kPrime_3 = 8)
  diagCheck <- list(ess = ess, rhat = NULL, minEss = 180, maxRhat = NA_real_)

  out <- capture.output(
    MkPrime:::.PrintProgressTable(
      diagCheck, nRuns = 1L, iter = 5000L, nSamples = 100L
    )
  )

  # kPrime row should still appear for informational display
  expect_true(any(grepl("kPrime \\(3\\)", out)))
  # But individual kPrime params should NOT appear as separate rows
  expect_false(any(grepl("kPrime_1", out)))
})
