# CONV-001: a convergence window in which every monitored scalar is constant
# must read as "cannot assess", not as "converged".
#
# R-hat and ESS are NA exactly when a chain is constant. Reducing an all-NA
# vector with `max(na.rm = TRUE)` yields -Inf and `min(na.rm = TRUE)` yields
# Inf, and both satisfy the stopping rules (`maxRhat <= threshold`,
# `minEss >= threshold`). So the worst mixing outcome there is -- a sampler
# that has not moved at all -- used to be reported as the best one, precisely
# where the user is relying on automatic stopping rather than reading traces.

test_that(".MinOrNA / .MaxOrNA report NA rather than a passing sentinel", {
  expect_identical(.MinOrNA(c(NA_real_, NA_real_)), NA_real_)
  expect_identical(.MaxOrNA(c(NA_real_, NA_real_)), NA_real_)
  expect_identical(.MinOrNA(numeric(0)), NA_real_)
  expect_identical(.MaxOrNA(numeric(0)), NA_real_)

  # Partially observed vectors still reduce normally, ignoring the NAs.
  expect_identical(.MinOrNA(c(NA, 5, 9)), 5)
  expect_identical(.MaxOrNA(c(NA, 5, 9)), 9)
})

test_that("an all-constant window does not report convergence", {
  paramNames <- c("log_posterior", "log_likelihood",
                  "tree_length", "rate_log_sd")
  n <- 60L

  # A frozen chain: every column constant across every retained draw.
  stuck <- matrix(rep(c(-100, -120, 2.5, 0.8), each = n), nrow = n,
                  dimnames = list(NULL, paramNames))
  runs <- list(list(samples = stuck, saved_idx = n),
               list(samples = stuck, saved_idx = n))

  res <- .CheckConvergence(runs, paramNames,
                           mcmc = list(minEss = 200, maxRhat = 1.01))

  expect_false(res$converged)
  expect_true(is.na(res$minEss))
  expect_true(is.na(res$maxRhat))
})

test_that("a mixing window is still assessed normally", {
  paramNames <- c("log_posterior", "log_likelihood",
                  "tree_length", "rate_log_sd")
  n <- 200L
  set.seed(1L)

  draws <- function() {
    matrix(c(rnorm(n, -100), rnorm(n, -120), rnorm(n, 2.5, 0.1),
             rnorm(n, 0.8, 0.05)), nrow = n,
           dimnames = list(NULL, paramNames))
  }
  runs <- list(list(samples = draws(), saved_idx = n),
               list(samples = draws(), saved_idx = n))

  res <- .CheckConvergence(runs, paramNames,
                           mcmc = list(minEss = 10, maxRhat = 1.5))

  # The point is that a real window still yields finite diagnostics -- the
  # guard must not turn every window into NA.
  expect_true(is.finite(res$minEss))
  expect_true(is.finite(res$maxRhat))
  expect_true(res$converged)
})


# --- #128: all-NA diagnostics must display as NA ---
#
# `.MinOrNA`/`.MaxOrNA` guard the computation. The print path reduces the
# k' rows separately and needs the same guards: a bare `min`/`max` over an
# all-NA vector warns, and yields `-Inf` for the display to render.

test_that("all-NA k' diagnostics print as NA, not Inf (#128)", {
  diag <- structure(
    list(ess = c(log_posterior = NA_real_, tree_length = NA_real_,
                 kPrime_1 = NA_real_, kPrime_2 = NA_real_),
         minEss = NA_real_,
         rhat = c(log_posterior = NA_real_,
                  kPrime_1 = NA_real_, kPrime_2 = NA_real_),
         maxRhat = NA_real_, treeEss = NULL,
         nRuns = 2L, nSamples = 0L, burnin = 0L),
    class = "MkpDiagnostics")

  expect_no_warning(out <- capture.output(print(diag)))
  expect_false(any(grepl("Inf", out, fixed = TRUE)))
  expect_true(any(grepl("kPrime (2)", out, fixed = TRUE)))
})


test_that("print on a zero-run MkPosterior is quiet and says so (#128)", {
  paramNames <- c("log_posterior", "log_likelihood", "tree_length")
  post <- MkPosterior(
    samples = matrix(numeric(0), 0L, length(paramNames),
                     dimnames = list(NULL, paramNames)),
    trees = list(), acceptance = numeric(0),
    model = NULL, data = NULL, mcmc = list(thin = 1L),
    warmup = 0L, tuning = list()
  )
  post$nRuns <- 0L
  post$nSamples <- 0L
  post$stop_reason <- "maxTime"

  # cli writes its alerts to the message stream, not to stdout.
  expect_no_warning(
    msgs <- capture.output(print(post), type = "message")
  )
  expect_true(any(grepl("No runs completed", cli::ansi_strip(msgs))))
})


# --- #197: an NA gate column blocks the verdict instead of leaving the gate ---

.GateRun <- function(rateLogSd, n = 200L) {
  samples <- cbind(log_posterior = rnorm(n, -100),
                   tree_length = rnorm(n, 2.5, 0.1),
                   rate_log_sd = rateLogSd, p = 0.5)
  list(samples = samples, saved_idx = n)
}
.gateNames <- c("log_posterior", "tree_length", "rate_log_sd", "p")

test_that("runs stuck at different values do not pass maxRhat (#197)", {
  set.seed(1971)
  res <- .CheckConvergence(list(.GateRun(0.3), .GateRun(0.7)), .gateNames,
                           mcmc = list(maxRhat = 1.01, fixedCols = "p"))
  expect_false(res$converged)
})

test_that("a frozen or non-finite gate column is unassessable (#197)", {
  set.seed(1972)
  frozen <- .CheckConvergence(list(.GateRun(0.3)), .gateNames,
                              mcmc = list(minEss = 10, fixedCols = "p"))
  expect_false(frozen$converged)
  expect_true(is.na(frozen$minEss))

  nonFinite <- .GateRun(rnorm(200, 0.8, 0.05))
  nonFinite$samples[5, "rate_log_sd"] <- Inf
  expect_false(.CheckConvergence(list(nonFinite), .gateNames,
                                 mcmc = list(minEss = 10,
                                             fixedCols = "p"))$converged)
})

test_that("a column no move updates stays exempt (#197)", {
  set.seed(1973)
  res <- .CheckConvergence(list(.GateRun(rnorm(200, 0.8, 0.05))), .gateNames,
                           mcmc = list(minEss = 10, fixedCols = "p"))
  expect_true(res$converged)
  expect_true(is.finite(res$minEss))
})

test_that("the tuning bandit scores a frozen gate column as unassessable (#197)", {
  set.seed(1974)
  mat <- cbind(log_posterior = rnorm(50), tree_length = 2.5, p = 0.5)
  expect_true(is.na(.MinEssRate(mat, 1)[["rate"]]))
  expect_true(is.finite(
    .MinEssRate(mat, 1, fixedCols = c("tree_length", "p"))[["rate"]]
  ))
})
