# Tests for .PrintProgressTable (M-097)

test_that(".PrintProgressTable returns line count", {
  diagCheck <- list(
    ess     = c(log_posterior = 150, tree_length = 80, rate_log_sd = 250, p = 200),
    psrf    = NULL,
    minEss  = 80,
    maxPsrf = NA_real_
  )

  nLines <- capture.output(
    result <- MkPrime:::.PrintProgressTable(
      diagCheck, nRuns = 1L, iter = 5000L, nSamples = 100L
    )
  )
  expect_type(result, "integer")
  expect_equal(result, length(nLines))
})


test_that(".PrintProgressTable returns correct count with PSRF", {
  diagCheck <- list(
    ess     = c(log_posterior = 150, tree_length = 80, rate_log_sd = 250),
    psrf    = c(log_posterior = 1.01, tree_length = 1.08, rate_log_sd = 1.00),
    minEss  = 80,
    maxPsrf = 1.08
  )

  nLines <- capture.output(
    result <- MkPrime:::.PrintProgressTable(
      diagCheck, nRuns = 2L, iter = 5000L, nSamples = 200L
    )
  )
  expect_equal(result, length(nLines))
  # PSRF legend adds one line
  expect_true(any(grepl("PSRF", nLines)))
})


test_that(".PrintProgressTable includes kPrime summary row", {
  ess <- c(log_posterior = 150, tree_length = 80, p = 200,
           kPrime_1 = 300, kPrime_2 = 100, kPrime_3 = 50)
  diagCheck <- list(ess = ess, psrf = NULL, minEss = 50, maxPsrf = NA_real_)

  out <- capture.output(
    result <- MkPrime:::.PrintProgressTable(
      diagCheck, nRuns = 1L, iter = 3000L, nSamples = 50L
    )
  )
  expect_true(any(grepl("kPrime \\(3\\)", out)))
  expect_equal(result, length(out))
})


test_that(".PrintProgressTable prevLines=0 does not emit ANSI codes", {
  diagCheck <- list(
    ess     = c(log_posterior = 150, tree_length = 80),
    psrf    = NULL,
    minEss  = 80,
    maxPsrf = NA_real_
  )

  # No ANSI cursor-up on first call (prevLines = 0)
  raw <- capture.output(
    MkPrime:::.PrintProgressTable(
      diagCheck, nRuns = 1L, iter = 1000L, nSamples = 10L, prevLines = 0L
    )
  )
  # \x1b[ is the ANSI escape prefix — should NOT appear when prevLines = 0
  expect_false(any(grepl("\x1b\\[", raw)))
})


test_that(".PrintProgressTable line count is stable across calls", {
  diagCheck <- list(
    ess     = c(log_posterior = 150, tree_length = 80, rate_log_sd = 250,
                p = 200, kPrime_1 = 300, kPrime_2 = 100),
    psrf    = NULL,
    minEss  = 80,
    maxPsrf = NA_real_
  )

  out1 <- capture.output(
    n1 <- MkPrime:::.PrintProgressTable(
      diagCheck, nRuns = 1L, iter = 1000L, nSamples = 10L
    )
  )

  # Change ESS values but keep same parameters — line count should be identical
  diagCheck$ess[] <- c(300, 200, 100, 50, 400, 250)
  out2 <- capture.output(
    n2 <- MkPrime:::.PrintProgressTable(
      diagCheck, nRuns = 1L, iter = 2000L, nSamples = 20L
    )
  )

  expect_equal(n1, n2)
})
