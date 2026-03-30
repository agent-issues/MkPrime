# Tests for progress ticker helpers
# (.BuildTickerPages, .TickerSummaryStr)

# --- .BuildTickerPages (M-139: simplified to summary only) ---

test_that(".BuildTickerPages returns single summary string", {
  diag <- list(
    ess = c(log_posterior = 150, tree_length = 80, rate_log_sd = 250),
    minEss = 80, maxRhat = NA_real_
  )
  pages <- MkPrime:::.BuildTickerPages(diag)
  expect_length(pages, 1L)
  expect_match(cli::ansi_strip(pages), "minESS: 80")
})

test_that(".BuildTickerPages includes Rhat in multi-run", {
  diag <- list(
    ess = c(log_posterior = 200), minEss = 200, maxRhat = 1.03
  )
  pages <- MkPrime:::.BuildTickerPages(diag)
  expect_length(pages, 1L)
  plain <- cli::ansi_strip(pages)
  expect_match(plain, "minESS: 200")
  expect_match(plain, "Rhat: 1.03")
})


# --- .TickerSummaryStr ---

test_that(".TickerSummaryStr shows minESS only for single-run", {
  diag <- list(minEss = 42, maxRhat = NA_real_)
  plain <- cli::ansi_strip(MkPrime:::.TickerSummaryStr(diag))
  expect_match(plain, "minESS: 42")
  expect_false(grepl("Rhat", plain))
})


test_that(".TickerSummaryStr includes Rhat when available", {
  diag <- list(minEss = 200, maxRhat = 1.03)
  plain <- cli::ansi_strip(MkPrime:::.TickerSummaryStr(diag))
  expect_match(plain, "minESS: 200")
  expect_match(plain, "Rhat: 1.03")
})


test_that(".TickerSummaryStr colours minESS red when < 100", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(minEss = 42, maxRhat = NA_real_)
  result <- MkPrime:::.TickerSummaryStr(diag)
  expect_true(cli::ansi_has_any(result))
  expect_match(cli::ansi_strip(result), "minESS: 42")
})


test_that(".TickerSummaryStr colours minESS yellow when 100-199", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(minEss = 150, maxRhat = NA_real_)
  result <- MkPrime:::.TickerSummaryStr(diag)
  expect_true(cli::ansi_has_any(result))
  expect_match(cli::ansi_strip(result), "minESS: 150")
})


test_that(".TickerSummaryStr colours minESS green when >= 200", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(minEss = 300, maxRhat = NA_real_)
  result <- MkPrime:::.TickerSummaryStr(diag)
  expect_true(cli::ansi_has_any(result))
  expect_match(cli::ansi_strip(result), "minESS: 300")
})


test_that(".TickerSummaryStr colours Rhat red when > 1.05", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(minEss = 300, maxRhat = 1.25)
  result <- MkPrime:::.TickerSummaryStr(diag)
  expect_true(cli::ansi_has_any(result))
  expect_match(cli::ansi_strip(result), "Rhat: 1.25")
})


test_that(".TickerSummaryStr colours Rhat yellow when 1.01-1.05", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(minEss = 300, maxRhat = 1.03)
  result <- MkPrime:::.TickerSummaryStr(diag)
  expect_true(cli::ansi_has_any(result))
  expect_match(cli::ansi_strip(result), "Rhat: 1.03")
})


test_that(".TickerSummaryStr colours Rhat green when <= 1.01", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(minEss = 300, maxRhat = 1.005)
  result <- MkPrime:::.TickerSummaryStr(diag)
  # Both minESS (green) and Rhat (green) produce ANSI
  expect_true(cli::ansi_has_any(result))
})


test_that(".TickerSummaryStr handles non-finite minESS", {
  diag <- list(minEss = Inf, maxRhat = NA_real_)
  plain <- cli::ansi_strip(MkPrime:::.TickerSummaryStr(diag))
  expect_match(plain, "minESS: \\?")
})
