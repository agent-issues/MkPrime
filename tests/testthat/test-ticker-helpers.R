# Tests for progress ticker helpers
# (.BuildTickerPages, .TickerSummaryStr)

# --- .BuildTickerPages (M-139: simplified to summary only) ---

test_that(".BuildTickerPages returns single summary string", {
  diag <- list(
    ess = c(log_posterior = 150, tree_length = 80, rate_log_sd = 250),
    minEss = 80, maxPsrf = NA_real_
  )
  pages <- MkPrime:::.BuildTickerPages(diag)
  expect_length(pages, 1L)
  expect_match(cli::ansi_strip(pages), "minESS: 80")
})

test_that(".BuildTickerPages includes PSRF in multi-run", {
  diag <- list(
    ess = c(log_posterior = 200), minEss = 200, maxPsrf = 1.03
  )
  pages <- MkPrime:::.BuildTickerPages(diag)
  expect_length(pages, 1L)
  plain <- cli::ansi_strip(pages)
  expect_match(plain, "minESS: 200")
  expect_match(plain, "PSRF: 1.03")
})


# --- .TickerSummaryStr ---

test_that(".TickerSummaryStr shows minESS only for single-run", {
  diag <- list(minEss = 42, maxPsrf = NA_real_)
  plain <- cli::ansi_strip(MkPrime:::.TickerSummaryStr(diag))
  expect_match(plain, "minESS: 42")
  expect_false(grepl("PSRF", plain))
})


test_that(".TickerSummaryStr includes PSRF when available", {
  diag <- list(minEss = 200, maxPsrf = 1.03)
  plain <- cli::ansi_strip(MkPrime:::.TickerSummaryStr(diag))
  expect_match(plain, "minESS: 200")
  expect_match(plain, "PSRF: 1.03")
})


test_that(".TickerSummaryStr colours minESS red when < 100", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(minEss = 42, maxPsrf = NA_real_)
  result <- MkPrime:::.TickerSummaryStr(diag)
  expect_true(cli::ansi_has_any(result))
  expect_match(cli::ansi_strip(result), "minESS: 42")
})


test_that(".TickerSummaryStr colours minESS yellow when 100-199", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(minEss = 150, maxPsrf = NA_real_)
  result <- MkPrime:::.TickerSummaryStr(diag)
  expect_true(cli::ansi_has_any(result))
  expect_match(cli::ansi_strip(result), "minESS: 150")
})


test_that(".TickerSummaryStr colours minESS green when >= 200", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(minEss = 300, maxPsrf = NA_real_)
  result <- MkPrime:::.TickerSummaryStr(diag)
  expect_true(cli::ansi_has_any(result))
  expect_match(cli::ansi_strip(result), "minESS: 300")
})


test_that(".TickerSummaryStr colours PSRF red when > 1.1", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(minEss = 300, maxPsrf = 1.25)
  result <- MkPrime:::.TickerSummaryStr(diag)
  expect_true(cli::ansi_has_any(result))
  expect_match(cli::ansi_strip(result), "PSRF: 1.25")
})


test_that(".TickerSummaryStr colours PSRF yellow when 1.05-1.1", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(minEss = 300, maxPsrf = 1.07)
  result <- MkPrime:::.TickerSummaryStr(diag)
  expect_true(cli::ansi_has_any(result))
  expect_match(cli::ansi_strip(result), "PSRF: 1.07")
})


test_that(".TickerSummaryStr colours PSRF green when <= 1.05", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(minEss = 300, maxPsrf = 1.02)
  result <- MkPrime:::.TickerSummaryStr(diag)
  # Both minESS (green) and PSRF (green) produce ANSI
  expect_true(cli::ansi_has_any(result))
})


test_that(".TickerSummaryStr handles non-finite minESS", {
  diag <- list(minEss = Inf, maxPsrf = NA_real_)
  plain <- cli::ansi_strip(MkPrime:::.TickerSummaryStr(diag))
  expect_match(plain, "minESS: \\?")
})
