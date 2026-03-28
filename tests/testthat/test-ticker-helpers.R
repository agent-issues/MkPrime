# Tests for M-097 rotating ticker helpers
# (.BuildTickerPages, .TickerSummaryStr)

# --- .BuildTickerPages ---

test_that(".BuildTickerPages returns summary-only when no scalar ESS", {
  diag <- list(
    ess = c(kPrime_1 = 300, kPrime_2 = 50),
    minEss = 50, maxPsrf = NA_real_
  )
  pages <- MkPrime:::.BuildTickerPages(diag)
  expect_length(pages, 1L)
  expect_match(cli::ansi_strip(pages), "minESS")
})


test_that(".BuildTickerPages interleaves summary and detail", {
  diag <- list(
    ess = c(log_posterior = 150, tree_length = 80, rate_log_sd = 250, p = 200),
    minEss = 80, maxPsrf = NA_real_
  )
  pages <- MkPrime:::.BuildTickerPages(diag)
  plain <- cli::ansi_strip(pages)
  # 4 params → 2 detail pages → [S, D1, S, D2] = 4 pages
  expect_length(pages, 4L)
  # Odd-numbered pages are summaries
  expect_match(plain[1], "minESS")
  expect_match(plain[3], "minESS")
  # Even-numbered pages are detail
  expect_match(plain[2], "^ESS ")
  expect_match(plain[4], "^ESS ")
})


test_that(".BuildTickerPages uses full parameter names", {
  diag <- list(
    ess = c(log_posterior = 150, tree_length = 80),
    minEss = 80, maxPsrf = NA_real_
  )
  pages <- MkPrime:::.BuildTickerPages(diag)
  plain <- cli::ansi_strip(pages)
  detail <- plain[grepl("^ESS", plain)]
  expect_match(detail, "log_posterior")
  expect_match(detail, "tree_length")
})


test_that(".BuildTickerPages puts <= 2 params per detail page", {
  diag <- list(
    ess = c(log_posterior = 100, tree_length = 200, rate_log_sd = 300,
            p = 400, rate_neo = 500),
    minEss = 100, maxPsrf = NA_real_
  )
  pages <- MkPrime:::.BuildTickerPages(diag)
  plain <- cli::ansi_strip(pages)
  details <- plain[grepl("^ESS", plain)]
  # 5 params → 3 detail pages
  expect_length(details, 3L)
  # Each detail page has at most 2 param names with ":"
  for (d in details) {
    colon_count <- length(gregexpr(":", d)[[1]])
    expect_lte(colon_count, 2L)
  }
})


test_that(".BuildTickerPages skips non-finite ESS params", {
  diag <- list(
    ess = c(log_posterior = 100, rate_loss = NA_real_, tree_length = 200),
    minEss = 100, maxPsrf = NA_real_
  )
  pages <- MkPrime:::.BuildTickerPages(diag)
  combined <- paste(cli::ansi_strip(pages), collapse = " ")
  expect_false(grepl("rate_loss", combined))
  expect_match(combined, "log_posterior")
  expect_match(combined, "tree_length")
})


test_that(".BuildTickerPages colours ESS values", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(
    ess = c(log_posterior = 50, tree_length = 150),
    minEss = 50, maxPsrf = NA_real_
  )
  pages <- MkPrime:::.BuildTickerPages(diag)
  detail <- pages[grepl("ESS", cli::ansi_strip(pages)) &
                    !grepl("minESS", cli::ansi_strip(pages))]
  # ESS 50 (red) and 150 (yellow) should produce ANSI codes

  expect_true(any(cli::ansi_has_any(detail)))
})


test_that(".BuildTickerPages handles single param", {
  diag <- list(
    ess = c(log_posterior = 200),
    minEss = 200, maxPsrf = NA_real_
  )
  pages <- MkPrime:::.BuildTickerPages(diag)
  # 1 param → 1 detail page → [S, D] = 2 pages
  expect_length(pages, 2L)
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


test_that(".TickerSummaryStr leaves minESS plain when >= 200", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(minEss = 300, maxPsrf = NA_real_)
  result <- MkPrime:::.TickerSummaryStr(diag)
  expect_false(cli::ansi_has_any(result))
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


test_that(".TickerSummaryStr leaves PSRF plain when <= 1.05", {
  withr::local_options(cli.num_colors = 256L)
  diag <- list(minEss = 300, maxPsrf = 1.02)
  result <- MkPrime:::.TickerSummaryStr(diag)
  expect_false(cli::ansi_has_any(result))
})


test_that(".TickerSummaryStr handles non-finite minESS", {
  diag <- list(minEss = Inf, maxPsrf = NA_real_)
  plain <- cli::ansi_strip(MkPrime:::.TickerSummaryStr(diag))
  expect_match(plain, "minESS: \\?")
})
