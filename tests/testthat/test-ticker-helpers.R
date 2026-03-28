# Tests for M-097 rotating ticker helpers
# (.CompactEssStr, .TickerSummaryStr)

# --- .CompactEssStr ---

test_that(".CompactEssStr produces abbreviated param:ESS string", {
  ess <- c(log_posterior = 138, tree_length = 5, rate_log_sd = 42, p = 139)
  result <- MkPrime:::.CompactEssStr(ess)
  plain <- cli::ansi_strip(result)
  expect_match(plain, "lP:138")
  expect_match(plain, "TL:5")
  expect_match(plain, "rsd:42")
  expect_match(plain, "p:139")
})


test_that(".CompactEssStr skips kPrime and br_ columns", {
  ess <- c(log_posterior = 100, kPrime_1 = 300, kPrime_2 = 50, br_1 = 80)
  plain <- cli::ansi_strip(MkPrime:::.CompactEssStr(ess))
  expect_match(plain, "lP:100")
  expect_false(grepl("kPrime", plain))
  expect_false(grepl("br_", plain))
})


test_that(".CompactEssStr returns '?' when no scalar params", {
  ess <- c(kPrime_1 = 300, kPrime_2 = 50)
  expect_equal(MkPrime:::.CompactEssStr(ess), "?")
})


test_that(".CompactEssStr handles non-finite ESS with '?'", {
  ess <- c(log_posterior = NA_real_, tree_length = Inf, rate_log_sd = 42)
  plain <- cli::ansi_strip(MkPrime:::.CompactEssStr(ess))
  expect_match(plain, "lP:\\?")
  expect_match(plain, "TL:\\?")
  expect_match(plain, "rsd:42")
})


test_that(".CompactEssStr includes optional params when present", {
  ess <- c(log_posterior = 100, tree_length = 50, rate_loss = 80,
           rate_log_sd = 42, p = 139, rate_neo = 60, beta_scale = 70)
  plain <- cli::ansi_strip(MkPrime:::.CompactEssStr(ess))
  expect_match(plain, "rL:80")
  expect_match(plain, "rN:60")
  expect_match(plain, "bS:70")
})


test_that(".CompactEssStr colours ESS values (red < 100, yellow 100-199)", {
  withr::local_options(cli.num_colors = 256L)
  ess <- c(log_posterior = 50, tree_length = 150, rate_log_sd = 250)
  result <- MkPrime:::.CompactEssStr(ess)

  # Red for ESS < 100 — ANSI codes present around "50"
  expect_true(cli::ansi_has_any(result))
  # Value 250 should appear uncoloured
  plain <- cli::ansi_strip(result)
  expect_match(plain, "rsd:250")
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
  # No ANSI codes when everything is fine
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
