# Tests for M-097 rotating ticker helpers (.CompactEssStr, .CompactKpStr)

test_that(".CompactEssStr produces abbreviated param:ESS string", {
  ess <- c(log_posterior = 138, tree_length = 5, rate_log_sd = 42, p = 139)
  result <- MkPrime:::.CompactEssStr(ess)
  expect_match(result, "lP:138")
  expect_match(result, "TL:5")
  expect_match(result, "rsd:42")
  expect_match(result, "p:139")
})


test_that(".CompactEssStr skips kPrime and br_ columns", {
  ess <- c(log_posterior = 100, kPrime_1 = 300, kPrime_2 = 50, br_1 = 80)
  result <- MkPrime:::.CompactEssStr(ess)
  expect_match(result, "lP:100")
  expect_false(grepl("kPrime", result))
  expect_false(grepl("br_", result))
})


test_that(".CompactEssStr returns '?' when no scalar params", {
  ess <- c(kPrime_1 = 300, kPrime_2 = 50)
  expect_equal(MkPrime:::.CompactEssStr(ess), "?")
})


test_that(".CompactEssStr handles non-finite ESS with '?'", {
  ess <- c(log_posterior = NA_real_, tree_length = Inf, rate_log_sd = 42)
  result <- MkPrime:::.CompactEssStr(ess)
  expect_match(result, "lP:\\?")
  expect_match(result, "TL:\\?")
  expect_match(result, "rsd:42")
})


test_that(".CompactEssStr includes optional params when present", {
  ess <- c(log_posterior = 100, tree_length = 50, rate_loss = 80,
           rate_log_sd = 42, p = 139, rate_neo = 60, beta_scale = 70)
  result <- MkPrime:::.CompactEssStr(ess)
  expect_match(result, "rL:80")
  expect_match(result, "rN:60")
  expect_match(result, "bS:70")
})


test_that(".CompactKpStr returns empty string when no kPrime", {
  ess <- c(log_posterior = 100, tree_length = 50)
  expect_equal(MkPrime:::.CompactKpStr(ess), "")
})


test_that(".CompactKpStr summarises kPrime as count: min/med/max", {
  ess <- c(log_posterior = 100, kPrime_1 = 300, kPrime_2 = 100, kPrime_3 = 50)
  result <- MkPrime:::.CompactKpStr(ess)
  expect_match(result, "kP\\(3\\)")
  expect_match(result, "50/100/300")
})


test_that(".CompactKpStr handles non-finite kPrime ESS", {
  ess <- c(kPrime_1 = NA_real_, kPrime_2 = NA_real_)
  result <- MkPrime:::.CompactKpStr(ess)
  expect_match(result, "kP\\(2\\): \\?")
})


test_that(".CompactKpStr handles single kPrime", {
  ess <- c(log_posterior = 100, kPrime_1 = 200)
  result <- MkPrime:::.CompactKpStr(ess)
  expect_match(result, "kP\\(1\\)")
  # min/med/max all the same
  expect_match(result, "200/200/200")
})
