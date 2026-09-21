# Unit tests for the tuning bandit's noise gate and payback freeze

test_that(".BeatsIncumbent adopts the first usable candidate", {
  expect_true(MkPrime:::.BeatsIncumbent(1, 100, -Inf, NA_real_))
  expect_true(MkPrime:::.BeatsIncumbent(1, 100, NA_real_, NA_real_))
})

test_that(".BeatsIncumbent rejects an unusable candidate", {
  expect_false(MkPrime:::.BeatsIncumbent(NA_real_, 100, 1, 100))
  expect_false(MkPrime:::.BeatsIncumbent(0, 100, 1, 100))
})

test_that(".BeatsIncumbent rejects a gain inside the estimator's noise", {
  # At ESS 50 apiece the relative standard error is 2 * sqrt(2 / 50) ~ 0.28,
  # so a 5% improvement is noise but a 50% one is not.
  expect_false(MkPrime:::.BeatsIncumbent(1.05, 50, 1, 50))
  expect_true(MkPrime:::.BeatsIncumbent(1.5, 50, 1, 50))
})

test_that(".BeatsIncumbent tightens as the estimates sharpen", {
  expect_false(MkPrime:::.BeatsIncumbent(1.1, 50, 1, 50))
  expect_true(MkPrime:::.BeatsIncumbent(1.1, 5000, 1, 5000))
})

test_that(".TuningPayback keeps tuning without an ESS target", {
  expect_equal(MkPrime:::.TuningPayback(0L, 1e6, 10, NULL),
               list(streak = 0L, freeze = FALSE))
})

test_that(".TuningPayback keeps tuning when the rate is unusable", {
  expect_false(MkPrime:::.TuningPayback(1L, 1e6, -Inf, 200)$freeze)
  expect_false(MkPrime:::.TuningPayback(1L, 1e6, NA_real_, 200)$freeze)
})

test_that(".TuningPayback keeps tuning while it is still cheap", {
  # 200 effective draws at 10/s is 20 s of sampling; 5 s spent is not yet 1:1.
  expect_equal(MkPrime:::.TuningPayback(1L, 5, 10, 200),
               list(streak = 0L, freeze = FALSE))
})

test_that(".TuningPayback needs two consecutive votes to freeze", {
  first <- MkPrime:::.TuningPayback(0L, 50, 10, 200)
  expect_equal(first, list(streak = 1L, freeze = FALSE))
  expect_true(MkPrime:::.TuningPayback(first$streak, 50, 10, 200)$freeze)
})

test_that(".TuningPayback shares the target between runs", {
  # Two runs each owe 100 draws, so 10 s of sampling: 12 s spent votes to
  # freeze where a single run owing 200 draws would not.
  expect_equal(MkPrime:::.TuningPayback(0L, 12, 10, 200, nRuns = 2L)$streak, 1L)
  expect_equal(MkPrime:::.TuningPayback(0L, 12, 10, 200, nRuns = 1L)$streak, 0L)
})
