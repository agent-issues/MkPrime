# Tests for mk_prime_relabel_log() — the falling-factorial relabelling
# correction:  log P(k', kObs) = log(k'!) - log((k'-kObs)!)

test_that("Relabelling correction equals log(kObs!) when kPrime == kObs", {
  # When k' = kObs, all states are observed.
  # P(k, k) = k! / 0! = k!, so log correction = log(k!).
  for (k in 2:6) {
    val <- MkPrime:::mk_prime_relabel_log(k, k)
    expected <- lgamma(k + 1)
    expect_equal(val, expected, tolerance = 1e-12,
                 label = sprintf("k=%d", k))
  }
})


test_that("Relabelling correction increases as kPrime increases", {
  # More model states → more equivalent label injections → larger correction
  kObs <- 3
  vals <- vapply(3:10, function(kp) MkPrime:::mk_prime_relabel_log(kp, kObs),
                 numeric(1))
  expect_true(all(diff(vals) > 0),
              label = "Correction should increase with increasing k'")
})


test_that("Relabelling correction errors when kPrime < kObs", {
  expect_error(MkPrime:::mk_prime_relabel_log(2, 3), "kPrime.*must be >= kObs")
})


test_that("Batch relabelling matches scalar version", {
  kObs <- 3
  kPrime_vec <- c(3L, 4L, 5L, 6L, 10L)
  batch_result <- MkPrime:::mk_prime_relabel_log_batch(kPrime_vec, kObs)

  scalar_results <- vapply(kPrime_vec, function(kp) {
    MkPrime:::mk_prime_relabel_log(kp, kObs)
  }, numeric(1))

  expect_equal(batch_result, scalar_results, tolerance = 1e-15)
})


test_that("Relabelling correction has correct analytical form", {
  # Falling factorial: P(k', kObs) = k'! / (k'-kObs)!
  #
  # k'=4, kObs=2: P(4,2) = 4!/2! = 12
  # log(12) = lgamma(5) - lgamma(3)
  kp <- 4
  ko <- 2
  expected <- lgamma(kp + 1) - lgamma(kp - ko + 1)
  expect_equal(MkPrime:::mk_prime_relabel_log(kp, ko), expected,
               tolerance = 1e-14)
})


test_that("Relabelling correction matches brute-force enumeration", {
  # Verified via exhaustive enumeration on a star tree (3 tips, t=0.3,
  # pattern (0,0,1), kObs=2).  See relabelling-correction-proof.md.
  #
  # P(k', 2) should equal k'*(k'-1) for each k':
  for (kp in 2:8) {
    val <- MkPrime:::mk_prime_relabel_log(kp, 2L)
    expected <- log(kp * (kp - 1))
    expect_equal(val, expected, tolerance = 1e-12,
                 label = sprintf("k'=%d, kObs=2", kp))
  }

  # P(k', 3) = k'*(k'-1)*(k'-2):
  for (kp in 3:8) {
    val <- MkPrime:::mk_prime_relabel_log(kp, 3L)
    expected <- log(kp * (kp - 1) * (kp - 2))
    expect_equal(val, expected, tolerance = 1e-12,
                 label = sprintf("k'=%d, kObs=3", kp))
  }
})


test_that("MH ratio delta favours higher k' (k'=3 vs k'=2, kObs=2)", {
  # This is the key property the old formula got wrong.
  # delta = log P(3,2) - log P(2,2) = log(6) - log(2) = log(3) ≈ 1.099
  c3 <- MkPrime:::mk_prime_relabel_log(3L, 2L)
  c2 <- MkPrime:::mk_prime_relabel_log(2L, 2L)
  expect_equal(c3 - c2, log(3), tolerance = 1e-12)
  expect_true(c3 > c2,
              label = "Correction at k'=3 must exceed correction at k'=2")
})
