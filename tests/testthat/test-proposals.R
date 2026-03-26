# Tests for proposal functions (M-015)

# --- Scale proposal ---

test_that("propose_scale returns positive values", {
  set.seed(4821)
  for (i in 1:50) {
    result <- MkPrime:::propose_scale(1.0, tuning = 2.0)
    expect_gt(result$value, 0)
  }
})

test_that("propose_scale Hastings ratio is log(x'/x)", {
  set.seed(7293)
  x <- 2.5
  result <- MkPrime:::propose_scale(x, tuning = 1.0)
  expect_equal(result$log_hastings, log(result$value / x), tolerance = 1e-12)
})

test_that("propose_scale: larger tuning gives wider proposals", {
  set.seed(1158)
  n <- 1000
  narrow <- replicate(n, MkPrime:::propose_scale(1.0, tuning = 0.1)$value)
  wide <- replicate(n, MkPrime:::propose_scale(1.0, tuning = 5.0)$value)
  expect_lt(sd(log(narrow)), sd(log(wide)))
})

test_that("propose_scale: tuning = 0 gives identity proposal", {
  set.seed(3082)
  result <- MkPrime:::propose_scale(3.7, tuning = 0)
  expect_equal(result$value, 3.7)
  expect_equal(result$log_hastings, 0)
})


# --- BetaSimplex proposal ---

test_that("propose_beta_simplex preserves simplex", {
  set.seed(6417)
  x <- c(0.3, 0.2, 0.15, 0.35)
  for (i in 1:20) {
    result <- MkPrime:::propose_beta_simplex(x, tuning = 10)
    expect_equal(sum(result$value), 1.0, tolerance = 1e-12)
    expect_true(all(result$value > 0))
  }
})

test_that("propose_beta_simplex modifies exactly two elements", {
  set.seed(2054)
  x <- c(0.25, 0.25, 0.25, 0.25)
  result <- MkPrime:::propose_beta_simplex(x, index = 1L, tuning = 10)
  changed <- which(result$value != x)
  expect_equal(length(changed), 2L)
  expect_true(1L %in% changed)
})

test_that("propose_beta_simplex Hastings ratio is finite", {
  set.seed(8541)
  x <- c(0.4, 0.3, 0.3)
  for (i in 1:20) {
    result <- MkPrime:::propose_beta_simplex(x, tuning = 10)
    expect_true(is.finite(result$log_hastings))
  }
})

test_that("propose_beta_simplex handles length-2 simplex", {
  set.seed(9876)
  x <- c(0.6, 0.4)
  result <- MkPrime:::propose_beta_simplex(x, index = 1L, tuning = 10)
  expect_equal(sum(result$value), 1.0, tolerance = 1e-12)
  expect_true(all(result$value > 0))
})

test_that("propose_beta_simplex: higher tuning = more conservative", {
  set.seed(5521)
  n <- 500
  x <- c(0.5, 0.3, 0.2)
  diffs_low <- replicate(n, {
    r <- MkPrime:::propose_beta_simplex(x, index = 1L, tuning = 5)
    abs(r$value[1] - x[1])
  })
  diffs_high <- replicate(n, {
    r <- MkPrime:::propose_beta_simplex(x, index = 1L, tuning = 100)
    abs(r$value[1] - x[1])
  })
  expect_gt(mean(diffs_low), mean(diffs_high))
})


# --- BoundedIntegerWalk proposal ---

test_that("propose_bounded_int_walk stays within bounds", {
  set.seed(3319)
  for (i in 1:100) {
    result <- MkPrime:::propose_bounded_int_walk(3L, lower = 2L, window = 2L)
    if (is.finite(result$log_hastings)) {
      expect_gte(result$value, 2L)
    }
  }
})

test_that("propose_bounded_int_walk rejects below lower bound", {
  set.seed(7842)
  # x = 2, lower = 2, window = 1: delta can be -1, 0, 1
  # If delta = -1: x' = 1 < 2 → reject (log_hastings = -Inf)
  rejected <- FALSE
  for (i in 1:50) {
    result <- MkPrime:::propose_bounded_int_walk(2L, lower = 2L, window = 1L)
    if (!is.finite(result$log_hastings)) {
      rejected <- TRUE
      expect_equal(result$value, 2L)
      break
    }
  }
  expect_true(rejected)
})

test_that("propose_bounded_int_walk is symmetric (Hastings = 0)", {
  set.seed(4556)
  result <- MkPrime:::propose_bounded_int_walk(5L, lower = 2L, window = 2L)
  if (is.finite(result$log_hastings)) {
    expect_equal(result$log_hastings, 0)
  }
})

test_that("propose_bounded_int_walk can propose delta = 0", {
  set.seed(2001)
  found_zero <- FALSE
  for (i in 1:100) {
    result <- MkPrime:::propose_bounded_int_walk(5L, lower = 2L, window = 1L)
    if (result$value == 5L && is.finite(result$log_hastings)) {
      found_zero <- TRUE
      break
    }
  }
  expect_true(found_zero)
})


# --- Reversibility check for Scale ---

test_that("Scale proposal satisfies detailed balance (statistical)", {
  set.seed(1673)
  # If we propose forward then backward, the product of Hastings ratios
  # should average to 1 (expectation over many trials)
  n <- 2000
  x0 <- 2.0
  hr_sum <- 0
  for (i in seq_len(n)) {
    fwd <- MkPrime:::propose_scale(x0, tuning = 1.0)
    bwd <- MkPrime:::propose_scale(fwd$value, tuning = 1.0)
    # The product exp(fwd$log_hastings) * exp(bwd$log_hastings) should
    # average to x0/x0 = 1 over many trials if reversible
    hr_sum <- hr_sum + exp(fwd$log_hastings + bwd$log_hastings)
  }
  # Mean should be close to 1
  expect_equal(hr_sum / n, 1.0, tolerance = 0.15)
})
