# Tests for proposal functions (M-015)

# --- Scale proposal ---

test_that("ProposeScale returns positive values", {
  set.seed(4821)
  for (i in 1:50) {
    result <- MkPrime:::ProposeScale(1.0, tuning = 2.0)
    expect_gt(result$value, 0)
  }
})

test_that("ProposeScale Hastings ratio is log(x'/x)", {
  set.seed(7293)
  x <- 2.5
  result <- MkPrime:::ProposeScale(x, tuning = 1.0)
  expect_equal(result$logHastings, log(result$value / x), tolerance = 1e-12)
})

test_that("ProposeScale: larger tuning gives wider proposals", {
  set.seed(1158)
  n <- 1000
  narrow <- replicate(n, MkPrime:::ProposeScale(1.0, tuning = 0.1)$value)
  wide <- replicate(n, MkPrime:::ProposeScale(1.0, tuning = 5.0)$value)
  expect_lt(sd(log(narrow)), sd(log(wide)))
})

test_that("ProposeScale: tuning = 0 gives identity proposal", {
  set.seed(3082)
  result <- MkPrime:::ProposeScale(3.7, tuning = 0)
  expect_equal(result$value, 3.7)
  expect_equal(result$logHastings, 0)
})


# --- BetaSimplex proposal ---

test_that("ProposeBetaSimplex preserves simplex", {
  set.seed(6417)
  x <- c(0.3, 0.2, 0.15, 0.35)
  for (i in 1:20) {
    result <- MkPrime:::ProposeBetaSimplex(x, tuning = 10)
    expect_equal(sum(result$value), 1.0, tolerance = 1e-12)
    expect_true(all(result$value > 0))
  }
})

test_that("ProposeBetaSimplex modifies exactly two elements", {
  set.seed(2054)
  x <- c(0.25, 0.25, 0.25, 0.25)
  result <- MkPrime:::ProposeBetaSimplex(x, index = 1L, tuning = 10)
  changed <- which(result$value != x)
  expect_equal(length(changed), 2L)
  expect_true(1L %in% changed)
})

test_that("ProposeBetaSimplex Hastings ratio is finite", {
  set.seed(8541)
  x <- c(0.4, 0.3, 0.3)
  for (i in 1:20) {
    result <- MkPrime:::ProposeBetaSimplex(x, tuning = 10)
    expect_true(is.finite(result$logHastings))
  }
})

test_that("ProposeBetaSimplex handles length-2 simplex", {
  set.seed(9876)
  x <- c(0.6, 0.4)
  result <- MkPrime:::ProposeBetaSimplex(x, index = 1L, tuning = 10)
  expect_equal(sum(result$value), 1.0, tolerance = 1e-12)
  expect_true(all(result$value > 0))
})

test_that("ProposeBetaSimplex: higher tuning = more conservative", {
  set.seed(5521)
  n <- 500
  x <- c(0.5, 0.3, 0.2)
  diffs_low <- replicate(n, {
    r <- MkPrime:::ProposeBetaSimplex(x, index = 1L, tuning = 5)
    abs(r$value[1] - x[1])
  })
  diffs_high <- replicate(n, {
    r <- MkPrime:::ProposeBetaSimplex(x, index = 1L, tuning = 100)
    abs(r$value[1] - x[1])
  })
  expect_gt(mean(diffs_low), mean(diffs_high))
})


# --- BoundedIntegerWalk proposal ---

test_that("ProposeBoundedIntWalk stays within bounds", {
  set.seed(3319)
  for (i in 1:100) {
    result <- MkPrime:::ProposeBoundedIntWalk(3L, lower = 2L, window = 2L)
    if (is.finite(result$logHastings)) {
      expect_gte(result$value, 2L)
    }
  }
})

test_that("ProposeBoundedIntWalk rejects below lower bound", {
  set.seed(7842)
  # x = 2, lower = 2, window = 1: delta can be -1, 0, 1
  # If delta = -1: x' = 1 < 2 → reject (logHastings = -Inf)
  rejected <- FALSE
  for (i in 1:50) {
    result <- MkPrime:::ProposeBoundedIntWalk(2L, lower = 2L, window = 1L)
    if (!is.finite(result$logHastings)) {
      rejected <- TRUE
      expect_equal(result$value, 2L)
      break
    }
  }
  expect_true(rejected)
})

test_that("ProposeBoundedIntWalk is symmetric (Hastings = 0)", {
  set.seed(4556)
  result <- MkPrime:::ProposeBoundedIntWalk(5L, lower = 2L, window = 2L)
  if (is.finite(result$logHastings)) {
    expect_equal(result$logHastings, 0)
  }
})

test_that("ProposeBoundedIntWalk can propose delta = 0", {
  set.seed(2001)
  found_zero <- FALSE
  for (i in 1:100) {
    result <- MkPrime:::ProposeBoundedIntWalk(5L, lower = 2L, window = 1L)
    if (result$value == 5L && is.finite(result$logHastings)) {
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
    fwd <- MkPrime:::ProposeScale(x0, tuning = 1.0)
    bwd <- MkPrime:::ProposeScale(fwd$value, tuning = 1.0)
    # The product exp(fwd$logHastings) * exp(bwd$logHastings) should
    # average to x0/x0 = 1 over many trials if reversible
    hr_sum <- hr_sum + exp(fwd$logHastings + bwd$logHastings)
  }
  # Mean should be close to 1
  expect_equal(hr_sum / n, 1.0, tolerance = 0.15)
})
