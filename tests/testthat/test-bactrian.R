# M-118: Bactrian proposal kernel (Yang & Rodríguez 2013)
#
# Tests that the Bactrian perturbation kernel has the expected statistical
# properties: symmetric, bimodal with modes at ±m, correct variance.

test_that("bactrian_draws() returns correct number of draws", {
  draws <- bactrian_draws(100L)
  expect_length(draws, 100)
})

test_that("Bactrian kernel is approximately symmetric around zero", {
  set.seed(4817)
  draws <- bactrian_draws(50000L)
  expect_lt(abs(mean(draws)), 0.02)
})

test_that("Bactrian kernel has unit variance", {
  # Var = m² + (1 - m²) = 1
  set.seed(7293)
  draws <- bactrian_draws(50000L)
  expect_gt(var(draws), 0.95)
  expect_lt(var(draws), 1.05)
})

test_that("Bactrian kernel is bimodal (density trough near zero)", {
  set.seed(3164)
  draws <- bactrian_draws(100000L)
  # Density near zero should be lower than density near ±0.95
  dens <- density(draws, bw = 0.1)
  idx_zero <- which.min(abs(dens$x))
  idx_pos  <- which.min(abs(dens$x - 0.95))
  idx_neg  <- which.min(abs(dens$x + 0.95))
  expect_lt(dens$y[idx_zero], dens$y[idx_pos])
  expect_lt(dens$y[idx_zero], dens$y[idx_neg])
})

test_that("Bactrian kernel rarely produces near-zero perturbations", {
  set.seed(5681)
  draws <- bactrian_draws(50000L)
  # Fewer than 5% of draws should be in (-0.3, 0.3)
  prop_near_zero <- mean(abs(draws) < 0.3)
  expect_lt(prop_near_zero, 0.05)
})
