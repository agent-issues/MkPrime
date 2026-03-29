# M-118: Bactrian proposal kernel (Yang & Rodríguez 2013)
#
# Tests that the Bactrian perturbation kernel has the expected statistical
# properties: symmetric, bimodal, variance-matched to Uniform(-0.5, 0.5).

# After normalization, draws have sd = 1/sqrt(12) ≈ 0.2887
# and modes at ±0.95/sqrt(12) ≈ ±0.274.
bactrian_mode <- 0.95 / sqrt(12)  # ≈ 0.274
target_var    <- 1 / 12            # ≈ 0.0833

test_that("bactrian_draws() returns correct number of draws", {
  draws <- bactrian_draws(100L)
  expect_length(draws, 100)
})

test_that("Bactrian kernel is approximately symmetric around zero", {
  set.seed(4817)
  draws <- bactrian_draws(50000L)
  expect_lt(abs(mean(draws)), 0.01)
})

test_that("Bactrian kernel variance matches Uniform(-0.5, 0.5)", {
  set.seed(7293)
  draws <- bactrian_draws(50000L)
  expect_gt(var(draws), target_var * 0.95)
  expect_lt(var(draws), target_var * 1.05)
})

test_that("Bactrian kernel is bimodal (density trough near zero)", {
  set.seed(3164)
  draws <- bactrian_draws(100000L)
  dens <- density(draws, bw = 0.03)
  idx_zero <- which.min(abs(dens$x))
  idx_pos  <- which.min(abs(dens$x - bactrian_mode))
  idx_neg  <- which.min(abs(dens$x + bactrian_mode))
  expect_lt(dens$y[idx_zero], dens$y[idx_pos])
  expect_lt(dens$y[idx_zero], dens$y[idx_neg])
})

test_that("Bactrian kernel rarely produces near-zero perturbations", {
  set.seed(5681)
  draws <- bactrian_draws(50000L)
  # Fewer than 5% of draws in (-0.08, 0.08) ≈ within ±0.3 * sd
  prop_near_zero <- mean(abs(draws) < 0.08)
  expect_lt(prop_near_zero, 0.05)
})
