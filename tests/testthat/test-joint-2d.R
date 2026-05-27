# M-120: 2D joint Bactrian proposals for correlated parameter pairs
#
# Tests for the bivariate Bactrian kernel and joint move integration.

# --- 2D Bactrian kernel properties ---

test_that("bactrian_2d_draws returns correct dimensions", {
  draws <- bactrian_2d_draws(100L, 0.5)
  expect_equal(dim(draws), c(100L, 2L))
})

test_that("2D Bactrian marginals match 1D Bactrian distribution", {
  set.seed(6142)
  draws_2d <- bactrian_2d_draws(50000L, 0.3)
  draws_1d <- bactrian_draws(50000L)

  target_var <- 1 / 12

  # Each marginal should have variance ≈ 1/12
  expect_gt(var(draws_2d[, 1]), target_var * 0.90)
  expect_lt(var(draws_2d[, 1]), target_var * 1.10)
  expect_gt(var(draws_2d[, 2]), target_var * 0.90)
  expect_lt(var(draws_2d[, 2]), target_var * 1.10)

  # Each marginal should be approximately symmetric
  expect_lt(abs(mean(draws_2d[, 1])), 0.01)
  expect_lt(abs(mean(draws_2d[, 2])), 0.01)
})

test_that("2D Bactrian correlation tracks rho parameter", {
  set.seed(8753)
  for (rho in c(0.0, 0.3, 0.7, -0.5)) {
    draws <- bactrian_2d_draws(50000L, rho)
    observed_cor <- cor(draws[, 1], draws[, 2])
    # Tolerance is wider because bimodal structure introduces additional
    # correlation beyond the Gaussian rho
    expect_lt(abs(observed_cor - rho), 0.10,
              label = sprintf("rho=%.1f, observed=%.3f", rho, observed_cor))
  }
})

test_that("2D Bactrian is symmetric: P(z1,z2) = P(-z1,-z2)", {
  set.seed(2914)
  draws <- bactrian_2d_draws(100000L, 0.5)
  # Mean of z1*z2 should be positive (correlated modes), but
  # mean of z1 and z2 should both be ≈ 0
  expect_lt(abs(mean(draws[, 1])), 0.01)
  expect_lt(abs(mean(draws[, 2])), 0.01)
})

test_that("2D Bactrian with rho=0 produces approximately uncorrelated draws", {
  set.seed(4501)
  draws <- bactrian_2d_draws(50000L, 0.0)
  observed_cor <- cor(draws[, 1], draws[, 2])
  # With rho=0, mode coupling p_same=0.5 (independent mode signs)
  # and noise rho=0 → nearly uncorrelated output
  expect_lt(abs(observed_cor), 0.05)
})

test_that("2D Bactrian is bimodal in both dimensions", {
  set.seed(5178)
  draws <- bactrian_2d_draws(100000L, 0.5)
  bactrian_mode <- 0.95 / sqrt(12)

  for (col in 1:2) {
    dens <- density(draws[, col], bw = 0.03)
    idx_zero <- which.min(abs(dens$x))
    idx_pos  <- which.min(abs(dens$x - bactrian_mode))
    expect_lt(dens$y[idx_zero], dens$y[idx_pos])
  }
})


# --- Rho estimation ---

test_that(".EstimateJointRhos returns zeros with too few samples", {
  rhos <- MkPrime:::.EstimateJointRhos(NULL, hasNeo = TRUE)
  expect_equal(rhos$rho_tl_rls, 0.0)
  expect_equal(rhos$rho_tl_rl, 0.0)

  small <- data.frame(tree_length = 1:20, rate_log_sd = 1:20)
  rhos2 <- MkPrime:::.EstimateJointRhos(as.matrix(small), hasNeo = TRUE)
  expect_equal(rhos2$rho_tl_rls, 0.0)
})

test_that("u.123: .EstimateJointRhos returns 0 for constant input (NaN cor)", {
  # If all values are identical, cor() returns NaN. The rho should stay at 0.
  n <- 100
  constant <- cbind(tree_length = rep(1.5, n), rate_log_sd = rep(0.3, n),
                    rate_loss = rep(0.5, n))
  rhos <- MkPrime:::.EstimateJointRhos(constant, hasNeo = TRUE)
  expect_equal(rhos$rho_tl_rls, 0.0)
  expect_equal(rhos$rho_tl_rl, 0.0)

  # One column constant, other varies → cor is NaN
  mixed <- cbind(tree_length = rep(2.0, n), rate_log_sd = rlnorm(n))
  rhos2 <- MkPrime:::.EstimateJointRhos(mixed, hasNeo = FALSE)
  expect_equal(rhos2$rho_tl_rls, 0.0)
})

test_that(".EstimateJointRhos recovers known correlation", {
  set.seed(9134)
  n <- 200
  x <- rlnorm(n, 0, 0.3)
  # y correlated with x in log space
  y <- exp(0.7 * log(x) + rnorm(n, 0, 0.2))
  samples <- cbind(tree_length = x, rate_log_sd = y)
  rhos <- MkPrime:::.EstimateJointRhos(samples, hasNeo = FALSE)

  # Should recover positive correlation
  expect_gt(rhos$rho_tl_rls, 0.3)
  # Capped at 0.95
  expect_lte(rhos$rho_tl_rls, 0.95)
  # No neo → rate_loss and rate_neo rhos stay 0
  expect_equal(rhos$rho_tl_rl, 0.0)
  expect_equal(rhos$rho_tl_rn, 0.0)
})

test_that(".EstimateJointRhos recovers (tree_length, rate_neo) correlation", {
  set.seed(2271)
  n <- 200
  tl <- rlnorm(n, 0, 0.3)
  rn <- exp(0.6 * log(tl) + rnorm(n, 0, 0.2))
  samples <- cbind(tree_length = tl, rate_neo = rn)
  rhos <- MkPrime:::.EstimateJointRhos(samples, hasNeo = TRUE)

  expect_gt(rhos$rho_tl_rn, 0.3)
  expect_lte(rhos$rho_tl_rn, 0.95)

  # Without hasNeo, rho_tl_rn stays 0 even if rate_neo column exists
  rhos_no_neo <- MkPrime:::.EstimateJointRhos(samples, hasNeo = FALSE)
  expect_equal(rhos_no_neo$rho_tl_rn, 0.0)
})

test_that(".BuildJointRhoMatrix fills correctly", {
  chainRhos <- list(
    list(rho_tl_rls = 0.5, rho_tl_rl = 0.3, rho_tl_rn = -0.2),
    list(rho_tl_rls = 0.5, rho_tl_rl = 0.3, rho_tl_rn = -0.2)
  )
  moves <- list(
    list(name = "tree_length", type = "scale"),
    list(name = "joint_tl_rls", type = "joint_2d"),
    list(name = "joint_tl_rl", type = "joint_2d"),
    list(name = "joint_tl_rn", type = "joint_2d"),
    list(name = "nni", type = "nni")
  )
  mat <- MkPrime:::.BuildJointRhoMatrix(chainRhos, moves, 2L)
  expect_equal(dim(mat), c(2, 5))
  expect_equal(mat[1, 2], 0.5)   # joint_tl_rls
  expect_equal(mat[1, 3], 0.3)   # joint_tl_rl
  expect_equal(mat[1, 4], -0.2)  # joint_tl_rn
  expect_equal(mat[1, 1], 0.0)   # non-joint
  expect_equal(mat[1, 5], 0.0)   # non-joint
})


# --- Integration test ---

test_that("Joint 2D moves produce valid posterior", {
  skip_on_cran()

  nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
  pd <- TreeTools::ReadAsPhyDat(nexFile)
  mkd <- MkPrimeData(pd)

  cfg <- MkPrimeMCMC(
    nIter = 3000L, nRuns = 1L, nChains = 1L,
    minWarmup = 500L, maxWarmup = 1000L,
    autoTune = FALSE, thin = 10L,
    joint2d = TRUE
  )

  set.seed(7382)
  res <- RunMkPrime(mkd, mcmc = cfg)

  expect_s3_class(res, "MkPosterior")
  expect_gt(nrow(res$samples), 0)
  expect_true("joint_tl_rls" %in% names(res$acceptance))
  expect_gt(res$acceptance[["joint_tl_rls"]], 0)
})

test_that("joint2d = FALSE excludes joint moves", {
  skip_on_cran()

  nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
  pd <- TreeTools::ReadAsPhyDat(nexFile)
  mkd <- MkPrimeData(pd)

  cfg <- MkPrimeMCMC(
    nIter = 2000L, nRuns = 1L, nChains = 1L,
    minWarmup = 500L, maxWarmup = 1000L,
    autoTune = FALSE, thin = 10L,
    joint2d = FALSE
  )

  set.seed(3896)
  res <- RunMkPrime(mkd, mcmc = cfg)

  expect_s3_class(res, "MkPosterior")
  expect_false("joint_tl_rls" %in% names(res$acceptance))
  expect_false("joint_tl_rl" %in% names(res$acceptance))
  expect_false("joint_tl_rn" %in% names(res$acceptance))
})

test_that("joint_tl_rn runs and produces acceptance on mixed data", {
  skip_on_cran()

  # Force a binary character to be neomorphic so hasNeo gates joint_tl_rn in.
  nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
  pd <- TreeTools::ReadAsPhyDat(nexFile)
  mkdBase <- MkPrimeData(pd)
  binaryChars <- which(mkdBase$kObs == 2L)
  skip_if(length(binaryChars) == 0L, "no binary characters to mark neomorphic")
  mkd <- MkPrimeData(pd, neomorphic = binaryChars[1])

  cfg <- MkPrimeMCMC(
    nIter = 3000L, nRuns = 1L, nChains = 1L,
    minWarmup = 500L, maxWarmup = 1000L,
    autoTune = FALSE, thin = 10L,
    joint2d = TRUE
  )

  set.seed(4719)
  res <- RunMkPrime(mkd, mcmc = cfg)

  expect_s3_class(res, "MkPosterior")
  expect_true("joint_tl_rn" %in% names(res$acceptance))
  expect_gt(res$acceptance[["joint_tl_rn"]], 0)
})
