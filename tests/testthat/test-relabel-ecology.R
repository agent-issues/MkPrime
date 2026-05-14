library("TreeTools")

# ---------------------------------------------------------------------------
# Helper: gamma_e and rate-factor formula mirroring mcmc_ecology.cpp:249-267.
# Used to verify that relabelling leaves per-(sample, char, ecology) rate
# factors unchanged.
# ---------------------------------------------------------------------------
.GammaE <- function(pi0, thetaE, phiE) {
  pi0 + (1 - pi0) * (thetaE * phiE + (1 - thetaE) / phiE)
}

.RateFactor <- function(z, phiE, gammaE) {
  mu <- if (z == 0L) 1.0 else if (z == 1L) phiE else 1.0 / phiE
  mu / gammaE
}

# Compute a scalar "fingerprint" that is invariant under the reflection
# for a single ecology column e (1-indexed) across all characters.
.ColRateFactors <- function(zCol, phiE, thetaE, pi0) {
  gE <- .GammaE(pi0, thetaE, phiE)
  vapply(zCol, .RateFactor, phiE = phiE, gammaE = gE, FUN.VALUE = numeric(1L))
}

# ---------------------------------------------------------------------------
# Minimal fake MkPosterior builder (no MCMC needed)
# ---------------------------------------------------------------------------
.FakePosterior <- function(samplesDf, zList, magnitudeMode = "global",
                            refEcology = 0L, kEcology = 2L) {
  samples <- as.matrix(samplesDf)
  model <- list(
    ecologyAware  = TRUE,
    magnitudeMode = magnitudeMode
  )
  data <- list(
    refEcology = refEcology,
    kEcology   = kEcology
  )
  result <- list(
    samples   = samples,
    z_samples = zList,
    model     = model,
    data      = data
  )
  class(result) <- "MkPosterior"
  result
}


# ---------------------------------------------------------------------------
# 1. Likelihood invariance (rate-factor check), global mode
# ---------------------------------------------------------------------------
test_that("global relabelling leaves per-cell rate factors unchanged", {
  pi0 <- 0.3

  # Four samples:
  #   s1: phi=0.25, theta=0.05, z has mix of 0/1/2
  #   s2: phi=0.50, theta=0.10, z has mix
  #   s3: phi=2.00, theta=0.90  (already canonical — should be no-op)
  #   s4: phi=4.00, theta=0.80
  phiVec   <- c(0.25, 0.50, 2.00, 4.00)
  thetaVec <- c(0.05, 0.10, 0.90, 0.80)

  # Fixed z matrices (nChar=6, kEco-1=1 column)
  zMats <- list(
    matrix(c(0L, 1L, 2L, 0L, 1L, 2L), nrow = 6L, ncol = 1L),
    matrix(c(2L, 2L, 1L, 0L, 0L, 1L), nrow = 6L, ncol = 1L),
    matrix(c(1L, 1L, 0L, 2L, 0L, 2L), nrow = 6L, ncol = 1L),
    matrix(c(0L, 2L, 2L, 1L, 1L, 0L), nrow = 6L, ncol = 1L)
  )

  # Compute rate factors before relabelling
  rfBefore <- lapply(seq_len(4L), function(i) {
    .ColRateFactors(as.integer(zMats[[i]][, 1L]), phiVec[i], thetaVec[i], pi0)
  })

  df <- data.frame(
    log_posterior  = rep(-10, 4L),
    log_likelihood = rep(-20, 4L),
    phi            = phiVec,
    pi0            = rep(pi0, 4L),
    theta_1        = thetaVec
  )
  res <- .FakePosterior(df, zMats, magnitudeMode = "global")
  res2 <- RelabelEcology(res, magnitudeMode = "global")

  # Compute rate factors after relabelling
  phiNew   <- res2$samples[, "phi"]
  thetaNew <- res2$samples[, "theta_1"]
  rfAfter <- lapply(seq_len(4L), function(i) {
    .ColRateFactors(
      as.integer(res2$z_samples[[i]][, 1L]),
      phiNew[i], thetaNew[i], pi0
    )
  })

  for (i in seq_len(4L)) {
    expect_equal(rfAfter[[i]], rfBefore[[i]], tolerance = 1e-12,
                 label = paste0("sample ", i, " rate factors"))
  }
})


# ---------------------------------------------------------------------------
# 2. Convention enforced: every phi >= 1 after relabelling (global mode)
# ---------------------------------------------------------------------------
test_that("global relabelling enforces phi >= 1 for all samples", {
  phiVec   <- c(0.1, 0.5, 1.0, 2.0, 5.0, 0.25)
  thetaVec <- c(0.9, 0.8, 0.5, 0.3, 0.7, 0.6)
  nS <- length(phiVec)

  zMats <- lapply(seq_len(nS), function(i) {
    matrix(sample(0L:2L, 10L, replace = TRUE), nrow = 10L, ncol = 1L)
  })

  df <- data.frame(
    phi     = phiVec,
    pi0     = rep(0.4, nS),
    theta_1 = thetaVec
  )
  res <- .FakePosterior(df, zMats)
  res2 <- RelabelEcology(res)

  expect_true(all(res2$samples[, "phi"] >= 1.0))
})


# ---------------------------------------------------------------------------
# 3. Theta flipped correctly: theta_new = 1 - theta_old when phi was < 1
# ---------------------------------------------------------------------------
test_that("theta is flipped exactly when phi < 1", {
  phiVec   <- c(0.2, 3.0, 0.8, 4.0)
  thetaVec <- c(0.1, 0.9, 0.3, 0.7)
  flip     <- phiVec < 1.0

  nS <- length(phiVec)
  zMats <- lapply(seq_len(nS), function(i) {
    matrix(c(0L, 1L, 2L), nrow = 3L, ncol = 1L)
  })

  df <- data.frame(phi = phiVec, pi0 = rep(0.3, nS), theta_1 = thetaVec)
  res <- .FakePosterior(df, zMats)
  res2 <- RelabelEcology(res)

  thetaNew <- res2$samples[, "theta_1"]
  expect_equal(thetaNew[!flip], thetaVec[!flip])
  expect_equal(thetaNew[flip],  1.0 - thetaVec[flip])
})


# ---------------------------------------------------------------------------
# 4. z codes swapped correctly: 0->0, 1->2, 2->1 when phi < 1
# ---------------------------------------------------------------------------
test_that("z codes 1 <-> 2 are swapped when phi < 1, 0 unchanged", {
  # Single sample with phi < 1; fixed z matrix
  zm <- matrix(c(0L, 1L, 2L, 0L, 2L, 1L), nrow = 6L, ncol = 1L)
  df <- data.frame(phi = 0.3, pi0 = 0.4, theta_1 = 0.05)
  res <- .FakePosterior(df, list(zm))
  res2 <- RelabelEcology(res)

  expected <- matrix(c(0L, 2L, 1L, 0L, 1L, 2L), nrow = 6L, ncol = 1L)
  expect_equal(res2$z_samples[[1L]], expected)

  # A sample with phi > 1 should have z unchanged
  df2 <- data.frame(phi = 3.0, pi0 = 0.4, theta_1 = 0.9)
  res3 <- .FakePosterior(df2, list(zm))
  res4 <- RelabelEcology(res3)
  expect_equal(res4$z_samples[[1L]], zm)
})


# ---------------------------------------------------------------------------
# 5. Idempotence: relabelling twice equals relabelling once
# ---------------------------------------------------------------------------
test_that("RelabelEcology is idempotent", {
  phiVec   <- c(0.2, 0.5, 3.0, 4.0)
  thetaVec <- c(0.1, 0.2, 0.8, 0.7)
  nS <- length(phiVec)

  zMats <- lapply(seq_len(nS), function(i) {
    matrix(c(0L, 1L, 2L, 1L, 0L, 2L), nrow = 6L, ncol = 1L)
  })

  df <- data.frame(phi = phiVec, pi0 = rep(0.3, nS), theta_1 = thetaVec)
  res <- .FakePosterior(df, zMats)

  res1  <- RelabelEcology(res)
  res2  <- RelabelEcology(res1)

  expect_equal(res2$samples, res1$samples)
  for (i in seq_len(nS)) {
    expect_equal(res2$z_samples[[i]], res1$z_samples[[i]])
  }
})


# ---------------------------------------------------------------------------
# 6. pi0 column is NOT modified
# ---------------------------------------------------------------------------
test_that("pi0 is left unchanged by relabelling", {
  pi0Vec <- c(0.3, 0.5, 0.7, 0.2)
  df <- data.frame(
    phi     = c(0.2, 0.5, 2.0, 4.0),
    pi0     = pi0Vec,
    theta_1 = c(0.1, 0.2, 0.8, 0.7)
  )
  nS <- nrow(df)
  zMats <- lapply(seq_len(nS), function(i) {
    matrix(c(0L, 1L, 2L), nrow = 3L, ncol = 1L)
  })
  res <- .FakePosterior(df, zMats)
  res2 <- RelabelEcology(res)
  expect_equal(res2$samples[, "pi0"], pi0Vec)
})


# ---------------------------------------------------------------------------
# 7. Non-ecology result returns unchanged with a message
# ---------------------------------------------------------------------------
test_that("non-ecology result returns unchanged with a message", {
  model <- list(ecologyAware = FALSE)
  data  <- list()
  res <- structure(
    list(model = model, data = data, samples = matrix(0, 1, 1),
         z_samples = NULL),
    class = "MkPosterior"
  )
  expect_message(out <- RelabelEcology(res), regexp = "No ecology model")
  expect_identical(out, res)
})


# ---------------------------------------------------------------------------
# 8. relabelled attribute is set to TRUE
# ---------------------------------------------------------------------------
test_that("relabelled attribute is set TRUE after relabelling", {
  df <- data.frame(phi = 0.5, pi0 = 0.3, theta_1 = 0.1)
  res <- .FakePosterior(df, list(matrix(c(0L, 1L, 2L), nrow = 3L, ncol = 1L)))
  res2 <- RelabelEcology(res)
  expect_true(isTRUE(attr(res2, "relabelled")))
})


# ---------------------------------------------------------------------------
# 9. Likelihood invariance: per_ecology mode (2 non-reference ecologies)
# ---------------------------------------------------------------------------
test_that("per-ecology relabelling leaves rate factors unchanged", {
  pi0 <- 0.25
  # kEco = 3, refEcology = 0 (0-indexed), so non-ref ecologies are 1 and 2 (0-idx)
  # phi_1 = phi[0] (ref, always 1, not flipped); phi_2 = phi[1]; phi_3 = phi[2]
  # theta_1, theta_2 correspond to non-ref ecologies 1 and 2 (0-indexed)
  # In R columns: phi_2 -> ecology 1, theta_1 -> ecology 1, zMat[,1]
  #               phi_3 -> ecology 2, theta_2 -> ecology 2, zMat[,2]

  # One sample with both non-ref phis < 1 (both should flip)
  df <- data.frame(
    phi_1   = 1.0,   # reference ecology phi (inert)
    phi_2   = 0.4,   # non-ref ecology 1 -> will flip
    phi_3   = 0.3,   # non-ref ecology 2 -> will flip
    pi0     = pi0,
    theta_1 = 0.05,
    theta_2 = 0.10
  )
  zm <- matrix(c(0L, 1L, 2L, 2L, 1L, 0L), nrow = 3L, ncol = 2L)

  # Rate factors before for each non-ref ecology column
  rfBefore1 <- .ColRateFactors(zm[, 1L], 0.4, 0.05, pi0)
  rfBefore2 <- .ColRateFactors(zm[, 2L], 0.3, 0.10, pi0)

  res <- .FakePosterior(
    df, list(zm),
    magnitudeMode = "per_ecology",
    refEcology = 0L, kEcology = 3L
  )
  res2 <- RelabelEcology(res, magnitudeMode = "per_ecology")

  phiNew   <- res2$samples[1L, ]
  zm2 <- res2$z_samples[[1L]]

  rfAfter1 <- .ColRateFactors(zm2[, 1L], phiNew["phi_2"], res2$samples[1L, "theta_1"], pi0)
  rfAfter2 <- .ColRateFactors(zm2[, 2L], phiNew["phi_3"], res2$samples[1L, "theta_2"], pi0)

  expect_equal(rfAfter1, rfBefore1, tolerance = 1e-12)
  expect_equal(rfAfter2, rfBefore2, tolerance = 1e-12)

  # Convention enforced
  expect_true(all(as.numeric(res2$samples[, c("phi_1", "phi_2", "phi_3")]) >= 1.0))
})


# ---------------------------------------------------------------------------
# 10. Round-trip on the real pilot result
# ---------------------------------------------------------------------------
test_that("round-trip on pilot result: phi > 1, theta_1 > 0.5 after relabelling", {
  # Pilot files are at the package root; resolve from the testthat wd.
  pkgRoot <- normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
  pilotRds <- file.path(
    pkgRoot,
    "dev/pilots/2026-05-13-step3-sim3-pilot/sim3-pilot-result.rds"
  )
  skip_if_not(file.exists(pilotRds), "Pilot RDS not found")

  x <- readRDS(pilotRds)
  res <- x$res

  # Log file path stored in result is relative to the project root; resolve it.
  logPath <- file.path(pkgRoot, res$logFile)
  skip_if_not(file.exists(logPath), "Pilot log file not found")

  preSamples <- ReadMkLog(logPath)
  expect_true(median(preSamples[, "phi"]) < 1.0)
  expect_true(median(preSamples[, "theta_1"]) < 0.5)

  # Override logFile with the resolved absolute path so RelabelEcology can
  # auto-load samples regardless of the calling working directory.
  res$logFile <- logPath
  res2 <- RelabelEcology(res)

  expect_true(median(res2$samples[, "phi"])     > 1.0,
              label = "median phi > 1 after relabelling")
  expect_true(median(res2$samples[, "theta_1"]) > 0.5,
              label = "median theta_1 > 0.5 after relabelling")

  # z code 1 (encouraged) should be the dominant non-zero code after relabelling
  # (before, code 2 dominated because the chain was in the reflected mode)
  allZ <- do.call(rbind, res2$z_samples)
  zTab <- table(allZ)
  expect_true(
    as.integer(zTab["1"]) > as.integer(zTab["2"]),
    label = "code 1 (encouraged) dominates code 2 after relabelling"
  )

  # Idempotence on the real result
  res3 <- RelabelEcology(res2)
  expect_equal(res3$samples, res2$samples)
})
