library("TreeTools")


# Build a small ecology fixture with kEco = 3 (two non-reference ecologies,
# so theta should have length 2 and produce theta_1 and theta_2 columns).
.LoggingFixture <- function() {
  set.seed(42L)
  tips <- paste0("t", 1:6)
  mat <- matrix(c(
    0, 1, 0, 1, 0, 1,   # char 1 (transformational)
    0, 0, 1, 1, 0, 1,   # char 2 (transformational)
    0, 1, 1, 0, 1, 0,   # char 3 (transformational)
    0, 0, 1, 1, 2, 2    # ecology (3 states → kEco = 3)
  ), nrow = 6, ncol = 4, byrow = FALSE,
     dimnames = list(tips, NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, ecology = 4L)
  tree <- Preorder(ape::rtree(6, tip.label = tips))
  list(mkd = mkd, tree = tree)
}


test_that("theta_1 and theta_2 columns appear in MCMC log when kEco = 3", {
  f <- .LoggingFixture()
  # kEco should be 3 → nTheta = 2
  expect_equal(f$mkd$kEcology, 3L)

  logPath <- tempfile(fileext = ".log")
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "none")
  mcmc  <- MkPrimeMCMC(nIter = 20L, nChains = 1L, nRuns = 1L,
                       thin = 10L, treeThin = 10L,
                       logFile = logPath,
                       checkpointFile = NULL,
                       maxWarmup = 0L, minWarmup = 0L)
  res <- RunMkPrime(f$mkd, tree = f$tree, model = model, mcmc = mcmc)

  # Read the log and confirm theta columns are present
  samp <- ReadMkLog(logPath)
  expect_true("theta_1" %in% colnames(samp),
              label = "theta_1 column present in log")
  expect_true("theta_2" %in% colnames(samp),
              label = "theta_2 column present in log")

  # Values should be in (0, 1) — theta is a Beta-distributed slab-shape param
  expect_true(all(samp[, "theta_1"] > 0 & samp[, "theta_1"] < 1),
              label = "theta_1 values in (0, 1)")
  expect_true(all(samp[, "theta_2"] > 0 & samp[, "theta_2"] < 1),
              label = "theta_2 values in (0, 1)")

  # The in-memory samples matrix should also carry the theta columns
  expect_true("theta_1" %in% colnames(res$samples),
              label = "theta_1 present in result$samples")
  expect_true("theta_2" %in% colnames(res$samples),
              label = "theta_2 present in result$samples")
})


test_that(".ParamNames emits theta_1 and theta_2 when nTheta = 2", {
  f <- .LoggingFixture()
  nms <- MkPrime:::.ParamNames(f$mkd, nEdge = 8L,
                               ecologyAware = TRUE,
                               nPhi = 1L,
                               nTheta = 2L)
  expect_true("pi0"     %in% nms)
  expect_true("theta_1" %in% nms)
  expect_true("theta_2" %in% nms)
  # theta columns must come after pi0
  idxPi0 <- which(nms == "pi0")
  expect_equal(which(nms == "theta_1"), idxPi0 + 1L)
  expect_equal(which(nms == "theta_2"), idxPi0 + 2L)
})
