library(TreeTools)

test_that("MkPrimeData classifies all characters as transformational by default", {
  mat <- matrix(c(0, 1, 0, 1, 2,
                  0, 0, 1, 1, 2,
                  1, 1, 0, 0, 1),
                nrow = 5, ncol = 3,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  expect_s3_class(mkd, "MkPrimeData")
  expect_equal(mkd$nTip, 5L)
  expect_equal(mkd$nChar, 3L)
  expect_equal(mkd$type, rep("transformational", 3))
  expect_equal(mkd$kObs, c(3L, 3L, 2L))
  expect_true(all(is.na(mkd$known_k)))
})


test_that("MkPrimeData classifies neomorphic characters", {
  mat <- matrix(c(0, 1, 0, 1, 0,
                  0, 0, 1, 1, 1,
                  0, 1, 0, 1, 0),
                nrow = 5, ncol = 3,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = c(1L, 3L))

  expect_equal(mkd$type, c("neomorphic", "transformational", "neomorphic"))
  expect_equal(mkd$kObs, c(2L, 2L, 2L))
})


test_that("MkPrimeData classifies known state space characters", {
  mat <- matrix(c(0, 1, 0, 1, 2,
                  0, 0, 1, 1, 2),
                nrow = 5, ncol = 2,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, knownStates = c("2" = 5L))

  expect_equal(mkd$type, c("transformational", "known"))
  expect_equal(mkd$known_k, c(NA_integer_, 5L))
})


test_that("MkPrimeData handles missing data correctly", {
  mat <- matrix(c("0", "1", "0", "1", "2",
                  "?", "0", "1", "1", "2"),
                nrow = 5, ncol = 2,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  expect_true(is.na(mkd$matrix[1, 2]))
  # kObs should exclude missing data
  expect_equal(mkd$kObs[2], 3L)
})


test_that("MkPrimeData rejects non-phyDat input", {
  expect_error(MkPrimeData(matrix(1:10, 5, 2)), "phyDat")
})


test_that("MkPrimeData rejects invalid neomorphic indices", {
  mat <- matrix(c(0, 1, 0, 1, 0), nrow = 5, ncol = 1,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  expect_error(MkPrimeData(pd, neomorphic = 2L), "between 1 and")
  expect_error(MkPrimeData(pd, neomorphic = 0L), "between 1 and")
})


test_that("MkPrimeData warns for non-binary neomorphic characters", {
  mat <- matrix(c(0, 1, 2, 1, 0), nrow = 5, ncol = 1,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  expect_warning(MkPrimeData(pd, neomorphic = 1L), "kObs != 2")
})


test_that("MkPrimeData rejects knownStates < kObs", {
  mat <- matrix(c(0, 1, 2, 1, 0), nrow = 5, ncol = 1,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  expect_error(MkPrimeData(pd, knownStates = c("1" = 2L)), "less than")
})


test_that("MkPrimeData rejects overlap between neomorphic and knownStates", {
  mat <- matrix(c(0, 1, 0, 1, 0), nrow = 5, ncol = 1,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  expect_error(
    MkPrimeData(pd, neomorphic = 1L, knownStates = c("1" = 2L)),
    "both"
  )
})


test_that("MkPrimeData integer matrix is 0-indexed", {
  mat <- matrix(c(0, 1, 2, 0, 1), nrow = 5, ncol = 1,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  expect_equal(range(mkd$matrix, na.rm = TRUE), c(0L, 2L))
})


# ===== AutoDetectNeomorphic =====

test_that("AutoDetectNeomorphic finds binary {0,1} characters", {
  tips <- paste0("t", 1:4)
  pd <- StringToPhyDat(setNames(c("01", "10", "01", "10"), tips),
                        tips = tips)
  neo <- AutoDetectNeomorphic(pd)
  expect_true(length(neo) > 0)
  # All detected characters should have exactly levels "0" and "1"
  mat <- MkPrime:::.PhyDatToIntMatrix(pd)
  lvls <- attr(pd, "levels")
  for (j in neo) {
    states <- unique(mat[, j][!is.na(mat[, j])])
    labels <- sort(lvls[states + 1L])
    expect_equal(labels, c("0", "1"))
  }
})

test_that("AutoDetectNeomorphic returns empty for non-binary data", {
  tips <- paste0("t", 1:4)
  pd <- StringToPhyDat(setNames(c("12", "23", "31", "12"), tips),
                        tips = tips)
  expect_length(AutoDetectNeomorphic(pd), 0L)
})

test_that("AutoDetectNeomorphic excludes chars with 0 and states > 1", {
  # MatrixToPhyDat gives cleaner control over states
  mat <- matrix(c(0, 1, 2, 0, 1), nrow = 5, ncol = 1,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  # This char has states {0, 1, 2} — not neomorphic
  expect_length(AutoDetectNeomorphic(pd), 0L)
})

test_that("AutoDetectNeomorphic works on real Nexus data", {
  skip_if_not(file.exists(
    system.file("datasets/Sun2018.nex", package = "TreeSearch")
  ))
  pd <- ReadAsPhyDat(
    system.file("datasets/Sun2018.nex", package = "TreeSearch")
  )
  neo <- AutoDetectNeomorphic(pd)
  # Should find a substantial number of neomorphic characters
  expect_true(length(neo) >= 100)
  # All detected chars should be binary {0,1} in levels
  mat  <- MkPrime:::.PhyDatToIntMatrix(pd)
  lvls <- attr(pd, "levels")
  lvls <- lvls[lvls != "-"]  # gap excluded, matching .PhyDatToIntMatrix indexing
  for (j in neo) {
    states <- unique(mat[, j][!is.na(mat[, j])])
    labels <- sort(lvls[states + 1L])
    expect_equal(labels, c("0", "1"))
  }
})

test_that("AutoDetectNeomorphic rejects non-phyDat input", {
  expect_error(AutoDetectNeomorphic("not_phyDat"), "phyDat")
})
