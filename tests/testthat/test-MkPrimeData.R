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
  mkd <- MkPrimeData(pd, known_states = c("2" = 5L))

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


test_that("MkPrimeData rejects known_states < kObs", {
  mat <- matrix(c(0, 1, 2, 1, 0), nrow = 5, ncol = 1,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  expect_error(MkPrimeData(pd, known_states = c("1" = 2L)), "less than")
})


test_that("MkPrimeData rejects overlap between neomorphic and known_states", {
  mat <- matrix(c(0, 1, 0, 1, 0), nrow = 5, ncol = 1,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  expect_error(
    MkPrimeData(pd, neomorphic = 1L, known_states = c("1" = 2L)),
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
