library(TreeTools)

test_that("Invariant characters are dropped with warning", {
  # Char 1: all state 0 (invariant); Char 2: variable
  mat <- matrix(c(0, 0, 0, 0, 0,
                  0, 1, 0, 1, 2),
                nrow = 5, ncol = 2,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)

  expect_warning(mkd <- MkPrimeData(pd), "invariant")
  expect_equal(mkd$nChar, 1L)
  expect_equal(mkd$kObs, 3L)
})


test_that("All-invariant dataset is an error", {
  mat <- matrix(c(0, 0, 0, 0, 0,
                  1, 1, 1, 1, 1),
                nrow = 5, ncol = 2,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  expect_error(
    suppressWarnings(MkPrimeData(pd)),
    "No variable characters"
  )
})


test_that("Character with only missing data is treated as invariant", {
  mat <- matrix(c("?", "?", "?", "?", "?",
                  "0", "1", "0", "1", "2"),
                nrow = 5, ncol = 2,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  expect_warning(mkd <- MkPrimeData(pd), "invariant")
  expect_equal(mkd$nChar, 1L)
})


test_that("Neomorphic index remapping after invariant drop", {
  # Char 1: invariant; Char 2: neomorphic (binary); Char 3: variable
  mat <- matrix(c(1, 1, 1, 1, 1,
                  0, 1, 0, 1, 0,
                  0, 1, 2, 1, 0),
                nrow = 5, ncol = 3,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)

  # neomorphic = 2 in original; becomes 1 after char 1 is dropped
  expect_warning(mkd <- MkPrimeData(pd, neomorphic = 2L), "invariant")
  expect_equal(mkd$nChar, 2L)
  expect_equal(mkd$type, c("neomorphic", "transformational"))
})


test_that("knownStates index remapping after invariant drop", {
  mat <- matrix(c(0, 0, 0, 0, 0,
                  0, 1, 2, 1, 0,
                  0, 1, 0, 1, 0),
                nrow = 5, ncol = 3,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)

  # knownStates on char 3 (original); becomes char 2 after drop
  expect_warning(
    mkd <- MkPrimeData(pd, knownStates = c("3" = 4L)),
    "invariant"
  )
  expect_equal(mkd$nChar, 2L)
  expect_equal(mkd$type, c("transformational", "known"))
  expect_equal(mkd$known_k, c(NA_integer_, 4L))
})


test_that("Missing data mixed with observed states", {
  mat <- matrix(c("0", "?", "1", "?", "2",
                  "?", "0", "?", "1", "?"),
                nrow = 5, ncol = 2,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  # Char 1: states 0,1,2 observed (kObs=3)
  # Char 2: states 0,1 observed (kObs=2)
  expect_equal(mkd$kObs, c(3L, 2L))
  expect_equal(sum(is.na(mkd$matrix[, 1])), 2L)  # 2 NAs in char 1
  expect_equal(sum(is.na(mkd$matrix[, 2])), 3L)  # 3 NAs in char 2
})


test_that("Single taxon with all states works", {
  mat <- matrix(c(0, 0, 1), nrow = 3, ncol = 1,
                dimnames = list(paste0("t", 1:3), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  expect_equal(mkd$nTip, 3L)
  expect_equal(mkd$kObs, 2L)
})


test_that("Dropped invariant neomorphic char does not trigger binary warning", {
  # Char 1: invariant (all 0), marked neomorphic; Char 2: variable
  mat <- matrix(c(0, 0, 0, 0, 0,
                  0, 1, 2, 1, 0),
                nrow = 5, ncol = 2,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)

  # Should warn about invariant drop, but NOT about kObs != 2
  expect_warning(
    mkd <- MkPrimeData(pd, neomorphic = 1L),
    "invariant"
  )
  expect_equal(mkd$nChar, 1L)
  expect_equal(mkd$type, "transformational")
})
