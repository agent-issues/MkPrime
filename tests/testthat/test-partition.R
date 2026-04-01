library("TreeTools")

test_that("Partitions group transformational characters by kObs", {
  # 3 chars: kObs = 3, 3, 2
  mat <- matrix(c(0, 1, 0, 1, 2,
                  0, 0, 1, 1, 2,
                  1, 1, 0, 0, 1),
                nrow = 5, ncol = 3,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  expect_length(mkd$partitions, 2L)  # kObs=2 and kObs=3

  # Sorted by kObs
  expect_equal(mkd$partitions[[1]]$kObs, 2L)
  expect_equal(mkd$partitions[[1]]$char_indices, 3L)
  expect_equal(mkd$partitions[[1]]$nChar, 1L)

  expect_equal(mkd$partitions[[2]]$kObs, 3L)
  expect_equal(mkd$partitions[[2]]$char_indices, c(1L, 2L))
  expect_equal(mkd$partitions[[2]]$nChar, 2L)
})


test_that("Neomorphic characters form a single partition", {
  mat <- matrix(c(0, 1, 0, 1, 0,
                  0, 0, 1, 1, 1,
                  0, 1, 0, 1, 0),
                nrow = 5, ncol = 3,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = c(1L, 3L))

  # 2 partitions: 1 neomorphic (chars 1,3), 1 transformational (char 2)
  expect_length(mkd$partitions, 2L)

  neo_part <- mkd$partitions[[1]]
  expect_equal(neo_part$type, "neomorphic")
  expect_equal(neo_part$kObs, 2L)
  expect_equal(neo_part$char_indices, c(1L, 3L))
  expect_equal(neo_part$nChar, 2L)
  expect_true(is.na(neo_part$k))
})


test_that("Known characters group by known_k", {
  mat <- matrix(c(0, 1, 2, 1, 0,
                  0, 1, 0, 1, 0,
                  0, 1, 2, 0, 1),
                nrow = 5, ncol = 3,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, knownStates = c("1" = 5L, "2" = 3L, "3" = 5L))

  # All known; group by known_k: k=3 (char 2) and k=5 (chars 1, 3)
  expect_length(mkd$partitions, 2L)

  expect_equal(mkd$partitions[[1]]$type, "known")
  expect_equal(mkd$partitions[[1]]$k, 3L)
  expect_equal(mkd$partitions[[1]]$char_indices, 2L)

  expect_equal(mkd$partitions[[2]]$type, "known")
  expect_equal(mkd$partitions[[2]]$k, 5L)
  expect_equal(mkd$partitions[[2]]$char_indices, c(1L, 3L))
})


test_that("Mixed types produce correct partition ordering", {
  mat <- matrix(c(0, 1, 0, 1, 0,   # binary, will be neomorphic
                  0, 1, 2, 1, 0,    # 3-state, transformational
                  0, 1, 0, 1, 0),   # binary, known k=4
                nrow = 5, ncol = 3,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = 1L, knownStates = c("3" = 4L))

  # Ordering: neomorphic, transformational, known
  expect_length(mkd$partitions, 3L)
  expect_equal(mkd$partitions[[1]]$type, "neomorphic")
  expect_equal(mkd$partitions[[2]]$type, "transformational")
  expect_equal(mkd$partitions[[3]]$type, "known")
})


test_that("Partition tip_states match original matrix subset", {
  mat <- matrix(c(0, 1, 0, 1, 2,
                  0, 0, 1, 1, 2,
                  1, 1, 0, 0, 1),
                nrow = 5, ncol = 3,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  for (part in mkd$partitions) {
    expect_equal(
      part$tip_states,
      mkd$matrix[, part$char_indices, drop = FALSE]
    )
  }
})


# ---------------------------------------------------------------------------
# M-172: unique_tip_states and pattern_index
# ---------------------------------------------------------------------------

test_that("pattern_index and unique_tip_states are correct for all-unique columns", {
  mat <- matrix(c(0, 1, 0, 1, 2,
                  0, 0, 1, 1, 2,
                  1, 1, 0, 0, 1),
                nrow = 5, ncol = 3,
                dimnames = list(paste0("t", 1:5), NULL))
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  for (part in mkd$partitions) {
    # When all columns are distinct, nUnique == nChar
    expect_equal(ncol(part$unique_tip_states), part$nChar)
    expect_equal(length(part$pattern_index), part$nChar)
    expect_true(all(part$pattern_index >= 0L))
    expect_true(all(part$pattern_index < part$nChar))
    # unique_tip_states columns are a subset of tip_states columns
    expect_equal(
      part$unique_tip_states,
      part$tip_states[, part$pattern_index + 1L, drop = FALSE]
    )
  }
})

test_that("pattern_index deduplicates identical columns", {
  # Columns 1 and 3 are identical; column 2 is distinct
  mat <- matrix(c(0, 1, 0, 1, 0,   # col 1
                  0, 0, 1, 1, 0,   # col 2
                  0, 1, 0, 1, 0),  # col 3 == col 1
                nrow = 5, ncol = 3,
                dimnames = list(paste0("t", 1:5), NULL))
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  # All three are binary transformational → one partition
  part <- mkd$partitions[[1]]
  expect_equal(part$nChar, 3L)
  expect_equal(ncol(part$unique_tip_states), 2L)  # 2 unique patterns
  # Columns 1 and 3 share the same pattern index; col 2 is different
  expect_equal(part$pattern_index[1], part$pattern_index[3])
  expect_false(part$pattern_index[1] == part$pattern_index[2])
})

test_that("pattern_index treats different NA placement as distinct patterns", {
  mat <- matrix(c("0", "1", "?", "1", "0",
                  "0", "1", "1", "?", "0"),
                nrow = 5, ncol = 2,
                dimnames = list(paste0("t", 1:5), NULL))
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  part <- mkd$partitions[[1]]
  # NA is in different positions → distinct patterns
  expect_equal(ncol(part$unique_tip_states), 2L)
  expect_false(part$pattern_index[1] == part$pattern_index[2])
})

test_that("unique_tip_states reconstruction matches tip_states via pattern_index", {
  set.seed(4419)
  mat <- matrix(sample(0:2, 10 * 8, replace = TRUE), 10, 8,
                dimnames = list(paste0("t", 1:10), NULL))
  pd  <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  for (part in mkd$partitions) {
    # Expanding unique_tip_states via pattern_index must recover tip_states
    reconstructed <- part$unique_tip_states[, part$pattern_index + 1L, drop = FALSE]
    expect_equal(reconstructed, part$tip_states)
  }
})

test_that("Partition tip_states preserve NA for ambiguous data", {
  mat <- matrix(c("0", "1", "?", "1", "2",
                  "0", "0", "1", "?", "2"),
                nrow = 5, ncol = 2,
                dimnames = list(paste0("t", 1:5), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  # Both chars have kObs=3, so one partition
  expect_length(mkd$partitions, 1L)
  tip_st <- mkd$partitions[[1]]$tip_states
  expect_true(is.na(tip_st[3, 1]))  # t3, char 1
  expect_true(is.na(tip_st[4, 2]))  # t4, char 2
})
