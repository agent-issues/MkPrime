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
