# Tests for .BuildPartitions(mkd, partition = ...) — the user-class sub-grouping
# added in Layer 1 plumbing. The existing tests in test-partition.R cover the
# partition = NULL path; these cover the partition-aware path before any of
# the consuming machinery (eta_neo, unlink, sibling C++) lands.

library("TreeTools")

# Helper: 5 tips × N chars, no invariants.
.make_mkd_nchar <- function(nChar, seed = 17L, neoIdx = integer(0)) {
  set.seed(seed)
  mat <- matrix(sample(0:1, 5L * nChar, replace = TRUE),
                nrow = 5L, ncol = nChar,
                dimnames = list(paste0("t", 1:5), NULL))
  # Force every char to be variable
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1L, j] <- 1L - mat[1L, j]
  }
  MkPrimeData(MatrixToPhyDat(mat), neomorphic = neoIdx)
}


test_that("partition = NULL gives every partition classIdx = 1L", {
  mkd <- .make_mkd_nchar(6L)
  parts <- mkd$partitions  # built by MkPrimeData via .BuildPartitions(mkd)
  for (p in parts) {
    expect_identical(p$classIdx, 1L)
  }
})


test_that("partition assignment splits a single (type, kObs) group across classes", {
  mkd <- .make_mkd_nchar(6L)
  # All chars are binary transformational → one base partition; split into 3 classes
  partition <- c(1L, 1L, 2L, 2L, 3L, 3L)
  parts <- .BuildPartitions(mkd, partition = partition)

  expect_length(parts, 3L)
  expect_identical(vapply(parts, `[[`, integer(1), "classIdx"), 1:3)
  expect_identical(vapply(parts, `[[`, integer(1), "nChar"), c(2L, 2L, 2L))

  # Union of char_indices == all chars; no overlap
  unioned <- sort(unlist(lapply(parts, `[[`, "char_indices")))
  expect_identical(unioned, seq_len(mkd$nChar))
})


test_that("a user class spanning multiple kObs values yields one PartInfo per kObs", {
  # Construct a matrix with mixed kObs in one class
  mat <- matrix(c(
    0, 1, 0, 1, 0,   # binary, transformational
    0, 1, 2, 1, 0,   # 3-state, transformational
    0, 1, 0, 1, 0    # binary, transformational
  ), nrow = 5L, ncol = 3L, dimnames = list(paste0("t", 1:5), NULL))
  mkd <- MkPrimeData(MatrixToPhyDat(mat))
  expect_identical(mkd$nChar, 3L)

  # All three in class 1
  parts <- .BuildPartitions(mkd, partition = rep(1L, 3L))

  # Expect 2 PartInfos (kObs = 2 and kObs = 3), both classIdx = 1L
  expect_length(parts, 2L)
  expect_setequal(vapply(parts, `[[`, integer(1), "kObs"), c(2L, 3L))
  expect_identical(unique(vapply(parts, `[[`, integer(1), "classIdx")), 1L)
})


test_that("mixed-type class (neomorphic + transformational) yields multiple PartInfos with same classIdx", {
  # 4 chars: chars 1 and 2 are neomorphic; chars 3 and 4 are transformational.
  # Put one of each into class 1, the other of each into class 2.
  mat <- matrix(c(
    0, 1, 0, 1, 0,   # char 1: neo, in class 1
    0, 1, 1, 1, 0,   # char 2: neo, in class 2
    0, 1, 2, 1, 0,   # char 3: trans, in class 1
    1, 0, 1, 0, 1    # char 4: trans, in class 2
  ), nrow = 5L, ncol = 4L, dimnames = list(paste0("t", 1:5), NULL))
  mkd <- MkPrimeData(MatrixToPhyDat(mat), neomorphic = c(1L, 2L))

  parts <- .BuildPartitions(mkd, partition = c(1L, 2L, 1L, 2L))

  # Expect 4 PartInfos: (class 1, neo), (class 1, trans), (class 2, neo), (class 2, trans)
  expect_length(parts, 4L)
  pf <- function(field) vapply(parts, `[[`, vector(typeof(parts[[1]][[field]]), 1), field)
  expect_identical(pf("classIdx"), c(1L, 1L, 2L, 2L))
  expect_identical(pf("type"),
                   c("neomorphic", "transformational", "neomorphic", "transformational"))
  expect_identical(pf("nChar"), c(1L, 1L, 1L, 1L))
})


test_that("an empty class is omitted (no PartInfo emitted for missing class IDs)", {
  # 4 chars, partition = c(1, 1, 1, 1) but caller passes class indices that
  # span 1:3 with class 2 empty? — actually the validator forbids this and
  # contiguity is enforced by .ValidatePartitionArgs. Here we test that
  # .BuildPartitions itself just emits no partition for empty classes when
  # called directly with a sparse partition (defensive).
  mkd <- .make_mkd_nchar(4L)
  # Skipping class 2 (only 1 and 3 populated). .BuildPartitions doesn't
  # validate contiguity itself; emits 2 non-empty groups.
  parts <- .BuildPartitions(mkd, partition = c(1L, 1L, 3L, 3L))
  expect_length(parts, 2L)
  expect_setequal(vapply(parts, `[[`, integer(1), "classIdx"), c(1L, 3L))
})


test_that("partition-aware tip_states subset matches the original matrix columns", {
  mkd <- .make_mkd_nchar(6L)
  partition <- c(1L, 2L, 1L, 2L, 1L, 2L)
  parts <- .BuildPartitions(mkd, partition = partition)
  for (p in parts) {
    expect_identical(p$tip_states, mkd$matrix[, p$char_indices, drop = FALSE])
    # All chars in a partition agree on classIdx
    expect_true(all(partition[p$char_indices] == p$classIdx))
  }
})


test_that("MkPrimeData still produces classIdx = 1L on every partition", {
  # Ensure existing MkPrimeData construction (which calls .BuildPartitions
  # with no partition arg) emits the new field everywhere.
  mat <- matrix(c(0, 1, 0, 1, 2,
                  0, 0, 1, 1, 2,
                  1, 1, 0, 0, 1),
                nrow = 5, ncol = 3,
                dimnames = list(paste0("t", 1:5), NULL))
  mkd <- MkPrimeData(MatrixToPhyDat(mat))
  for (p in mkd$partitions) {
    expect_true("classIdx" %in% names(p))
    expect_identical(p$classIdx, 1L)
  }
})
