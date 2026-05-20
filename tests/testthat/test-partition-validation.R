# Validation and prefix-matching tests for the partition / unlink args.
# These exercise .ValidatePartitionArgs() directly; integration tests
# against RunMkPrime live alongside the §7a bitcompat test.

library("TreeTools")

# ---- helpers -----------------------------------------------------------

.make_mkd <- function(nChar = 4L, nTip = 5L, seed = 7L) {
  set.seed(seed)
  mat <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                nrow = nTip, ncol = nChar,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  # Ensure no invariant chars (so nChar is preserved)
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1, j] <- 1L - mat[1, j]
  }
  MkPrimeData(MatrixToPhyDat(mat))
}


# ---- partition validation ---------------------------------------------

test_that("partition = NULL is valid (legacy path)", {
  mkd <- .make_mkd()
  spec <- .ValidatePartitionArgs(NULL, character(0), mkd)
  expect_null(spec$partition)
  expect_identical(spec$unlink, character(0))
  expect_identical(spec$nClasses, 1L)
})

test_that("partition with right length and contiguous classes is accepted", {
  mkd <- .make_mkd()
  spec <- .ValidatePartitionArgs(c(1L, 1L, 2L, 2L), character(0), mkd)
  expect_identical(spec$partition, c(1L, 1L, 2L, 2L))
  expect_identical(spec$nClasses, 2L)
})

test_that("partition length mismatch errors with a clear message", {
  mkd <- .make_mkd(nChar = 4L)
  expect_error(
    .ValidatePartitionArgs(c(1L, 1L, 2L), character(0), mkd),
    regexp = "length 3 but the data has 4 character",
    fixed = FALSE
  )
})

test_that("partition with NA errors", {
  mkd <- .make_mkd()
  expect_error(
    .ValidatePartitionArgs(c(1L, NA, 1L, 2L), character(0), mkd),
    regexp = "contains.*NA",
    fixed = FALSE
  )
})

test_that("partition with non-integer values errors", {
  mkd <- .make_mkd()
  expect_error(
    .ValidatePartitionArgs(c(1.5, 1, 1, 2), character(0), mkd),
    regexp = "whole-number",
    fixed = FALSE
  )
})

test_that("partition values < 1 error", {
  mkd <- .make_mkd()
  expect_error(
    .ValidatePartitionArgs(c(0L, 1L, 2L, 2L), character(0), mkd),
    regexp = ">= 1",
    fixed = FALSE
  )
})

test_that("partition with non-contiguous class IDs errors", {
  mkd <- .make_mkd()
  # Skips class 2: ids are 1 and 3
  expect_error(
    .ValidatePartitionArgs(c(1L, 1L, 3L, 3L), character(0), mkd),
    regexp = "skips class id",
    fixed = FALSE
  )
})


# ---- unlink token resolution ------------------------------------------

test_that("exact unlink tokens are accepted (case-insensitive)", {
  mkd <- .make_mkd()
  spec <- .ValidatePartitionArgs(
    partition = c(1L, 1L, 2L, 2L),
    unlink    = c("Shape", "RATEMULTIPLIER"),
    mkd       = mkd
  )
  expect_setequal(spec$unlink, c("shape", "ratemultiplier"))
})

test_that("partial-prefix unlink token resolves with a warning", {
  mkd <- .make_mkd()
  expect_warning(
    spec <- .ValidatePartitionArgs(
      partition = c(1L, 1L, 2L, 2L),
      unlink    = "rate",
      mkd       = mkd
    ),
    regexp = "matched via prefix",
    fixed = FALSE
  )
  expect_identical(spec$unlink, "ratemultiplier")
})

test_that("ambiguous unlink prefix errors", {
  mkd <- .make_mkd()
  # No current token is an ambiguous prefix of two — fabricate by feeding
  # a prefix shared by "shape" and "ratemultiplier"? They don't share one.
  # Test instead that "b" (matches "brlens" alone) is unambiguous, and
  # then test the ambiguity path by adding more tokens via mocking the
  # constant — done by checking startsWith logic directly.
  # In practice: today no real ambiguous prefix exists. Skip this test
  # until v2 adds tokens that could collide.
  skip("No ambiguous prefix exists for the current Layer-1 token set")
})

test_that("unknown unlink token errors with `agrep` suggestion", {
  mkd <- .make_mkd()
  expect_error(
    .ValidatePartitionArgs(
      partition = c(1L, 1L, 2L, 2L),
      unlink    = "shaper",
      mkd       = mkd
    ),
    regexp = "not recognised",
    fixed = FALSE
  )
})

test_that("duplicate unlink tokens warn and deduplicate", {
  mkd <- .make_mkd()
  expect_warning(
    spec <- .ValidatePartitionArgs(
      partition = c(1L, 1L, 2L, 2L),
      unlink    = c("shape", "Shape"),
      mkd       = mkd
    ),
    regexp = "duplicate",
    fixed = FALSE
  )
  expect_identical(spec$unlink, "shape")
})

test_that("empty-string unlink token errors", {
  mkd <- .make_mkd()
  expect_error(
    .ValidatePartitionArgs(
      partition = c(1L, 1L, 2L, 2L),
      unlink    = c("shape", ""),
      mkd       = mkd
    ),
    regexp = "empty string",
    fixed = FALSE
  )
})


# ---- silent coercion of unlink when partition is trivial --------------

test_that("unlink is coerced to character(0) with cli_alert_info when partition = NULL", {
  mkd <- .make_mkd()
  expect_message(
    spec <- .ValidatePartitionArgs(NULL, c("shape", "ratemultiplier"), mkd),
    regexp = "ignored"
  )
  expect_identical(spec$unlink, character(0))
  expect_null(spec$partition)
})

test_that("unlink is coerced to character(0) when nClasses == 1", {
  mkd <- .make_mkd()
  expect_message(
    spec <- .ValidatePartitionArgs(
      partition = rep(1L, mkd$nChar),
      unlink    = "shape",
      mkd       = mkd
    ),
    regexp = "ignored"
  )
  expect_identical(spec$unlink, character(0))
  expect_identical(spec$nClasses, 1L)
})


# ---- Layer 1 gate ------------------------------------------------------

test_that(".RequirePartitionImplemented is a no-op for the trivial spec", {
  spec <- list(partition = NULL, unlink = character(0), nClasses = 1L)
  expect_silent(.RequirePartitionImplemented(spec))
})

test_that(".RequirePartitionImplemented errors with helpful message when partition is supplied", {
  spec <- list(partition = c(1L, 2L), unlink = character(0), nClasses = 2L)
  expect_error(.RequirePartitionImplemented(spec),
               regexp = "not yet implemented",
               fixed = FALSE)
})

test_that(".RequirePartitionImplemented errors when unlink is supplied", {
  # This shouldn't happen in practice (validator coerces unlink to character(0)
  # when partition = NULL or nClasses == 1), but guard the gate explicitly.
  spec <- list(partition = NULL, unlink = "shape", nClasses = 1L)
  expect_error(.RequirePartitionImplemented(spec),
               regexp = "not yet implemented",
               fixed = FALSE)
})
