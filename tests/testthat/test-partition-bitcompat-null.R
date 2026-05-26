# §7a: bit-identity guarantee for the legacy partition = NULL code path.
#
# When `partition = NULL`, RunMkPrime must route through the unchanged legacy
# code path (no per-class state, no eta_neo). After every edit on this branch
# the locked reference must reproduce exactly.
#
# If this test fails on your branch:
#   (1) Confirm the failure is intentional (e.g. you introduced a justified
#       behavioural change to the legacy path); discuss with reviewers.
#   (2) If intentional and approved, regenerate the reference by running
#       `Rscript tests/testthat/_reference/generate-partition-bitcompat-null.R`
#       and add a NEWS.md entry recording the change.
#   (3) If unintentional, the failure indicates drift in the legacy code
#       path that the partition API was contracted not to introduce — fix.
#
# The reference (samples matrix) was generated on 2026-05-20 from the
# pre-partition-API state of feature/partition-api (commit 6f7d36f).

test_that("§7a bit-identity: partition = NULL reproduces stored reference", {
  ref_path <- .PartitionBitcompatReferencePath()
  testthat::skip_if_not(file.exists(ref_path),
    paste("Reference RDS not found at", ref_path,
          "— regenerate via tests/testthat/_reference/generate-partition-bitcompat-null.R."))

  ref <- readRDS(ref_path)

  result <- .RunPartitionBitcompatReference()

  # Schema first: a column-name or dimension diff is easier to read than a
  # value diff. Catch schema drift in .ParamNames() / .StateToRow() early.
  expect_identical(colnames(result$samples), ref$samples_cols)
  expect_identical(dim(result$samples),       ref$samples_dim)
  expect_identical(result$nSamples,           ref$nSamples)
  expect_identical(result$stop_reason,        ref$stop_reason)
  expect_identical(result$actual_iter,        ref$actual_iter)
  expect_identical(result$treeThin,           ref$treeThin)

  # The load-bearing assertion: every sampled value must match bit-for-bit.
  expect_identical(result$samples, ref$samples)
})
