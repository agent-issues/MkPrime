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
# Before concluding (2), check that the difference really is the package and
# not the environment: build the commit that generated the reference and run
# the fixture under it. If that also fails to reproduce the stored values,
# the fixture has an unpinned input and regenerating would hide nothing —
# fix the input instead. That is how the 2026-09-18 breakage was diagnosed.
#
# The reference (samples matrix) was regenerated on 2026-09-18 against a
# pinned starting tree (_reference/partition-bitcompat-null-start.nwk). The
# previous reference dated from 2b3c054 and delegated its starting topology
# to TreeSearch::AdditionTree(), so it silently depended on the installed
# version of a Suggests package; see helper-partition-ref.R.

test_that("§7a bit-identity: partition = NULL reproduces stored reference", {
  ref_path <- .PartitionBitcompatReferencePath()
  testthat::skip_if_not(file.exists(ref_path),
    paste("Reference RDS not found at", ref_path,
          "— regenerate via tests/testthat/_reference/generate-partition-bitcompat-null.R."))

  tree_path <- .PartitionBitcompatStartTreePath()
  testthat::skip_if_not(file.exists(tree_path),
    paste("Pinned starting tree not found at", tree_path,
          "— the fixture cannot be reproduced without it."))

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
