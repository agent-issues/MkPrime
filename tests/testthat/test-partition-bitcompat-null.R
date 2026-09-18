# §7a: bit-identity guarantee for the legacy partition = NULL code path.
#
# When `partition = NULL`, RunMkPrime must route through the unchanged legacy
# code path (no per-class state, no eta_neo). After every edit on this branch
# the locked reference must reproduce exactly.
#
# If this test fails on your branch:
#   (1) Confirm the failure is intentional (e.g. you introduced a justified
#       behavioural change to the legacy path); discuss with reviewers.
#   (2) Before concluding it is intentional, check that the difference is
#       really the package and not the environment: build the commit that
#       generated the reference and run the fixture under it. If that build
#       also fails to reproduce the stored values, the fixture has an
#       unpinned input and regenerating would hide nothing — fix the input
#       instead. That is how the 2026-09-18 breakage was diagnosed.
#   (3) If it is a genuine, approved behavioural change, prefer PINNING the
#       offending move off in helper-partition-ref.R to regenerating. The
#       reference guards the legacy path as a whole, and every regeneration
#       is indistinguishable from evading the guard. Regenerate only after
#       showing, on a build of the merge base and a build of your branch,
#       that the pinned schedule gives bit-identical payloads; then run
#       `Rscript tests/testthat/_reference/generate-partition-bitcompat-null.R`
#       once and add a NEWS.md entry recording the change.
#   (4) If unintentional, the failure indicates drift in the legacy code
#       path that the partition API was contracted not to introduce — fix.
#
# The reference (samples matrix) was first generated on 2026-05-20 from the
# pre-partition-API state of feature/partition-api (commit 6f7d36f), and
# regenerated once on 2026-09-18 when the fixture was made hermetic: it now
# pins its starting tree (_reference/partition-bitcompat-null-start.nwk) and
# pins `gibbsSpr = FALSE`. The previous reference delegated its starting
# topology to TreeSearch::AdditionTree(), so it silently depended on the
# installed version of a Suggests package and no CI job could ever reproduce
# it; see helper-partition-ref.R.
#
# What the pinned schedule does NOT cover: `gibbs_spr` (pinned off),
# `gibbs_subtree_swap`, `joint2d`, and the default-off weighted / block-Gibbs
# moves. Everything else in the legacy path — the starting tree, the four
# unweighted topology moves (nni / spr / tbr / pspr), the branch-length and
# Dirichlet moves, the k-prime and rate moves, the likelihood, and the sample
# schema — is still asserted bit-for-bit.
#
# The comparator is bound to the build that generated it: a last-bit
# difference early in the chain flips one accept/reject and the trajectory
# diverges wholesale. Generate and check it with the same installed build
# (build-agent.sh / test-agent.sh), never devtools::load_all().
#
# It is also bound to the CPU architecture that generated it, which is why the
# value assertion is bit-exact only on a matching arch and tolerant elsewhere.
# See .PartitionBitcompatReferenceArch() in helper-partition-ref.R.

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

  # The load-bearing assertion. Bit-identity is asserted only on the
  # architecture that generated the reference; see
  # .PartitionBitcompatReferenceArch() in helper-partition-ref.R for why the
  # two cannot both hold, and why the tolerant branch still detects drift.
  if (identical(R.version$arch, .PartitionBitcompatReferenceArch())) {
    expect_identical(result$samples, ref$samples)
  } else {
    expect_equal(result$samples, ref$samples, tolerance = 1e-9)
  }
})
