# Test-suite wall-clock discipline
#
# Design rule: every test_that() block must complete in < 30 s on a typical
# developer machine.  Enforce this by keeping MCMC iteration counts low:
#
#   - nIter  <= 1000  for single-run tests
#   - nIter  <=  500  for multi-run (nRuns >= 2) tests
#   - nIter  <=  200  for tests that call treess (O(n^2) in tree count)
#
# If a test genuinely needs more iterations (e.g. testing stopping rules),
# use maxTime = 0.5 instead of a large nIter, so it self-limits.
#
# For any subprocess Rscript call made by agents during development, always
# add setTimeLimit(elapsed = 60, transient = FALSE) at the top of the -e
# block, or wrap with: timeout 60 Rscript -e "..."
# (Rtools45 ships coreutils 'timeout' on PATH).

# Shared tiny dataset re-used across many tests to avoid repeated construction.
.mkp_test_tree <- function() {
  ape::read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
}
.mkp_test_pd <- function() {
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  TreeTools::MatrixToPhyDat(mat)
}
