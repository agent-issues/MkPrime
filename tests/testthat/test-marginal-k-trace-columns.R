# Tests for Plan §7.4: no-op kPrime_* columns absent from trace under
# marginal_k mode (feat/marginal-k).
#
# Three subtests:
#   1. marginal_k: colnames(samples) has NO ^kPrime_ entries.
#   2. sampled_k:  colnames(samples) HAS  ^kPrime_ entries (regression guard).
#   3. ConvergenceDiagnostics on a marginal_k result returns without error
#      and does not include a kPrime summary row (empty group is silent).

library("ape")
library("TreeTools")

# ---------------------------------------------------------------------------
# Tiny fixture: 8 tips × 4 transformational chars (all binary, kObs = 2).
# Shared by all three subtests.
# ---------------------------------------------------------------------------

.mktc_tree <- function() {
  read.tree(text = paste0(
    "(((t1:0.05,t2:0.07):0.04,(t3:0.06,t4:0.05):0.03):0.05,",
    "((t5:0.04,t6:0.06):0.05,(t7:0.05,t8:0.04):0.06):0.04);"
  ))
}

.mktc_mkd <- function() {
  tips <- paste0("t", 1:8)
  mat <- matrix(
    c(0, 0, 1, 1, 0, 1, 0, 1,
      0, 1, 1, 0, 1, 0, 0, 1,
      0, 1, 0, 0, 1, 1, 1, 0,
      1, 1, 0, 1, 0, 0, 1, 0),
    nrow = 8, ncol = 4,
    dimnames = list(tips, NULL)
  )
  MatrixToPhyDat(mat)
}

# Short MCMC settings that complete quickly in tests.
.mktc_mcmc <- function() {
  MkPrimeMCMC(nIter = 200L, thin = 10L,
              maxWarmup = 100L, minWarmup = 100L,
              autoTune = FALSE)
}


# ---------------------------------------------------------------------------
# Subtest 1: marginal_k produces no kPrime_* columns
# ---------------------------------------------------------------------------

test_that("marginal-k mode produces a trace with no kPrime_* columns", {
  set.seed(7421L)
  result <- RunMkPrime(
    .mktc_mkd(), .mktc_tree(),
    model = MkPrimeModel(kPrimePrior = "geometric",
                         likelihoodMode = "marginal_k",
                         coding = "none",
                         relabel = FALSE),
    mcmc = .mktc_mcmc()
  )

  cn <- colnames(result$samples)
  kp_cols <- grep("^kPrime_", cn, value = TRUE)
  expect_length(kp_cols, 0L)
})


# ---------------------------------------------------------------------------
# Subtest 2: sampled_k still emits kPrime_* columns (regression guard)
# ---------------------------------------------------------------------------

test_that("sampled-k mode still emits kPrime_* columns (regression guard)", {
  set.seed(7422L)
  result <- RunMkPrime(
    .mktc_mkd(), .mktc_tree(),
    model = MkPrimeModel(kPrimePrior = "geometric",
                         likelihoodMode = "sampled_k",
                         coding = "none",
                         relabel = FALSE),
    mcmc = .mktc_mcmc()
  )

  cn <- colnames(result$samples)
  kp_cols <- grep("^kPrime_", cn, value = TRUE)
  # 4 transformational characters → expect kPrime_1 .. kPrime_4
  expect_true(length(kp_cols) > 0L)
})


# ---------------------------------------------------------------------------
# Subtest 3: ConvergenceDiagnostics on marginal_k result is error-free
# (empty kPrime group should be silently absent, not crash or print "kPrime 0")
# ---------------------------------------------------------------------------

test_that("ConvergenceDiagnostics handles empty kPrime group under marginal-k", {
  set.seed(7423L)
  result <- RunMkPrime(
    .mktc_mkd(), .mktc_tree(),
    model = MkPrimeModel(kPrimePrior = "geometric",
                         likelihoodMode = "marginal_k",
                         coding = "none",
                         relabel = FALSE),
    mcmc = .mktc_mcmc()
  )

  # Must complete without error
  diag <- expect_no_error(ConvergenceDiagnostics(result))

  # No kPrime entries in the ESS vector
  kp_names <- grep("^kPrime_", names(diag$ess), value = TRUE)
  expect_length(kp_names, 0L)

  # minEss is finite (not poisoned by a missing kPrime group)
  expect_true(is.finite(diag$minEss))
})
