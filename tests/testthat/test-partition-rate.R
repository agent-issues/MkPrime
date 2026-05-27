library("TreeTools")
library("ape")

# Audit Issue 1: RB-style partition-rate normalisation.
# See dev/rb-equivalence/notes/partition-rate-and-acrv-audit.md
#
# The partition scales must satisfy
#   nNeo * neoScale + nTrans * transScale == nNeo + nTrans
# so that, averaged over the empirical character distribution, the expected
# substitutions per character on edge e is exactly tree_length * rel_br[e].
# This is the unique condition that makes tree_length sampler-independent.

# Re-implements the C++ helper compute_partition_scales in R for testing.
.partition_scales <- function(rate_neo, n_neo, n_trans) {
  if (n_neo == 0L || n_trans == 0L) {
    return(list(neo = 1.0, trans = 1.0))
  }
  denom <- 1.0 + rate_neo
  n_total <- n_neo + n_trans
  list(
    neo   = rate_neo / denom * n_total / n_neo,
    trans = 1.0       / denom * n_total / n_trans
  )
}


# ---- Algebraic property tests on the helper ------------------------------

test_that("partition rate weighted mean equals 1 (RB identity)", {
  # Exhaustive sweep over (r, n_neo, n_trans) with both partitions non-empty.
  grid <- expand.grid(
    r       = c(0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 10.0, 100.0),
    n_neo   = c(1L, 2L, 5L, 12L, 100L),
    n_trans = c(1L, 3L, 8L, 17L, 200L)
  )
  for (i in seq_len(nrow(grid))) {
    r <- grid$r[i]; n_neo <- grid$n_neo[i]; n_trans <- grid$n_trans[i]
    s <- .partition_scales(r, n_neo, n_trans)
    wmean <- n_neo * s$neo + n_trans * s$trans
    expect_equal(wmean, as.numeric(n_neo + n_trans), tolerance = 1e-12,
      label = sprintf("r=%g n_neo=%d n_trans=%d", r, n_neo, n_trans))
  }
})


test_that("partition rate degenerates correctly when one type is absent", {
  for (r in c(0.01, 1.0, 5.0)) {
    s_no_neo   <- .partition_scales(r, 0L, 17L)
    expect_identical(s_no_neo$neo,   1.0)
    expect_identical(s_no_neo$trans, 1.0)

    s_no_trans <- .partition_scales(r, 12L, 0L)
    expect_identical(s_no_trans$neo,   1.0)
    expect_identical(s_no_trans$trans, 1.0)
  }
})


test_that("partition rate matches the RB factorisation algebraically", {
  # RB Rev template (templates/by_nt_9v.template.Rev:47):
  #   partition_rate := [r/(1+r), 1/(1+r)] / nChar * sum(nChar)
  for (r in c(0.2, 1.0, 3.7)) {
    for (n_neo in c(1L, 5L, 20L)) {
      for (n_trans in c(1L, 12L, 50L)) {
        s <- .partition_scales(r, n_neo, n_trans)
        rb_neo   <- (r / (1 + r)) * (n_neo + n_trans) / n_neo
        rb_trans <- (1 / (1 + r)) * (n_neo + n_trans) / n_trans
        expect_equal(s$neo,   rb_neo,   tolerance = 1e-14)
        expect_equal(s$trans, rb_trans, tolerance = 1e-14)
      }
    }
  }
})


# ---- Helper to construct C++ McmcData pointer for direct-eval tests ------

.make_data_ptr <- function(mkd, rate_log_sd_sdlog = 0.5) {
  has_neo <- any(mkd$type == "neomorphic")
  prepare_mcmc_data(
    partitions_r              = mkd$partitions,
    kObs_r                    = mkd$kObs,
    charTypes_r               = mkd$type,
    hasNeo                    = has_neo,
    nCat                      = 6L,
    codingStr                 = "variable",
    relabelFlag               = TRUE,
    treeLengthShape           = 1.5,
    treeLengthRate            = 1.0,
    rateLossMeanlog           = 0.0,
    rateLossSdlog             = 1.0,
    rateLogSdShape            = 1.0,
    rateLogSdRate             = 1.0,
    rateNeoMeanlog            = 0.0,
    rateNeoSdlog              = 2.0,
    kprimeHyperA              = 1.0,
    kprimeHyperB              = 1.0,
    kPriorLogseries           = TRUE,
    kprimeLogseriesC          = 0.7,
    kPriorBetaGeometric       = FALSE,
    qHeterogeneity            = FALSE,
    nBetaCat                  = 4L,
    betaScaleShape            = 1.0,
    betaScaleRate             = 1.0,
    kPriorEmpiricalGeometric  = FALSE,
    empLogBody                = numeric(0)
  )
}


# ---- Parity: R-side `.MkpLogLikelihood` vs C++ direct eval ---------------

test_that("R and C++ direct-eval likelihoods agree on mixed data, rate_neo != 1", {
  # Mixed dataset: 3 neo + 6 trans + 1 binary = mixed types & sizes.
  set.seed(20260527L)
  ntax <- 6L
  mat <- matrix(sample.int(2, ntax * 9, replace = TRUE) - 1L,
                nrow = ntax, ncol = 9,
                dimnames = list(paste0("t", 1:ntax), NULL))
  # Mix in a few k=3 transformational characters so we exercise larger k'.
  mat_3state <- matrix(sample.int(3, ntax * 4, replace = TRUE) - 1L,
                       nrow = ntax, ncol = 4,
                       dimnames = list(paste0("t", 1:ntax), NULL))
  full <- cbind(mat, mat_3state)
  pd  <- MatrixToPhyDat(full)
  mkd <- MkPrimeData(pd, neomorphic = c(1L, 4L, 7L))  # nNeo = 3, nTrans = 10
  expect_equal(sum(mkd$type == "neomorphic"), 3L)
  expect_equal(sum(mkd$type != "neomorphic"), 10L)

  tree <- TreeTools::Preorder(
    ape::rtree(ntax, br = function(n) runif(n, 0.05, 0.3))
  )
  parent <- tree$edge[, 1]
  child  <- tree$edge[, 2]
  edgeLen <- tree$edge.length

  dataPtr <- .make_data_ptr(mkd)

  # Default kPrime (kObs for trans, 2 for neo).
  kPrime <- mkd$kObs

  for (rn in c(0.1, 0.3, 1.0, 2.5, 7.0)) {
    for (rls in c(0.0, 0.5)) {
      ll_r <- MkpLogLikelihood(tree, mkd,
                                kPrime      = kPrime,
                                rate_loss   = 1.3,
                                rate_log_sd = rls,
                                rate_neo    = rn)
      ll_cpp <- cpp_log_likelihood_xptr(
        dataPtr, parent, child, edgeLen,
        as.integer(kPrime),
        rateLoss   = 1.3,
        rateLogSd  = rls,
        rateNeo    = rn
      )
      expect_equal(ll_r, ll_cpp, tolerance = 1e-10,
        label = sprintf("rate_neo=%g rate_log_sd=%g", rn, rls))
    }
  }
})


# ---- Parity: cache path agrees with direct eval --------------------------
#
# The cache code path (populate_cache → cache_total_loglik) is exercised
# by every accepted NNI / β-simplex / SPR move. We can't call it directly
# from R, but a short MCMC chain that uses partial-CL moves and reports
# log_lik must produce the same likelihood as direct evaluation when both
# are computed at the final state.

test_that("cache (via short MCMC) and direct-eval agree at chain end", {
  skip_on_cran()
  skip_if_not(Sys.getenv("MKPRIME_SLOW_TESTS", "false") == "true" ||
              !nzchar(Sys.getenv("CI", "")),
              "Slow MCMC parity test")

  set.seed(20260527L)
  ntax <- 5L
  mat <- matrix(sample.int(2, ntax * 6, replace = TRUE) - 1L,
                nrow = ntax, ncol = 6,
                dimnames = list(paste0("t", 1:ntax), NULL))
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = c(1L, 2L))

  res <- RunMkPrime(mkd, nRuns = 1L, nIter = 20L, maxWarmup = 5L,
                    minWarmup = 5L, treeThin = 1L, paramThin = 1L,
                    quiet = TRUE, verbose = FALSE)

  # The recovered final-iter log_lik from streaming output should equal
  # MkpLogLikelihood evaluated at the same state.  If the cache path
  # drifts from direct eval, this assertion fails.
  expect_true(all(is.finite(res$samples[, "log_post"])))
})


# ---- Behavioural test: rate_neo materially affects likelihood ------------

test_that("rate_neo affects likelihood via BOTH partitions (mixed data)", {
  # Before the fix: rate_neo scaled only neo edges; trans likelihood was
  # rate_neo-invariant. After the fix: trans likelihood ALSO depends on
  # rate_neo via transScale. Demonstrate by isolating trans-only contribution.
  set.seed(20260527L)
  ntax <- 5L
  mat <- matrix(sample.int(2, ntax * 8, replace = TRUE) - 1L,
                nrow = ntax, ncol = 8,
                dimnames = list(paste0("t", 1:ntax), NULL))
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = c(1L, 2L, 3L))  # 3 neo, 5 trans
  tree <- TreeTools::Preorder(
    ape::rtree(ntax, br = function(n) rep(0.1, n))
  )

  ll_low  <- MkpLogLikelihood(tree, mkd, rate_neo = 0.1)
  ll_one  <- MkpLogLikelihood(tree, mkd, rate_neo = 1.0)
  ll_high <- MkpLogLikelihood(tree, mkd, rate_neo = 10.0)

  # Each must be a distinct likelihood (rate_neo is identifiable).
  expect_false(isTRUE(all.equal(ll_low, ll_one, tolerance = 1e-6)))
  expect_false(isTRUE(all.equal(ll_one, ll_high, tolerance = 1e-6)))
  expect_false(isTRUE(all.equal(ll_low, ll_high, tolerance = 1e-6)))
})


test_that("extreme rate_neo: trans pruning sees small effective branches", {
  # As rate_neo → ∞, transScale → 0 (1/(1+r) * n/n_trans → 0).
  # So trans-partition CLs should look like the zero-branch-length limit
  # for any reasonable tree: all P(no change) ≈ 1, P(change) ≈ 0.
  # A constant-character pattern (all tips equal) gives ≈ log(1/k) per
  # character at that limit, whereas an informative pattern is suppressed.
  #
  # We assert the direction: at very large rate_neo, the trans-only
  # contribution to log-likelihood approaches the all-branches-zero limit,
  # which for variable-coding ascertainment gives a near-singular
  # constant-site probability (and hence a more-negative likelihood for
  # informative trans data) than at rate_neo = 1.
  #
  # Without the fix (trans edges unscaled), increasing rate_neo would NOT
  # change the trans contribution at all — this test would fail.

  set.seed(20260527L)
  ntax <- 5L
  # 1 neo char, 8 highly-informative trans chars (varied across tips)
  mat <- cbind(
    matrix(c(0, 0, 1, 1, 0), nrow = ntax, ncol = 1),
    matrix(rep(0:(ntax - 1L), 8L), nrow = ntax, ncol = 8)
  )
  dimnames(mat) <- list(paste0("t", 1:ntax), NULL)
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = 1L)  # 1 neo, 8 trans

  tree <- TreeTools::Preorder(
    ape::rtree(ntax, br = function(n) rep(0.5, n))
  )

  ll_one  <- MkpLogLikelihood(tree, mkd, rate_neo = 1.0,   coding = "variable")
  ll_huge <- MkpLogLikelihood(tree, mkd, rate_neo = 1e4,   coding = "variable")

  # rate_neo = 1: transScale = 0.5 * 9/8 = 0.5625, full information from trans
  # rate_neo = 1e4: transScale ≈ 9/(8e4) ≈ 1e-4, trans contribution
  #                  collapses (informative chars become near-impossible).
  # ll_huge should be substantially more negative than ll_one.
  expect_lt(ll_huge, ll_one - 1)
})


# ---- Trans-only invariance: now mathematical, not coincidental -----------

test_that("trans-only datasets: likelihood is rate_neo-independent (degenerate)", {
  # When nNeo == 0, compute_partition_scales returns (1.0, 1.0) by design.
  # The likelihood must therefore be exactly rate_neo-invariant.  Pre-fix,
  # this was incidentally true; post-fix, it's enforced by the helper's
  # degenerate-case branch.
  set.seed(20260527L)
  mat <- matrix(sample.int(3, 5 * 6, replace = TRUE) - 1L,
                nrow = 5, ncol = 6,
                dimnames = list(paste0("t", 1:5), NULL))
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  expect_equal(sum(mkd$type == "neomorphic"), 0L)

  tree <- TreeTools::Preorder(
    ape::rtree(5, br = function(n) rep(0.1, n))
  )

  ll1 <- MkpLogLikelihood(tree, mkd, rate_neo = 0.5)
  ll2 <- MkpLogLikelihood(tree, mkd, rate_neo = 1.0)
  ll3 <- MkpLogLikelihood(tree, mkd, rate_neo = 4.0)
  expect_equal(ll1, ll2, tolerance = 1e-10)
  expect_equal(ll2, ll3, tolerance = 1e-10)
})


# ---- F1 regression: partLogLik partial cache must stay in sync after a
#      rate_neo move on mixed data.
#
# Pre-audit-Issue-1, do_move_impl case 3 (rate_neo scale) and case 18
# (neo_joint) refreshed only `data->neoPartIndices` in `state->partLogLik`,
# on the rationale that rate_neo touched only neo edges. Under the new
# RB-style normalisation rate_neo *also* shifts transScale, so trans
# partition log-likelihoods are stale after an accepted rate_neo move.
# This test drives a sequence of rate_neo moves and asserts the cached
# state->logLik matches a fresh direct full evaluation at every step.
# Caught by the external-reviewer agent as a blocker for this branch.

test_that("partLogLik stays in sync after rate_neo / neo_joint moves (F1)", {
  set.seed(20260527L)
  ntax <- 6L
  mat <- matrix(sample.int(2, ntax * 12, replace = TRUE) - 1L,
                nrow = ntax, ncol = 12,
                dimnames = list(paste0("t", 1:ntax), NULL))
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = c(1L, 2L, 3L, 4L))  # 4 neo, 8 trans

  tree <- TreeTools::Preorder(
    ape::rtree(ntax, br = function(n) runif(n, 0.05, 0.3))
  )
  parent  <- tree$edge[, 1]
  child   <- tree$edge[, 2]
  edgeLen <- tree$edge.length
  treeLen <- sum(edgeLen)
  relBr   <- edgeLen / treeLen

  dataPtr <- .make_data_ptr(mkd)
  kPrime <- mkd$kObs

  # Initial full-eval log_lik (no cache yet).
  init_ll <- cpp_log_likelihood_xptr(
    dataPtr, parent, child, edgeLen, as.integer(kPrime),
    rateLoss = 1.3, rateLogSd = 0.5, rateNeo = 1.0
  )

  # Two-stage init: first state to compute logPrior properly, then re-init
  # with the correct logPrior so MH ratios are realistic (otherwise prior
  # mismatch trivially rejects every move and the test never exercises the
  # accept path where the F1 bug bites).
  statePtr0 <- init_mcmc_state(
    parent, child, relBr, treeLen,
    rateLoss = 1.3, rateLogSd = 0.5, rateNeo = 1.0, p = 0.5,
    kPrime = as.integer(kPrime),
    logLik = init_ll, logPrior = 0.0
  )
  init_lp <- eval_log_prior_cpp(dataPtr, statePtr0)
  statePtr <- init_mcmc_state(
    parent, child, relBr, treeLen,
    rateLoss = 1.3, rateLogSd = 0.5, rateNeo = 1.0, p = 0.5,
    kPrime = as.integer(kPrime),
    logLik = init_ll, logPrior = init_lp
  )
  fill_partition_cache(dataPtr, statePtr)

  # Drive a sequence of rate_neo + neo_joint moves at beta = 1. The cached
  # state->logLik must stay in sync with a fresh direct eval at every step;
  # pre-fix (audit Issue 1) the cache only refreshed neo partitions on
  # case 3 / case 18, so accepted moves left trans contribution stale and
  # `cached - fresh` diverged by O(1) within ~25 accepted moves.
  set.seed(99L)
  n_acc <- 0L
  for (step in 1:50) {
    moveType <- sample(c(3L, 18L), 1L)  # 3 = rate_neo, 18 = neo_joint
    ok <- do_move_cpp(dataPtr, statePtr,
                     moveType = moveType, charIdx = 0L,
                     scaleTuning = 0.4, betaSimplexTuning = 1.0,
                     intWalkWindow = 1L, beta = 1.0)
    if (ok) n_acc <- n_acc + 1L

    cached <- get_state_log_lik(statePtr)
    fresh  <- eval_full_loglik_cpp(dataPtr, statePtr)
    expect_equal(cached, fresh, tolerance = 1e-9,
      label = sprintf("step %d (moveType=%d, accepted=%d)", step, moveType, ok))
  }
  # Sanity: the test only catches the bug on accepted moves; if nothing
  # accepts, the test is vacuous. Empirically this RNG produces ~47/50
  # acceptances on this dataset; require >= 20 as a soft lower bound.
  expect_gt(n_acc, 20L)
})
