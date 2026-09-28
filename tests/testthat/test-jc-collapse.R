# Likelihood-identity tests for the state-collapsing optimisation.
#
# When a JC(kFull) character has only kObs < kFull distinct observed states,
# the collapsed pruning kernels evaluate it on kEff = kObs + 1 CL columns and
# must return EXACTLY the same likelihood as the uncollapsed kernel called on
# the same tip data with kStates = kFull (the tips already use states
# 0..kObs-1, so the uncollapsed kernel simply never populates the upper
# columns — both computations describe the same Markov chain).

skip_if_not_installed("ape")
skip_if_not_installed("TreeTools")

mk_test_tree <- function(nTip = 7L, seed = 1L) {
  set.seed(seed)
  tr <- ape::rtree(nTip, br = function(n) runif(n, 0.05, 0.6))
  Preorder(tr)
}

mk_test_tipstates <- function(nTip, kObs, seed = 2L,
                              missing_frac = 0, n_char = 4L) {
  set.seed(seed)
  m <- matrix(sample.int(kObs, nTip * n_char, replace = TRUE) - 1L,
              nrow = nTip, ncol = n_char)
  if (missing_frac > 0) {
    n_na <- ceiling(missing_frac * length(m))
    m[sample.int(length(m), n_na)] <- -1L
  }
  storage.mode(m) <- "integer"
  m
}


test_that("pruning_jc_collapsed matches pruning_jc across kFull >= kObs+1", {
  tr <- mk_test_tree(7L)
  parent <- tr$edge[, 1L]
  child  <- tr$edge[, 2L]
  el     <- tr$edge.length
  nTip   <- length(tr$tip.label)

  for (kObs in c(1L, 2L, 3L, 4L)) {
    tip_states <- mk_test_tipstates(nTip, kObs, seed = 10L + kObs)
    for (kFull in c(kObs + 1L, kObs + 3L, max(kObs + 5L, 13L))) {
      rf_full <- rep(1.0 / kFull, kFull)
      ll_full <- pruning_jc(parent, child, el, tip_states, kFull, rf_full)
      ll_coll <- pruning_jc_collapsed(parent, child, el, tip_states,
                                      kFull, kObs)
      expect_equal(ll_coll, ll_full, tolerance = 1e-10,
                   info = sprintf("kObs=%d, kFull=%d", kObs, kFull))
    }
  }
})


test_that("pruning_jc_collapsed handles missing tip states", {
  tr <- mk_test_tree(6L)
  parent <- tr$edge[, 1L]
  child  <- tr$edge[, 2L]
  el     <- tr$edge.length
  nTip   <- length(tr$tip.label)

  for (kObs in c(2L, 3L)) {
    tip_states <- mk_test_tipstates(nTip, kObs, seed = 20L + kObs,
                                    missing_frac = 0.3)
    for (kFull in c(kObs + 1L, kObs + 5L)) {
      rf_full <- rep(1.0 / kFull, kFull)
      ll_full <- pruning_jc(parent, child, el, tip_states, kFull, rf_full)
      ll_coll <- pruning_jc_collapsed(parent, child, el, tip_states,
                                      kFull, kObs)
      expect_equal(ll_coll, ll_full, tolerance = 1e-10,
                   info = sprintf("kObs=%d kFull=%d (missing)", kObs, kFull))
    }
  }
})


test_that("pruning_jc_acrv_collapsed matches pruning_jc_acrv", {
  tr <- mk_test_tree(6L, seed = 3L)
  parent <- tr$edge[, 1L]
  child  <- tr$edge[, 2L]
  el     <- tr$edge.length
  nTip   <- length(tr$tip.label)

  rates <- DiscreteLognormalRates(0.5, 6L)

  for (kObs in c(2L, 3L)) {
    tip_states <- mk_test_tipstates(nTip, kObs, seed = 30L + kObs)
    for (kFull in c(kObs + 1L, kObs + 4L)) {
      rf_full <- rep(1.0 / kFull, kFull)
      ll_full <- pruning_jc_acrv(parent, child, el, tip_states,
                                 kFull, rf_full, rates)
      ll_coll <- pruning_jc_acrv_collapsed(parent, child, el, tip_states,
                                           kFull, kObs, rates)
      expect_equal(ll_coll, ll_full, tolerance = 1e-10,
                   info = sprintf("ACRV kObs=%d kFull=%d", kObs, kFull))
    }
  }
})


test_that("constant_site_prob_jc_collapsed matches constant_site_prob_jc", {
  tr <- mk_test_tree(5L, seed = 4L)
  parent <- tr$edge[, 1L]
  child  <- tr$edge[, 2L]
  el     <- tr$edge.length
  nTip   <- length(tr$tip.label)

  rates_flat <- 1.0
  rates_acrv <- DiscreteLognormalRates(0.4, 6L)

  for (kObs in c(1L, 2L, 3L)) {
    for (kFull in c(kObs + 1L, kObs + 5L)) {
      rf_full <- rep(1.0 / kFull, kFull)
      for (rates in list(rates_flat, rates_acrv)) {
        p_full <- constant_site_prob_jc(parent, child, el, nTip,
                                        kFull, rf_full, rates)
        p_coll <- constant_site_prob_jc_collapsed(parent, child, el, nTip,
                                                  kFull, kObs, rates)
        expect_equal(p_coll, p_full, tolerance = 1e-10,
                     info = sprintf("const kObs=%d kFull=%d nCat=%d",
                                    kObs, kFull, length(rates)))
      }
    }
  }
})


test_that("singleton_site_prob_jc_collapsed matches singleton_site_prob_jc", {
  tr <- mk_test_tree(5L, seed = 5L)
  parent <- tr$edge[, 1L]
  child  <- tr$edge[, 2L]
  el     <- tr$edge.length
  nTip   <- length(tr$tip.label)

  rates_flat <- 1.0
  rates_acrv <- DiscreteLognormalRates(0.4, 4L)

  for (kObs in c(1L, 2L, 3L)) {
    for (kFull in c(kObs + 1L, kObs + 4L)) {
      rf_full <- rep(1.0 / kFull, kFull)
      for (rates in list(rates_flat, rates_acrv)) {
        p_full <- singleton_site_prob_jc(parent, child, el, nTip,
                                         kFull, rf_full, rates)
        p_coll <- singleton_site_prob_jc_collapsed(parent, child, el, nTip,
                                                   kFull, kObs, rates)
        expect_equal(p_coll, p_full, tolerance = 1e-10,
                     info = sprintf("singleton kObs=%d kFull=%d nCat=%d",
                                    kObs, kFull, length(rates)))
      }
    }
  }
})


test_that("MkpLogLikelihood (auto-dispatched) matches uncollapsed reference", {
  # End-to-end: MkpLogLikelihood now silently routes kPrime > kObs through
  # the collapsed kernels. Compare against the value obtained by calling
  # the uncollapsed C++ pruning directly on the same partition tip states.
  tr <- mk_test_tree(8L, seed = 7L)
  nTip <- length(tr$tip.label)

  # Transformational character with kObs = 3
  mat <- matrix(c(0, 1, 2, 0, 1, 2, 1, 0), nrow = nTip, ncol = 1,
                dimnames = list(tr$tip.label, NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  # Reference: manually call uncollapsed kernels
  manual_ll <- function(kp) {
    tip_states <- mkd$matrix
    tip_states[is.na(tip_states)] <- -1L
    storage.mode(tip_states) <- "integer"
    parent <- tr$edge[, 1L]
    child  <- tr$edge[, 2L]
    el     <- tr$edge.length
    rf <- rep(1.0 / kp, kp)
    ll <- pruning_jc(parent, child, el, tip_states, kp, rf)
    # ascertainment: variable coding
    puninf <- constant_site_prob_jc(parent, child, el, nTip, kp, rf, 1.0)
    ll - log(1 - puninf)
  }

  for (kp in c(3L, 5L, 8L, 13L)) {
    ll_disp <- MkpLogLikelihood(tr, mkd, kPrime = kp,
                                rate_log_sd = 0, coding = "variable",
                                relabel = FALSE)
    expect_equal(ll_disp, manual_ll(kp), tolerance = 1e-10,
                 info = sprintf("auto-dispatched kPrime=%d", kp))
  }
})


test_that("Mk (known partition) mixed-kObs: collapse matches uncollapsed", {
  # Mk use case: fixed k known per character. Within a single known partition
  # different chars may have different kObs values. The batched collapse with
  # kEff = max(kObs) + 1 must still give identical likelihood — by JC
  # lumpability, chars with small kObs simply never light up the upper
  # observed columns, and that is still a valid (finer) state partition.
  tr <- mk_test_tree(7L, seed = 11L)
  nTip <- length(tr$tip.label)
  parent <- tr$edge[, 1L]
  child  <- tr$edge[, 2L]
  el     <- tr$edge.length

  kFull <- 9L

  # Build a tip matrix with 4 chars of mixed kObs: 2, 3, 4, 1
  set.seed(99L)
  mk_col <- function(kObs) {
    sample.int(kObs, nTip, replace = TRUE) - 1L
  }
  mat <- cbind(mk_col(2L), mk_col(3L), mk_col(4L), mk_col(1L))
  storage.mode(mat) <- "integer"
  rownames(mat) <- tr$tip.label

  # Reference: uncollapsed pruning_jc on the same data at kStates = kFull
  rf_full <- rep(1.0 / kFull, kFull)
  ll_full <- pruning_jc(parent, child, el, mat, kFull, rf_full)

  # Batched collapse with kObsMax = 4 (max over the four chars)
  ll_coll <- pruning_jc_collapsed(parent, child, el, mat, kFull, 4L)
  expect_equal(ll_coll, ll_full, tolerance = 1e-10,
               info = "mixed-kObs Mk batched collapse")

  # Same check via constant_site_prob (ascertainment in Mk variable coding)
  p_full <- constant_site_prob_jc(parent, child, el, nTip, kFull, rf_full, 1.0)
  p_coll <- constant_site_prob_jc_collapsed(parent, child, el, nTip,
                                            kFull, 4L, 1.0)
  expect_equal(p_coll, p_full, tolerance = 1e-10,
               info = "mixed-kObs Mk constant-site")
})


test_that("MkpLogLikelihood collapse-path informative coding matches uncollapsed", {
  tr <- mk_test_tree(6L, seed = 9L)
  nTip <- length(tr$tip.label)

  mat <- matrix(c(0, 1, 2, 0, 1, 2), nrow = nTip, ncol = 1,
                dimnames = list(tr$tip.label, NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  manual_ll <- function(kp, rate_sd, coding) {
    tip_states <- mkd$matrix
    tip_states[is.na(tip_states)] <- -1L
    storage.mode(tip_states) <- "integer"
    parent <- tr$edge[, 1L]
    child  <- tr$edge[, 2L]
    el     <- tr$edge.length
    rf <- rep(1.0 / kp, kp)
    rates <- DiscreteLognormalRates(rate_sd, 6L)
    if (rate_sd > 0) {
      ll <- pruning_jc_acrv(parent, child, el, tip_states, kp, rf, rates)
    } else {
      ll <- pruning_jc(parent, child, el, tip_states, kp, rf)
    }
    puninf <- constant_site_prob_jc(parent, child, el, nTip, kp, rf, rates)
    if (coding == "informative") {
      puninf <- puninf +
        uninf_nonconst_prob_jc(parent, child, el, nTip, kp, rates)
    }
    ll - log(1 - puninf)
  }

  for (kp in c(4L, 7L)) {
    for (rs in c(0, 0.5)) {
      for (cd in c("variable", "informative")) {
        ll_disp <- MkpLogLikelihood(tr, mkd, kPrime = kp,
                                    rate_log_sd = rs, coding = cd,
                                    relabel = FALSE)
        expect_equal(ll_disp, manual_ll(kp, rs, cd), tolerance = 1e-10,
                     info = sprintf("kp=%d rs=%g cd=%s", kp, rs, cd))
      }
    }
  }
})
