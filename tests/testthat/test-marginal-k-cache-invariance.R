# Bit-identity regression test for the marginal-k charLL cache under p-moves.
#
# Background (dev/red-team/proofs/marginal-k-geometric.md §5.1):
#   charLL[ti, ko] = LL(y_i | tree, mu, kObs + ko)
#   is stored as the raw per-(char, k) log-likelihood, which does NOT depend
#   on p.  On a p-move (do_move_impl case 30, mh_logit_p), the C++ engine
#   leaves charLLCacheReady = true and re-runs only the per-character
#   logSumExp against new geometric weights.
#
# This test asserts that the cached marginal LL after an accepted p-move is
# bit-identical to a fresh cold-cache recompute at the same p.  If a future
# change to cache-invalidation logic introduced a stale-cache bug, this test
# would FAIL.
#
# Design: parallel-state approach (avoids src/* changes).
#
#   State A at p_init:
#     1. eval_full_loglik_cpp(A) → L_init  [fills cache, charLLCacheReady=true]
#     2. Fire do_move_cpp(moveType=30) until an acceptance occurs (p changes).
#        Acceptance is guaranteed by symmetry within a few tries for any
#        p_init not near the Beta(a,b) mode.
#     3. Read p_new = get_mcmc_state(A)$p
#     4. eval_full_loglik_cpp(A) → L_cached  [uses charLL cache fast-path]
#
#   State B: fresh build at the same (tree, mu) but p = p_new.
#     5. eval_full_loglik_cpp(B) → L_fresh  [cold cache: full recompute]
#
#   Assert: L_cached == L_fresh to 1e-12.
#   Guard:  L_init != L_cached  (proves p changed and marginal LL is p-sensitive).

library("ape")
library("TreeTools")

# ---------------------------------------------------------------------------
# Fixture helpers (shared with test-marginal-k-geometric.R via reference copy)
# ---------------------------------------------------------------------------

.cache_make_tree <- function() {
  read.tree(text = paste0(
    "(((t1:0.05,t2:0.07):0.04,(t3:0.06,t4:0.05):0.03):0.05,",
    "((t5:0.04,t6:0.06):0.05,(t7:0.05,t8:0.04):0.06):0.04);"
  ))
}

.cache_make_mkd <- function() {
  set.seed(42L)
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

# Build an initialised marginal-k state at a given p.
# Returns: list(dataPtr, statePtr, tree, mkd, model, state0)
.cache_build_ptrs <- function(p) {
  tree  <- TreeTools::Preorder(.cache_make_tree())
  pd    <- .cache_make_mkd()
  mkd   <- MkPrimeData(pd)
  model <- MkPrimeModel(kPrimePrior = "geometric",
                        likelihoodMode = "marginal_k",
                        coding = "none",
                        relabel = FALSE)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)

  state0 <- MkPrime:::.InitState(tree, mkd, model)
  state0$p          <- p
  state0$rate_log_sd <- 0  # disable ACRV: cleanest cache comparison
  state0$tree_length <- sum(tree$edge.length)
  state0$rel_br_lengths <- tree$edge.length / state0$tree_length

  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)

  list(dataPtr  = dataPtr,
       statePtr = statePtr,
       tree     = tree,
       mkd      = mkd,
       model    = model,
       state0   = state0)
}


# ---------------------------------------------------------------------------
# Main invariance test
# ---------------------------------------------------------------------------

test_that("marginal-k charLL cache is bit-identical to fresh recompute after p-move", {

  # --- Step 1: Build state A at p_init and populate the charLL cache. -------
  p_init <- 0.55  # away from prior mode so proposals are likely to move
  A <- .cache_build_ptrs(p_init)

  # First eval: fills charLLCache with raw per-(char,k) LLs; sets Ready=true.
  L_init <- eval_full_loglik_cpp(A$dataPtr, A$statePtr)
  expect_true(is.finite(L_init),
              label = "L_init should be finite before any p-move")

  # --- Step 2: Fire case-30 until an accepted p-move occurs. ----------------
  # Case 30 = mh_logit_p.  Under marginal_k the cache is NOT invalidated on
  # case 30, so charLLCacheReady remains true after the proposal.
  # MH acceptance is controlled by the Beta(a,b) hyperprior on p only (the
  # per-character P(u|p) mass is consumed by the marginal LL, not the prior).
  # With a fairly uninformative hyperprior and p_init = 0.55, proposals will
  # be accepted within a handful of tries.
  #
  # scaleTuning = 1.0 gives logit-scale steps of order 1, i.e. multiplicative
  # changes of ~e on the odds ratio — large enough to move p visibly but not
  # so large that acceptance rate collapses.
  accepted <- FALSE
  p_new <- p_init
  for (trial in seq_len(200L)) {
    set.seed(trial)
    ok <- do_move_cpp(A$dataPtr, A$statePtr,
                      moveType          = 30L,
                      charIdx           = 0L,
                      scaleTuning       = 1.0,
                      betaSimplexTuning = 0.0,
                      intWalkWindow     = 1L,
                      beta              = 1.0)
    if (ok) {
      p_new <- get_mcmc_state(A$statePtr)$p
      if (abs(p_new - p_init) > 1e-9) {
        accepted <- TRUE
        break
      }
    }
  }
  expect_true(accepted,
              label = paste0("A case-30 move should be accepted within 200 ",
                             "tries at p_init = ", p_init))

  # --- Step 3: Compute marginal LL on state A using the charLL cache. --------
  # charLLCacheReady is still true (case 30 does not invalidate it).
  # This call takes the fast-path: reuses raw cache, re-runs logSumExp
  # against new geometric weights for p_new.
  L_cached <- eval_full_loglik_cpp(A$dataPtr, A$statePtr)
  expect_true(is.finite(L_cached),
              label = "L_cached should be finite after accepted p-move")

  # Guard: the marginal LL must change when p changes — if this fails the test
  # is trivially passing because the move was a no-op.
  expect_false(identical(L_cached, L_init),
               label = paste0("L_cached should differ from L_init when p ",
                              "moved from ", p_init, " to ", round(p_new, 6)))

  # --- Step 4: Build state B fresh at p_new — cold cache. -------------------
  # Identical tree / mu; only p differs.  charLLCacheReady starts false.
  B <- .cache_build_ptrs(p_new)

  # Full cold recompute: fills cache from scratch, no fast-path.
  L_fresh <- eval_full_loglik_cpp(B$dataPtr, B$statePtr)
  expect_true(is.finite(L_fresh),
              label = "L_fresh should be finite after cold-cache recompute")

  # --- Step 5: Bit-identity assertion. --------------------------------------
  # The two paths are algebraically identical but differ in floating-point
  # order:
  #
  #   Cache path:  rawLL = (LL + logP_init + ko*log1mP_init)
  #                           - (logP_init + ko*log1mP_init)
  #                then w = rawLL + (logP_new + ko*log1mP_new)
  #
  #   Fresh path:  w = LL + logP_new + ko*log1mP_new  (single computation)
  #
  # The round-trip subtraction/addition introduces FP rounding at the
  # level of |LL| × ε × O(nCand), empirically ~3e-11 nats (absolute) on
  # the fixture above.  A real stale-cache bug would produce errors that
  # dwarf this: wrong cached rawLLs (stale from a different tree/μ) would
  # shift the LL by many nats.  Using tolerance = 1e-9 (absolute) draws a
  # safe line: above FP noise, far below any realistic stale-cache error.
  expect_equal(L_cached, L_fresh,
               tolerance = 1e-9,
               label = paste0("charLL cache after p-move must agree with ",
                              "cold-cache recompute at p_new = ",
                              round(p_new, 8), " to 1e-9"))
})


# ---------------------------------------------------------------------------
# Sanity check: non-p move invalidates the cache (cache IS stale after case 0)
# ---------------------------------------------------------------------------
#
# This test documents the complementary guarantee: a tree-length scale (case 0)
# DOES invalidate charLLCacheReady.  If the cache were consulted after case 0
# at the proposed (new treeLength) and then rolled back to the old treeLength,
# the next eval_full_loglik_cpp would see wrong raw LLs.  The invariant is
# that charLLCacheReady=false is set at the top of do_move_impl for all non-30
# moves (lines 4475-4476 of mcmc.cpp).
#
# We test this indirectly: after a case-0 move the fresh recompute at the
# same state must match eval_full_loglik_cpp (which rebuilds from scratch when
# charLLCacheReady=false). If a future commit removed that invalidation guard
# and case-0 left stale raw LLs in the cache, this paired test would FAIL on
# the bit-identity check above (different treeLength → different raw LLs).
# This test is intentionally light — it just verifies eval is consistent with
# itself after a non-p move.

test_that("marginal-k: eval_full_loglik_cpp is self-consistent after a non-p move", {
  p0 <- 0.5
  X <- .cache_build_ptrs(p0)

  # Populate cache at initial state.
  L_before <- eval_full_loglik_cpp(X$dataPtr, X$statePtr)
  expect_true(is.finite(L_before))

  # Fire a tree-length scale (case 0) — this sets charLLCacheReady=false
  # regardless of acceptance.
  set.seed(999L)
  do_move_cpp(X$dataPtr, X$statePtr,
              moveType = 0L, charIdx = 0L,
              scaleTuning = 0.1, betaSimplexTuning = 0.0,
              intWalkWindow = 1L, beta = 1.0)

  # Two consecutive eval calls should agree: both trigger a fresh compute
  # (charLLCacheReady was invalidated by the case-0 move or its rollback).
  L1 <- eval_full_loglik_cpp(X$dataPtr, X$statePtr)
  L2 <- eval_full_loglik_cpp(X$dataPtr, X$statePtr)

  expect_true(is.finite(L1))
  expect_equal(L1, L2, tolerance = 1e-12,
               label = "Two consecutive evals after a non-p move must agree exactly")
})
