# Regression test for FU-5: p-moves change the marginal LL under marginal-k.
#
# Background
# ----------
# Under sampled-k semantics, p enters only the prior on kPrime, and per-
# character LL = L(y_i | tree, mu, kPrime[i]). A p-move therefore leaves
# state->logLik untouched, and the historical code at src/mcmc.cpp:5142
# carved out cases 8 / 30 from `likChanges`.
#
# Under MARGINAL-K, kPrime is marginalised out and p enters the per-
# character marginal LL via the geometric weights P(u | p) = p (1 - p)^u
# consumed inside the per-character logSumExp. A p-move therefore DOES
# change state->logLik even though it changes no other state.
#
# The buggy `likChanges = (moveType != 8 && moveType != 30)` skipped the
# recompute in the MH acceptance ratio, collapsing logAlpha to
#   prior_ratio + logHastings (Jacobian)
# Under the default Beta(1, 1) hyperprior the prior ratio is zero, so the
# chain random-walked on p ignoring the data. This test guards against a
# regression of that line.
#
# Design
# ------
# Single-step bit-identity check. Build a marginal-k state at known
# (tree, mu, p_init). Fire one accepted case-30 move via do_move_cpp.
# Read state->logLik (must reflect the new p under the fix) and compare
# to an independent fresh recompute via eval_full_loglik_cpp.
#
# Pre-fix outcome: state->logLik == L(p_init) (stale, because case-30
#   takes the `!likChanges` branch and reuses state->logLik in MH). The
#   fresh recompute at p_new differs by O(|nTrans| * |Delta log p|), of
#   order 1 nat for the fixture below. The expect_equal would FAIL.
#
# Post-fix outcome: state->logLik == L(p_new). They agree to FP noise.
#
# We also verify the move was actually accepted (p changed) and that
# state->logLik moved by a non-trivial amount, so the test is not
# trivially passing on a no-op move.

library("ape")
library("TreeTools")

# ---------------------------------------------------------------------------
# Fixture (mirrors test-marginal-k-cache-invariance.R for consistency)
# ---------------------------------------------------------------------------

.pmove_make_tree <- function() {
  read.tree(text = paste0(
    "(((t1:0.05,t2:0.07):0.04,(t3:0.06,t4:0.05):0.03):0.05,",
    "((t5:0.04,t6:0.06):0.05,(t7:0.05,t8:0.04):0.06):0.04);"
  ))
}

.pmove_make_mkd <- function() {
  tips <- paste0("t", 1:8)
  mat <- matrix(
    c(0, 0, 1, 1, 0, 1, 0, 1,
      0, 1, 1, 0, 1, 0, 0, 1,
      0, 1, 0, 0, 1, 1, 1, 0,
      1, 1, 0, 1, 0, 0, 1, 0,
      0, 1, 1, 1, 0, 0, 1, 0,
      1, 0, 1, 0, 1, 0, 1, 0,
      0, 0, 1, 1, 1, 0, 0, 1,
      1, 1, 0, 0, 0, 1, 1, 0),
    nrow = 8, ncol = 8,
    dimnames = list(tips, NULL)
  )
  MatrixToPhyDat(mat)
}

.pmove_build_ptrs <- function(p) {
  tree  <- TreeTools::Preorder(.pmove_make_tree())
  pd    <- .pmove_make_mkd()
  mkd   <- MkPrimeData(pd)
  model <- MkPrimeModel(kPrimePrior = "geometric",
                        likelihoodMode = "marginal_k",
                        coding = "none",
                        relabel = FALSE)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)

  state0 <- MkPrime:::.InitState(tree, mkd, model)
  state0$p          <- p
  state0$rate_log_sd <- 0  # disable ACRV
  state0$tree_length <- sum(tree$edge.length)
  state0$rel_br_lengths <- tree$edge.length / state0$tree_length

  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)

  list(dataPtr = dataPtr, statePtr = statePtr)
}

# ---------------------------------------------------------------------------
# FU-5 fix guard: state->logLik tracks the marginal LL after a p-move
# ---------------------------------------------------------------------------

test_that("marginal-k: case-30 p-move writes the new marginal LL to state->logLik", {

  # Direct bit-identity check: after an accepted case-30 move under
  # marginal-k, state->logLik must equal a fresh independent recompute
  # (eval_full_loglik_cpp) at the new p.
  #
  # Pre-fix:  case 30 takes the `!likChanges` branch, sets
  #   newLogLik = state->logLik, and state->logLik is therefore
  #   unchanged on acceptance. It remains stuck at its prior value
  #   (the sampled-k init from .MkpLogLikelihood, or the post-move
  #   value from a previous accepted likelihood-changing move).
  #   The fresh marginal recompute at the new p disagrees by O(few nats).
  #
  # Post-fix: do_move_cpp routes through compute_full_loglik_at under
  #   marginal-k and writes newLogLik (= L_marg(new p)) to state->logLik
  #   on acceptance. The two values agree to FP noise.

  p_init <- 0.55
  X <- .pmove_build_ptrs(p_init)

  # Fire case-30 until acceptance with a noticeable p shift.
  accepted <- FALSE
  p_new <- p_init
  for (trial in seq_len(500L)) {
    set.seed(trial)
    ok <- do_move_cpp(X$dataPtr, X$statePtr,
                      moveType          = 30L,
                      charIdx           = 0L,
                      scaleTuning       = 1.0,
                      betaSimplexTuning = 0.0,
                      intWalkWindow     = 1L,
                      beta              = 1.0)
    if (ok) {
      p_new <- get_mcmc_state(X$statePtr)$p
      if (abs(p_new - p_init) > 0.02) {
        accepted <- TRUE
        break
      }
    }
  }
  expect_true(accepted,
              label = "A case-30 move with |delta p| > 0.02 should be accepted within 500 tries")

  L_stored <- get_mcmc_state(X$statePtr)$logLik
  L_fresh  <- eval_full_loglik_cpp(X$dataPtr, X$statePtr)
  expect_true(is.finite(L_stored))
  expect_true(is.finite(L_fresh))

  # Sanity: a fresh recompute at p_init differs from the recompute at
  # p_new (proves marginal LL depends on p — guards a vacuous pass).
  Y <- .pmove_build_ptrs(p_init)
  L_fresh_p_init <- eval_full_loglik_cpp(Y$dataPtr, Y$statePtr)
  expect_false(isTRUE(all.equal(L_fresh, L_fresh_p_init, tolerance = 1e-6)),
               label = paste0("Marginal LL should depend on p — ",
                              "L_fresh(p_new = ", round(p_new, 4),
                              ") = ", round(L_fresh, 4),
                              " vs L_fresh(p_init) = ",
                              round(L_fresh_p_init, 4)))

  # CORE REGRESSION ASSERTION.
  expect_equal(L_stored, L_fresh, tolerance = 1e-9,
               label = paste0("After case-30 (p: ", round(p_init, 6),
                              " -> ", round(p_new, 6), "), ",
                              "state->logLik = ", round(L_stored, 6),
                              " must equal fresh marginal recompute = ",
                              round(L_fresh, 6),
                              ". If this fails, src/mcmc.cpp:5142 is excluding ",
                              "case 30 from `likChanges` under marginal-k -- the ",
                              "MH ratio is missing the LL contribution and the ",
                              "chain is random-walking on p."))
})


# ---------------------------------------------------------------------------
# Belt-and-braces: chain-level check that LL varies with p across iterations
# ---------------------------------------------------------------------------
#
# A short marginal-k chain run with only the case-30 move active should
# show LL varying across recorded samples whenever p has changed. Under
# the bug the LL in get_mcmc_state would lag p (only updated on a non-30
# move acceptance, of which there are zero in this scenario), so any
# pair (p_i, LL_i) with differing p_i but identical LL_i is a smoking
# gun. We turn that into a positive assertion: across a handful of
# accepted p-moves, the stored LL must take more than one distinct value.

test_that("marginal-k: stored logLik varies across a series of accepted p-moves", {

  X <- .pmove_build_ptrs(0.5)
  L0 <- eval_full_loglik_cpp(X$dataPtr, X$statePtr)
  expect_true(is.finite(L0))

  ll_samples <- numeric(0)
  p_samples  <- numeric(0)
  for (trial in seq_len(500L)) {
    set.seed(1000L + trial)
    ok <- do_move_cpp(X$dataPtr, X$statePtr,
                      moveType = 30L, charIdx = 0L,
                      scaleTuning = 1.0, betaSimplexTuning = 0.0,
                      intWalkWindow = 1L, beta = 1.0)
    if (ok) {
      st <- get_mcmc_state(X$statePtr)
      ll_samples <- c(ll_samples, st$logLik)
      p_samples  <- c(p_samples, st$p)
      if (length(ll_samples) >= 20L) break
    }
  }

  expect_gte(length(ll_samples), 20L)
  expect_gt(length(unique(p_samples)),  1L)
  expect_gt(length(unique(ll_samples)), 1L)
})
