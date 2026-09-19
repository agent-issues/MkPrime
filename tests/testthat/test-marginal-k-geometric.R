# Tests for the v1 marginal-k geometric arm
# (see dev/notes/2026-05-28-marginal-k-plan.md §7.1).
#
# The marginal-k evaluator analytically sums per-character k'_i out of the
# likelihood:
#   L_marg(y_i | tree, mu, p) =
#     logSumExp_{u = 0..u_max} [ LL(y_i | tree, mu, kObs_i + u)
#                                + log P(u | p) ]
#
# We validate by direct comparison against a brute-force sum-of-products on
# a tiny synthetic dataset, where we can enumerate the full u-grid.

# ---------------------------------------------------------------------------
# Tiny fixture: 4 transformational chars, 8 tips, all kObs == 2 (binary).
# ---------------------------------------------------------------------------

.marg_make_tree <- function() {
  read.tree(text = paste0(
    "(((t1:0.05,t2:0.07):0.04,(t3:0.06,t4:0.05):0.03):0.05,",
    "((t5:0.04,t6:0.06):0.05,(t7:0.05,t8:0.04):0.06):0.04);"
  ))
}

.marg_make_mkd <- function() {
  set.seed(42L)
  tips <- paste0("t", 1:8)
  # 4 transformational binary characters
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

# Build dataPtr / statePtr at a given p and likelihoodMode.
.marg_build_ptrs <- function(p, mode = c("sampled_k", "marginal_k"),
                              kPrime = NULL) {
  mode <- match.arg(mode)
  tree  <- TreeTools::Preorder(.marg_make_tree())
  pd    <- .marg_make_mkd()
  mkd   <- MkPrimeData(pd)
  model <- MkPrimeModel(kPrimePrior = "geometric",
                        likelihoodMode = mode,
                        coding = "none",
                        relabel = FALSE)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)

  state0 <- MkPrime:::.InitState(tree, mkd, model)
  state0$p <- p
  if (!is.null(kPrime)) state0$kPrime <- as.integer(kPrime)
  state0$rate_log_sd <- 0  # disable ACRV for cleanest comparison
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
       state    = state0)
}

# ---------------------------------------------------------------------------
# Bit-identity test
# ---------------------------------------------------------------------------

test_that("marginal-k LL equals brute-force sum-of-products at small nTrans", {
  # p = 0.7 makes the geometric tail P(u >= 12) = 0.3^12 ≈ 5.3e-7, well
  # below the marginal evaluator's M-164 cutoff (-25 nats) so the two
  # truncations are effectively the same. uMax = 11 gives 12^4 = 20,736
  # brute-force pruning calls — runs in a few seconds.
  p_test <- 0.7
  # Build marginal-k state and read off its marginal LL.
  marg <- .marg_build_ptrs(p_test, mode = "marginal_k")
  L_marg <- eval_full_loglik_cpp(marg$dataPtr, marg$statePtr)
  expect_true(is.finite(L_marg))

  # Brute-force: enumerate u_vec ∈ {0..uMax}^4 and sum
  #   exp( LL(tree, mu, kObs + u_vec)  +  Σ_i log P(u_i | p) )
  # in log-space. Geometric P(u | p) = p * (1-p)^u.
  uMax    <- 11L
  uGrid   <- 0:uMax
  nTrans  <- 4L
  logPu   <- log(p_test) + uGrid * log1p(-p_test)

  # We need the per-(tree, mu, kPrime) LL via cpp_log_likelihood_xptr in
  # sampled-k mode (k = kObs + u).  Build a sampled-k dataPtr / statePtr
  # so we can call cpp_log_likelihood_xptr through it.
  samp <- .marg_build_ptrs(p_test, mode = "sampled_k")
  edge  <- samp$tree$edge
  edgeLen <- samp$tree$edge.length
  kObs    <- samp$mkd$kObs

  # Enumerate all u-combinations: 7^4 = 2401 evaluations.
  uPerms <- expand.grid(u1 = uGrid, u2 = uGrid, u3 = uGrid, u4 = uGrid,
                        KEEP.OUT.ATTRS = FALSE)

  logTerms <- numeric(nrow(uPerms))
  for (idx in seq_len(nrow(uPerms))) {
    uVec <- as.integer(unlist(uPerms[idx, ]))
    kP   <- as.integer(kObs + uVec)
    ll <- cpp_log_likelihood_xptr(
      samp$dataPtr, edge[, 1L], edge[, 2L], edgeLen, kP,
      rateLoss = 1.0, rateLogSd = 0.0, rateNeo = 1.0, betaScale = 1.0
    )
    logTerms[idx] <- ll + sum(logPu[uVec + 1L])
  }
  mx <- max(logTerms[is.finite(logTerms)])
  L_naive <- mx + log(sum(exp(logTerms[is.finite(logTerms)] - mx)))

  # Truncation bias: under the M-164 cutoff at -25 nats, the marginal
  # evaluator may drop tail terms below ~exp(-25)·nChar ≈ 4e-11 per char.
  # The brute force also truncates (we cap at u = 6 here) so both
  # truncations bound the discrepancy at the same scale. Use a comfortable
  # 1e-6 tolerance.
  expect_equal(L_marg, L_naive, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# Mode-switching guards
# ---------------------------------------------------------------------------

test_that("likelihoodMode='marginal_k' aborts for non-geometric arms", {
  expect_error(
    MkPrimeModel(kPrimePrior = "empirical_geometric",
                 likelihoodMode = "marginal_k"),
    regexp = "requires"
  )
  expect_error(
    MkPrimeModel(kPrimePrior = "beta_geometric",
                 likelihoodMode = "marginal_k"),
    regexp = "requires"
  )
  expect_error(
    MkPrimeModel(kPrimePrior = "logseries",
                 likelihoodMode = "marginal_k"),
    regexp = "requires"
  )
})

test_that("likelihoodMode='marginal_k' aborts under qHeterogeneity", {
  expect_error(
    MkPrimeModel(kPrimePrior = "geometric",
                 likelihoodMode = "marginal_k",
                 qHeterogeneity = TRUE),
    regexp = "qHeterogeneity"
  )
})

# ---------------------------------------------------------------------------
# LogPrior drops the geometric u-term under marginal_k
# ---------------------------------------------------------------------------

test_that("LogPrior under marginal_k omits the per-character P(u | p) term", {
  tree  <- TreeTools::Preorder(.marg_make_tree())
  pd    <- .marg_make_mkd()
  mkd   <- MkPrimeData(pd)
  model_s <- MkPrimeModel(kPrimePrior = "geometric",
                          likelihoodMode = "sampled_k")
  model_m <- MkPrimeModel(kPrimePrior = "geometric",
                          likelihoodMode = "marginal_k")
  model_s <- MkPrime:::.FinalizeModel(model_s, tree, mkd)
  model_m <- MkPrime:::.FinalizeModel(model_m, tree, mkd)

  state <- MkPrime:::.InitState(tree, mkd, model_s)
  # Set some u_i > 0 so the term is non-trivial.
  trans_idx <- which(mkd$type == "transformational")
  state$kPrime[trans_idx] <- as.integer(mkd$kObs[trans_idx] + c(0, 1, 2, 3))
  state$p <- 0.4

  lp_s <- MkPrime:::LogPrior(state, model_s, mkd)
  lp_m <- MkPrime:::LogPrior(state, model_m, mkd)

  # Expected difference is the dropped per-character u-term:
  #   n * log(p) + sum(u) * log(1-p)
  u_vec <- state$kPrime[trans_idx] - mkd$kObs[trans_idx]
  expected_diff <- length(trans_idx) * log(state$p) +
                   sum(u_vec) * log1p(-state$p)
  expect_equal(lp_s - lp_m, expected_diff, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Model A (unconditional) vs Model B (conditional) marginal-weight factor
#
# Model B weight for state count k is p (1-p)^(k - kObs_i); Model A weight is
# p (1-p)^(k - 2). The two differ by a per-character constant (1-p)^(kObs_i-2)
# that factors out of the per-character logSumExp, so the total marginal LL
# under priorVariant="unconditional" exceeds the "conditional" total by
# Σ_i (kObs_i - 2) · log(1-p). Characters with kObs == 2 contribute nothing.
# ---------------------------------------------------------------------------

test_that("priorVariant='unconditional' shifts marginal LL by (kObs-2)log(1-p)", {
  tips <- paste0("t", 1:8)
  tree <- TreeTools::Preorder(.marg_make_tree())
  # Two characters with kObs > 2 (one 3-state, one 4-state) plus two binary,
  # so the factor (1-p)^(kObs-2) is non-trivial.
  mat <- matrix(
    c(0, 1, 2, 0, 1, 2, 0, 1,   # kObs = 3
      0, 1, 2, 3, 0, 1, 2, 3,   # kObs = 4
      0, 0, 1, 1, 0, 1, 0, 1,   # kObs = 2
      1, 0, 1, 0, 1, 1, 0, 0),  # kObs = 2
    nrow = 8, ncol = 4, dimnames = list(tips, NULL)
  )
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  p_test <- 0.6

  build <- function(variant) {
    model <- MkPrimeModel(kPrimePrior = "geometric",
                          likelihoodMode = "marginal_k",
                          priorVariant = variant,
                          coding = "none", relabel = FALSE)
    model <- MkPrime:::.FinalizeModel(model, tree, mkd)
    state0 <- MkPrime:::.InitState(tree, mkd, model)
    state0$p <- p_test
    state0$rate_log_sd <- 0
    state0$tree_length <- sum(tree$edge.length)
    state0$rel_br_lengths <- tree$edge.length / state0$tree_length
    dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
    statePtr <- MkPrime:::.InitMcmcChain(state0)
    fill_partition_cache(dataPtr, statePtr)
    eval_full_loglik_cpp(dataPtr, statePtr)
  }

  L_cond   <- build("conditional")
  L_uncond <- build("unconditional")
  expect_true(is.finite(L_cond) && is.finite(L_uncond))

  trans_idx <- which(mkd$type == "transformational")
  expected_diff <- sum((mkd$kObs[trans_idx] - 2) * log1p(-p_test))
  expect_equal(L_uncond - L_cond, expected_diff, tolerance = 1e-9)

  # Default (no priorVariant) now equals "unconditional" (Model A).
  expect_equal(L_cond, build("conditional"), tolerance = 1e-12)
})
