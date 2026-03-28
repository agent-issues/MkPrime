# Tests for M-052: Q-matrix heterogeneity (Dirichlet-marginal Het)
#
# Tests exercise the Het code path via the C++ MCMC engine, verifying:
#   1. JC recovery: large beta_scale → Het likelihood ≈ non-Het likelihood
#   2. Rotation symmetry: relabelling states produces the same likelihood
#   3. Hand-computed 3-taxon likelihood under Het
#   4. Binary backward compatibility: Het with k=2 matches expected behaviour
#   5. beta_scale prior density correctness

library(ape)
library(TreeTools)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Build C++ pointers for a given tree, phyDat, and model.
# Returns list(dataPtr, statePtr, parent, child, relBr, treeLen).
.het_ptrs <- function(tree, pd, model,
                      neomorphic = integer(0),
                      knownStates = NULL,
                      rate_loss = 1.0,
                      beta_scale = 1.0,
                      rate_log_sd = 0.0) {
  mkd   <- MkPrimeData(pd, neomorphic = neomorphic,
                        knownStates = knownStates)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  tree  <- Preorder(tree)

  state0 <- MkPrime:::.InitState(tree, mkd, model)
  state0$rate_loss   <- rate_loss
  state0$beta_scale  <- beta_scale
  state0$rate_log_sd <- rate_log_sd
  # Recompute log_prior with updated state
  state0$log_prior <- MkPrime:::LogPrior(state0, model, mkd)

  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)

  s <- get_mcmc_state(statePtr)
  list(dataPtr  = dataPtr,
       statePtr = statePtr,
       parent   = s$edge[, 1L],
       child    = s$edge[, 2L],
       relBr    = s$relBrLengths,
       treeLen  = s$treeLength,
       mkd      = mkd,
       model    = model)
}

# Simple 3-taxon unrooted tree: ((t1:b1,t2:b2):0,t3:b3)
# (unrooted, so the internal edge has zero length — star tree for simplicity)
.three_taxon_tree <- function(b1 = 0.1, b2 = 0.2, b3 = 0.15) {
  tr <- read.tree(text = sprintf("((t1:%g,t2:%g):0.0,t3:%g);", b1, b2, b3))
  tr
}

# 4-taxon tree with informative topology
.four_taxon_tree <- function() {
  read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
}


# ===========================================================================
# 1. JC recovery: Het with large alpha ≈ non-Het
# ===========================================================================

test_that("Het with large beta_scale recovers JC likelihood (binary chars)", {
  tree <- .four_taxon_tree()
  mat <- matrix(c(0, 1, 0, 1,
                  1, 0, 1, 0),
                nrow = 4, ncol = 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  # Non-Het baseline
  modelJC <- MkPrimeModel(nCat = 1L, coding = "none")
  ptsJC   <- .het_ptrs(tree, pd, modelJC)
  llJC    <- eval_full_loglik_cpp(ptsJC$dataPtr, ptsJC$statePtr)

  # Het with very large alpha (beta_scale = 1000) → bins concentrate at 1/k
  modelHet <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L,
                            nCat = 1L, coding = "none")
  ptsHet   <- .het_ptrs(tree, pd, modelHet, beta_scale = 1000)
  llHet    <- eval_full_loglik_cpp(ptsHet$dataPtr, ptsHet$statePtr)

  expect_equal(llHet, llJC, tolerance = 1e-4,
               label = "Het with large alpha should recover JC likelihood")
})


test_that("Het with large beta_scale recovers JC likelihood (multistate chars)", {
  tree <- .four_taxon_tree()
  mat <- matrix(c(0, 1, 2, 0,
                  1, 2, 0, 1),
                nrow = 4, ncol = 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  # Non-Het with known k=3
  modelJC <- MkPrimeModel(nCat = 1L, coding = "none")
  ptsJC   <- .het_ptrs(tree, pd, modelJC, knownStates = c("1" = 3L, "2" = 3L))
  llJC    <- eval_full_loglik_cpp(ptsJC$dataPtr, ptsJC$statePtr)

  # Het with large alpha
  modelHet <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 8L,
                            nCat = 1L, coding = "none")
  ptsHet   <- .het_ptrs(tree, pd, modelHet, beta_scale = 5000,
                         knownStates = c("1" = 3L, "2" = 3L))
  llHet    <- eval_full_loglik_cpp(ptsHet$dataPtr, ptsHet$statePtr)

  # Tolerance is wider for multistate because k rotations each introduce
  # small discretisation artifacts even at large alpha.
  expect_equal(llHet, llJC, tolerance = 1e-3,
               label = "Het with large alpha should recover JC(3) likelihood")
})


# ===========================================================================
# 2. Rotation symmetry: relabelling states → same likelihood
# ===========================================================================

test_that("Relabelling k=3 states gives identical Het likelihood", {
  tree <- .four_taxon_tree()
  # Original: states {0, 1, 2}
  mat_orig <- matrix(c(0, 1, 2, 0),
                     nrow = 4, ncol = 1,
                     dimnames = list(paste0("t", 1:4), NULL))

  # Relabelled: 0→1, 1→2, 2→0
  mat_relabel <- matrix(c(1, 2, 0, 1),
                        nrow = 4, ncol = 1,
                        dimnames = list(paste0("t", 1:4), NULL))

  model <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L,
                         nCat = 1L, coding = "none")

  pts_orig    <- .het_ptrs(tree, MatrixToPhyDat(mat_orig), model,
                           beta_scale = 2.0,
                           knownStates = c("1" = 3L))
  pts_relabel <- .het_ptrs(tree, MatrixToPhyDat(mat_relabel), model,
                           beta_scale = 2.0,
                           knownStates = c("1" = 3L))

  ll_orig    <- eval_full_loglik_cpp(pts_orig$dataPtr, pts_orig$statePtr)
  ll_relabel <- eval_full_loglik_cpp(pts_relabel$dataPtr, pts_relabel$statePtr)

  expect_equal(ll_orig, ll_relabel, tolerance = 1e-10,
               label = "Rotation symmetry: state relabelling must not change likelihood")
})


test_that("Relabelling binary states gives identical Het likelihood", {
  tree <- .four_taxon_tree()
  mat_a <- matrix(c(0, 1, 0, 1), nrow = 4, ncol = 1,
                  dimnames = list(paste0("t", 1:4), NULL))
  mat_b <- matrix(c(1, 0, 1, 0), nrow = 4, ncol = 1,
                  dimnames = list(paste0("t", 1:4), NULL))

  model <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L,
                         nCat = 1L, coding = "none")

  pts_a <- .het_ptrs(tree, MatrixToPhyDat(mat_a), model, beta_scale = 1.5)
  pts_b <- .het_ptrs(tree, MatrixToPhyDat(mat_b), model, beta_scale = 1.5)

  ll_a <- eval_full_loglik_cpp(pts_a$dataPtr, pts_a$statePtr)
  ll_b <- eval_full_loglik_cpp(pts_b$dataPtr, pts_b$statePtr)

  expect_equal(ll_a, ll_b, tolerance = 1e-10,
               label = "Binary state swap must not change Het likelihood")
})


# ===========================================================================
# 3. Hand-computed 3-taxon Het likelihood (k=2, 1 character)
# ===========================================================================

test_that("3-taxon Het likelihood matches hand computation (k=2)", {
  # 3-taxon tree: ((t1:0.1,t2:0.2):0,t3:0.15)
  tree <- .three_taxon_tree(0.1, 0.2, 0.15)
  # Single binary character: t1=0, t2=1, t3=0
  mat <- matrix(c(0, 1, 0), nrow = 3, ncol = 1,
                dimnames = list(paste0("t", 1:3), NULL))
  pd <- MatrixToPhyDat(mat)

  alpha <- 3.0
  nBins <- 4L
  # relabel = FALSE so the raw pruning likelihood matches the hand computation
  model <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = nBins,
                         nCat = 1L, coding = "none", relabel = FALSE)

  pts <- .het_ptrs(tree, pd, model, beta_scale = alpha)
  ll_cpp <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)

  # Hand computation in R:
  # Discretize Beta(alpha, alpha) into nBins equal-probability bins
  # (conditional mean within each bin)
  a <- alpha; b <- alpha  # k=2 → Beta(α, α)
  bins <- numeric(nBins)
  for (i in seq_len(nBins)) {
    lo <- qbeta((i - 1) / nBins, a, b)
    hi <- qbeta(i / nBins, a, b)
    p_lo <- pbeta(lo, a + 1, b)
    p_hi <- pbeta(hi, a + 1, b)
    denom <- pbeta(hi, a, b) - pbeta(lo, a, b)
    bins[i] <- (a / (a + b)) * (p_hi - p_lo) / denom
  }

  # Get tree structure from the C++ state (preordered)
  ptree <- Preorder(tree)
  parent <- ptree$edge[, 1]
  child  <- ptree$edge[, 2]
  el     <- pts$treeLen * pts$relBr
  nTip   <- 3L
  k      <- 2L
  nNode  <- ptree$Nnode
  root   <- nTip + 1L

  # Tip states: t1=0, t2=1, t3=0 (indexed by tip number in ptree)
  states <- c(0L, 1L, 0L)

  # For each bin, run full Felsenstein pruning on the actual tree
  char_lik <- 0.0
  for (bi in seq_along(bins)) {
    beta_val <- bins[bi]
    pi_vec <- c(beta_val, 1 - beta_val)
    mu <- 1 / (1 - sum(pi_vec^2))

    # F81 transition matrix
    f81_P <- function(t) {
      P <- matrix(0, k, k)
      d <- exp(-mu * t)
      for (ii in seq_len(k)) for (jj in seq_len(k)) {
        delta <- if (ii == jj) 1 else 0
        P[ii, jj] <- pi_vec[jj] * (1 - d) + delta * d
      }
      P
    }

    # Conditional likelihoods per node
    cl <- matrix(0, nrow = nTip + nNode, ncol = k)
    for (tip in seq_len(nTip)) cl[tip, states[tip] + 1L] <- 1

    # Bottom-up traversal (reverse preorder)
    for (e in rev(seq_along(parent))) {
      p_node <- parent[e]; ch_node <- child[e]
      P <- f81_P(el[e])
      contrib <- as.numeric(P %*% cl[ch_node, ])
      if (all(cl[p_node, ] == 0)) {
        cl[p_node, ] <- contrib
      } else {
        cl[p_node, ] <- cl[p_node, ] * contrib
      }
    }

    root_lik <- sum(pi_vec * cl[root, ])
    char_lik <- char_lik + root_lik
  }
  # Average over nBins (1 rotation for k=2)
  char_lik <- char_lik / nBins

  ll_hand <- log(char_lik)

  expect_equal(ll_cpp, ll_hand, tolerance = 1e-10,
               label = "3-taxon Het likelihood must match hand computation")
})


# ===========================================================================
# 4. Het affects likelihood (non-trivial heterogeneity)
# ===========================================================================

test_that("Het produces different likelihood than JC at moderate alpha", {
  tree <- .four_taxon_tree()
  mat <- matrix(c(0, 1, 0, 1,
                  1, 0, 1, 0),
                nrow = 4, ncol = 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  modelJC  <- MkPrimeModel(nCat = 1L, coding = "none")
  modelHet <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L,
                            nCat = 1L, coding = "none")

  ptsJC  <- .het_ptrs(tree, pd, modelJC)
  ptsHet <- .het_ptrs(tree, pd, modelHet, beta_scale = 1.0)

  llJC  <- eval_full_loglik_cpp(ptsJC$dataPtr, ptsJC$statePtr)
  llHet <- eval_full_loglik_cpp(ptsHet$dataPtr, ptsHet$statePtr)

  expect_true(is.finite(llJC))
  expect_true(is.finite(llHet))
  # At alpha=1, bins spread significantly away from 0.5 → measurable difference
  expect_false(isTRUE(all.equal(llJC, llHet, tolerance = 1e-6)),
               label = "Het at moderate alpha should differ from JC")
})


# ===========================================================================
# 5. Het with ACRV (composition)
# ===========================================================================

test_that("Het + ACRV produces finite likelihood and differs from Het alone", {
  tree <- .four_taxon_tree()
  mat <- matrix(c(0, 1, 0, 1, 1, 0, 1, 0),
                nrow = 4, ncol = 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  modelHet <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L,
                            nCat = 1L, coding = "none")
  modelHA  <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L,
                            nCat = 4L, coding = "none")

  ptsHet <- .het_ptrs(tree, pd, modelHet, beta_scale = 2.0)
  ptsHA  <- .het_ptrs(tree, pd, modelHA, beta_scale = 2.0,
                       rate_log_sd = 0.5)

  llHet <- eval_full_loglik_cpp(ptsHet$dataPtr, ptsHet$statePtr)
  llHA  <- eval_full_loglik_cpp(ptsHA$dataPtr, ptsHA$statePtr)

  expect_true(is.finite(llHet))
  expect_true(is.finite(llHA))
  expect_false(isTRUE(all.equal(llHet, llHA, tolerance = 1e-6)),
               label = "Het+ACRV should differ from Het without ACRV")
})


# ===========================================================================
# 6. Het with neomorphic rate_loss composition
# ===========================================================================

test_that("Het + rate_loss produces finite likelihood for neomorphic chars", {
  tree <- .four_taxon_tree()
  mat <- matrix(c(0, 1, 0, 1, 1, 0, 0, 1),
                nrow = 4, ncol = 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  model <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L,
                         nCat = 1L, coding = "none")

  pts <- .het_ptrs(tree, pd, model, neomorphic = 1:2,
                    rate_loss = 2.0, beta_scale = 1.5)
  ll <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)
  expect_true(is.finite(ll))
})

test_that("Het + rate_loss=1 on neomorphic ≈ Het on transformational (k=2)", {
  # When rate_loss = 1 (symmetric), neomorphic and transformational k=2
  # characters should yield the same Het likelihood (both are symmetric F81).
  tree <- .four_taxon_tree()
  mat <- matrix(c(0, 1, 0, 1), nrow = 4, ncol = 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  model <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L,
                         nCat = 1L, coding = "none", relabel = FALSE)

  # Neomorphic with rate_loss = 1 (symmetric)
  pts_neo <- .het_ptrs(tree, pd, model, neomorphic = 1L,
                        rate_loss = 1.0, beta_scale = 2.0)
  # Transformational (default, also k=2 symmetric)
  pts_trans <- .het_ptrs(tree, pd, model, beta_scale = 2.0)

  ll_neo   <- eval_full_loglik_cpp(pts_neo$dataPtr, pts_neo$statePtr)
  ll_trans <- eval_full_loglik_cpp(pts_trans$dataPtr, pts_trans$statePtr)

  expect_equal(ll_neo, ll_trans, tolerance = 1e-10,
               label = "Symmetric neomorphic Het should equal transformational Het")
})


# ===========================================================================
# 7. Het with ascertainment correction
# ===========================================================================

test_that("Het with variable coding gives different (larger abs) likelihood", {
  tree <- .four_taxon_tree()
  mat <- matrix(c(0, 1, 0, 1, 1, 0, 1, 0),
                nrow = 4, ncol = 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  modelNone <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L,
                              nCat = 1L, coding = "none")
  modelVar  <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L,
                              nCat = 1L, coding = "variable")

  ptsNone <- .het_ptrs(tree, pd, modelNone, beta_scale = 2.0)
  ptsVar  <- .het_ptrs(tree, pd, modelVar,  beta_scale = 2.0)

  llNone <- eval_full_loglik_cpp(ptsNone$dataPtr, ptsNone$statePtr)
  llVar  <- eval_full_loglik_cpp(ptsVar$dataPtr,  ptsVar$statePtr)

  expect_true(is.finite(llNone))
  expect_true(is.finite(llVar))
  # Variable coding conditions out constant sites → log-lik should increase
  expect_gt(llVar, llNone)
})


# ===========================================================================
# 8. beta_scale prior density
# ===========================================================================

test_that("LogPrior includes Gamma prior for beta_scale when Het is enabled", {
  tree <- .four_taxon_tree()
  mat <- matrix(c(0, 1, 0, 1), nrow = 4, ncol = 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  model <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L,
                         betaScaleShape = 2, betaScaleRate = 0.5,
                         expSteps = 10)

  state <- list(
    tree_length = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss = 1.0,
    rate_log_sd = 0.3,
    kPrime = 2L,
    p = 0.5,
    beta_scale = 3.0
  )

  lp <- MkPrime:::LogPrior(state, model, mkd)

  # Manual: same as non-Het prior + Gamma(2, 0.5) density on beta_scale
  lp_base <- MkPrime:::LogPrior(
    state[names(state) != "beta_scale"],
    MkPrimeModel(expSteps = 10), mkd
  )
  lp_het_term <- dgamma(3.0, shape = 2, rate = 0.5, log = TRUE)

  expect_equal(lp, lp_base + lp_het_term, tolerance = 1e-12,
               label = "LogPrior should add Gamma(shape,rate) for beta_scale")
})


test_that("LogPrior returns -Inf for beta_scale <= 0", {
  tree <- .four_taxon_tree()
  mat <- matrix(c(0, 1, 0, 1), nrow = 4, ncol = 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(qHeterogeneity = TRUE, expSteps = 10)

  state <- list(
    tree_length = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss = 1.0,
    rate_log_sd = 0.3,
    kPrime = 2L,
    p = 0.5,
    beta_scale = -1.0
  )

  expect_equal(MkPrime:::LogPrior(state, model, mkd), -Inf)

  state$beta_scale <- 0.0
  expect_equal(MkPrime:::LogPrior(state, model, mkd), -Inf)
})


# ===========================================================================
# 9. Determinism: repeated Het evaluations are identical
# ===========================================================================

test_that("Het likelihood evaluation is deterministic", {
  tree <- .four_taxon_tree()
  mat <- matrix(c(0, 1, 2, 0, 1, 2, 0, 1),
                nrow = 4, ncol = 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  model <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L,
                         nCat = 4L, coding = "variable")
  pts <- .het_ptrs(tree, pd, model, beta_scale = 2.0,
                    rate_log_sd = 0.5,
                    knownStates = c("1" = 3L, "2" = 3L))

  ll1 <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)
  ll2 <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)

  expect_equal(ll1, ll2, tolerance = 1e-14)
})


# ===========================================================================
# 10. Regression: checkpoint round-trip preserves beta_scale
# ===========================================================================

test_that("Checkpoint serialization includes beta_scale", {
  tree <- .four_taxon_tree()
  mat <- matrix(c(0, 1, 0, 1, 1, 0, 1, 0),
                nrow = 4, ncol = 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  model <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L,
                         nCat = 1L, coding = "none")

  pts <- .het_ptrs(tree, pd, model, beta_scale = 7.5)

  # Simulate what .SaveCheckpoint does: get_mcmc_state → serialize → reconstruct
  s <- get_mcmc_state(pts$statePtr)

  # The serialized list now includes beta_scale
  serialized <- list(
    log_lik        = s$logLik,
    log_prior      = s$logPrior,
    log_post       = s$logPost,
    tree_length    = s$treeLength,
    rel_br_lengths = s$relBrLengths,
    rate_loss      = s$rateLoss,
    rate_log_sd    = s$rateLogSd,
    rate_neo       = s$rateNeo,
    p              = s$p,
    kPrime         = s$kPrime,
    edge           = s$edge,
    beta_scale     = s$betaScale
  )
  expect_equal(serialized$beta_scale, 7.5)

  # Reconstruct C++ state from serialized list (simulates resume path)
  newPtr <- init_mcmc_state(
    serialized$edge[, 1L], serialized$edge[, 2L],
    serialized$rel_br_lengths, serialized$tree_length,
    serialized$rate_loss, serialized$rate_log_sd,
    serialized$rate_neo %||% 1.0, serialized$p %||% 0.5,
    as.integer(serialized$kPrime),
    serialized$log_lik, serialized$log_prior,
    serialized$beta_scale %||% 1.0
  )
  s2 <- get_mcmc_state(newPtr)
  expect_equal(s2$betaScale, 7.5)
})


# ===========================================================================
# 11. Regression: nBetaCat validation prevents buffer overrun
# ===========================================================================

test_that("MkPrimeModel rejects nBetaCat > 16", {
  expect_error(
    MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 20L),
    "nBetaCat"
  )
  expect_error(
    MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 0L),
    "nBetaCat"
  )
  # nBetaCat = 1 and 16 are valid
  expect_s3_class(
    MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 1L),
    "MkPrimeModel"
  )
  expect_s3_class(
    MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 16L),
    "MkPrimeModel"
  )
})


# ===========================================================================
# 12. Regression: small alpha produces finite Het likelihood (no NaN bins)
# ===========================================================================

test_that("Het with very small beta_scale produces finite likelihood", {
  tree <- .four_taxon_tree()
  mat <- matrix(c(0, 1, 0, 1, 1, 0, 1, 0),
                nrow = 4, ncol = 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  model <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L,
                         nCat = 1L, coding = "none")

  # alpha = 0.01 previously caused NaN in compute_het_bins for the last bin
  pts <- .het_ptrs(tree, pd, model, beta_scale = 0.01)
  ll <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)

  expect_true(is.finite(ll),
              label = "Small beta_scale should not produce NaN likelihood")
})


# ===========================================================================
# 13. Het with nBetaCat = 1 (single bin, edge case)
# ===========================================================================

test_that("Het with nBetaCat = 1 produces finite likelihood", {
  tree <- .four_taxon_tree()
  mat <- matrix(c(0, 1, 0, 1), nrow = 4, ncol = 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)

  model <- MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 1L,
                         nCat = 1L, coding = "none")
  pts <- .het_ptrs(tree, pd, model, beta_scale = 2.0)
  ll <- eval_full_loglik_cpp(pts$dataPtr, pts$statePtr)
  expect_true(is.finite(ll))
})


# ===========================================================================
# 14. M-100: Het + informative coding rejected at model construction time
# ===========================================================================

test_that("MkPrimeModel rejects qHeterogeneity + informative coding", {
  expect_error(
    MkPrimeModel(qHeterogeneity = TRUE, coding = "informative"),
    "qHeterogeneity.*informative|informative.*qHeterogeneity"
  )
  # variable and none are fine with Het
  expect_s3_class(
    MkPrimeModel(qHeterogeneity = TRUE, coding = "variable"),
    "MkPrimeModel"
  )
  expect_s3_class(
    MkPrimeModel(qHeterogeneity = TRUE, coding = "none"),
    "MkPrimeModel"
  )
})
