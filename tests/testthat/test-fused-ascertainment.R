test_that("fused JC ascertainment matches standalone", {
  set.seed(7284)
  tree <- ape::rtree(8, br = runif)
  tree$edge.length <- abs(tree$edge.length) + 0.01
  tree <- ape::unroot(tree)
  nTip <- ape::Ntip(tree)
  parent <- tree$edge[, 1]
  child  <- tree$edge[, 2]

  kStates <- 3L
  nChar <- 5
  tipStates <- matrix(sample(0:(kStates - 1), nTip * nChar, replace = TRUE),
                      nrow = nTip, ncol = nChar,
                      dimnames = list(tree$tip.label, NULL))
  pd <- MatrixToPhyDat(tipStates)

  mkd <- MkPrimeData(pd, knownStates = rep(kStates, nChar))

  ll_none <- MkpLogLikelihood(tree, mkd, coding = "none", relabel = FALSE)
  ll_var  <- MkpLogLikelihood(tree, mkd, coding = "variable", relabel = FALSE)

  rootFreqs <- rep(1 / kStates, kStates)
  p_standalone <- constant_site_prob_jc(parent, child, tree$edge.length,
                                         nTip, kStates, rootFreqs, 1.0)

  expected_correction <- -nChar * log(1 - p_standalone)
  expect_equal(ll_var - ll_none, expected_correction, tolerance = 1e-10)
})


test_that("fused MkN ascertainment matches standalone", {
  set.seed(3419)
  tree <- ape::rtree(10, br = runif)
  tree$edge.length <- abs(tree$edge.length) + 0.01
  tree <- ape::unroot(tree)
  nTip <- ape::Ntip(tree)
  parent <- tree$edge[, 1]
  child  <- tree$edge[, 2]

  nChar <- 8
  tipStates <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                      nrow = nTip, ncol = nChar,
                      dimnames = list(tree$tip.label, NULL))
  pd <- MatrixToPhyDat(tipStates)

  mkd <- MkPrimeData(pd, neomorphic = seq_len(nChar))

  rateLoss <- 1.5

  ll_none <- MkpLogLikelihood(tree, mkd, coding = "none", relabel = FALSE,
                               rate_loss = rateLoss)
  ll_var  <- MkpLogLikelihood(tree, mkd, coding = "variable", relabel = FALSE,
                               rate_loss = rateLoss)

  rootFreqs <- c(rateLoss / (1 + rateLoss), 1 / (1 + rateLoss))
  p_standalone <- constant_site_prob_mkn(parent, child, tree$edge.length,
                                          nTip, rateLoss, rootFreqs, 1.0)

  expected_correction <- -nChar * log(1 - p_standalone)
  expect_equal(ll_var - ll_none, expected_correction, tolerance = 1e-10)
})


test_that("fused JC ACRV ascertainment matches standalone", {
  set.seed(8176)
  tree <- ape::rtree(12, br = runif)
  tree$edge.length <- abs(tree$edge.length) + 0.01
  tree <- ape::unroot(tree)
  nTip <- ape::Ntip(tree)
  parent <- tree$edge[, 1]
  child  <- tree$edge[, 2]

  kStates <- 2L
  nChar <- 6
  tipStates <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                      nrow = nTip, ncol = nChar,
                      dimnames = list(tree$tip.label, NULL))
  pd <- MatrixToPhyDat(tipStates)

  mkd <- MkPrimeData(pd, knownStates = rep(kStates, nChar))

  rateLogSd <- 0.5
  nCat <- 4L

  ll_none <- MkpLogLikelihood(tree, mkd, coding = "none", relabel = FALSE,
                               nCat = nCat, rate_log_sd = rateLogSd)
  ll_var  <- MkpLogLikelihood(tree, mkd, coding = "variable", relabel = FALSE,
                               nCat = nCat, rate_log_sd = rateLogSd)

  rootFreqs <- rep(1 / kStates, kStates)
  acrv_rates <- MkPrime:::DiscreteLognormalRates(rateLogSd, nCat)
  p_standalone <- constant_site_prob_jc(parent, child, tree$edge.length,
                                         nTip, kStates, rootFreqs,
                                         acrv_rates)

  expected_correction <- -nChar * log(1 - p_standalone)
  expect_equal(ll_var - ll_none, expected_correction, tolerance = 1e-10)
})


test_that("fused Het ascertainment matches R reference (Felsenstein)", {
  # M-169: verify F81-Het constant-site probability computed via
  # Felsenstein pruning (k pseudo-characters) against a pure R reference.
  # This catches the pre-M-169 bug where the standalone function used a
  # per-edge-product formula ignoring tree topology.
  set.seed(5982)
  tree <- ape::rtree(8, br = runif)
  tree$edge.length <- abs(tree$edge.length) + 0.01
  tree <- ape::unroot(tree)
  nTip <- ape::Ntip(tree)

  nChar <- 4
  kStates <- 3L
  # Ensure all characters observe all k states (force kObs == kStates)
  tipStates <- matrix(sample(0:(kStates - 1), nTip * nChar, replace = TRUE),
                      nrow = nTip, ncol = nChar,
                      dimnames = list(tree$tip.label, NULL))
  for (j in seq_len(nChar)) {
    for (s in 0:(kStates - 1)) {
      if (!s %in% tipStates[, j]) tipStates[s + 1L, j] <- s
    }
  }
  pd <- MatrixToPhyDat(tipStates)
  # knownStates must be a named vector (names = char indices, values = k)
  knownStates <- setNames(rep(kStates, nChar), seq_len(nChar))
  mkd <- MkPrimeData(pd, knownStates = knownStates)

  betaScale <- 2.0
  rateLogSd <- 0.3
  nCat <- 2L
  nBetaCat <- 3L

  # --- R reference: F81 Felsenstein pruning for P(constant site) ---
  # Compute discretized Beta bins (same as C++)
  a <- betaScale
  b <- (kStates - 1) * betaScale
  betaBins <- vapply(seq_len(nBetaCat), function(i) {
    lo <- qbeta((i - 1) / nBetaCat, a, b)
    hi <- qbeta(i / nBetaCat, a, b)
    if (hi - lo < 1e-15) return(0.5 * (lo + hi))
    p_lo <- pbeta(lo, a + 1, b)
    p_hi <- pbeta(hi, a + 1, b)
    denom <- pbeta(hi, a, b) - pbeta(lo, a, b)
    if (denom < 1e-300) return(0.5 * (lo + hi))
    (a / (a + b)) * (p_hi - p_lo) / denom
  }, double(1))

  acrvRates <- MkPrime:::DiscreteLognormalRates(rateLogSd, nCat)

  tree0 <- Preorder(tree)
  parent <- tree0$edge[, 1]
  child  <- tree0$edge[, 2]
  edgeLen <- tree0$edge.length
  nEdge <- length(parent)
  maxNode <- 2 * nTip - 1
  root <- nTip + 1L

  nRot <- kStates  # k >= 3
  totalComp <- nCat * nBetaCat * nRot

  totalP <- 0.0
  for (cat in seq_len(nCat)) {
    acrvRate <- acrvRates[cat]
    for (bi in seq_len(nBetaCat)) {
      beta_val <- betaBins[bi]
      for (rot in seq_len(nRot)) {
        r <- (1 - beta_val) / (kStates - 1)
        pi_vec <- rep(r, kStates)
        pi_vec[rot] <- beta_val
        mu <- 1 / (1 - sum(pi_vec^2))

        # CL buffer: pseudo-chars (rows) × states (cols)
        cl <- matrix(0, nrow = maxNode, ncol = kStates * kStates)
        init <- logical(maxNode)
        # Init tips: pseudo-char s has CL = e_s
        for (tip in seq_len(nTip)) {
          for (s in seq_len(kStates)) {
            off <- (s - 1) * kStates
            cl[tip, off + s] <- 1.0
          }
          init[tip] <- TRUE
        }

        # Postorder traversal
        for (e in rev(seq_len(nEdge))) {
          par <- parent[e]
          ch  <- child[e]
          t <- edgeLen[e] * acrvRate
          d <- exp(-mu * t)
          one_minus_d <- 1 - d

          for (cc in seq_len(kStates)) {
            off <- (cc - 1) * kStates
            idx <- off + seq_len(kStates)
            sum_pi_a <- sum(pi_vec * cl[ch, idx])
            aBase <- one_minus_d * sum_pi_a
            trans <- aBase + d * cl[ch, idx]
            if (!init[par]) {
              cl[par, idx] <- trans
            } else {
              cl[par, idx] <- cl[par, idx] * trans
            }
          }
          init[par] <- TRUE
        }

        # Root: sum over pseudo-chars of π-weighted root CL
        for (cc in seq_len(kStates)) {
          off <- (cc - 1) * kStates
          idx <- off + seq_len(kStates)
          totalP <- totalP + sum(pi_vec * cl[root, idx])
        }
      }
    }
  }
  p_const_R <- totalP / totalComp

  # --- C++ pipeline: get correction from ll_var - ll_none ---
  # All chars in one known-k partition, so P(const) is uniform across chars
  # and ll_var - ll_none = -nChar * log(1 - P_const).
  model_var  <- MkPrimeModel(coding = "variable", nCat = nCat,
                              qHeterogeneity = TRUE, nBetaCat = nBetaCat)
  model_none <- MkPrimeModel(coding = "none", nCat = nCat,
                              qHeterogeneity = TRUE, nBetaCat = nBetaCat)

  .build <- function(model) {
    model <- MkPrime:::.FinalizeModel(model, tree, mkd)
    tree0 <- Preorder(tree)
    state0 <- MkPrime:::.InitState(tree0, mkd, model)
    state0$beta_scale <- betaScale
    state0$rate_log_sd <- rateLogSd
    state0$log_prior <- MkPrime:::LogPrior(state0, model, mkd)
    dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
    statePtr <- MkPrime:::.InitMcmcChain(state0)
    fill_partition_cache(dataPtr, statePtr)
    allocate_cl_workspace(dataPtr, statePtr)
    eval_full_loglik_cpp(dataPtr, statePtr)
  }

  ll_none <- .build(model_none)
  ll_var  <- .build(model_var)
  # Extract P(const) from the correction: ll_var - ll_none = -nChar*log(1-p)
  p_const_cpp <- 1 - exp(-(ll_var - ll_none) / nChar)

  expect_true(is.finite(ll_none))
  expect_true(is.finite(ll_var))
  expect_gt(ll_var, ll_none)
  expect_equal(p_const_cpp, p_const_R, tolerance = 1e-10)
})
