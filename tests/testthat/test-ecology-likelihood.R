library("TreeTools")


# JC-K transition probabilities
.JC_P <- function(k, t) {
  e <- exp(-k * t / (k - 1))
  ps <- 1 / k + (1 - 1 / k) * e
  pd <- 1 / k - 1 / k * e
  list(same = ps, diff = pd)
}


# Brute-force per-node marginals by enumerating all internal-node state
# assignments. Returns nNode-by-kStates matrix indexed as node n -> row n.
.BruteForceMarginals <- function(parent, child, edgeLen, tipStates, kStates) {
  nTip <- length(tipStates)
  maxNode <- 2L * nTip - 1L
  nInternal <- nTip - 1L
  root <- nTip + 1L

  # State of node n given a vector of internal-node states (length nInternal,
  # indexed by node - root + 1).
  state_of <- function(n, internal) {
    if (n <= nTip) tipStates[n] else internal[n - root + 1L]
  }

  # All internal-state assignments
  grid <- expand.grid(
    rep(list(seq_len(kStates) - 1L), nInternal),
    KEEP.OUT.ATTRS = FALSE
  )

  # Compute joint probability for each assignment
  joint <- numeric(nrow(grid))
  for (g in seq_len(nrow(grid))) {
    internal <- as.integer(grid[g, ])
    # Root prior: uniform (1/k each)
    logp <- log(1 / kStates)
    for (e in seq_along(parent)) {
      p <- .JC_P(kStates, edgeLen[e])
      pa <- state_of(parent[e], internal)
      ch <- state_of(child[e], internal)
      logp <- logp + log(if (pa == ch) p$same else p$diff)
    }
    joint[g] <- exp(logp)
  }
  z <- sum(joint)

  # Marginals
  out <- matrix(0, nrow = maxNode, ncol = kStates)
  for (n in seq_len(maxNode)) {
    if (n <= nTip) {
      out[n, tipStates[n] + 1L] <- 1
    } else {
      for (s in seq_len(kStates) - 1L) {
        idx <- grid[[n - root + 1L]] == s
        out[n, s + 1L] <- sum(joint[idx]) / z
      }
    }
  }
  out
}


test_that("EcologyNodeMarginals: tips return one-hot for observed states", {
  parent <- c(4L, 5L, 5L, 4L)
  child  <- c(5L, 1L, 2L, 3L)
  edgeLen <- c(1, 1, 1, 1)
  tipStates <- c(0L, 0L, 1L)

  marg <- MkPrime:::.EcologyNodeMarginals(parent, child, edgeLen,
                                          tipStates, kStates = 2L)
  expect_equal(marg[1, ], c(1, 0))
  expect_equal(marg[2, ], c(1, 0))
  expect_equal(marg[3, ], c(0, 1))
})


test_that("EcologyNodeMarginals: rows sum to 1", {
  parent <- c(4L, 5L, 5L, 4L)
  child  <- c(5L, 1L, 2L, 3L)
  edgeLen <- c(0.5, 1.2, 0.7, 1.5)
  tipStates <- c(0L, 1L, 1L)

  marg <- MkPrime:::.EcologyNodeMarginals(parent, child, edgeLen,
                                          tipStates, kStates = 3L)
  expect_equal(rowSums(marg), rep(1, nrow(marg)), tolerance = 1e-12)
})


test_that("EcologyNodeMarginals matches hand-derived JC-2 marginals on 3 tips", {
  # Tree ((t1, t2), t3); ecology: 0, 0, 1; all edges length 1.
  parent <- c(4L, 5L, 5L, 4L)
  child  <- c(5L, 1L, 2L, 3L)
  edgeLen <- c(1, 1, 1, 1)
  tipStates <- c(0L, 0L, 1L)

  marg <- MkPrime:::.EcologyNodeMarginals(parent, child, edgeLen,
                                          tipStates, kStates = 2L)

  e2 <- exp(-2)
  ps <- (1 + e2) / 2
  pd <- (1 - e2) / 2
  z  <- 2 * ps^3 * pd + pd^4 + ps^2 * pd^2

  # Node 5 = inner (parent of t1, t2)
  expect_equal(marg[5, 1], 2 * ps^3 * pd / z, tolerance = 1e-12)
  expect_equal(marg[5, 2], (pd^4 + ps^2 * pd^2) / z, tolerance = 1e-12)
  # Node 4 = root (parent of inner and t3)
  expect_equal(marg[4, 1], (ps^3 * pd + pd^4) / z, tolerance = 1e-12)
  expect_equal(marg[4, 2], (ps^3 * pd + ps^2 * pd^2) / z, tolerance = 1e-12)
})


test_that("EcologyNodeMarginals matches brute-force on a 4-tip K=2 case", {
  # Tree (((t1,t2),t3),t4)
  # Internal: 5 = parent(1,2); 6 = parent(5,3); 7 = root, parent(6,4)
  # Wait: ape convention: tips 1..4, root = 5, then 6, 7
  # Tree: (((t1,t2),t3),t4) — let's lay out:
  #   root (5) → (inner6, t4)
  #   inner6 → (inner7, t3)
  #   inner7 → (t1, t2)
  parent <- c(5L, 6L, 7L, 7L, 6L, 5L)
  child  <- c(6L, 7L, 1L, 2L, 3L, 4L)
  edgeLen <- c(0.5, 0.3, 0.7, 0.4, 0.6, 0.8)
  tipStates <- c(0L, 1L, 0L, 1L)
  kStates <- 2L

  marg <- MkPrime:::.EcologyNodeMarginals(parent, child, edgeLen,
                                          tipStates, kStates)
  bf <- .BruteForceMarginals(parent, child, edgeLen, tipStates, kStates)
  expect_equal(marg, bf, tolerance = 1e-10)
})


test_that("EcologyNodeMarginals matches brute-force on a 4-tip K=3 case", {
  parent <- c(5L, 6L, 7L, 7L, 6L, 5L)
  child  <- c(6L, 7L, 1L, 2L, 3L, 4L)
  edgeLen <- c(0.6, 0.4, 0.9, 0.2, 1.1, 0.5)
  tipStates <- c(0L, 2L, 1L, 0L)
  kStates <- 3L

  marg <- MkPrime:::.EcologyNodeMarginals(parent, child, edgeLen,
                                          tipStates, kStates)
  bf <- .BruteForceMarginals(parent, child, edgeLen, tipStates, kStates)
  expect_equal(marg, bf, tolerance = 1e-10)
})


test_that("EcologyNodeMarginals handles missing tip states", {
  # Missing tip should get posterior predictive from the rest
  parent <- c(4L, 5L, 5L, 4L)
  child  <- c(5L, 1L, 2L, 3L)
  edgeLen <- c(1, 1, 1, 1)
  tipStates <- c(0L, 0L, -1L)  # t3 missing
  kStates <- 2L

  marg <- MkPrime:::.EcologyNodeMarginals(parent, child, edgeLen,
                                          tipStates, kStates)
  # Missing tip's marginal should be a proper distribution (not one-hot)
  expect_equal(sum(marg[3, ]), 1, tolerance = 1e-12)
  expect_true(all(marg[3, ] > 0))
  # Other tips remain one-hot
  expect_equal(marg[1, ], c(1, 0))
  expect_equal(marg[2, ], c(1, 0))
})


test_that("EcologyNodeMarginals: very long branches drive root toward uniform", {
  # All branches very long → tips give almost no information at the root
  parent <- c(4L, 5L, 5L, 4L)
  child  <- c(5L, 1L, 2L, 3L)
  edgeLen <- c(100, 100, 100, 100)
  tipStates <- c(0L, 0L, 1L)

  marg <- MkPrime:::.EcologyNodeMarginals(parent, child, edgeLen,
                                          tipStates, kStates = 2L)
  expect_equal(marg[4, ], c(0.5, 0.5), tolerance = 1e-6)
  expect_equal(marg[5, ], c(0.5, 0.5), tolerance = 1e-6)
})


# ===== EcologyEdgeWeights =====

test_that("EcologyEdgeWeights: mean of parent and child marginals", {
  marg <- rbind(
    c(1.0, 0.0),
    c(0.0, 1.0),
    c(0.5, 0.5),
    c(0.4, 0.6),
    c(0.7, 0.3)
  )
  parent <- c(4L, 4L, 5L, 5L)
  child  <- c(5L, 3L, 1L, 2L)

  w <- MkPrime:::.EcologyEdgeWeights(marg, parent, child)
  expect_equal(w[1, ], c(0.5 * (0.4 + 0.7), 0.5 * (0.6 + 0.3)))
  expect_equal(w[2, ], c(0.5 * (0.4 + 0.5), 0.5 * (0.6 + 0.5)))
  expect_equal(w[3, ], c(0.5 * (0.7 + 1.0), 0.5 * (0.3 + 0.0)))
  expect_equal(w[4, ], c(0.5 * (0.7 + 0.0), 0.5 * (0.3 + 1.0)))
  expect_equal(rowSums(w), rep(1, 4), tolerance = 1e-14)
})


test_that("EcologyEdgeWeights rows are proper probability distributions", {
  parent <- c(5L, 6L, 7L, 7L, 6L, 5L)
  child  <- c(6L, 7L, 1L, 2L, 3L, 4L)
  edgeLen <- c(0.5, 0.3, 0.7, 0.4, 0.6, 0.8)
  tipStates <- c(0L, 1L, 0L, 1L)

  marg <- MkPrime:::.EcologyNodeMarginals(parent, child, edgeLen,
                                          tipStates, kStates = 2L)
  w <- MkPrime:::.EcologyEdgeWeights(marg, parent, child)

  expect_equal(rowSums(w), rep(1, length(parent)), tolerance = 1e-12)
  expect_true(all(w >= 0))
  expect_true(all(w <= 1))
})


# ===== PruningJcEcology: transformational ecology mixture =====


# Reference implementation: explicit matrix-mult mixture pruning.
# Same algorithm as the C++ but written with no JC shortcuts, so a shared
# bug between C++ and R is unlikely.
# v2: zMat has (kEco - 1) columns (one per non-reference ecology). For
# ecology state s, the rate factor is 1 at refEcology and mu_z / gamma_e
# elsewhere, where mu_z = 1, phi, 1/phi for z in {0, 1, 2}.
.EcologyJcLogLik_R <- function(parent, child, edgeLen, tipStates,
                                kStates, rootFreqs, rateMultipliers,
                                wEdge, zMat, phi, mode,
                                refEcology = 0L, theta = NULL,
                                pi0 = 0.0) {
  nTip  <- nrow(tipStates)
  nChar <- ncol(tipStates)
  nCat  <- length(rateMultipliers)
  kEco  <- ncol(wEdge)
  maxNode <- 2L * nTip - 1L
  root <- nTip + 1L
  if (is.null(theta)) theta <- rep(0.5, max(0L, kEco - 1L))

  # gamma_e per ecology state (1 at refEcology). theta/zMat columns are
  # 0-indexed (s < refEcology) ? s : (s - 1).
  gammaE <- rep(1.0, kEco)
  for (s in seq_len(kEco) - 1L) {
    if (s == refEcology) next
    j1 <- if (s < refEcology) s + 1L else s   # 1-based theta/zMat col
    phi_s <- if (mode == 0L) phi[1] else phi[s + 1L]
    gammaE[s + 1L] <- pi0 + (1 - pi0) *
      (theta[j1] * phi_s + (1 - theta[j1]) / phi_s)
  }

  zLookup <- function(c, s) {
    if (s == refEcology) return(0L)
    zCol <- if (s < refEcology) s + 1L else s
    zMat[c, zCol]
  }

  rateFactor <- function(z, s) {
    if (s == refEcology) return(1)
    p <- if (mode == 0L) phi[1] else phi[s + 1L]
    mu <- if (z == 0L) 1 else if (z == 1L) p else 1 / p
    mu / gammaE[s + 1L]
  }

  siteLik <- numeric(nChar)
  for (cat_idx in seq_len(nCat)) {
    rate <- rateMultipliers[cat_idx]
    cl <- array(0, dim = c(maxNode, nChar, kStates))
    init <- rep(FALSE, maxNode)

    for (i in seq_len(nTip)) {
      for (c in seq_len(nChar)) {
        st <- tipStates[i, c]
        if (st < 0) cl[i, c, ] <- 1
        else        cl[i, c, st + 1L] <- 1
      }
      init[i] <- TRUE
    }

    for (e in rev(seq_along(parent))) {
      pa <- parent[e]; ch <- child[e]
      t_base <- edgeLen[e] * rate
      for (c in seq_len(nChar)) {
        Pmix <- matrix(0, kStates, kStates)
        for (s in seq_len(kEco) - 1L) {
          lambda <- rateFactor(zLookup(c, s), s)
          t_eff <- t_base * lambda
          e_term <- exp(-kStates * t_eff / (kStates - 1))
          ps <- 1 / kStates + (1 - 1 / kStates) * e_term
          pd <- 1 / kStates - 1 / kStates * e_term
          Ps <- matrix(pd, kStates, kStates)
          diag(Ps) <- ps
          Pmix <- Pmix + wEdge[e, s + 1L] * Ps
        }
        msg <- Pmix %*% cl[ch, c, ]
        if (!init[pa]) cl[pa, c, ] <- as.vector(msg)
        else            cl[pa, c, ] <- cl[pa, c, ] * as.vector(msg)
      }
      init[pa] <- TRUE
    }

    for (c in seq_len(nChar))
      siteLik[c] <- siteLik[c] + sum(rootFreqs * cl[root, c, ])
  }

  total <- 0
  for (c in seq_len(nChar)) {
    avg <- siteLik[c] / nCat
    if (avg <= 0) return(-Inf)
    total <- total + log(avg)
  }
  total
}


# Build a small test setup (4-tip tree, 1 ecology, varying chars/states).
.MakeJcEcoSetup <- function(kStates, nChar, kEco,
                             ecoTipStates = c(0L, 0L, 1L, 1L)) {
  parent <- c(5L, 6L, 7L, 7L, 6L, 5L)
  child  <- c(6L, 7L, 1L, 2L, 3L, 4L)
  edgeLen <- c(0.4, 0.5, 0.6, 0.3, 0.7, 0.2)
  marg <- MkPrime:::.EcologyNodeMarginals(parent, child, edgeLen,
                                          ecoTipStates, kEco)
  wEdge <- MkPrime:::.EcologyEdgeWeights(marg, parent, child)
  list(parent = parent, child = child, edgeLen = edgeLen,
       wEdge = wEdge, kEco = kEco, kStates = kStates,
       rootFreqs = rep(1 / kStates, kStates))
}


test_that("PruningJcEcology with z = 0 and phi = 1 matches reference R", {
  s <- .MakeJcEcoSetup(kStates = 3L, nChar = 2L, kEco = 2L)
  set.seed(7)
  tipStates <- matrix(sample.int(s$kStates, 4 * 2, replace = TRUE) - 1L,
                      nrow = 4, ncol = 2)
  zMat <- matrix(0L, nrow = 2, ncol = s$kEco - 1L)
  ll_cpp <- MkPrime:::.PruningJcEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, rateMultipliers = 1,
    s$wEdge, zMat, phi = 1, mode = 0L, refEcology = 0L
  )
  ll_r <- .EcologyJcLogLik_R(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 1, mode = 0L
  )
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


# v2: with z = 0 and non-reference ecology, the rate factor is 1/gamma_e
# which depends on phi, so the v1 phi-invariance no longer holds. The test
# is preserved (phi = 1 in both calls) to keep the smoke check that z = 0
# is a stable baseline.
test_that("PruningJcEcology: z = 0 is stable across calls (phi = 1)", {
  s <- .MakeJcEcoSetup(kStates = 2L, nChar = 3L, kEco = 2L)
  set.seed(11)
  tipStates <- matrix(sample.int(s$kStates, 4 * 3, replace = TRUE) - 1L,
                      nrow = 4, ncol = 3)
  zMat <- matrix(0L, nrow = 3, ncol = s$kEco - 1L)
  ll1 <- MkPrime:::.PruningJcEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 1, mode = 0L, refEcology = 0L)
  ll2 <- MkPrime:::.PruningJcEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 1, mode = 0L, refEcology = 0L)
  expect_equal(ll1, ll2, tolerance = 1e-12)
})


test_that("PruningJcEcology matches reference for mixed z and phi != 1", {
  skip("retired (EBE): pins C++ == legacy rate-model R reference (.EcologyJcLogLik_R), which EBE replaces. Under EBE transformational chars carry no ecology effect (R8); neomorphic correctness is covered by test-ebe-likelihood.R.")
  s <- .MakeJcEcoSetup(kStates = 3L, nChar = 4L, kEco = 2L)
  set.seed(31)
  tipStates <- matrix(sample.int(s$kStates, 4 * 4, replace = TRUE) - 1L,
                      nrow = 4, ncol = 4)
  # v2: zMat has kEco - 1 columns (one for the non-reference ecology).
  zMat <- matrix(c(1L, 0L, 2L, 2L), nrow = 4, ncol = 1L)
  ll_cpp <- MkPrime:::.PruningJcEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 2.0, mode = 0L, refEcology = 0L)
  ll_r <- .EcologyJcLogLik_R(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 2.0, mode = 0L)
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("PruningJcEcology matches reference with per_ecology phi", {
  skip("retired (EBE): pins C++ == legacy rate-model R reference (.EcologyJcLogLik_R), which EBE replaces. Under EBE transformational chars carry no ecology effect (R8); neomorphic correctness is covered by test-ebe-likelihood.R.")
  s <- .MakeJcEcoSetup(kStates = 2L, nChar = 3L, kEco = 3L,
                       ecoTipStates = c(0L, 1L, 2L, 0L))
  set.seed(42)
  tipStates <- matrix(sample.int(s$kStates, 4 * 3, replace = TRUE) - 1L,
                      nrow = 4, ncol = 3)
  # v2: zMat is nChar x (kEco - 1); columns map to non-ref ecologies 1, 2.
  zMat <- rbind(c(1L, 2L),
                c(2L, 0L),
                c(0L, 1L))
  phi <- c(1.5, 2.0, 0.7)
  ll_cpp <- MkPrime:::.PruningJcEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, 1,
    s$wEdge, zMat, phi = phi, mode = 1L, refEcology = 0L)
  ll_r <- .EcologyJcLogLik_R(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, 1,
    s$wEdge, zMat, phi = phi, mode = 1L)
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("PruningJcEcology matches reference with ACRV (nCat = 4)", {
  skip("retired (EBE): pins C++ == legacy rate-model R reference (.EcologyJcLogLik_R), which EBE replaces. Under EBE transformational chars carry no ecology effect (R8); neomorphic correctness is covered by test-ebe-likelihood.R.")
  s <- .MakeJcEcoSetup(kStates = 3L, nChar = 5L, kEco = 2L)
  set.seed(99)
  tipStates <- matrix(sample.int(s$kStates, 4 * 5, replace = TRUE) - 1L,
                      nrow = 4, ncol = 5)
  zMat <- matrix(as.integer(sample(0:2, 5 * 1, replace = TRUE)),
                 nrow = 5, ncol = 1L)
  rateMults <- c(0.5, 0.8, 1.2, 1.5)
  ll_cpp <- MkPrime:::.PruningJcEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, rateMults,
    s$wEdge, zMat, phi = 1.8, mode = 0L, refEcology = 0L)
  ll_r <- .EcologyJcLogLik_R(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, rateMults,
    s$wEdge, zMat, phi = 1.8, mode = 0L)
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


# ===== PruningMknEcology: neomorphic ecology mixture =====


# Reference implementation for the binary asymmetric ecology mixture.
# v2: zMat has (kEco - 1) columns; rates for non-ref ecology are scaled by
# (mu_z / gamma_e); refEcology keeps base rates.
.EcologyMknLogLik_R <- function(parent, child, edgeLen, tipStates,
                                 rateLoss, rootFreqs, rateMultipliers,
                                 wEdge, zMat, phi, mode,
                                 refEcology = 0L, theta = NULL,
                                 pi0 = 0.0) {
  nTip  <- nrow(tipStates)
  nChar <- ncol(tipStates)
  nCat  <- length(rateMultipliers)
  kEco  <- ncol(wEdge)
  maxNode <- 2L * nTip - 1L
  root <- nTip + 1L
  if (is.null(theta)) theta <- rep(0.5, max(0L, kEco - 1L))

  r01_base <- 2 / (1 + rateLoss)
  r10_base <- 2 * rateLoss / (1 + rateLoss)

  gammaE <- rep(1.0, kEco)
  for (s in seq_len(kEco) - 1L) {
    if (s == refEcology) next
    j1 <- if (s < refEcology) s + 1L else s
    phi_s <- if (mode == 0L) phi[1] else phi[s + 1L]
    gammaE[s + 1L] <- pi0 + (1 - pi0) *
      (theta[j1] * phi_s + (1 - theta[j1]) / phi_s)
  }

  zLookup <- function(c, s) {
    if (s == refEcology) return(0L)
    zCol <- if (s < refEcology) s + 1L else s
    zMat[c, zCol]
  }

  ratesForState <- function(z, s) {
    if (s == refEcology) return(c(r01_base, r10_base))
    p <- if (mode == 0L) phi[1] else phi[s + 1L]
    gE <- gammaE[s + 1L]
    if (z == 0L) c(r01_base / gE, r10_base / gE)
    else if (z == 1L) c(r01_base * p / gE, r10_base / p / gE)
    else c(r01_base / p / gE, r10_base * p / gE)
  }

  siteLik <- numeric(nChar)
  for (cat_idx in seq_len(nCat)) {
    rate <- rateMultipliers[cat_idx]
    cl <- array(0, dim = c(maxNode, nChar, 2))
    init <- rep(FALSE, maxNode)

    for (i in seq_len(nTip)) {
      for (c in seq_len(nChar)) {
        st <- tipStates[i, c]
        if (st < 0) cl[i, c, ] <- 1
        else        cl[i, c, st + 1L] <- 1
      }
      init[i] <- TRUE
    }

    for (e in rev(seq_along(parent))) {
      pa <- parent[e]; ch <- child[e]
      t <- edgeLen[e] * rate
      for (c in seq_len(nChar)) {
        Pmix <- matrix(0, 2, 2)
        for (s in seq_len(kEco) - 1L) {
          r <- ratesForState(zLookup(c, s), s)
          r01 <- r[1]; r10 <- r[2]
          lam <- r01 + r10
          pi0_s <- r10 / lam; pi1_s <- r01 / lam
          ex  <- exp(-lam * t)
          Ps <- matrix(c(pi0_s + pi1_s * ex, pi1_s - pi1_s * ex,
                         pi0_s - pi0_s * ex, pi1_s + pi0_s * ex),
                        nrow = 2, byrow = TRUE)
          Pmix <- Pmix + wEdge[e, s + 1L] * Ps
        }
        msg <- Pmix %*% cl[ch, c, ]
        if (!init[pa]) cl[pa, c, ] <- as.vector(msg)
        else            cl[pa, c, ] <- cl[pa, c, ] * as.vector(msg)
      }
      init[pa] <- TRUE
    }

    for (c in seq_len(nChar))
      siteLik[c] <- siteLik[c] + sum(rootFreqs * cl[root, c, ])
  }

  total <- 0
  for (c in seq_len(nChar)) {
    avg <- siteLik[c] / nCat
    if (avg <= 0) return(-Inf)
    total <- total + log(avg)
  }
  total
}


.MakeMknEcoSetup <- function(nChar, kEco, rateLoss = 1.0,
                              ecoTipStates = c(0L, 0L, 1L, 1L)) {
  parent <- c(5L, 6L, 7L, 7L, 6L, 5L)
  child  <- c(6L, 7L, 1L, 2L, 3L, 4L)
  edgeLen <- c(0.4, 0.5, 0.6, 0.3, 0.7, 0.2)
  marg <- MkPrime:::.EcologyNodeMarginals(parent, child, edgeLen,
                                          ecoTipStates, kEco)
  wEdge <- MkPrime:::.EcologyEdgeWeights(marg, parent, child)
  pi0 <- rateLoss / (1 + rateLoss)
  pi1 <- 1 / (1 + rateLoss)
  list(parent = parent, child = child, edgeLen = edgeLen,
       wEdge = wEdge, kEco = kEco,
       rateLoss = rateLoss,
       rootFreqs = c(pi0, pi1))
}


test_that("PruningMknEcology with z = 0 and phi = 1 matches reference R", {
  s <- .MakeMknEcoSetup(nChar = 3L, kEco = 2L)
  tipStates <- matrix(c(0L, 1L, 0L, 1L,
                        1L, 0L, 1L, 0L,
                        0L, 0L, 1L, 1L), nrow = 4, byrow = FALSE)
  zMat <- matrix(0L, nrow = 3, ncol = s$kEco - 1L)
  ll_cpp <- MkPrime:::.PruningMknEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$rateLoss, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 1, mode = 0L, refEcology = 0L)
  ll_r <- .EcologyMknLogLik_R(
    s$parent, s$child, s$edgeLen, tipStates,
    s$rateLoss, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 1, mode = 0L)
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


# v2: gamma-normalisation breaks the v1 phi-invariance under z = 0.
# Test preserved as a stability smoke check at phi = 1.
test_that("PruningMknEcology: z = 0 is stable across calls (phi = 1)", {
  s <- .MakeMknEcoSetup(nChar = 2L, kEco = 2L, rateLoss = 0.8)
  tipStates <- matrix(c(0L, 1L, 0L, 1L,
                        1L, 1L, 0L, 0L), nrow = 4, byrow = FALSE)
  zMat <- matrix(0L, nrow = 2, ncol = s$kEco - 1L)
  ll1 <- MkPrime:::.PruningMknEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$rateLoss, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 1, mode = 0L, refEcology = 0L)
  ll2 <- MkPrime:::.PruningMknEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$rateLoss, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 1, mode = 0L, refEcology = 0L)
  expect_equal(ll1, ll2, tolerance = 1e-12)
})


test_that("PruningMknEcology matches reference for mixed z and phi != 1", {
  skip("retired (EBE): pins C++ == legacy rate-model R reference (.EcologyMknLogLik_R), which EBE replaces. EBE neomorphic equilibrium-tilt correctness is covered by test-ebe-likelihood.R (T2 oracle).")
  s <- .MakeMknEcoSetup(nChar = 4L, kEco = 2L, rateLoss = 1.5)
  set.seed(13)
  tipStates <- matrix(sample(0:1, 4 * 4, replace = TRUE), nrow = 4, ncol = 4)
  tipStates <- matrix(as.integer(tipStates), nrow = 4, ncol = 4)
  # v2: zMat is nChar x (kEco - 1).
  zMat <- matrix(c(1L, 2L, 0L, 1L), nrow = 4, ncol = 1L)
  ll_cpp <- MkPrime:::.PruningMknEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$rateLoss, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 2.5, mode = 0L, refEcology = 0L)
  ll_r <- .EcologyMknLogLik_R(
    s$parent, s$child, s$edgeLen, tipStates,
    s$rateLoss, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 2.5, mode = 0L)
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("PruningMknEcology matches reference with per_ecology phi", {
  skip("retired (EBE): pins C++ == legacy rate-model R reference (.EcologyMknLogLik_R), which EBE replaces. EBE neomorphic equilibrium-tilt correctness is covered by test-ebe-likelihood.R (T2 oracle).")
  s <- .MakeMknEcoSetup(nChar = 3L, kEco = 3L, rateLoss = 0.5,
                        ecoTipStates = c(0L, 1L, 2L, 0L))
  set.seed(101)
  tipStates <- matrix(sample(0:1, 4 * 3, replace = TRUE), nrow = 4, ncol = 3)
  tipStates <- matrix(as.integer(tipStates), nrow = 4, ncol = 3)
  # v2: zMat is nChar x (kEco - 1) — columns map to non-ref ecologies 1, 2.
  zMat <- rbind(c(1L, 2L),
                c(2L, 0L),
                c(0L, 1L))
  phi <- c(1.2, 2.5, 0.6)
  ll_cpp <- MkPrime:::.PruningMknEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$rateLoss, s$rootFreqs, 1,
    s$wEdge, zMat, phi = phi, mode = 1L, refEcology = 0L)
  ll_r <- .EcologyMknLogLik_R(
    s$parent, s$child, s$edgeLen, tipStates,
    s$rateLoss, s$rootFreqs, 1,
    s$wEdge, zMat, phi = phi, mode = 1L)
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("PruningMknEcology matches reference with ACRV (nCat = 4)", {
  skip("retired (EBE): pins C++ == legacy rate-model R reference (.EcologyMknLogLik_R), which EBE replaces. EBE neomorphic equilibrium-tilt correctness is covered by test-ebe-likelihood.R (T2 oracle).")
  s <- .MakeMknEcoSetup(nChar = 5L, kEco = 2L, rateLoss = 1.2)
  set.seed(77)
  tipStates <- matrix(sample(0:1, 4 * 5, replace = TRUE), nrow = 4, ncol = 5)
  tipStates <- matrix(as.integer(tipStates), nrow = 4, ncol = 5)
  zMat <- matrix(as.integer(sample(0:2, 5 * 1, replace = TRUE)),
                  nrow = 5, ncol = 1L)
  rateMults <- c(0.4, 0.9, 1.1, 1.6)
  ll_cpp <- MkPrime:::.PruningMknEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$rateLoss, s$rootFreqs, rateMults,
    s$wEdge, zMat, phi = 1.7, mode = 0L, refEcology = 0L)
  ll_r <- .EcologyMknLogLik_R(
    s$parent, s$child, s$edgeLen, tipStates,
    s$rateLoss, s$rootFreqs, rateMults,
    s$wEdge, zMat, phi = 1.7, mode = 0L)
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("PruningMknEcology validates argument shapes", {
  s <- .MakeMknEcoSetup(nChar = 2L, kEco = 2L)
  tipStates <- matrix(c(0L, 1L, 0L, 1L,
                        1L, 0L, 1L, 0L), nrow = 4, ncol = 2)
  zMat <- matrix(0L, nrow = 2, ncol = s$kEco - 1L)
  # v2 allows phi length 1 or kEco in global mode; use kEco + 1 to trigger.
  expect_error(
    MkPrime:::.PruningMknEcology(
      s$parent, s$child, s$edgeLen, tipStates,
      s$rateLoss, s$rootFreqs, 1, s$wEdge, zMat,
      phi = c(1, 1, 1), mode = 0L, refEcology = 0L),
    "global mode requires"
  )
  expect_error(
    MkPrime:::.PruningMknEcology(
      s$parent, s$child, s$edgeLen, tipStates,
      s$rateLoss, c(0.5, 0.5), 1, s$wEdge, zMat,
      phi = 1, mode = 1L, refEcology = 0L),
    "per_ecology mode requires"
  )
  expect_error(
    MkPrime:::.PruningMknEcology(
      s$parent, s$child, s$edgeLen, tipStates,
      -1, s$rootFreqs, 1, s$wEdge, zMat,
      phi = 1, mode = 0L, refEcology = 0L),
    "rateLoss must be positive"
  )
})


# ===== .MkpEcologyLogLikelihood: full orchestrator integration =====


test_that("T-010: orchestrator output reproducible across calls (cache-safe)", {
  # The T-010 refactor splits cpp_log_likelihood_ecology into a per-partition
  # function plus a thin outer loop. The orchestrator-total = sum(per-char)
  # invariant is exercised by test-ecology-plumbing.R::"per-char ecology log-liks
  # sum to total (...)" across multiple coding/mode combinations.
  #
  # This additional probe asserts that two back-to-back orchestrator calls with
  # the same inputs return bit-identical results — a regression guard against
  # any cache-state-leakage bug introduced by the T-010 refactor.
  set.seed(20260521)
  tips <- paste0("t", 1:8)
  nChar <- 12L
  mat <- matrix(0L, nrow = 8L, ncol = nChar + 1L,
                dimnames = list(tips, NULL))
  for (c in seq_len(nChar)) mat[, c] <- sample.int(3L, 8L, replace = TRUE) - 1L
  mat[, nChar + 1L] <- sample.int(4L, 8L, replace = TRUE) - 1L  # ecology kEco=4
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, ecology = nChar + 1L)

  tree <- TreeTools::Preorder(ape::rtree(8, tip.label = tips))

  zCols <- mkd$kEcology - 1L
  zMat <- matrix(sample.int(3L, mkd$nChar * zCols, replace = TRUE) - 1L,
                 nrow = mkd$nChar, ncol = zCols)
  storage.mode(zMat) <- "integer"

  ll1 <- MkPrime:::.MkpEcologyLogLikelihood(
    tree, mkd, kPrime = mkd$kObs,
    rate_loss = 1.0, rate_log_sd = 0, nCat = 1L,
    rate_neo = 1.0, relabel = TRUE,
    phi = 1.6, zMat = zMat, magnitudeMode = "global",
    refEcology = 0L, theta = rep(0.5, zCols), pi0 = 0.5
  )
  ll2 <- MkPrime:::.MkpEcologyLogLikelihood(
    tree, mkd, kPrime = mkd$kObs,
    rate_loss = 1.0, rate_log_sd = 0, nCat = 1L,
    rate_neo = 1.0, relabel = TRUE,
    phi = 1.6, zMat = zMat, magnitudeMode = "global",
    refEcology = 0L, theta = rep(0.5, zCols), pi0 = 0.5
  )
  expect_equal(ll1, ll2, tolerance = 0)
  expect_true(is.finite(ll1))
})


test_that(".MkpEcologyLogLikelihood with z = 0 matches non-ecology likelihood", {
  # Construct a small MkPrimeData with one ecology column.
  # Tree: 6-tip random; chars: 4 transformational + 2 neomorphic + 1 ecology.
  set.seed(1234)
  tips <- paste0("t", 1:6)
  mat <- matrix(c(
    # transformational chars (3 states each)
    0, 1, 2, 0, 1, 2,
    1, 0, 1, 2, 2, 0,
    2, 1, 0, 1, 0, 2,
    0, 0, 1, 1, 2, 2,
    # neomorphic chars (binary)
    0, 1, 0, 1, 0, 1,
    1, 0, 1, 0, 1, 0,
    # ecology (3 states)
    0, 0, 1, 1, 2, 2
  ), nrow = 6, ncol = 7, byrow = FALSE,
     dimnames = list(tips, NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = c(5L, 6L), ecology = 7L)

  tree <- TreeTools::Preorder(ape::rtree(6, tip.label = tips))

  # Ecology-blind likelihood (coding = "none" to match the eco path)
  ll_baseline <- MkPrime:::.MkpLogLikelihood(
    tree, mkd,
    kPrime = mkd$kObs,
    rate_loss = 1.0, rate_log_sd = 0, nCat = 1L,
    coding = "none", rate_neo = 1.0, relabel = TRUE
  )

  # Ecology-aware with z = 0 everywhere: rate factors all 1 at phi = 1
  # (v2: gamma-normalisation makes non-ref factors phi-dependent otherwise).
  zMat <- matrix(0L, nrow = mkd$nChar, ncol = mkd$kEcology - 1L)
  ll_eco <- MkPrime:::.MkpEcologyLogLikelihood(
    tree, mkd, kPrime = mkd$kObs,
    rate_loss = 1.0, rate_log_sd = 0, nCat = 1L,
    rate_neo = 1.0, relabel = TRUE,
    phi = 1, zMat = zMat, magnitudeMode = "global"
  )

  expect_equal(ll_eco, ll_baseline, tolerance = 1e-10)
})


# v2: phi-invariance under z = 0 was a v1 property removed by gamma-
# normalisation. Test preserved as a degenerate stability check at phi = 1.
test_that(".MkpEcologyLogLikelihood: z = 0 is stable across calls (phi = 1)", {
  set.seed(5678)
  tips <- paste0("t", 1:5)
  mat <- matrix(c(
    0, 1, 0, 1, 2,
    1, 0, 1, 2, 0,
    0, 0, 1, 1, 1,    # neomorphic
    0, 1, 1, 2, 2     # ecology
  ), nrow = 5, ncol = 4, byrow = FALSE,
     dimnames = list(tips, NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = 3L, ecology = 4L)
  tree <- TreeTools::Preorder(ape::rtree(5, tip.label = tips))

  zMat <- matrix(0L, nrow = mkd$nChar, ncol = mkd$kEcology - 1L)
  ll1 <- MkPrime:::.MkpEcologyLogLikelihood(
    tree, mkd, kPrime = mkd$kObs,
    rate_loss = 0.8, rate_log_sd = 0, nCat = 1L,
    rate_neo = 1.2, relabel = TRUE,
    phi = 1, zMat = zMat, magnitudeMode = "global"
  )
  ll2 <- MkPrime:::.MkpEcologyLogLikelihood(
    tree, mkd, kPrime = mkd$kObs,
    rate_loss = 0.8, rate_log_sd = 0, nCat = 1L,
    rate_neo = 1.2, relabel = TRUE,
    phi = 1, zMat = zMat, magnitudeMode = "global"
  )
  expect_equal(ll1, ll2, tolerance = 1e-12)
})


test_that(".MkpEcologyLogLikelihood: phi > 1 with non-zero z changes likelihood", {
  set.seed(91011)
  tips <- paste0("t", 1:5)
  mat <- matrix(c(
    0, 1, 1, 0, 1,
    1, 0, 0, 1, 0,
    0, 0, 1, 1, 1,
    0, 1, 1, 2, 2
  ), nrow = 5, ncol = 4, byrow = FALSE,
     dimnames = list(tips, NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = 3L, ecology = 4L)
  tree <- TreeTools::Preorder(ape::rtree(5, tip.label = tips))

  # v2: zMat is nChar x (kEco - 1) (kEcology is 3 → 2 non-ref columns).
  zMat <- matrix(c(1L, 2L,
                   2L, 0L,
                   0L, 1L), nrow = 3, ncol = 2, byrow = TRUE)
  ll1 <- MkPrime:::.MkpEcologyLogLikelihood(
    tree, mkd, kPrime = mkd$kObs,
    rate_loss = 1.0, rate_log_sd = 0, nCat = 1L,
    rate_neo = 1.0, relabel = TRUE,
    phi = 1, zMat = zMat, magnitudeMode = "global"
  )
  ll2 <- MkPrime:::.MkpEcologyLogLikelihood(
    tree, mkd, kPrime = mkd$kObs,
    rate_loss = 1.0, rate_log_sd = 0, nCat = 1L,
    rate_neo = 1.0, relabel = TRUE,
    phi = 2.5, zMat = zMat, magnitudeMode = "global"
  )
  # With non-zero z and phi != 1, the likelihood should change.
  expect_false(isTRUE(all.equal(ll1, ll2, tolerance = 1e-6)))
  expect_true(is.finite(ll1))
  expect_true(is.finite(ll2))
})


test_that(".MkpEcologyLogLikelihood with z = 0 + coding=variable matches standard", {
  set.seed(8211)
  tips <- paste0("t", 1:6)
  mat <- matrix(c(
    0, 1, 2, 0, 1, 2,
    1, 0, 1, 2, 2, 0,
    2, 1, 0, 1, 0, 2,
    0, 1, 0, 1, 0, 1,
    1, 0, 1, 0, 1, 0,
    0, 0, 1, 1, 2, 2
  ), nrow = 6, ncol = 6, byrow = FALSE, dimnames = list(tips, NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = c(4L, 5L), ecology = 6L)
  tree <- TreeTools::Preorder(ape::rtree(6, tip.label = tips))

  ll_baseline <- MkPrime:::.MkpLogLikelihood(
    tree, mkd, kPrime = mkd$kObs,
    rate_loss = 1.0, rate_log_sd = 0, nCat = 1L,
    coding = "variable", rate_neo = 1.0, relabel = TRUE
  )

  # v2: phi must be 1 for z = 0 to match the non-ecology baseline; with phi != 1
  # the non-reference ecology rate factor 1/gamma_e drifts away from 1.
  zMat <- matrix(0L, nrow = mkd$nChar, ncol = mkd$kEcology - 1L)
  ll_eco <- MkPrime:::.MkpEcologyLogLikelihood(
    tree, mkd, kPrime = mkd$kObs,
    rate_loss = 1.0, rate_log_sd = 0, nCat = 1L,
    rate_neo = 1.0, relabel = TRUE,
    phi = 1, zMat = zMat, magnitudeMode = "global",
    coding = "variable"
  )
  expect_equal(ll_eco, ll_baseline, tolerance = 1e-10)
})


test_that(".MkpEcologyLogLikelihood with ACRV + coding=variable + z = 0 matches standard", {
  set.seed(91)
  tips <- paste0("t", 1:6)
  mat <- matrix(c(
    0, 1, 2, 0, 1, 2,
    1, 0, 1, 2, 2, 0,
    0, 1, 0, 1, 0, 1,
    1, 0, 1, 0, 1, 0,
    0, 0, 1, 1, 2, 2
  ), nrow = 6, ncol = 5, byrow = FALSE, dimnames = list(tips, NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = c(3L, 4L), ecology = 5L)
  tree <- TreeTools::Preorder(ape::rtree(6, tip.label = tips))

  ll_baseline <- MkPrime:::.MkpLogLikelihood(
    tree, mkd, kPrime = mkd$kObs,
    rate_loss = 1.0, rate_log_sd = 0.4, nCat = 4L,
    coding = "variable", rate_neo = 1.0, relabel = TRUE
  )
  zMat <- matrix(0L, nrow = mkd$nChar, ncol = mkd$kEcology - 1L)
  ll_eco <- MkPrime:::.MkpEcologyLogLikelihood(
    tree, mkd, kPrime = mkd$kObs,
    rate_loss = 1.0, rate_log_sd = 0.4, nCat = 4L,
    rate_neo = 1.0, relabel = TRUE,
    phi = 1.0, zMat = zMat, magnitudeMode = "global",
    coding = "variable"
  )
  expect_equal(ll_eco, ll_baseline, tolerance = 1e-10)
})


test_that(".MkpEcologyLogLikelihood errors on unsupported coding (informative)", {
  set.seed(771)
  tips <- paste0("t", 1:5)
  mat <- matrix(c(0, 1, 0, 1, 2,
                  1, 0, 1, 2, 0,
                  0, 1, 1, 2, 2),
                nrow = 5, ncol = 3,
                dimnames = list(tips, NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, ecology = 3L)
  tree <- TreeTools::Preorder(ape::rtree(5, tip.label = tips))
  zMat <- matrix(0L, nrow = mkd$nChar, ncol = mkd$kEcology - 1L)
  expect_error(
    MkPrime:::.MkpEcologyLogLikelihood(
      tree, mkd, kPrime = mkd$kObs,
      rate_loss = 1.0, rate_log_sd = 0, nCat = 1L,
      rate_neo = 1.0, relabel = TRUE,
      phi = 1.0, zMat = zMat, magnitudeMode = "global",
      coding = "informative"
    ),
    "should be one of"
  )
})


test_that(".MkpEcologyLogLikelihood errors when mkd has no ecology", {
  tips <- paste0("t", 1:5)
  mat <- matrix(c(0, 1, 0, 1, 2,
                  1, 0, 1, 2, 0), nrow = 5, ncol = 2,
                dimnames = list(tips, NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)  # no ecology
  tree <- TreeTools::Preorder(ape::rtree(5, tip.label = tips))
  zMat <- matrix(0L, nrow = mkd$nChar, ncol = 2)
  expect_error(
    MkPrime:::.MkpEcologyLogLikelihood(
      tree, mkd, kPrime = mkd$kObs,
      rate_loss = 1.0, rate_log_sd = 0, nCat = 1L,
      rate_neo = 1.0, relabel = TRUE,
      phi = 1, zMat = zMat, magnitudeMode = "global"
    ),
    "requires.*ecology"
  )
})


test_that("PruningJcEcology validates argument shapes", {
  s <- .MakeJcEcoSetup(kStates = 2L, nChar = 2L, kEco = 2L)
  tipStates <- matrix(c(0L, 0L, 1L, 1L, 1L, 0L, 0L, 1L),
                      nrow = 4, ncol = 2)
  zMat <- matrix(0L, nrow = 2, ncol = s$kEco - 1L)
  # v2 allows phi length 1 or kEco in global mode; use kEco + 1 to trigger.
  expect_error(
    MkPrime:::.PruningJcEcology(
      s$parent, s$child, s$edgeLen, tipStates,
      s$kStates, s$rootFreqs, 1, s$wEdge, zMat,
      phi = c(1, 1, 1), mode = 0L, refEcology = 0L),
    "global mode requires"
  )
  # mode = 1 requires phi length kEco
  expect_error(
    MkPrime:::.PruningJcEcology(
      s$parent, s$child, s$edgeLen, tipStates,
      s$kStates, s$rootFreqs, 1, s$wEdge, zMat,
      phi = 1, mode = 1L, refEcology = 0L),
    "per_ecology mode requires"
  )
})
