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
.EcologyJcLogLik_R <- function(parent, child, edgeLen, tipStates,
                                kStates, rootFreqs, rateMultipliers,
                                wEdge, zMat, phi, mode) {
  nTip  <- nrow(tipStates)
  nChar <- ncol(tipStates)
  nCat  <- length(rateMultipliers)
  kEco  <- ncol(wEdge)
  maxNode <- 2L * nTip - 1L
  root <- nTip + 1L

  rateFactor <- function(z, s) {
    if (z == 0L) return(1)
    p <- if (mode == 0L) phi[1] else phi[s + 1L]
    if (z == 1L) p else 1 / p
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
          lambda <- rateFactor(zMat[c, s + 1L], s)
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
  zMat <- matrix(0L, nrow = 2, ncol = s$kEco)
  ll_cpp <- MkPrime:::.PruningJcEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, rateMultipliers = 1,
    s$wEdge, zMat, phi = 1, mode = 0L
  )
  ll_r <- .EcologyJcLogLik_R(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 1, mode = 0L
  )
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("PruningJcEcology: z = 0 implies result is independent of phi", {
  s <- .MakeJcEcoSetup(kStates = 2L, nChar = 3L, kEco = 2L)
  set.seed(11)
  tipStates <- matrix(sample.int(s$kStates, 4 * 3, replace = TRUE) - 1L,
                      nrow = 4, ncol = 3)
  zMat <- matrix(0L, nrow = 3, ncol = s$kEco)
  ll1 <- MkPrime:::.PruningJcEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 1, mode = 0L)
  ll2 <- MkPrime:::.PruningJcEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 3.7, mode = 0L)
  expect_equal(ll1, ll2, tolerance = 1e-12)
})


test_that("PruningJcEcology matches reference for mixed z and phi != 1", {
  s <- .MakeJcEcoSetup(kStates = 3L, nChar = 4L, kEco = 2L)
  set.seed(31)
  tipStates <- matrix(sample.int(s$kStates, 4 * 4, replace = TRUE) - 1L,
                      nrow = 4, ncol = 4)
  zMat <- rbind(c(0L, 1L),
                c(1L, 0L),
                c(2L, 2L),
                c(0L, 2L))
  ll_cpp <- MkPrime:::.PruningJcEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 2.0, mode = 0L)
  ll_r <- .EcologyJcLogLik_R(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, 1,
    s$wEdge, zMat, phi = 2.0, mode = 0L)
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("PruningJcEcology matches reference with per_ecology phi", {
  s <- .MakeJcEcoSetup(kStates = 2L, nChar = 3L, kEco = 3L,
                       ecoTipStates = c(0L, 1L, 2L, 0L))
  set.seed(42)
  tipStates <- matrix(sample.int(s$kStates, 4 * 3, replace = TRUE) - 1L,
                      nrow = 4, ncol = 3)
  zMat <- rbind(c(0L, 1L, 2L),
                c(1L, 2L, 0L),
                c(2L, 0L, 1L))
  phi <- c(1.5, 2.0, 0.7)
  ll_cpp <- MkPrime:::.PruningJcEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, 1,
    s$wEdge, zMat, phi = phi, mode = 1L)
  ll_r <- .EcologyJcLogLik_R(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, 1,
    s$wEdge, zMat, phi = phi, mode = 1L)
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("PruningJcEcology matches reference with ACRV (nCat = 4)", {
  s <- .MakeJcEcoSetup(kStates = 3L, nChar = 5L, kEco = 2L)
  set.seed(99)
  tipStates <- matrix(sample.int(s$kStates, 4 * 5, replace = TRUE) - 1L,
                      nrow = 4, ncol = 5)
  zMat <- matrix(sample(0:2, 5 * 2, replace = TRUE), nrow = 5, ncol = 2L)
  zMat <- matrix(as.integer(zMat), nrow = 5, ncol = 2L)
  rateMults <- c(0.5, 0.8, 1.2, 1.5)
  ll_cpp <- MkPrime:::.PruningJcEcology(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, rateMults,
    s$wEdge, zMat, phi = 1.8, mode = 0L)
  ll_r <- .EcologyJcLogLik_R(
    s$parent, s$child, s$edgeLen, tipStates,
    s$kStates, s$rootFreqs, rateMults,
    s$wEdge, zMat, phi = 1.8, mode = 0L)
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("PruningJcEcology validates argument shapes", {
  s <- .MakeJcEcoSetup(kStates = 2L, nChar = 2L, kEco = 2L)
  tipStates <- matrix(c(0L, 0L, 1L, 1L, 1L, 0L, 0L, 1L),
                      nrow = 4, ncol = 2)
  zMat <- matrix(0L, nrow = 2, ncol = 2)
  # mode = 0 requires phi length 1
  expect_error(
    MkPrime:::.PruningJcEcology(
      s$parent, s$child, s$edgeLen, tipStates,
      s$kStates, s$rootFreqs, 1, s$wEdge, zMat,
      phi = c(1, 1), mode = 0L),
    "global mode requires"
  )
  # mode = 1 requires phi length kEco
  expect_error(
    MkPrime:::.PruningJcEcology(
      s$parent, s$child, s$edgeLen, tipStates,
      s$kStates, s$rootFreqs, 1, s$wEdge, zMat,
      phi = 1, mode = 1L),
    "per_ecology mode requires"
  )
})
