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
