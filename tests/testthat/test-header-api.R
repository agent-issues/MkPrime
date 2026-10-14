# Public C++ API in inst/include/MkPrime, exercised through the hooks in
# src/header_api.cpp.

.HeaderFixture <- function(nTip = 10, nChar = 25, k = 3, seed = 6107,
                           missing = 0.15) {
  set.seed(seed)
  tree <- Preorder(RandomTree(nTip, root = FALSE))
  tree$edge.length <- rexp(nrow(tree$edge), 6)
  states <- matrix(sample.int(k, nTip * nChar, replace = TRUE) - 1L, nTip)
  states[sample.int(length(states), round(missing * length(states)))] <- -1L
  list(tree = tree, states = states, k = k)
}

.TreeList <- function(tree) {
  list(parent = tree$edge[, 1], child = tree$edge[, 2],
       edge_length = tree$edge.length)
}

# Felsenstein pruning written independently, with explicit P(t) matrices.
.RefSiteLik <- function(tree, states, k, PMatrix, rootFreqs, rates = 1) {
  tree <- Postorder(tree)
  edge <- tree$edge
  nTip <- NTip(tree)
  nChar <- ncol(states)
  rootNode <- edge[nrow(edge), 1]
  total <- numeric(nChar)
  for (rate in rates) {
    cl <- array(1, c(max(edge), k, nChar))
    for (tip in seq_len(nTip)) {
      for (ch in seq_len(nChar)) {
        if (states[tip, ch] >= 0) {
          cl[tip, , ch] <- 0
          cl[tip, states[tip, ch] + 1, ch] <- 1
        }
      }
    }
    for (e in seq_len(nrow(edge))) {
      P <- PMatrix(tree$edge.length[e] * rate)
      cl[edge[e, 1], , ] <- cl[edge[e, 1], , , drop = FALSE] *
        array(P %*% matrix(cl[edge[e, 2], , ], k), c(1, k, nChar))
    }
    total <- total + colSums(rootFreqs * matrix(cl[rootNode, , ], k))
  }
  # Return:
  total / length(rates)
}

.PathToRoot <- function(tree, node) {
  n <- 0L
  while (length(node)) {
    n <- n + 1L
    node <- tree$edge[tree$edge[, 2] == node, 1]
  }
  # Return:
  n
}

.F81P <- function(pi, mu) {
  function(t) {
    d <- exp(-mu * t)
    d * diag(length(pi)) + (1 - d) * matrix(pi, length(pi), length(pi),
                                            byrow = TRUE)
  }
}

test_that("header API version is exposed", {
  expect_identical(header_api_version(), 2L)
})

test_that("transition kernels match the exported P(t) matrices", {
  tt <- c(0, 1e-12, 1e-6, 0.01, 0.3, 2, 40)
  for (k in c(2L, 3L, 7L)) {
    p <- header_jc_probs(k, tt)
    for (i in seq_along(tt)) {
      P <- jc_transition_probs(k, tt[i])
      expect_identical(p[i, 1:2], c(P[1, 1], P[1, 2]))
    }
    expect_equal(p[, 3], p[, 4], tolerance = 1e-15)
  }
  p <- header_mkn_probs(0.6, tt)
  for (i in seq_along(tt)) {
    expect_identical(p[i, ], as.vector(t(mkn_transition_probs(0.6, tt[i]))))
  }
})

test_that("lognormal ACRV rates have mean 1", {
  z <- qnorm((seq_len(5) - 0.5) / 5)
  r <- header_lognormal_rates(0.7, z)
  expect_equal(mean(r), 1, tolerance = 1e-15)
  expect_equal(r, exp(0.7 * z) / mean(exp(0.7 * z)), tolerance = 1e-14)
  expect_identical(header_lognormal_rates(0, z), rep(1, 5))
})

test_that("Pruner matches independent pruning for each model", {
  fx <- .HeaderFixture()
  tr <- fx$tree
  edge <- tr$edge
  rates <- header_lognormal_rates(0.5, qnorm((1:3 - 0.5) / 3))
  jc <- header_site_likelihoods("jc", edge[, 1], edge[, 2], tr$edge.length,
                                fx$states, 3L, 0, numeric(0), rates)
  expect_equal(jc, .RefSiteLik(tr, fx$states, 3, function(t) {
    jc_transition_probs(3L, t)
  }, rep(1 / 3, 3), rates), tolerance = 1e-13)

  binary <- pmin(fx$states, 1L)
  mkn <- header_site_likelihoods("mkn", edge[, 1], edge[, 2], tr$edge.length,
                                 binary, 2L, 0.4, numeric(0), 1)
  expect_equal(mkn, .RefSiteLik(tr, binary, 2, function(t) {
    mkn_transition_probs(0.4, t)
  }, mkn_stationary_freqs(0.4)), tolerance = 1e-13)

  pi <- c(0.5, 0.3, 0.2)
  mu <- 1 / (1 - sum(pi ^ 2))
  f81 <- header_site_likelihoods("f81", edge[, 1], edge[, 2], tr$edge.length,
                                 fx$states, 3L, mu, pi, 1)
  expect_equal(f81, .RefSiteLik(tr, fx$states, 3, .F81P(pi, mu), pi),
               tolerance = 1e-13)
})

test_that("rooting does not change reversible likelihoods", {
  fx <- .HeaderFixture(nTip = 9)
  unrooted <- fx$tree
  rooted <- Preorder(ape::root(unrooted, 3, resolve.root = TRUE))
  expect_true(TreeIsRooted(rooted))
  expect_equal(sum(rooted$edge.length), sum(unrooted$edge.length))
  ll <- function(tr) {
    pruning_jc(tr$edge[, 1], tr$edge[, 2], tr$edge.length, fx$states, 3L,
               rep(1 / 3, 3))
  }
  expect_equal(ll(rooted), ll(unrooted), tolerance = 1e-12)
})

test_that("lumped JC(k') equals full JC(k')", {
  fx <- .HeaderFixture(k = 2)
  tr <- fx$tree
  full <- pruning_jc(tr$edge[, 1], tr$edge[, 2], tr$edge.length, fx$states,
                     5L, rep(0.2, 5))
  lumped <- pruning_jc_collapsed(tr$edge[, 1], tr$edge[, 2], tr$edge.length,
                                 fx$states, 5L, 2L)
  expect_equal(lumped, full, tolerance = 1e-12)
})

test_that("constant patterns give the ascertainment probability", {
  fx <- .HeaderFixture()
  tr <- fx$tree
  edge <- tr$edge
  nTip <- NTip(tr)
  rates <- c(0.6, 1.4)
  pats <- header_constant_patterns(nTip, 3L, logical(0))
  expect_identical(pats, matrix(rep(0:2, each = nTip), nTip))
  pConst <- sum(header_site_likelihoods(
    "jc", edge[, 1], edge[, 2], tr$edge.length, pats, 3L, 0, numeric(0), rates
  ))
  expect_equal(pConst, constant_site_prob_jc(
    edge[, 1], edge[, 2], tr$edge.length, nTip, 3L, rep(1 / 3, 3), rates
  ), tolerance = 1e-13)

  miss <- seq_len(nTip) %in% c(2, 5)
  pats <- header_constant_patterns(nTip, 1L, miss)
  expect_identical(pats[, 1], ifelse(miss, -1L, 0L))
  pMasked <- 3 * sum(header_site_likelihoods(
    "jc", edge[, 1], edge[, 2], tr$edge.length, pats, 3L, 0, numeric(0), rates
  ))
  expect_equal(pMasked, asc_site_prob_missing(
    edge[, 1], edge[, 2], tr$edge.length, nTip, 3L, FALSE, 1, rates, miss,
    FALSE
  ), tolerance = 1e-13)
})

test_that("cached likelihoods equal full recomputation; reject restores", {
  fx <- .HeaderFixture(nTip = 12)
  tr <- fx$tree
  nEdge <- nrow(tr$edge)
  set.seed(2290)
  other <- lapply(1:3, function(i) {
    t2 <- Preorder(RandomTree(tr$tip.label, root = FALSE))
    t2$edge.length <- rexp(nEdge, 6)
    t2 <- RenumberTips(t2, tr$tip.label)
    .TreeList(t2)
  })
  steps <- list(
    list(edge = 4L, length = 0.31, accept = TRUE),
    list(edge = c(1L, nEdge), length = c(0.02, 0.5), accept = FALSE),
    c(other[[1]], accept = TRUE),
    list(edge = 7L, length = 0.11, accept = FALSE),
    c(other[[2]], accept = FALSE),
    list(rates = c(0.5, 1.5), accept = FALSE),
    list(rates = c(0.5, 1.5), accept = TRUE),
    c(other[[3]], list(edge = 2L, length = 0.05), accept = TRUE),
    list(edge = 3L, length = 0.2, accept = TRUE)
  )
  accepted <- vapply(steps, `[[`, logical(1), "accept")
  for (model in c("jc", "jc_lumped", "mkn", "f81")) {
    states <- if (model == "mkn") pmin(fx$states, 1L) else fx$states
    k <- switch(model, mkn = 2L, jc_lumped = 6L, 3L)
    param <- switch(model, mkn = 0.7, f81 = 1.6, 0)
    mSteps <- if (model == "mkn") {
      c(steps, list(list(param = 2.5, accept = FALSE),
                    list(param = 2.5, accept = TRUE)))
    } else {
      steps
    }
    mAccepted <- vapply(mSteps, `[[`, logical(1), "accept")
    trace <- header_pruner_trace(model, states, k, 3L, param,
                                 c(0.4, 0.4, 0.2), 1, .TreeList(tr), mSteps)
    expect_identical(trace$cached, trace$fresh)
    expect_identical(trace$after_reject[!mAccepted], trace$before[!mAccepted])
    # A single-edge change recomputes only its parent's path to the root
    expect_identical(trace$recomputed[1], .PathToRoot(tr, tr$edge[4, 1]))
  }
  expect_true(any(!accepted))
})

test_that("Pruner rejects malformed input", {
  fx <- .HeaderFixture(nTip = 6, nChar = 3)
  fx$states[1, 1] <- 2L
  tr <- fx$tree
  edge <- tr$edge
  expect_error(pruning_jc(edge[, 1], edge[, 2], tr$edge.length,
                          fx$states, 2L, c(0.5, 0.5)),
               "exceeds n_states")
  twoParents <- edge[, 2]
  twoParents[2] <- twoParents[1]
  expect_error(pruning_jc(edge[, 1], twoParents, tr$edge.length,
                          fx$states, 3L, rep(1 / 3, 3)),
               "two parents|not in the tree")
  expect_error(pruning_jc_acrv(edge[, 1], edge[, 2], tr$edge.length,
                               fx$states, 3L, rep(1 / 3, 3), numeric(0)),
               "must not be empty")
  expect_identical(pruning_jc(edge[, 1], edge[, 2], tr$edge.length,
                              fx$states[, 0, drop = FALSE], 3L, rep(1 / 3, 3)),
                   0)
})
