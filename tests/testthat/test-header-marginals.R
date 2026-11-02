# Header API v2: transposed kernels, edge-aware models and the upward pass
# (Pruner::update_marginals / node_marginals), via src/header_api.cpp.

.F81P <- function(pi, mu) {
  function(t) {
    e <- exp(-mu * t)
    e * diag(length(pi)) + (1 - e) * matrix(pi, length(pi), length(pi),
                                             byrow = TRUE)
  }
}

# Posterior state probabilities at every node by enumerating the states of
# all internal nodes: k x nChar x nNode, averaged over rate categories in
# proportion to their likelihoods.
.BruteMarginals <- function(tree, states, k, PMatrix, rootFreqs, rates = 1) {
  edge <- tree$edge
  nTip <- NTip(tree)
  nNode <- max(edge)
  internal <- (nTip + 1):nNode
  root <- setdiff(edge[, 1], edge[, 2])
  configs <- as.matrix(expand.grid(rep(list(seq_len(k)), length(internal))))
  out <- array(0, c(k, ncol(states), nNode))
  for (ch in seq_len(ncol(states))) {
    acc <- matrix(0, k, nNode)
    for (rate in rates) {
      Ps <- lapply(tree$edge.length * rate, PMatrix)
      for (i in seq_len(nrow(configs))) {
        s <- integer(nNode)
        s[internal] <- configs[i, ]
        w <- rootFreqs[s[root]]
        tipVec <- vector("list", nTip)
        for (e in seq_len(nrow(edge))) {
          p <- edge[e, 1]
          ch2 <- edge[e, 2]
          if (ch2 > nTip) {
            w <- w * Ps[[e]][s[p], s[ch2]]
          } else {
            obs <- states[ch2, ch]
            cl <- if (obs < 0) rep(1, k) else as.numeric(seq_len(k) == obs + 1)
            v <- Ps[[e]][s[p], ] * cl
            tipVec[[ch2]] <- v / sum(v)
            w <- w * sum(v)
          }
        }
        for (n in internal) acc[s[n], n] <- acc[s[n], n] + w
        for (tip in seq_len(nTip)) acc[, tip] <- acc[, tip] + w * tipVec[[tip]]
      }
    }
    out[, ch, ] <- sweep(acc, 2, colSums(acc), "/")
  }
  # Return:
  out
}

.MargFixture <- function(seed = 4418) {
  set.seed(seed)
  tree <- Preorder(RootTree(RandomTree(5, root = FALSE), 1))
  tree$edge.length <- rexp(nrow(tree$edge), 3)
  states <- matrix(sample.int(3, 5 * 6, replace = TRUE) - 1L, 5)
  states[c(2, 9, 17)] <- -1L
  list(tree = tree, states = states)
}

test_that("transition_transpose() applies P(t)'", {
  expect_identical(header_api_version(), 3L)
  pi <- c(0.2, 0.5, 0.3)
  for (spec in list(list("jc", 4L, 0), list("mkn", 2L, 0.7),
                    list("f81", 3L, 1.3))) {
    m <- header_transition_matrices(spec[[1]], spec[[2]], spec[[3]], pi, 0.37)
    expect_equal(m$Pt, t(m$P), tolerance = 1e-14, label = spec[[1]])
  }
  m <- header_transition_matrices("f81", 3L, 1.3, pi, 0.37)
  expect_equal(m$P, .F81P(pi, 1.3)(0.37), tolerance = 1e-14)
})

test_that("node marginals match brute-force enumeration", {
  fx <- .MargFixture()
  pi <- c(0.2, 0.5, 0.3)
  mu <- 1.3
  rates <- c(0.4, 1.6)
  tr <- fx$tree
  got <- header_marginals("f81", tr$edge[, 1], tr$edge[, 2], tr$edge.length,
                          fx$states, 3L, mu, pi, rates, numeric(0))
  ref <- .BruteMarginals(tr, fx$states, 3L, .F81P(pi, mu), pi, rates)
  expect_equal(got$marginals, ref, tolerance = 1e-12)
  # down(v) . up(v) is the site likelihood at every node.
  expect_equal(got$node_lik,
               matrix(got$site_lik, ncol(fx$states), max(tr$edge)),
               tolerance = 1e-12)
  # Observed tips are certain.
  obs <- which(fx$states[1, ] >= 0)
  expect_equal(got$marginals[cbind(fx$states[1, obs] + 1L, obs, 1L)],
               rep(1, length(obs)))
})

test_that("node marginals hold for unrooted trees and MkN", {
  set.seed(4419)
  tree <- Preorder(RandomTree(6, root = FALSE))
  tree$edge.length <- rexp(nrow(tree$edge), 4)
  states <- matrix(sample.int(2, 6 * 8, replace = TRUE) - 1L, 6)
  got <- header_marginals("mkn", tree$edge[, 1], tree$edge[, 2],
                          tree$edge.length, states, 2L, 0.6, numeric(0), 1,
                          numeric(0))
  expect_equal(got$node_lik, matrix(got$site_lik, 8, max(tree$edge)),
               tolerance = 1e-12)
  expect_equal(apply(got$marginals, 2:3, sum),
               matrix(1, 8, max(tree$edge)), tolerance = 1e-12)
})

test_that("edge-aware models see their edge", {
  fx <- .MargFixture(4420)
  tr <- fx$tree
  pi <- c(0.2, 0.5, 0.3)
  scale <- runif(nrow(tr$edge), 0.5, 2)
  aware <- header_marginals("f81", tr$edge[, 1], tr$edge[, 2], tr$edge.length,
                            fx$states, 3L, 1.3, pi, c(0.7, 1.3), scale)
  plain <- header_marginals("f81", tr$edge[, 1], tr$edge[, 2],
                            tr$edge.length * scale, fx$states, 3L, 1.3, pi,
                            c(0.7, 1.3), numeric(0))
  expect_equal(aware$site_lik, plain$site_lik, tolerance = 1e-12)
  expect_equal(aware$marginals, plain$marginals, tolerance = 1e-12)

  # Renumbering the edges of the same tree moves each kernel, so the cache
  # must recompute.
  perm <- rev(seq_len(nrow(tr$edge)))
  ll <- header_edge_renumber(tr$edge[, 1], tr$edge[, 2], tr$edge.length,
                             fx$states, 3L, scale, perm)
  expect_equal(ll[[2]], ll[[3]], tolerance = 1e-12)
  expect_false(isTRUE(all.equal(ll[[1]], ll[[3]])))
})
