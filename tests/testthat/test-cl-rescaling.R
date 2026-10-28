# Conditional-likelihood rescaling (#346). A near-random character with many
# states on a few hundred tips has a log-likelihood below -745, where an
# unscaled CL vector underflows to zero and the whole state scores -Inf.
# Every kernel is checked against an independent, renormalising R pruner.

# Root CLs (nChar x k) by pruning with every node renormalised, and their
# log scale; `Pt(t)` gives the transition matrix.
.ScaledRootCl <- function(tree, states, Pt, k, rate) {
  edge <- tree$edge
  nTip <- NTip(tree)
  nNode <- max(edge)
  cl <- lapply(seq_len(nNode), function(n) matrix(1, ncol(states), k))
  for (i in seq_len(nTip)) {
    known <- states[i, ] >= 0
    cl[[i]][known, ] <- 0
    cl[[i]][cbind(which(known), states[i, known] + 1)] <- 1
  }
  logS <- matrix(0, nNode, ncol(states))
  for (e in rev(seq_len(nrow(edge)))) {
    p <- edge[e, 1]
    ch <- edge[e, 2]
    cl[[p]] <- cl[[p]] * (cl[[ch]] %*% t(Pt(tree$edge.length[e] * rate)))
    m <- apply(cl[[p]], 1, max)
    cl[[p]] <- cl[[p]] / m
    logS[p, ] <- logS[p, ] + logS[ch, ] + log(m)
  }
  # Return:
  list(cl = cl[[nTip + 1]], logS = logS[nTip + 1, ])
}

# Per-character log-likelihoods, averaged over ACRV multipliers `rates`.
.LogPruneRef <- function(tree, states, Pt, rootFreq, rates = 1) {
  catLL <- vapply(rates, function(rate) {
    root <- .ScaledRootCl(tree, states, Pt, length(rootFreq), rate)
    log(root$cl %*% rootFreq)[, 1] + root$logS
  }, numeric(ncol(states)))
  catLL <- matrix(catLL, ncol = length(rates))
  mx <- apply(catLL, 1, max)
  # Return:
  mx + log(rowMeans(exp(catLL - mx)))
}

# Posterior state probabilities at the root of the first character.
.RootPosterior <- function(tree, states, Pt, rootFreq, rates) {
  roots <- lapply(rates, function(rate) {
    .ScaledRootCl(tree, states, Pt, length(rootFreq), rate)
  })
  logS <- vapply(roots, function(r) r$logS[1], double(1))
  post <- Reduce(`+`, lapply(seq_along(roots), function(i) {
    roots[[i]]$cl[1, ] * rootFreq * exp(logS[i] - max(logS))
  }))
  # Return:
  post / sum(post)
}

.JcP <- function(k) {
  function(t) {
    e <- exp(-k * t / (k - 1))
    matrix((1 - e) / k, k, k) + diag(e, k)
  }
}

.MknP <- function(rateLoss) {
  pi0 <- rateLoss / (1 + rateLoss)
  lambda <- (1 + rateLoss)^2 / (2 * rateLoss)
  function(t) {
    e <- exp(-lambda * t)
    matrix(c(pi0, pi0, 1 - pi0, 1 - pi0), 2) +
      e * matrix(c(1 - pi0, -pi0, pi0 - 1, pi0), 2)
  }
}

# The issue's case: one near-random 8-state character on 200 tips, total
# tree length 20, followed by ordinary characters.
.UnderflowCase <- function(nTip = 200L, nState = 8L, nOrdinary = 3L,
                           seed = 346L, treeLength = 20,
                           shape = RandomTree) {
  set.seed(seed)
  tree <- Preorder(shape(nTip))
  tree$edge.length <- rexp(nrow(tree$edge))
  tree$edge.length <- tree$edge.length * treeLength / sum(tree$edge.length)
  ordinary <- matrix(rep(0:1, each = nTip / 2), nTip, nOrdinary)
  states <- cbind(sample(0:(nState - 1), nTip, replace = TRUE), ordinary)
  states[sample(length(states), 10)] <- -1L
  storage.mode(states) <- "integer"
  dimnames(states) <- list(tree$tip.label, NULL)
  list(tree = tree, states = states, p = tree$edge[, 1],
       ch = tree$edge[, 2], el = tree$edge.length)
}

test_that("a near-random character no longer underflows to -Inf (#346)", {
  x <- .UnderflowCase()
  rates <- c(0.4, 1.6)
  for (k in c(16L, 64L)) {
    ref <- .LogPruneRef(x$tree, x$states, .JcP(k), rep(1 / k, k))
    expect_lt(ref[1], -745)
    refAcrv <- .LogPruneRef(x$tree, x$states, .JcP(k), rep(1 / k, k), rates)

    # Public header (MkPrime::Pruner)
    expect_equal(pruning_jc(x$p, x$ch, x$el, x$states, k, rep(1 / k, k)),
                 sum(ref), tolerance = 1e-12)
    expect_equal(pruning_jc_collapsed(x$p, x$ch, x$el, x$states, k, 8L),
                 sum(ref), tolerance = 1e-12)
    expect_equal(pruning_jc_acrv(x$p, x$ch, x$el, x$states, k, rep(1 / k, k),
                                 rates), sum(refAcrv), tolerance = 1e-12)
    expect_equal(pruning_jc_acrv_collapsed(x$p, x$ch, x$el, x$states, k, 8L,
                                           rates),
                 sum(refAcrv), tolerance = 1e-12)
    expect_equal(
      header_site_likelihoods("jc", x$p, x$ch, x$el, x$states, k, 0, 0, rates,
                              log = TRUE),
      refAcrv, tolerance = 1e-12
    )

    # Flat-buffer MCMC kernels
    expect_equal(test_flat_jc_uncollapsed(x$p, x$ch, x$el, x$states, k, FALSE,
                                          1), sum(ref), tolerance = 1e-12)
    expect_equal(test_flat_jc_uncollapsed(x$p, x$ch, x$el, x$states, k, TRUE,
                                          rates), sum(refAcrv),
                 tolerance = 1e-12)
    expect_equal(test_flat_jc_collapsed(x$p, x$ch, x$el, x$states, k, 8L,
                                        FALSE, 1), sum(ref), tolerance = 1e-12)
    expect_equal(test_flat_jc_collapsed(x$p, x$ch, x$el, x$states, k, 8L,
                                        TRUE, rates), sum(refAcrv),
                 tolerance = 1e-12)

    # Per-site kernels of the k' sweep and marginal_k
    expect_equal(test_persite_uncollapsed(x$p, x$ch, x$el, x$states, k,
                                          rates), refAcrv, tolerance = 1e-12)
    expect_equal(test_persite_collapsed(x$p, x$ch, x$el, x$states, k, 8L,
                                        rates), refAcrv, tolerance = 1e-12)
  }
})

test_that("MkN and F81 kernels rescale", {
  set.seed(2041)
  nTip <- 1200L
  tree <- Preorder(RandomTree(nTip, root = FALSE))
  tree$edge.length <- rexp(nrow(tree$edge)) * 2
  states <- matrix(sample(0:1, nTip * 2, replace = TRUE), nTip,
                   dimnames = list(tree$tip.label, NULL))
  p <- tree$edge[, 1]
  ch <- tree$edge[, 2]
  el <- tree$edge.length
  rates <- c(0.5, 1.5)
  rl <- 1.7
  freqs <- c(rl, 1) / (1 + rl)
  ref <- .LogPruneRef(tree, states, .MknP(rl), freqs)
  expect_lt(ref[1], -745)
  expect_equal(pruning_mkn(p, ch, el, states, rl, freqs), sum(ref),
               tolerance = 1e-12)
  expect_equal(
    pruning_mkn_acrv(p, ch, el, states, rl, freqs, rates),
    sum(.LogPruneRef(tree, states, .MknP(rl), freqs, rates)),
    tolerance = 1e-12
  )

  freqF81 <- c(0.7, 0.3)
  mu <- 1 / (1 - sum(freqF81^2))
  F81P <- function(t) {
    e <- exp(-mu * t)
    e * diag(2) + (1 - e) * matrix(freqF81, 2, 2, byrow = TRUE)
  }
  expect_equal(
    header_site_likelihoods("f81", p, ch, el, states, 2L, mu, freqF81, rates,
                            log = TRUE),
    .LogPruneRef(tree, states, F81P, freqF81, rates), tolerance = 1e-12
  )

  # Flat kernels behind MkpLogLikelihood()
  mkd <- MkPrimeData(MatrixToPhyDat(states), neomorphic = 1:2)
  expect_equal(
    MkpLogLikelihood(tree, mkd, rate_loss = rl, coding = "none"),
    sum(ref), tolerance = 1e-12
  )
  acrvRates <- header_lognormal_rates(0.8, qnorm((1:4 - 0.5) / 4))
  expect_equal(
    MkpLogLikelihood(tree, mkd, rate_loss = rl, rate_log_sd = 0.8, nCat = 4L,
                     coding = "none"),
    sum(.LogPruneRef(tree, states, .MknP(rl), freqs, acrvRates)),
    tolerance = 1e-12
  )
})

# Pins every character's k, so that the likelihood has no relabelling or
# ascertainment term and the partial-CL evaluators can be compared with a
# plain full evaluation.
.UnderflowChain <- function(model = MkPrimeModel(coding = "none"), ...) {
  x <- .UnderflowCase(...)
  mkd <- MkPrimeData(MatrixToPhyDat(x$states),
                     knownStates = setNames(c(16L, 2L, 2L, 2L), 1:4))
  model <- MkPrime:::.FinalizeModel(model, x$tree, mkd)
  state0 <- MkPrime:::.InitState(x$tree, mkd, model)
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  list(x = x, mkd = mkd, dataPtr = dataPtr, statePtr = statePtr)
}

test_that("the node CL cache keeps its shifts across partial updates", {
  for (model in list(MkPrimeModel(coding = "none", nCat = 1L),
                     MkPrimeModel(coding = "none", nCat = 4L))) {
    chain <- .UnderflowChain(model)
    cold <- .ColdLogLik(chain$dataPtr, chain$statePtr)
    expect_true(is.finite(cold))
    expect_lt(cold, -745)
    expect_equal(get_state_log_lik(chain$statePtr), cold, tolerance = 1e-12)
    set.seed(931)
    nMoved <- 0L
    for (i in seq_len(40)) {
      before <- get_state_log_lik(chain$statePtr)
      move <- if (i %% 2) 5L else 4L  # NNI, beta_simplex
      do_move_cpp(chain$dataPtr, chain$statePtr, move, 0L, 0.5, 2, 1L, 1.0)
      expect_equal(get_state_log_lik(chain$statePtr),
                   .ColdLogLik(chain$dataPtr, chain$statePtr),
                   tolerance = 1e-12)
      nMoved <- nMoved + (get_state_log_lik(chain$statePtr) != before)
    }
    expect_gt(nMoved, 0L)
    counts <- get_mcmc_state(chain$statePtr)
    expect_gt(counts$diagNniPartial, 0L)
    expect_gt(counts$diagBsPartial, 0L)
  }
})

test_that("Gibbs partial-CL evaluators rescale", {
  chain <- .UnderflowChain(MkPrimeModel(coding = "none", nCat = 2L))
  st <- get_mcmc_state(chain$statePtr)
  edge <- st$edge
  len <- st$treeLength * st$relBrLengths
  nTip <- NTip(chain$x$tree)

  pruneRow <- which(edge[, 1] != nTip + 1L & edge[, 2] > nTip)[3]
  en <- gibbs_spr_enumerate_cpp(chain$dataPtr, chain$statePtr, pruneRow)
  expect_true(all(is.finite(en$logLik)))
  expect_lt(max(en$logLik), -745)
  # The merged edge is the subtree's own position, split in half
  parentRow <- which(edge[, 2] == en$u)
  sibRow <- which(edge[, 1] == en$u & edge[, 2] != en$v)
  merged <- len
  merged[c(parentRow, sibRow)] <- en$lMerge / 2
  expect_equal(
    en$logLik[length(en$logLik)],
    cpp_log_likelihood_xptr(chain$dataPtr, edge[, 1], edge[, 2], merged,
                            st$kPrime, st$rateLoss, st$rateLogSd, st$rateNeo,
                            st$p, st$kprimeAlpha, st$kprimeBeta,
                            st$betaScale),
    tolerance = 1e-12
  )

  # Ascertainment pseudo-characters take the same partial-CL path. Their
  # cached CLs are rescaled, so only a long path from the regraft point to
  # the root underflows: prune the deepest clade of a pectinate tree
  variable <- .UnderflowChain(MkPrimeModel(coding = "variable", nCat = 1L),
                              nTip = 400L, treeLength = 800,
                              shape = PectinateTree)
  vst <- get_mcmc_state(variable$statePtr)
  vEdge <- vst$edge
  vRow <- max(which(vEdge[, 2] > 400L))
  ven <- gibbs_spr_enumerate_cpp(variable$dataPtr, variable$statePtr, vRow)
  expect_true(all(is.finite(ven$logLik)))
  vMerged <- vst$treeLength * vst$relBrLengths
  vMerged[c(which(vEdge[, 2] == ven$u),
            which(vEdge[, 1] == ven$u & vEdge[, 2] != ven$v))] <- ven$lMerge / 2
  expect_equal(
    ven$logLik[length(ven$logLik)],
    cpp_log_likelihood_xptr(variable$dataPtr, vEdge[, 1], vEdge[, 2], vMerged,
                            vst$kPrime, vst$rateLoss, vst$rateLogSd,
                            vst$rateNeo, vst$p, vst$kprimeAlpha,
                            vst$kprimeBeta, vst$betaScale),
    tolerance = 1e-12
  )

  for (nodeA in edge[c(5, 60, 150), 2]) {
    res <- validate_swap_partial_cl(chain$dataPtr, chain$statePtr, nodeA)
    expect_gt(nrow(res), 0)
    expect_true(all(is.finite(res$ll_partial)))
    expect_equal(res$ll_partial, res$ll_full, tolerance = 1e-12)
  }
})

test_that("Q-heterogeneity kernel rescales", {
  x <- .UnderflowCase()
  mkd <- MkPrimeData(MatrixToPhyDat(x$states),
                     knownStates = setNames(c(16L, 2L, 2L, 2L), 1:4))
  model <- MkPrime:::.FinalizeModel(
    MkPrimeModel(coding = "none", qHeterogeneity = TRUE), x$tree, mkd
  )
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)
  k <- c(16L, 2L, 2L, 2L)
  # As beta_scale grows, every F81 component tends to JC
  ref <- sum(vapply(seq_along(k), function(i) {
    .LogPruneRef(x$tree, x$states[, i, drop = FALSE], .JcP(k[i]),
                 rep(1 / k[i], k[i]))
  }, double(1)))
  expect_equal(
    cpp_log_likelihood_xptr(dataPtr, x$p, x$ch, x$el, k, 1, 0, 1, 0.5, 1, 1,
                            betaScale = 1e9),
    ref, tolerance = 1e-6
  )
})

test_that("Q-het Gibbs topology moves stop rather than drop a character", {
  chain <- .UnderflowChain(
    MkPrimeModel(coding = "none", nCat = 1L, qHeterogeneity = TRUE)
  )
  expect_true(is.finite(get_state_log_lik(chain$statePtr)))
  for (move in c(10L, 11L)) {
    expect_error(
      do_move_cpp(chain$dataPtr, chain$statePtr, move, 0L, 0.5, 2, 1L, 1.0),
      "Q-heterogeneity Gibbs topology moves"
    )
  }
})

test_that("Pruner restores shifts on reject() and keeps marginals finite", {
  x <- .UnderflowCase(nTip = 120L)
  k <- 32L
  tr <- list(parent = x$p, child = x$ch, edge_length = x$el)
  steps <- list(
    list(edge = c(3L, 50L), length = c(0.4, 0.01), accept = FALSE),
    list(edge = 7L, length = 2, accept = TRUE),
    list(rates = c(0.3, 1, 1.7), accept = FALSE),
    list(rates = c(0.5, 1.5), accept = TRUE),
    list(edge = 100L, length = 0.2, accept = FALSE)
  )
  trace <- header_pruner_trace("jc", x$states, k, 0L, 0, 0, 1, tr, steps)
  expect_lt(trace$fresh[1], -745)
  expect_equal(trace$cached, trace$fresh, tolerance = 1e-12)
  rejected <- !is.na(trace$after_reject)
  expect_equal(trace$after_reject[rejected], trace$before[rejected],
               tolerance = 1e-12)

  rates <- c(0.5, 1.5)
  marg <- header_marginals("jc", x$p, x$ch, x$el, x$states, k, 0, 0, rates,
                           numeric(0))$marginals
  expect_true(all(is.finite(marg)))
  expect_equal(apply(marg, 2:3, sum), matrix(1, ncol(x$states), dim(marg)[3]),
               tolerance = 1e-12)
  # At the root the upward vector is the root frequencies
  root <- NTip(x$tree) + 1L
  expect_equal(
    marg[, 1, root],
    .RootPosterior(x$tree, x$states[, 1, drop = FALSE], .JcP(k),
                   rep(1 / k, k), rates),
    tolerance = 1e-10
  )
})

test_that("Pruner rescales within a polytomy", {
  set.seed(3)
  nSub <- 6L
  k <- 64L
  subs <- lapply(seq_len(nSub), function(i) {
    sub <- ape::rtree(40, tip.label = paste0("s", i, "_", 1:40))
    sub$edge.length <- rep(5, nrow(sub$edge))
    sub$root.edge <- 0.01
    sub
  })
  tree <- ape::stree(nSub)
  tree$edge.length <- rep(0.01, nSub)
  tree$tip.label <- paste0("X", seq_len(nSub))
  for (i in seq_len(nSub)) {
    tree <- ape::bind.tree(tree, subs[[i]],
                           where = which(tree$tip.label == paste0("X", i)))
  }
  tree <- Preorder(tree)
  states <- matrix(sample(0:(k - 1), NTip(tree), replace = TRUE),
                   dimnames = list(tree$tip.label, NULL))
  storage.mode(states) <- "integer"
  expect_equal(sum(tree$edge[, 1] == NTip(tree) + 1L), nSub)
  expect_equal(
    header_site_likelihoods("jc", tree$edge[, 1], tree$edge[, 2],
                            tree$edge.length, states, k, 0, 0, 1, log = TRUE),
    .LogPruneRef(tree, states, .JcP(k), rep(1 / k, k)),
    tolerance = 1e-12
  )
})

test_that("MCMC workspace kernels rescale MkN characters", {
  set.seed(2041)
  nTip <- 1200L
  # A pectinate tree visits an internal child first at most nodes; the small
  # tree needs no rescaling. The "variable" block below reuses the last tree
  for (tree in list(RandomTree(20L, root = FALSE),
                    RandomTree(nTip, root = FALSE), PectinateTree(nTip))) {
    tree <- Preorder(tree)
    tree$edge.length <- rexp(nrow(tree$edge)) * 2
    states <- matrix(sample(0:1, NTip(tree), replace = TRUE), NTip(tree),
                     dimnames = list(tree$tip.label, NULL))
    mkd <- MkPrimeData(MatrixToPhyDat(states), neomorphic = 1L)
    for (nCat in c(1L, 4L)) {
      model <- MkPrime:::.FinalizeModel(
        MkPrimeModel(coding = "none", nCat = nCat), tree, mkd
      )
      state0 <- MkPrime:::.InitState(tree, mkd, model)
      if (nCat == 1L) state0$rate_log_sd <- 0
      dataPtr <- MkPrime:::.InitMcmcData(mkd, model)
      statePtr <- MkPrime:::.InitMcmcChain(state0)
      fill_partition_cache(dataPtr, statePtr)
      allocate_cl_workspace(dataPtr, statePtr)
      st <- get_mcmc_state(statePtr)
      edged <- tree
      edged$edge <- st$edge
      edged$edge.length <- st$treeLength * st$relBrLengths
      freqs <- c(st$rateLoss, 1) / (1 + st$rateLoss)
      rates <- if (nCat == 1L) 1 else {
        header_lognormal_rates(st$rateLogSd, qnorm((1:nCat - 0.5) / nCat))
      }
      ref <- .LogPruneRef(edged, states, .MknP(st$rateLoss), freqs, rates)
      expect_equal(ref < -745, NTip(tree) == nTip)
      expect_equal(eval_full_loglik_cpp(dataPtr, statePtr), ref,
                   tolerance = 1e-12)
    }
  }

  # Under "variable" coding the constant-pattern pseudo-characters underflow
  # too, so the Gibbs evaluators take their rescaling pass on them
  model <- MkPrime:::.FinalizeModel(
    MkPrimeModel(coding = "variable", nCat = 1L), tree, mkd
  )
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(MkPrime:::.InitState(tree, mkd, model))
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  edge <- get_mcmc_state(statePtr)$edge
  pruneRow <- which(edge[, 1] != nTip + 1L & edge[, 2] > nTip)[2]
  expect_true(all(is.finite(
    gibbs_spr_enumerate_cpp(dataPtr, statePtr, pruneRow)$logLik
  )))
})

test_that("Pruner handles unary nodes, zero likelihoods and empty categories", {
  x <- .UnderflowCase()
  k <- 16L
  ref <- .LogPruneRef(x$tree, x$states, .JcP(k), rep(1 / k, k), 1.5)
  # A zero-rate category scores every variable character 0
  expect_equal(
    header_site_likelihoods("jc", x$p, x$ch, x$el, x$states, k, 0, 0,
                            c(1.5, 0), log = TRUE)[1],
    ref[1] + log(0.5), tolerance = 1e-12
  )
  # Node 5 has a single child; with lengths summed over it, the tree is
  # the same as one without it. Tips 1 and 2, on zero-length edges, differ
  # in the second character, whose likelihood is therefore 0.
  tips <- matrix(c(0L, 0L, 1L, 0L, 1L, 1L), 3)
  for (model in c("jc", "mkn", "f81")) {
    param <- if (model == "jc") 0 else 2
    pi <- if (model == "f81") c(0.3, 0.7) else 0
    unary <- header_site_likelihoods(model, c(4L, 5L, 6L, 6L, 4L),
                                     c(5L, 6L, 1L, 2L, 3L),
                                     c(0.3, 0.2, 0, 0, 0.4), tips, 2L,
                                     param, pi, 1, log = TRUE)
    merged <- header_site_likelihoods(model, c(4L, 5L, 5L, 4L),
                                      c(5L, 1L, 2L, 3L), c(0.5, 0, 0, 0.4),
                                      tips, 2L, param, pi, 1, log = TRUE)
    expect_equal(unary, merged, tolerance = 1e-14)
    expect_identical(unary[2], -Inf)
  }
})

test_that("kernels rescale a first-visited child and many observed states", {
  set.seed(1311)
  nTip <- 200L
  tree <- Preorder(PectinateTree(nTip))
  tree$edge.length <- rexp(nrow(tree$edge))
  tree$edge.length <- tree$edge.length * 20 / sum(tree$edge.length)
  states <- matrix(sample(0:12, nTip, replace = TRUE), nTip,
                   dimnames = list(tree$tip.label, NULL))
  p <- tree$edge[, 1]
  ch <- tree$edge[, 2]
  el <- tree$edge.length
  k <- 64L
  rates <- c(0.4, 1.6)
  ref <- .LogPruneRef(tree, states, .JcP(k), rep(1 / k, k), rates)
  expect_lt(ref, -745)
  expect_equal(test_persite_uncollapsed(p, ch, el, states, k, rates), ref,
               tolerance = 1e-12)
  expect_equal(test_persite_collapsed(p, ch, el, states, k, 13L, rates), ref,
               tolerance = 1e-12)
  expect_equal(test_persite_uncollapsed(p, ch, el, states, 16L, rates),
               .LogPruneRef(tree, states, .JcP(16L), rep(1 / 16, 16), rates),
               tolerance = 1e-12)
})
