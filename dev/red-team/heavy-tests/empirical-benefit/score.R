# Scoring for the empirical-vs-geometric k' prior study (#179).
#
# .ExactKPosterior() gives each character's posterior on k' on a known tree
# with no sampler: per-k' likelihoods are enumerated over [2, K] and p is
# integrated on a grid. It isolates what the prior alone does to k', and is the
# yardstick the MCMC arm is read against. Needs MkPrime loaded.

.Mkp <- function(name) get(name, envir = asNamespace("MkPrime"))

# log Pr(y_j | k) for each column of `states` and k = 2..K (column k - 1);
# -Inf where k < kObs. Rows of `states` follow the tips' node numbers.
.PerCharLogLik <- function(tree, states, kObs, K) {
  tree <- TreeTools::Preorder(tree)
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  el <- tree$edge.length
  storage.mode(states) <- "integer"
  out <- matrix(-Inf, ncol(states), K - 1L)
  for (j in seq_len(ncol(states))) {
    y <- states[, j, drop = FALSE]
    kj <- kObs[[j]]
    out[j, kj - 1L] <- .Mkp("pruning_jc")(parent, child, el, y, kj,
                                          rep(1 / kj, kj))
    for (k in seq.int(kj + 1L, length.out = K - kj)) {
      out[j, k - 1L] <- .Mkp("pruning_jc_collapsed")(parent, child, el, y,
                                                     k, kj)
    }
    out[j, ] <- out[j, ] +
      c(rep(0, kj - 2L), .Mkp("mk_prime_relabel_log_batch")(kj:K, kj))
  }
  # Return:
  out
}

# log Pr(V | k), k = 2..K, for a character with no missing data.
.LogPVariable <- function(tree, K) {
  tree <- TreeTools::Preorder(tree)
  nTip <- length(tree$tip.label)
  const <- .Mkp("lumped_const_probs")(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length, nTip, 1,
    matrix(0L, nTip, 1L), as.numeric(2:K)
  )
  # Return:
  log1p(-drop(const))
}

.LogSumExpRows <- function(x) {
  mx <- apply(x, 1, max)
  # Return:
  mx + log(rowSums(exp(x - mx)))
}

# Posterior on k' (rows: characters; columns: k = 2..K) under `arm`, with p
# integrated over a logit-spaced grid under its Beta(a, b) hyperprior. Exact
# for a known tree, up to the grid. pPost is the posterior on that grid.
.ExactKPosterior <- function(logLik, logPV, arm, emp, K, a = 1, b = 1,
                             pGrid = stats::plogis(seq(-8, 8, by = 0.05))) {
  model <- list(kPrimePrior = arm, kprimeTruncK = K, empiricalNObs = emp)
  nChar <- nrow(logLik)
  logW <- stats::dbeta(pGrid, a, b, log = TRUE) + log(pGrid * (1 - pGrid))
  perP <- lapply(pGrid, function(p) {
    logPrior <- .Mkp(".KprimeLogPriorTable")(model, p, 1, 1)
    joint <- sweep(logLik, 2, logPrior, `+`)
    lse <- .LogSumExpRows(joint)
    logD <- log(sum(exp(logPrior + logPV)))
    list(post = exp(joint - lse), logMarg = sum(lse) - nChar * logD)
  })
  logPP <- logW + vapply(perP, `[[`, double(1), "logMarg")
  pPost <- exp(logPP - max(logPP))
  pPost <- pPost / sum(pPost)
  kPost <- Reduce(`+`, Map(function(x, w) w * x$post, perP, pPost))
  # Return:
  list(kPost = kPost, pPost = pPost, pGrid = pGrid)
}

# Posterior on k' from MCMC draws; `draws` has one column per character.
.SampledKPosterior <- function(draws, K) {
  counts <- t(apply(draws, 2, function(x) tabulate(x, nbins = K)[-1L]))
  # Return:
  counts / rowSums(counts)
}

# Per-character scores of a k' posterior against the truth:
# logScore = log Pr(k' = kTrue), cover = kTrue in the 90% HPD set,
# err = posterior mean - kTrue. `floor` keeps a kTrue that no draw visited
# finite; pass half a draw's mass for sampled posteriors.
.ScoreK <- function(kPost, kTrue, level = 0.9, floor = 0) {
  ks <- seq_len(ncol(kPost)) + 1L
  cover <- vapply(seq_len(nrow(kPost)), function(i) {
    o <- order(kPost[i, ], decreasing = TRUE)
    inSet <- o[seq_len(which(cumsum(kPost[i, o]) >= level - 1e-12)[1])]
    (kTrue[[i]] - 1L) %in% inSet
  }, logical(1))
  # Return:
  data.frame(
    logScore = log(pmax(kPost[cbind(seq_along(kTrue), kTrue - 1L)], floor)),
    cover = cover,
    err = drop(kPost %*% ks) - kTrue
  )
}

# Mean distances of posterior trees to the truth. Both are split-based, so
# blind to rooting, unlike ape::prop.part.
.ScoreTrees <- function(trees, trueTree) {
  class(trees) <- "multiPhylo"
  # Return:
  c(
    cid = mean(TreeDist::ClusteringInfoDistance(trees, trueTree,
                                                normalize = TRUE)),
    rf = mean(TreeDist::RobinsonFoulds(trees, trueTree, normalize = TRUE))
  )
}
