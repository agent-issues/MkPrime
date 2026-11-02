# beta = 1 invariance harness for tree and branch-length moves.
#
# At beta = 0 every selection or bin weight is exp(0) = 1, so the prior gates
# in helper-prior-invariance.R cannot see how the weights enter a Hastings
# ratio.  At beta = 1 the target is a posterior, which is not i.i.d.-
# sampleable, but on five tips it can be integrated: each of the fifteen
# unrooted topologies carries mass E[L(t, f)] over flat Dirichlet edge
# fractions f, estimated by plain Monte Carlo.  The likelihood is invariant to
# the storage root, so that mass is the chain's target however the root
# wanders.
#
# Chain frequencies are compared with the target by batch-means z, never by a
# chi-squared test on a thinned chain, which measures autocorrelation rather
# than bias.  The z denominator adds the target's own Monte Carlo error.

.PosteriorContext <- function(model = MkPrimeModel(), nTip = 5L, nChar = 12L,
                              states = 0:1) {
  tipLabel <- paste0("t", seq_len(nTip))
  mat <- matrix(sample(states, nTip * nChar, replace = TRUE), nrow = nTip,
                dimnames = list(tipLabel, NULL))
  mkd <- suppressWarnings(MkPrimeData(TreeTools::MatrixToPhyDat(mat)))
  tree <- TreeTools::Preorder(TreeTools::RandomTree(tipLabel, root = FALSE))
  nEdge <- nrow(tree$edge)
  tree$edge.length <- rep(1 / nEdge, nEdge)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(MkPrime:::.InitState(tree, mkd, model))
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  treeLength <- get_mcmc_state(statePtr)$treeLength
  list(
    dataPtr = dataPtr, statePtr = statePtr, tipLabel = tipLabel,
    nEdge = nEdge,
    LogLik = function(edge, f) {
      eval_full_loglik_at_cpp(dataPtr, statePtr, edge[, 1], edge[, 2],
                              treeLength * f)
    }
  )
}

# Posterior mass of each unrooted topology, with its Monte Carlo SE.
.TopologyTarget <- function(ctx, nDraw = 4000L) {
  nTip <- length(ctx$tipLabel)
  topologies <- lapply(
    ape::as.phylo(0:14, nTip, ctx$tipLabel),
    function(t) TreeTools::Preorder(TreeTools::UnrootTree(t))
  )
  keys <- vapply(topologies, function(t) .TopologyKey(t$edge, ctx$tipLabel),
                 "")
  moments <- vapply(topologies, function(t) {
    ll <- vapply(seq_len(nDraw), function(i) {
      f <- rexp(ctx$nEdge)
      ctx$LogLik(t$edge, f / sum(f))
    }, numeric(1))
    w <- exp(ll - max(ll))
    c(max(ll) + log(mean(w)), sd(w) / mean(w) / sqrt(nDraw))
  }, numeric(2))
  target <- exp(moments[1, ] - max(moments[1, ]))
  target <- target / sum(target)
  list(keys = keys, target = target, se = target * moments[2, ])
}

# Batch-means z of each topology's visit frequency over `nIter` applications
# of `moveType` at beta = 1.
.TopologyZ <- function(ctx, moveType, nIter, nBurn = 1000L, nBatch = 50L) {
  tt <- .TopologyTarget(ctx)
  for (i in seq_len(nBurn)) {
    do_move_cpp(ctx$dataPtr, ctx$statePtr, moveType, 0L, 0.5, 0.5, 1L, 1.0)
  }
  visited <- vapply(seq_len(nIter), function(i) {
    do_move_cpp(ctx$dataPtr, ctx$statePtr, moveType, 0L, 0.5, 0.5, 1L, 1.0)
    .TopologyKey(get_mcmc_state(ctx$statePtr)$edge, ctx$tipLabel)
  }, "")
  visited <- factor(visited, tt$keys)
  expect_false(anyNA(visited))
  .BatchZ(vapply(levels(visited), function(k) visited == k,
                 logical(nIter)), tt$target, tt$se, nBatch)
}

# Batch-means z of each column mean of `x` against `target`, whose own
# standard error is `targetSE`.
.BatchZ <- function(x, target, targetSE, nBatch) {
  batch <- rep(seq_len(nBatch), each = nrow(x) %/% nBatch)
  x <- x[seq_along(batch), , drop = FALSE]
  batchMeans <- apply(x, 2, function(col) tapply(col, batch, mean))
  chainSE <- apply(batchMeans, 2, sd) / sqrt(nBatch)
  # Return:
  (colMeans(x) - target) / sqrt(chainSE^2 + targetSE^2)
}

# Batch-means z of each edge's mean fraction on a fixed topology over `nIter`
# applications of `moveType` at beta = 1, against a self-normalized importance
# sample from the flat Dirichlet prior.
.EdgeFractionZ <- function(ctx, moveType, nIter, nBurn = 1000L, nBatch = 50L,
                           nDraw = 40000L) {
  edge <- get_mcmc_state(ctx$statePtr)$edge
  f <- matrix(rexp(nDraw * ctx$nEdge), nDraw)
  f <- f / rowSums(f)
  ll <- apply(f, 1, function(fi) ctx$LogLik(edge, fi))
  w <- exp(ll - max(ll))
  w <- w / sum(w)
  target <- colSums(w * f)
  # Delta-method SE of a self-normalized importance-sampling mean
  targetSE <- sqrt(colSums(w^2 * (f - rep(target, each = nDraw))^2))

  for (i in seq_len(nBurn)) {
    do_move_cpp(ctx$dataPtr, ctx$statePtr, moveType, 0L, 0.5, 0.5, 1L, 1.0)
  }
  sampled <- t(vapply(seq_len(nIter), function(i) {
    do_move_cpp(ctx$dataPtr, ctx$statePtr, moveType, 0L, 0.5, 0.5, 1L, 1.0)
    st <- get_mcmc_state(ctx$statePtr)
    # Fixed-topology moves never renumber, so edge rows keep their identity
    if (!identical(st$edge, edge)) stop("Move ", moveType, " renumbered edges")
    st$relBrLengths
  }, numeric(ctx$nEdge)))
  .BatchZ(sampled, target, targetSE, nBatch)
}
