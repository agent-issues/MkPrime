# beta = 0 invariance harness for tree and branch-length moves.
#
# At beta = 0 the target is the prior, which is exactly i.i.d.-sampleable:
# a uniform labelled unrooted topology, flat Dirichlet edge fractions, and --
# because the moves see the storage root -- a uniformly chosen internal node as
# the trifurcation.  Starting every replicate from an exact draw tests
# pi K^n = pi with no burn-in or mixing assumption.

.PriorDraw <- function(tipLabel) {
  nTip <- length(tipLabel)
  tree <- TreeTools::RandomTree(tipLabel, root = FALSE)
  tree <- TreeTools::Preorder(ape::root(
    tree, node = nTip + sample.int(nTip - 2L, 1L), resolve.root = FALSE))
  f <- rexp(nrow(tree$edge))
  tree$edge.length <- f / sum(f)
  tree
}

# Smallest, largest and summed-internal edge fraction; pendant fraction of t1;
# whether t1 hangs from the storage root; number of cherries.
.PriorStats <- function(edge, f) {
  nTip <- min(edge[, 1]) - 1L
  c(min(f), max(f), sum(f[edge[, 2] > nTip]), f[edge[, 2] == 1L],
    edge[edge[, 2] == 1L, 1] == nTip + 1L,
    sum(tabulate(edge[edge[, 2] <= nTip, 1], 2L * nTip) == 2L))
}

# p-values comparing .PriorStats after `nMove` applications of `moveType` at
# beta = 0 with ten times as many direct draws: KS for the four continuous
# statistics, chi-squared for the two discrete ones.
.PriorInvarianceP <- function(moveType, nMove, nReps, nTip = 6L) {
  tipLabel <- paste0("t", seq_len(nTip))
  mat <- matrix(sample(0:1, nTip * 4L, replace = TRUE), nrow = nTip,
                dimnames = list(tipLabel, NULL))
  mkd <- suppressWarnings(MkPrimeData(TreeTools::MatrixToPhyDat(mat)))
  model <- MkPrimeModel()

  sampled <- t(vapply(seq_len(nReps), function(r) {
    tree <- .PriorDraw(tipLabel)
    finalModel <- MkPrime:::.FinalizeModel(model, tree, mkd)
    dataPtr  <- MkPrime:::.InitMcmcData(mkd, finalModel)
    statePtr <- MkPrime:::.InitMcmcChain(
      MkPrime:::.InitState(tree, mkd, finalModel))
    fill_partition_cache(dataPtr, statePtr)
    allocate_cl_workspace(dataPtr, statePtr)
    for (i in seq_len(nMove))
      do_move_cpp(dataPtr, statePtr, moveType, 0L, 0.5, 0.5, 1L, 0.0)
    st <- get_mcmc_state(statePtr)
    .PriorStats(st$edge, st$relBrLengths)
  }, numeric(6L)))

  reference <- t(vapply(seq_len(nReps * 10L), function(i) {
    tree <- .PriorDraw(tipLabel)
    .PriorStats(tree$edge, tree$edge.length)
  }, numeric(6L)))

  c(
    vapply(1:4, function(j) suppressWarnings(
      ks.test(sampled[, j], reference[, j])$p.value), numeric(1)),
    vapply(5:6, function(j) {
      lev <- sort(unique(c(sampled[, j], reference[, j])))
      if (length(lev) < 2L) return(1)
      chisq.test(rbind(table(factor(sampled[, j], lev)),
                       table(factor(reference[, j], lev))))$p.value
    }, numeric(1))
  )
}
