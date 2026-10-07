# beta = 0 invariance harness for tree and branch-length moves.
#
# At beta = 0 the target is the prior, which is exactly i.i.d.-sampleable:
# a uniform labelled unrooted topology, flat Dirichlet edge fractions, and --
# because the moves see the storage root -- a uniformly chosen internal node as
# the trifurcation.  Starting every replicate from an exact draw tests
# pi K^n = pi with no burn-in or mixing assumption.
#
# `rooted = TRUE` stores the bifurcating rooted trees a production run starts
# from: a uniform labelled rooted topology with flat Dirichlet fractions on all
# 2n - 2 stored edges.  Rooting a uniform unrooted topology on a uniformly
# chosen edge draws it.

.PriorDraw <- function(tipLabel, rooted = FALSE) {
  nTip <- length(tipLabel)
  tree <- TreeTools::RandomTree(tipLabel, root = FALSE)
  tree <- if (rooted) {
    edge <- tree$edge
    clade <- edge[sample.int(nrow(edge), 1L), 2]
    repeat {
      extra <- setdiff(edge[edge[, 1] %in% clade, 2], clade)
      if (!length(extra)) break
      clade <- c(clade, extra)
    }
    TreeTools::Preorder(ape::root(tree, outgroup = clade[clade <= nTip],
                                  resolve.root = TRUE))
  } else {
    TreeTools::Preorder(ape::root(
      tree, node = nTip + sample.int(nTip - 2L, 1L), resolve.root = FALSE))
  }
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
.PriorInvarianceP <- function(moveType, nMove, nReps, nTip = 6L,
                              rooted = FALSE) {
  tipLabel <- paste0("t", seq_len(nTip))
  mat <- matrix(sample(0:1, nTip * 4L, replace = TRUE), nrow = nTip,
                dimnames = list(tipLabel, NULL))
  mkd <- suppressWarnings(MkPrimeData(TreeTools::MatrixToPhyDat(mat)))
  model <- MkPrimeModel()

  sampled <- t(vapply(seq_len(nReps), function(r) {
    tree <- .PriorDraw(tipLabel, rooted)
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
    tree <- .PriorDraw(tipLabel, rooted)
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
