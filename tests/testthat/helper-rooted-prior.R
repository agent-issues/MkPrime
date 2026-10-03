# beta = 0 invariance harness on rooted trees (#338).
#
# The default start trees have a degree-2 root (2n - 2 edges), and no move
# changes the root's degree, so such a chain walks labelled rooted topologies.
# At beta = 0 its target is uniform over those, with flat Dirichlet edge
# fractions.  A uniform rooted tree on n tips is a uniform unrooted tree on
# n + 1 tips rooted on the extra tip and shorn of it.

.RootedPriorDraw <- function(tipLabel) {
  tree <- TreeTools::RandomTree(c(tipLabel, ".root"), root = ".root")
  tree <- TreeTools::Preorder(TreeTools::DropTip(tree, ".root"))
  f <- rexp(nrow(tree$edge))
  tree$edge.length <- f / sum(f)
  tree
}

# Smallest, largest and summed-internal edge fraction; pendant fraction of t1;
# tips on the smaller side of the root; number of cherries.
.RootedPriorStats <- function(edge, f) {
  nTip <- min(edge[, 1]) - 1L
  rootKids <- edge[edge[, 1] == nTip + 1L, 2]
  below <- vapply(rootKids, function(node) {
    nodes <- node
    frontier <- node
    while (length(frontier)) {
      frontier <- edge[edge[, 1] %in% frontier, 2]
      nodes <- c(nodes, frontier)
    }
    sum(nodes <= nTip)
  }, numeric(1))
  c(min(f), max(f), sum(f[edge[, 2] > nTip]), f[edge[, 2] == 1L],
    min(below), sum(tabulate(edge[edge[, 2] <= nTip, 1], 2L * nTip) == 2L))
}

# p-values comparing .RootedPriorStats after `nMove` applications of
# `moveType` at beta = 0 with ten times as many direct draws: KS for the four
# continuous statistics, chi-squared for the two discrete ones.
.RootedInvarianceP <- function(moveType, nMove, nReps, nTip = 6L) {
  tipLabel <- paste0("t", seq_len(nTip))
  mat <- matrix(sample(0:1, nTip * 4L, replace = TRUE), nrow = nTip,
                dimnames = list(tipLabel, NULL))
  mkd <- suppressWarnings(MkPrimeData(TreeTools::MatrixToPhyDat(mat)))
  model <- MkPrimeModel()

  sampled <- t(vapply(seq_len(nReps), function(r) {
    tree <- .RootedPriorDraw(tipLabel)
    finalModel <- MkPrime:::.FinalizeModel(model, tree, mkd)
    dataPtr  <- MkPrime:::.InitMcmcData(mkd, finalModel)
    statePtr <- MkPrime:::.InitMcmcChain(
      MkPrime:::.InitState(tree, mkd, finalModel))
    fill_partition_cache(dataPtr, statePtr)
    allocate_cl_workspace(dataPtr, statePtr)
    for (i in seq_len(nMove))
      do_move_cpp(dataPtr, statePtr, moveType, 0L, 0.5, 0.5, 1L, 0.0)
    st <- get_mcmc_state(statePtr)
    stopifnot(nrow(st$edge) == 2L * nTip - 2L)
    .RootedPriorStats(st$edge, st$relBrLengths)
  }, numeric(6L)))

  reference <- t(vapply(seq_len(nReps * 10L), function(i) {
    tree <- .RootedPriorDraw(tipLabel)
    .RootedPriorStats(tree$edge, tree$edge.length)
  }, numeric(6L)))

  c(
    vapply(1:4, function(j) suppressWarnings(
      ks.test(sampled[, j], reference[, j])$p.value), numeric(1)),
    vapply(5:6, function(j) {
      lev <- sort(unique(c(sampled[, j], reference[, j])))
      if (length(lev) < 2L) return(1)
      suppressWarnings(chisq.test(rbind(
        table(factor(sampled[, j], lev)),
        table(factor(reference[, j], lev))))$p.value)
    }, numeric(1))
  )
}
