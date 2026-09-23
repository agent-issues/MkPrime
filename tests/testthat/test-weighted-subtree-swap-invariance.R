# WSWAP-001: weighted_subtree_swap (moveType 14) leaves pi invariant.
#
# The move anchors on a node, enumerates (partner, bin) configurations, then
# draws the committed fraction from a Beta centred on the chosen bin.  Two
# terms fail to cancel: the bin mixture, which the reverse move scores on the
# reverse topology, and the selection normaliser of the anchored
# neighbourhood, which differs between the endpoints as it does for
# gibbs_subtree_swap.  Derivation:
# dev/red-team/proofs/weighted-subtree-swap-hastings.md.

test_that("weighted_subtree_swap holds the prior at beta = 0", {
  skip_under_memcheck()
  set.seed(1404L)
  p <- .PriorInvarianceP(moveType = 14L, nMove = 25L, nReps = 1000L)
  expect_true(all(p > 1e-4), info = paste(signif(p, 2), collapse = ", "))
})



# Replicates of a few moves each are weak on topology.  A long chain is not:
# under the prior, 1/7 of six-tip trees have three cherries, and t1 hangs from
# the storage root 1/4 of the time.  Enumerating the reverse neighbourhood on
# the canonicalised proposal, whose internal nodes are renumbered, visits
# three-cherry trees 18% of the time (z = 5.9 here) yet passes the gate above.
test_that("weighted_subtree_swap holds the topology prior at beta = 0", {
  skip_under_memcheck()
  set.seed(1405L)
  nTip <- 6L
  nIter <- 40000L
  tipLabel <- paste0("t", seq_len(nTip))
  mat <- matrix(sample(0:1, nTip * 4L, replace = TRUE), nrow = nTip,
                dimnames = list(tipLabel, NULL))
  mkd <- suppressWarnings(MkPrimeData(TreeTools::MatrixToPhyDat(mat)))
  tree <- .PriorDraw(tipLabel)
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(), tree, mkd)
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(MkPrime:::.InitState(tree, mkd, model))
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)

  shape <- t(vapply(seq_len(nIter), function(i) {
    do_move_cpp(dataPtr, statePtr, 14L, 0L, 0.5, 0.5, 1L, 0.0)
    st <- get_mcmc_state(statePtr)
    .PriorStats(st$edge, st$relBrLengths)[5:6]
  }, numeric(2)))
  observed <- c(mean(shape[, 1]), mean(shape[, 2] == 3))
  batches <- rep(1:50, each = nIter / 50L)
  se <- c(sd(tapply(shape[, 1], batches, mean)),
          sd(tapply(shape[, 2] == 3, batches, mean))) / sqrt(50)
  z <- (observed - c(1 / 4, 1 / 7)) / se
  expect_true(all(abs(z) < 4), info = paste(signif(z, 3), collapse = ", "))
})

# At beta = 0 every weight is 1, so the normaliser ratio is only a ratio of
# neighbourhood sizes, and a kernel that drops it is barely distinguishable
# there.  At beta = 1 it is not.  Five tips give fifteen unrooted topologies,
# whose posterior mass integrates over the Dirichlet edge fractions by plain
# Monte Carlo; the likelihood is invariant to the storage root, so that mass is
# the chain's target however the root wanders.
test_that("weighted_subtree_swap samples its topology posterior at beta = 1", {
  skip_under_memcheck()
  set.seed(2L)
  nTip <- 5L
  nEdge <- 2L * nTip - 3L
  nIter <- 20000L
  tipLabel <- paste0("t", seq_len(nTip))
  mat <- matrix(sample(0:1, nTip * 12L, replace = TRUE), nrow = nTip,
                dimnames = list(tipLabel, NULL))
  mkd <- suppressWarnings(MkPrimeData(TreeTools::MatrixToPhyDat(mat)))
  tree <- TreeTools::Preorder(TreeTools::RandomTree(tipLabel, root = FALSE))
  tree$edge.length <- rep(1 / nEdge, nEdge)
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(), tree, mkd)
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(MkPrime:::.InitState(tree, mkd, model))
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  st <- get_mcmc_state(statePtr)

  LogLik <- function(edge, f) {
    cpp_log_likelihood_xptr(dataPtr, edge[, 1], edge[, 2], st$treeLength * f,
                            st$kPrime, st$rateLoss, st$rateLogSd, st$rateNeo)
  }
  expect_equal(LogLik(st$edge, st$relBrLengths), st$logLik)

  TopologyKey <- function(edge) {
    tr <- structure(list(edge = edge, tip.label = tipLabel, Nnode = nTip - 2L),
                    class = "phylo")
    splits <- as.logical(TreeTools::as.Splits(tr))
    paste(sort(apply(splits, 1, function(s) {
      paste(which(if (s[1]) !s else s), collapse = "")
    })), collapse = "|")
  }
  topologies <- lapply(ape::as.phylo(0:14, nTip, tipLabel),
                       function(t) TreeTools::Preorder(TreeTools::UnrootTree(t)))
  keys <- vapply(topologies, function(t) TopologyKey(t$edge), "")
  logMass <- vapply(topologies, function(t) {
    ll <- vapply(seq_len(2000L), function(i) {
      f <- rexp(nEdge)
      LogLik(t$edge, f / sum(f))
    }, numeric(1))
    max(ll) + log(mean(exp(ll - max(ll))))
  }, numeric(1))
  target <- exp(logMass - max(logMass))
  target <- target / sum(target)

  visited <- vapply(seq_len(nIter), function(i) {
    do_move_cpp(dataPtr, statePtr, 14L, 0L, 0.5, 0.5, 1L, 1.0)
    TopologyKey(get_mcmc_state(statePtr)$edge)
  }, "")
  visited <- factor(visited, keys)
  empirical <- as.numeric(table(visited)) / nIter
  batchFreq <- vapply(split(visited, rep(1:50, each = nIter / 50L)),
                      function(b) as.numeric(table(b)) / length(b),
                      numeric(length(keys)))
  z <- (empirical - target) / pmax(apply(batchFreq, 1, sd) / sqrt(50), 1e-3)

  # Corrected: TV 0.023, max |z| 2.3.  Uncorrected: TV 0.21, max |z| 25.  The
  # bin mixture alone, without the normaliser ratio: TV 0.13, max |z| 18.
  expect_lt(max(abs(z)), 4.5)
  expect_lt(sum(abs(empirical - target)) / 2, 0.05)
})
