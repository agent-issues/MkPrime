# Detailed balance for the topology proposals as the MCMC engine runs them:
# through do_move_cpp, on the bifurcating rooted trees TreeSearch::AdditionTree
# supplies by default, one move type at a time.
#
# beta = 0 leaves the prior as the target.  It puts no mass on topology and is
# Dirichlet(1, ..., 1) on the branch-length simplex, so a correct proposal
# samples labelled topologies uniformly.  The pSPR candidate weights come from
# the parsimony scores of the data, not from beta, so they stay in play.

# --- Shared fixtures -------------------------------------------------------

# Rooted bifurcating chain state: nEdge = 2 * nTip - 2, root = nTip + 1.
.DbState <- function(seed, nTip = 5L, nChar = 8L, qHeterogeneity = FALSE) {
  set.seed(seed)
  tree <- Preorder(RandomTree(nTip, root = TRUE))
  tree$edge.length <- runif(nrow(tree$edge), 0.05, 0.5)
  mat <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                nrow = nTip, ncol = nChar,
                dimnames = list(tree$tip.label, NULL))
  mkd <- allow_warning(MkPrimeData(MatrixToPhyDat(mat)),
                       "Dropping [0-9]+ invariant characters")
  model <- .FinalizeModel(MkPrimeModel(qHeterogeneity = qHeterogeneity),
                          tree, mkd)
  statePtr <- .InitMcmcChain(.InitState(tree, mkd, model))
  dataPtr <- .InitMcmcData(mkd, model)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  # Return:
  list(dataPtr = dataPtr, statePtr = statePtr, tips = tree$tip.label,
       nTip = nTip, tipStates = mkd$partitions[[1]]$tip_states)
}

# Flat-target chain driven through the production move dispatcher.
.DbFlatChain <- function(moveType, fixtureSeed, chainSeed,
                         nSample = 4000L, thin = 250L,
                         qHeterogeneity = FALSE) {
  pts <- .DbState(fixtureSeed, qHeterogeneity = qHeterogeneity)
  keys <- character(nSample)
  set.seed(chainSeed)
  for (i in seq_len(nSample * thin)) {
    do_move_cpp(pts$dataPtr, pts$statePtr, moveType, 0L, 0.5, 0.5, 1L, 0.0)
    if (i %% thin == 0L) {
      keys[i %/% thin] <-
        .TopologyKey(get_mcmc_state(pts$statePtr)$edge, pts$tips)
    }
  }
  # Return:
  table(keys)
}

.DbExpectUniform <- function(freq) {
  expect_equal(length(freq), 15L)
  expect_gt(chisq.test(as.numeric(freq))$p.value, 0.005)
  expect_lt(max(freq) / min(freq), 2)
}


# --- pSPR proposal density is reversible -----------------------------------

# The pSPR Hastings ratio assumes a candidate weight is a function of the
# topology the regraft produces, so that the reverse move's normalizer is
# sumW + wOrig - w[chosen] and its chosen weight is wOrig.  Checking that
# identity against a reverse enumeration is sharper than any chain: it holds
# exactly or not at all.
test_that("pSPR proposal density is reversible", {
  alpha <- 0.1  # PSPR_ALPHA
  nTip <- 5L
  root <- nTip + 1L
  fixture <- .DbState(504L)
  tipStates <- fixture$tipStates[
    match(fixture$tips, rownames(fixture$tipStates)), , drop = FALSE]
  # Mirrors kEff in pspr_proposal_impl.  An NA here would make every state set
  # full and every weight equal, which the scoreSpread guard below also catches.
  kStates <- max(2L, max(tipStates, na.rm = TRUE) + 1L)

  Canonical <- function(edge) RenumberTree(edge[, 1], edge[, 2])
  TipsBelow <- function(edge, node) {
    below <- node
    repeat {
      extra <- setdiff(edge[edge[, 1] %in% below, 2], below)
      if (!length(extra)) break
      below <- c(below, extra)
    }
    sort(below[below <= nTip])
  }
  # Candidate weights and the original position's weight, sharing the minScore
  # offset that cancels out of the ratio.
  Enumerate <- function(edge, pruneRow) {
    sites <- .PsprSites(edge, nTip, pruneRow)
    if (is.null(sites)) {
      return(NULL)
    }
    scores <- fitch_score_candidates_r(
      edge[, 1], edge[, 2], tipStates, nTip, kStates,
      pruneRow, sites$parentRow, sites$sibRow,
      sites$u, sites$v, sites$sibNode, sites$candidates)
    canonical <- Canonical(edge)
    scoreOrig <- fitch_score_r(canonical[, 1], canonical[, 2],
                               tipStates, nTip, kStates)
    minScore <- min(c(scores, scoreOrig))
    # Return:
    list(sites = sites,
         w = exp(-alpha * (scores - minScore)),
         wOrig = exp(-alpha * (scoreOrig - minScore)))
  }

  nChecked <- 0L
  scoreSpread <- integer(0)
  set.seed(1L)
  for (trial in seq_len(20L)) {
    edge <- Preorder(RandomTree(nTip, root = TRUE))$edge
    for (pruneRow in which(edge[, 1] != root)) {
      fwd <- Enumerate(edge, pruneRow)
      if (is.null(fwd)) next
      sumW <- sum(fwd$w)
      scoreSpread <- c(scoreSpread, fwd$w)
      prunedTips <- TipsBelow(edge, fwd$sites$v)
      backKey <- paste(Canonical(edge), collapse = ",")
      for (ci in seq_along(fwd$sites$candidates)) {
        moved <- Canonical(.PsprRegraft(edge, fwd$sites,
                                        fwd$sites$candidates[ci]))
        # The same pruned subtree, renumbered by the canonical reordering.
        movedV <- Find(function(nd) identical(TipsBelow(moved, nd), prunedTips),
                       unique(c(moved[, 2], moved[, 1])))
        bwd <- Enumerate(moved, which(moved[, 2] == movedV))
        expect_false(is.null(bwd))
        back <- which(vapply(bwd$sites$candidates, function(rr) {
          identical(paste(Canonical(.PsprRegraft(moved, bwd$sites, rr)),
                          collapse = ","), backKey)
        }, logical(1)))
        expect_length(back, 1L)

        trueRatio <- log(bwd$w[back]) - log(sum(bwd$w)) -
          log(fwd$w[ci]) + log(sumW)
        sumWRev <- sumW + fwd$wOrig - fwd$w[ci]
        codeRatio <- log(fwd$wOrig) - log(fwd$w[ci]) +
          log(sumW) - log(sumWRev)
        expect_equal(codeRatio, trueRatio, tolerance = 1e-9)
        nChecked <- nChecked + 1L
      }
    }
  }
  expect_gt(nChecked, 400L)
  expect_gt(length(unique(scoreSpread)), 1L)
})


# --- The dispatcher reaches the branches under test ------------------------

test_that("case 5 takes both the in-place and the reordering NNI branch", {
  pts <- .DbState(211L)
  nMoved <- 0L
  for (i in seq_len(400L)) {
    nMoved <- nMoved +
      do_move_cpp(pts$dataPtr, pts$statePtr, 5L, 0L, 0.5, 0.5, 1L, 0.0)
  }
  inPlace <- get_mcmc_state(pts$statePtr)$diagNniPartial
  expect_gt(inPlace, 0L)        # wRow > edgeRow: swap applied in place
  expect_lt(inPlace, nMoved)    # wRow <= edgeRow: canonicalized via reorder
})

test_that("case 6 takes the TreeNav SPR branch", {
  pts <- .DbState(212L)
  do_move_cpp(pts$dataPtr, pts$statePtr, 6L, 0L, 0.5, 0.5, 1L, 0.0)
  before <- get_mcmc_state(pts$statePtr)$edge
  everCold <- !node_cl_ready(pts$statePtr)
  moved <- FALSE
  for (i in seq_len(50L)) {
    do_move_cpp(pts$dataPtr, pts$statePtr, 6L, 0L, 0.5, 0.5, 1L, 0.0)
    everCold <- everCold || !node_cl_ready(pts$statePtr)
    moved <- moved || !identical(before, get_mcmc_state(pts$statePtr)$edge)
  }
  expect_false(everCold)
  expect_true(moved)
})


# --- Flat-target uniformity, one move type at a time -----------------------

# These confirm the whole committed move — proposal, Hastings ratio, canonical
# reordering, cache invalidation — leaves a flat target flat.  They cannot
# resolve a bias of a few per cent, so they complement rather than replace the
# reversibility test above.

test_that("NNI satisfies detailed balance on rooted bifurcating trees", {
  skip_slow_tests()
  skip_under_memcheck()
  .DbExpectUniform(.DbFlatChain(5L, fixtureSeed = 4729L, chainSeed = 13L))
})

test_that("TreeNav SPR satisfies detailed balance on rooted bifurcating trees", {
  skip_slow_tests()
  skip_under_memcheck()
  .DbExpectUniform(.DbFlatChain(6L, fixtureSeed = 4729L, chainSeed = 13L))
})

# qHeterogeneity routes case 6 to spr_proposal_impl whatever the cache holds.
test_that("legacy SPR satisfies detailed balance on rooted bifurcating trees", {
  skip_slow_tests()
  skip_under_memcheck()
  .DbExpectUniform(.DbFlatChain(6L, fixtureSeed = 4729L, chainSeed = 13L,
                                qHeterogeneity = TRUE))
})

test_that("TBR satisfies detailed balance on rooted bifurcating trees", {
  skip_slow_tests()
  skip_under_memcheck()
  .DbExpectUniform(.DbFlatChain(17L, fixtureSeed = 4729L, chainSeed = 13L))
})

test_that("pSPR satisfies detailed balance on rooted bifurcating trees", {
  skip_slow_tests()
  skip_under_memcheck()
  .DbExpectUniform(.DbFlatChain(20L, fixtureSeed = 4729L, chainSeed = 13L))
})
