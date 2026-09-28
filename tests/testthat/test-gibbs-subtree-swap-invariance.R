# GSWAP-001: gibbs_subtree_swap (moveType 11) leaves pi invariant.
#
# The move anchors on a uniformly chosen node, enumerates its valid swap
# partners, and draws one proportional to exp(beta * logLik).  The anchored
# neighbourhood is not the same in both directions, so the selection
# normaliser does not cancel and the min(1, Z_x / Z_y) accept step is load
# bearing: dev/red-team/proofs/gibbs-subtree-swap-hastings.md.
#
# Giving every edge the same length makes the reachable state space finite --
# the move only permutes slot contents, so one shared length value is a fixed
# point -- and the target is then exp(logLik) normalised over the states the
# chain visits, every state carrying the same exchangeable Dirichlet prior.


.SwapChain <- function(nTip = 5L, nIter = 60000L, seed = 7L) {
  set.seed(seed)
  tree <- TreeTools::Preorder(ape::rtree(nTip, rooted = FALSE))
  tree$edge.length <- rep(1 / nrow(tree$edge), nrow(tree$edge))
  mat <- matrix(sample(0:1, nTip * 8L, replace = TRUE), nrow = nTip,
                dimnames = list(tree$tip.label, NULL))
  mkd   <- suppressWarnings(MkPrimeData(MatrixToPhyDat(mat)))
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(), tree, mkd)
  statePtr <- MkPrime:::.InitMcmcChain(MkPrime:::.InitState(tree, mkd, model))
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)

  keys <- character(nIter)
  lls  <- numeric(nIter)
  for (i in seq_len(nIter)) {
    do_move_cpp(dataPtr, statePtr, 11L, 0L, 0.5, 0.5, 1L, 1.0)
    st <- get_mcmc_state(statePtr)
    keys[i] <- paste(st$edge, collapse = ",")
    lls[i]  <- st$logLik
  }

  tab <- table(keys)
  uk  <- names(tab)
  tgt <- exp(lls[match(uk, keys)] - max(lls))
  list(empirical = as.numeric(tab) / nIter, target = tgt / sum(tgt))
}


test_that("gibbs_subtree_swap samples its exact target", {
  chain <- .SwapChain()

  expect_gt(length(chain$target), 20L)

  # A kernel targeting pi * Z over-visits high-probability states, which shows
  # up as a positive slope of log(empirical / target) on log(target): t > 7
  # without the accept step, however long the run.
  fit <- summary(lm(log(chain$empirical / chain$target) ~ log(chain$target)))
  expect_lt(abs(fit$coefficients[2L, 3L]), 4)

  # Total variation is only a loose sanity bound at this sample size: its
  # corrected and uncorrected values (0.016 and 0.029 here) are a single RNG
  # stream apart.  The slope above is the discriminating statistic.
  expect_lt(sum(abs(chain$empirical - chain$target)) / 2, 0.05)
})


# Each of the three evaluation paths -- partial CL, Q-heterogeneity, and the
# full evaluator that coding = "informative" forces -- must score the tree the
# commit then writes.  The evaluator attaches each subtree with the stem of
# the slot it moves INTO; a commit giving each subtree the stem it arrived
# WITH leaves state$logLik holding the likelihood of a different tree,
# several log-units adrift per accepted swap.
.SwapLogLikDrift <- function(model, seed = 5719L, nTip = 8L, nMove = 400L,
                             partitioned = FALSE) {
  set.seed(seed)
  tree <- TreeTools::Preorder(ape::rtree(nTip, rooted = FALSE))
  mat <- matrix(sample(0:2, nTip * 10L, replace = TRUE), nrow = nTip,
                dimnames = list(tree$tip.label, paste0("c", seq_len(10L))))
  mkd <- suppressWarnings(MkPrimeData(MatrixToPhyDat(mat)))
  if (partitioned) {
    chain <- .PartitionedChain(tree, mkd, model)
    dataPtr <- chain$dataPtr
    statePtr <- chain$statePtr
    ColdLogLik <- function(st) .PartitionedLogLik(chain)
  } else {
    finalModel <- MkPrime:::.FinalizeModel(model, tree, mkd)
    statePtr <- MkPrime:::.InitMcmcChain(
      MkPrime:::.InitState(tree, mkd, finalModel))
    dataPtr <- MkPrime:::.InitMcmcData(mkd, finalModel)
    fill_partition_cache(dataPtr, statePtr)
    allocate_cl_workspace(dataPtr, statePtr)
    ColdLogLik <- function(st) {
      cpp_log_likelihood_xptr(
        dataPtr, st$edge[, 1], st$edge[, 2],
        st$treeLength * st$relBrLengths, st$kPrime,
        st$rateLoss, st$rateLogSd, st$rateNeo, st$betaScale)
    }
  }

  set.seed(11L)
  accepted <- 0L
  worst <- 0
  for (i in seq_len(nMove)) {
    if (do_move_cpp(dataPtr, statePtr, 11L, 0L, 0.5, 0.5, 1L, 1.0)) {
      accepted <- accepted + 1L
      st <- get_mcmc_state(statePtr)
      worst <- max(worst, abs(st$logLik - ColdLogLik(st)))
    }
  }
  list(accepted = accepted, drift = worst)
}


test_that("gibbs_subtree_swap commits the tree it evaluated", {
  for (model in list(MkPrimeModel(),
                     MkPrimeModel(qHeterogeneity = TRUE),
                     MkPrimeModel(coding = "informative"))) {
    result <- .SwapLogLikDrift(model)
    expect_gt(result$accepted, 50L)
    expect_lt(result$drift, 1e-8)
  }
})


# Under the partition API the neighbourhood must be scored with each class's
# own shape and rate multiplier, as state$logLik is: a class-blind score
# commits the likelihood of a different model (#163).
test_that("gibbs_subtree_swap scores the partitioned likelihood", {
  for (model in list(MkPrimeModel(),
                     MkPrimeModel(qHeterogeneity = TRUE),
                     MkPrimeModel(coding = "informative"))) {
    result <- .SwapLogLikDrift(model, partitioned = TRUE)
    expect_gt(result$accepted, 50L)
    expect_lt(result$drift, 1e-8)
  }
})
