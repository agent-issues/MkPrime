# gibbs_spr (moveType 10) leaves pi invariant at beta = 1.
#
# The move draws a regraft edge with probability proportional to
# exp(beta * logLik) at a reference split, so its MH ratio must carry the
# selection-probability ratio beta * (candLL[merged] - candLL[chosen]) (#19,
# GSPR-004).  Without that term the kernel targets pi * Z and over-visits
# high-likelihood topologies; the beta = 0 tests cannot see this because the
# term vanishes there.
#
# Tree length and every model parameter stay fixed, so the target over the
# stored rooted topology r is proportional to the likelihood integrated over
# the flat Dirichlet prior on edge fractions, Z_r = E[L(r, b)], which Monte
# Carlo over exact Dirichlet draws estimates directly.

# Rooted topology key: the tip sets below each internal node.
.RootedKey <- function(edge, nTip) {
  Below <- function(node) {
    below <- node
    repeat {
      extra <- setdiff(edge[edge[, 1] %in% below, 2], below)
      if (!length(extra)) break
      below <- c(below, extra)
    }
    paste(sort(below[below <= nTip]), collapse = "")
  }
  # Return:
  paste(sort(vapply(unique(edge[, 1]), Below, character(1))), collapse = "/")
}

test_that("gibbs_spr samples its exact target at beta = 1", {
  skip_under_memcheck()
  # t1+t2 and t1+t4 each have two supporting characters, t1+t3 one.
  mat <- matrix(c(1, 1, 0, 0,  1, 1, 0, 0,  1, 0, 1, 0,  1, 0, 0, 1,
                  0, 1, 1, 0), nrow = 4L,
                dimnames = list(paste0("t", 1:4), NULL))
  mkd <- suppressWarnings(MkPrimeData(MatrixToPhyDat(mat)))
  tips <- rownames(mkd$matrix)
  nTip <- length(tips)
  tree <- Preorder(RenumberTips(
    ape::read.tree(text = "((t1:0.1,t2:0.1):0.2,(t3:0.1,t4:0.1):0.2);"),
    tips))
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(), tree, mkd)
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(MkPrime:::.InitState(tree, mkd, model))
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  st0 <- get_mcmc_state(statePtr)

  set.seed(3231L)
  nIter <- 20000L
  keys <- character(nIter)
  edges <- list()
  for (i in seq_len(nIter)) {
    do_move_cpp(dataPtr, statePtr, 10L, 0L, 0.5, 0.5, 1L, 1.0)
    edge <- get_mcmc_state(statePtr)$edge
    keys[[i]] <- .RootedKey(edge, nTip)
    if (is.null(edges[[keys[[i]]]])) edges[[keys[[i]]]] <- edge + 0L
  }
  # All 15 rooted topologies of four tips.
  expect_equal(length(edges), 15L)

  nDraw <- 4000L
  fractions <- matrix(rexp(nDraw * nrow(st0$edge)), nDraw)
  fractions <- fractions / rowSums(fractions)
  logZ <- vapply(edges, function(edge) {
    logLik <- apply(fractions, 1, function(f) cpp_log_likelihood_xptr(
      dataPtr, edge[, 1], edge[, 2], st0$treeLength * f, st0$kPrime,
      st0$rateLoss, st0$rateLogSd, st0$rateNeo, st0$betaScale))
    max(logLik) + log(mean(exp(logLik - max(logLik))))
  }, double(1))
  target <- exp(logZ - max(logZ))
  target <- target / sum(target)
  empirical <- as.numeric(table(factor(keys, names(edges)))) / nIter

  # A kernel targeting pi * Z over-visits likely states: log(empirical /
  # target) rises with log(target), t = 7 to 12 without the selection
  # ratio against |t| <= 2 with it.
  fit <- summary(lm(log(empirical / target) ~ log(target)))
  expect_lt(abs(fit$coefficients[2L, 3L]), 4)

  # The minority unrooted topology, t1+t3: target 0.14; 0.057 without the
  # selection ratio.
  minorityKey <- .TopologyKey(Preorder(RenumberTips(
    ape::read.tree(text = "((t1,t3),(t2,t4));"), tips))$edge, tips)
  minority <- vapply(edges, function(edge) {
    .TopologyKey(edge, tips) == minorityKey
  }, logical(1))
  expect_lt(abs(sum(empirical[minority]) - sum(target[minority])), 0.03)
})
