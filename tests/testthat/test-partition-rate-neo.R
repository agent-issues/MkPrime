# rate_neo under the partition API (issue #241).
#
# rate_neo splits the rate between neomorphic and other characters exactly as
# on the legacy path; class rates multiply that split, renormalised so the
# char-weighted mean rate stays 1. Before #241 every partitioned evaluator
# scored rate_neo = 1 whatever the state held.

# 10 tips; characters 1-8 neomorphic, 9-18 transformational.
.NeoTransData <- function(seed = 21L) {
  set.seed(seed)
  nTip <- 10L
  nChar <- 18L
  mat <- matrix(sample(0:2, nTip * nChar, TRUE), nTip,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  mat[, 1:8] <- sample(0:1, nTip * 8L, TRUE)
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1L, j] <- 1L - (mat[1L, j] > 0)
  }
  pd <- MatrixToPhyDat(mat)
  tree <- Preorder(NJTree(pd, edgeLengths = TRUE))
  tree$edge.length <- pmax(tree$edge.length, 0.05)
  list(mkd = MkPrimeData(pd, neomorphic = 1:8), tree = tree)
}

# Full log-likelihood of the initial state, with or without a partition.
.RateNeoLogLik <- function(d, rateNeo, partition = NULL) {
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(), d$tree, d$mkd)
  mkd <- d$mkd
  if (is.null(partition)) {
    s <- MkPrime:::.InitState(d$tree, mkd, model)
  } else {
    mkd$partitions <- MkPrime:::.BuildPartitions(mkd, partition)
    K <- max(partition)
    s <- MkPrime:::.InitStatePartitioned(d$tree, mkd, model, list(
      partition = partition, nClasses = K,
      unlink = if (K > 1L) "ratemultiplier" else character(0)))
  }
  ch <- list(edge = s$tree$edge, rel_br_lengths = s$rel_br_lengths,
             tree_length = s$tree_length, rate_loss = s$rate_loss,
             rate_log_sd = s$rate_log_sd, rate_neo = rateNeo, p = s$p,
             kPrime = as.integer(s$kPrime),
             log_lik = 0, log_prior = 0, log_post = 0,
             class_rate_log_sd = s$class_rate_log_sd, class_w = s$class_w,
             class_rate = s$class_rate, nChar_c = s$nChar_c,
             eta_neo = s$eta_neo)
  eval_full_loglik_cpp(MkPrime:::.InitMcmcData(mkd, model),
                       MkPrime:::.DeserialiseChain(ch))
}

test_that("a partition scores rate_neo as the legacy path does", {
  d <- .NeoTransData()
  expect_equal(table(d$mkd$type)[["neomorphic"]], 8L)
  one <- rep(1L, d$mkd$nChar)
  two <- rep(1:2, length.out = d$mkd$nChar)
  for (rateNeo in c(0.1, 10)) {
    legacy <- .RateNeoLogLik(d, rateNeo)
    expect_equal(.RateNeoLogLik(d, rateNeo, one), legacy, tolerance = 1e-10)
    expect_equal(.RateNeoLogLik(d, rateNeo, two), legacy, tolerance = 1e-10)
  }
  expect_gt(abs(.RateNeoLogLik(d, 0.1, one) - .RateNeoLogLik(d, 10, one)), 1)
})

# With one class of neomorphic and one of transformational characters, class
# rates (c1, c2) and rate_neo r set the neo:trans rate ratio to r * c1 / c2,
# which the legacy path reaches with rate_neo = r * c1 / c2. Both paths hold
# the char-weighted mean rate at 1, so the likelihoods must agree.
test_that("class rates compose with rate_neo at a char-weighted mean of 1", {
  d <- .NeoTransData()
  partition <- ifelse(d$mkd$type == "neomorphic", 1L, 2L)
  classW <- c(0.3, 0.7)
  classRate <- classW * d$mkd$nChar / tabulate(partition, 2L)
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(), d$tree, d$mkd)
  legacyData <- MkPrime:::.InitMcmcData(d$mkd, model)
  for (rateNeo in c(0.2, 1, 5)) {
    chain <- .PartitionedChain(d$tree, d$mkd, classW = classW,
                               classRateLogSd = c(0.5, 0.5),
                               partition = partition, rateNeo = rateNeo)
    st <- get_mcmc_state(chain$statePtr)
    legacy <- cpp_log_likelihood_xptr(
      legacyData, st$edge[, 1], st$edge[, 2],
      st$treeLength * st$relBrLengths, st$kPrime, st$rateLoss, 0.5,
      rateNeo * classRate[[1]] / classRate[[2]], st$betaScale)
    expect_equal(eval_full_loglik_cpp(chain$dataPtr, chain$statePtr), legacy,
                 tolerance = 1e-10, label = paste("rate_neo", rateNeo))
  }
})

# Every evaluator (node-CL cache, k' sweep slots, CL groups) must score the
# rate_neo the state holds.
test_that("every move commits the partitioned likelihood at rate_neo != 1", {
  d <- .NeoTransData()
  moves <- c(scale_tree_length = 0L, rate_neo = 3L, beta_simplex = 4L,
             nni = 5L, spr = 6L, int_walk = 7L, gibbs_spr = 10L,
             gibbs_subtree_swap = 11L, tbr = 17L, slice_rate_neo = 19L,
             pspr = 20L, dirichlet_branch = 23L, local_dirichlet = 24L,
             gibbs_kprime_sweep = 25L, block_kprime_shift = 26L,
             dirichlet_simplex_class_w = 32L, joint_tl_rn = 33L)
  partitions <- list(one = rep(1L, d$mkd$nChar),
                     two = rep(1:2, length.out = d$mkd$nChar))
  for (pName in names(partitions)) {
    part <- partitions[[pName]]
    oneClass <- pName == "one"
    for (move in names(moves)) {
      if (oneClass && move == "dirichlet_simplex_class_w") next
      chain <- .PartitionedChain(
        d$tree, d$mkd, partition = part, rateNeo = 3,
        classW = if (oneClass) 1 else c(0.3, 0.7),
        classRateLogSd = if (oneClass) 0.5 else c(0.5, 0.5))
      set.seed(1L)
      accepted <- 0L
      drift <- 0
      for (i in seq_len(40L)) {
        charIdx <- switch(move,
                          int_walk = sample(which(d$mkd$type != "neomorphic"),
                                            1L) - 1L,
                          slice_rate_neo = 3L,
                          0L)
        if (do_move_cpp(chain$dataPtr, chain$statePtr, moves[[move]],
                        charIdx, 0.5, 0.5, 3L, 1.0)) {
          accepted <- accepted + 1L
          drift <- max(drift, abs(get_state_log_lik(chain$statePtr) -
                                    .PartitionedLogLik(chain)))
        }
      }
      label <- paste(pName, move)
      expect_gt(accepted, 0L, label = paste(label, "acceptances"))
      expect_lt(drift, 1e-9, label = paste(label, "logLik drift"))
      if (oneClass && move %in% c("nni", "beta_simplex", "dirichlet_branch")) {
        st <- get_mcmc_state(chain$statePtr)
        expect_gt(st$diagNniPartial + st$diagBsPartial + st$diagDirPartial, 0,
                  label = paste(label, "node-CL cache evaluations"))
      }
    }
  }
})

test_that("rate_neo is sampled only when neomorphic and other characters mix", {
  d <- .NeoTransData()
  mkd <- d$mkd
  allNeo <- mkd
  allNeo$type[] <- "neomorphic"
  noNeo <- mkd
  noNeo$type[] <- "transformational"
  expect_true(MkPrime:::.RateNeoLive(mkd))
  expect_false(MkPrime:::.RateNeoLive(allNeo))
  expect_false(MkPrime:::.RateNeoLive(noNeo))

  rateNeoMoves <- c("rate_neo", "slice_rate_neo", "joint_tl_rn")
  mcmc <- MkPrimeMCMC(joint2d = TRUE)
  spec <- list(partition = rep(1:2, length.out = mkd$nChar), nClasses = 2L,
               unlink = "ratemultiplier")
  for (live in c(TRUE, FALSE)) {
    legacy <- MkPrime:::.BuildMoves(20L, 10L, TRUE, mcmc, rateNeoLive = live)
    part <- MkPrime:::.BuildMovesPartitioned(20L, 10L, TRUE, mcmc, spec,
                                             rateNeoLive = live)
    for (moves in list(legacy, part)) {
      nms <- vapply(moves, `[[`, character(1), "name")
      expect_equal(rateNeoMoves %in% nms, rep(live, 3L))
      expect_true(all(c("rate_loss", "slice_rate_loss", "joint_tl_rl") %in%
                        nms))
    }
  }
  expect_true("rate_neo" %in% MkPrime:::.FixedCols("rate_neo", list()))
})

test_that("RunMkPrime holds rate_neo at 1 on all-neomorphic data", {
  d <- .NeoTransData()
  set.seed(3L)
  mat <- matrix(sample(0:1, 10L * 12L, TRUE), 10L,
                dimnames = list(paste0("t", 1:10), NULL))
  mat[1, ] <- 0L
  mat[2, ] <- 1L
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = seq_len(ncol(mat)))
  expect_true(all(mkd$type == "neomorphic"))
  set.seed(4L)
  res <- allow_warning(RunMkPrime(
    mkd, d$tree, nIter = 100L, maxWarmup = 50L, minWarmup = 50L,
    nChains = 1L, nRuns = 1L, thin = 1L, maxTime = 30,
    partition = rep(1:2, length.out = mkd$nChar), unlink = "ratemultiplier"
  ), "")
  expect_false(any(c("rate_neo", "slice_rate_neo", "joint_tl_rn") %in%
                     names(res$moveWeights)))
  samples <- as.matrix(res$samples %||% res$runSamples[[1]])
  expect_true(all(samples[, "rate_neo"] == 1))
})
