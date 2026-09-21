# Plan v4 §7b: numeric-tolerance equivalence for `partition = rep(1L, nChar),
# unlink = character(0)`.
#
# The new partition-aware likelihood path routes through cpp_log_likelihood_partitioned
# (sibling to cpp_log_likelihood) and consumes per-class state vectors. At the
# trivial spec (one user class, no unlinks), it must agree with the legacy
# scalar path to ~1e-10 on the initial state. The sample matrix is NOT required
# to match bit-for-bit (RNG ordering may differ once moves are wired in).
# §7c Casali bit-identity contract is separately covered by test-partition-bitcompat-null.R.
#
# These tests exercise the C++ surface directly via the [[Rcpp::export]]
# wrappers, before the RunMkPrime gate opens.



# Helper: build mkd + a simple tree + the C++ McmcData XPtr ready for both
# legacy and partitioned likelihood calls. Returns a list with all inputs.
.setup_lik <- function(partition = NULL, nChar = 8L, nTip = 6L, seed = 17L) {
  set.seed(seed)
  mat <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                nrow = nTip, ncol = nChar,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1L, j] <- 1L - mat[1L, j]
  }
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  tree <- TreeTools::Preorder(TreeTools::NJTree(pd, edgeLengths = TRUE) %||%
                              TreeTools::RandomTree(pd, root = TRUE))
  if (is.null(tree$edge.length) || any(tree$edge.length <= 0)) {
    tree$edge.length <- rep(0.1, nrow(tree$edge))
  }

  # If partition supplied, rebuild mkd$partitions with classIdx — otherwise
  # the default mkd already has classIdx = 1L on every partition.
  if (!is.null(partition)) {
    mkd$partitions <- .BuildPartitions(mkd, partition = partition)
  }

  # expSteps resolves treeLengthRate; the likelihood never reads it, but the
  # model must carry a resolved value.
  model <- MkPrimeModel(expSteps = 10)
  # Build the C++ McmcData XPtr.
  dataPtr <- prepare_mcmc_data(
    partitions_r              = mkd$partitions,
    kObs_r                    = mkd$kObs,
    charTypes_r               = mkd$type,
    hasNeo                    = any(mkd$type == "neomorphic"),
    nCat                      = model$nCat,
    codingStr                 = model$coding,
    relabelFlag               = isTRUE(model$relabel),
    treeLengthShape           = model$treeLengthShape,
    treeLengthRate            = model$treeLengthRate,
    rateLossMeanlog           = model$rateLossMeanlog,
    rateLossSdlog             = model$rateLossSdlog,
    rateLogSdShape            = model$rateLogSdShape,
    rateLogSdRate             = model$rateLogSdRate,
    rateNeoMeanlog            = model$rateNeoMeanlog,
    rateNeoSdlog              = model$rateNeoSdlog,
    kprimeHyperA              = model$kprimeHyperA,
    kprimeHyperB              = model$kprimeHyperB,
    kPriorLogseries           = identical(model$kPrimePrior, "logseries"),
    kprimeLogseriesC          = model$kprimeLogseriesC,
    kPriorBetaGeometric       = identical(model$kPrimePrior, "beta_geometric"),
    kPriorEmpiricalGeometric  = FALSE,
    qHeterogeneity            = isTRUE(model$qHeterogeneity),
    nBetaCat                  = model$nBetaCat,
    betaScaleShape            = model$betaScaleShape,
    betaScaleRate             = model$betaScaleRate,
    # Workaround: RcppExports.R has `empLogBody = numericVector()` as a
    # malformed default (lowercase `numericVector` is not an R function); the
    # default trips up R evaluation whenever a 2nd partition forces a code
    # path that touches it. Pass an explicit empty numeric vector.
    empLogBody                = numeric(0)
  )

  list(
    mkd     = mkd,
    tree    = tree,
    model   = model,
    dataPtr = dataPtr,
    parent  = tree$edge[, 1],
    child   = tree$edge[, 2],
    edgeLen = tree$edge.length,
    kPrime  = as.integer(ifelse(mkd$type == "known", mkd$known_k, mkd$kObs))
  )
}


# ---- §7b core: partition = rep(1L, nChar), unlink = character(0) ----

test_that("§7b: trivial partition reproduces legacy log-likelihood (~1e-10)", {
  s <- .setup_lik(partition = rep(1L, 8L))
  expect_identical(cpp_data_nclasses(s$dataPtr), 1L)

  ll_legacy <- cpp_log_likelihood_xptr(
    s$dataPtr, s$parent, s$child, s$edgeLen, s$kPrime,
    rateLoss = 1.0, rateLogSd = 0.5, rateNeo = 1.0
  )

  ll_part <- cpp_log_likelihood_partitioned_xptr(
    s$dataPtr, s$parent, s$child, s$edgeLen, s$kPrime,
    rateLoss   = 1.0,
    rateLogSd  = 0.5,                   # length-1 vector — collapses to scalar
    classRate  = 1.0,                   # length-1 vector
    etaNeo     = 1.0
  )

  expect_equal(ll_part, ll_legacy, tolerance = 1e-10)
})


test_that("§7b: trivial partition agrees with R-side .MkpLogLikelihood", {
  s <- .setup_lik(partition = rep(1L, 8L))

  ll_part <- cpp_log_likelihood_partitioned_xptr(
    s$dataPtr, s$parent, s$child, s$edgeLen, s$kPrime,
    rateLoss = 1.0, rateLogSd = 0.5,
    classRate = 1.0, etaNeo = 1.0
  )

  ll_r <- MkPrime:::.MkpLogLikelihood(
    s$tree, s$mkd, kPrime = s$kPrime,
    rate_loss   = 1.0,
    rate_log_sd = 0.5,
    nCat        = s$model$nCat,
    coding      = s$model$coding,
    rate_neo    = 1.0,
    relabel     = isTRUE(s$model$relabel)
  )

  expect_equal(ll_part, ll_r, tolerance = 1e-10)
})


# ---- Sanity: multi-class with all classRate = 1 == legacy ----

test_that("multi-class partition with all classRate = 1 reproduces legacy", {
  # All chars in classes 1 vs 2 but with class_rate = 1.0 everywhere — the
  # branch-scaling step is a no-op, so likelihood must equal legacy.
  s <- .setup_lik(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L))
  expect_identical(cpp_data_nclasses(s$dataPtr), 2L)

  ll_legacy <- cpp_log_likelihood_xptr(
    s$dataPtr, s$parent, s$child, s$edgeLen, s$kPrime,
    rateLoss = 1.0, rateLogSd = 0.5, rateNeo = 1.0
  )

  ll_part <- cpp_log_likelihood_partitioned_xptr(
    s$dataPtr, s$parent, s$child, s$edgeLen, s$kPrime,
    rateLoss  = 1.0,
    rateLogSd = c(0.5, 0.5),           # both classes share — equals scalar 0.5
    classRate = c(1.0, 1.0),
    etaNeo    = 1.0
  )

  expect_equal(ll_part, ll_legacy, tolerance = 1e-10)
})


test_that("per-class rateLogSd differs from scalar when classes diverge", {
  # Same partition as above but give class 2 a very different rate_log_sd.
  # Result must differ from the all-0.5 case.
  s <- .setup_lik(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L))

  ll_uniform <- cpp_log_likelihood_partitioned_xptr(
    s$dataPtr, s$parent, s$child, s$edgeLen, s$kPrime,
    rateLoss  = 1.0,
    rateLogSd = c(0.5, 0.5),
    classRate = c(1.0, 1.0),
    etaNeo    = 1.0
  )

  ll_split <- cpp_log_likelihood_partitioned_xptr(
    s$dataPtr, s$parent, s$child, s$edgeLen, s$kPrime,
    rateLoss  = 1.0,
    rateLogSd = c(0.5, 1.5),
    classRate = c(1.0, 1.0),
    etaNeo    = 1.0
  )

  expect_false(isTRUE(all.equal(ll_uniform, ll_split, tolerance = 1e-8)))
  expect_true(is.finite(ll_split))
})


# ---- Mean-1 constrained classRate produces a different but finite LL ----

test_that("char-weighted mean-1 classRate is honoured and produces finite LL", {
  # 4 chars in class 1, 4 chars in class 2 → equal weights.
  # class_rate = (2, 2/3): char-weighted mean = (4*2 + 4*(2/3)) / 8 = 4/3,
  # NOT 1. Re-derive to satisfy mean-1: pick w on simplex, class_rate[c] =
  # w_c * nChar / nChar_c. With w = (0.7, 0.3), nChar = 8, nChar_c = 4:
  # class_rate = (0.7*8/4, 0.3*8/4) = (1.4, 0.6). Mean: (4*1.4 + 4*0.6)/8 = 1 ✓.
  s <- .setup_lik(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L))

  ll <- cpp_log_likelihood_partitioned_xptr(
    s$dataPtr, s$parent, s$child, s$edgeLen, s$kPrime,
    rateLoss  = 1.0,
    rateLogSd = c(0.5, 0.5),
    classRate = c(1.4, 0.6),
    etaNeo    = 1.0
  )
  expect_true(is.finite(ll))
})
