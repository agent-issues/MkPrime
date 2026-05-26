# Tests for Layer 1 per-class MCMC moves:
#   case 31: scale_class_rate_log_sd (Bactrian MH on classRateLogSd[c])
#   case 32: dirichlet_simplex_class_w (Dirichlet simplex on class_w)
#
# IMPORTANT: cpp_log_prior does not yet include the per-class Gamma (shape)
# or Dirichlet (ratemultiplier) prior contributions — that is the parallel
# agent's work. Acceptance ratios are therefore driven by LL-ratio only.
# These tests assert:
#   (a) chains run to completion with all-finite LLs (moves fire without error),
#   (b) per-class columns are present in result$samples,
#   (c) per-class values evolve across samples (moves are accepted sometimes),
#   (d) simplex constraint: classW columns sum to ~1 every sample.
#
# Precise posterior correctness is NOT asserted here; the integration commit
# (third agent) will add the prior terms and a numeric-accuracy test.

library("TreeTools")

# ---- shared data setup -------------------------------------------------------

.setup_2class_data <- function(seed = 77L, nChar = 10L, nTip = 6L) {
  set.seed(seed)
  mat <- matrix(
    sample(0:1, nTip * nChar, replace = TRUE),
    nrow = nTip, ncol = nChar,
    dimnames = list(paste0("t", seq_len(nTip)), NULL)
  )
  # Ensure every character is variable
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1L, j] <- 1L - mat[1L, j]
  }
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  tree <- TreeTools::Preorder(
    TreeTools::NJTree(pd, edgeLengths = TRUE)
  )
  if (is.null(tree$edge.length) || any(tree$edge.length <= 0)) {
    tree$edge.length <- pmax(
      tree$edge.length %||% rep(0.1, nrow(tree$edge)), 1e-8
    )
  }
  # 2-class partition: first half class 1, second half class 2
  n  <- mkd$nChar
  part <- c(rep(1L, floor(n / 2)), rep(2L, n - floor(n / 2)))
  list(mkd = mkd, tree = tree, part = part)
}

.short_mcmc <- function() {
  MkPrimeMCMC(
    nIter       = 100L,
    maxWarmup   = 50L,
    minWarmup   = 50L,
    nChains     = 1L,
    thin        = 1L,
    autoTune    = FALSE,
    gibbsSubtreeSwap = FALSE   # avoid non-determinism from extra topology moves
  )
}


# ==============================================================================
# 1. unlink = "shape"
# ==============================================================================

test_that("unlink=shape chain runs, per-class shape columns present and vary", {
  d <- .setup_2class_data()
  set.seed(42L)
  result <- RunMkPrime(
    data      = d$mkd,
    tree      = d$tree,
    mcmc      = .short_mcmc(),
    partition = d$part,
    unlink    = "shape"
  )

  samp <- result$samples

  # (a) All LLs finite
  expect_true(all(is.finite(samp[, "log_likelihood"])))
  expect_gt(result$nSamples, 0L)

  # (b) Per-class columns present
  expect_true("class1_rate_log_sd" %in% colnames(samp))
  expect_true("class2_rate_log_sd" %in% colnames(samp))

  # (c) Values are positive (rateLogSd is always positive)
  expect_true(all(samp[, "class1_rate_log_sd"] > 0))
  expect_true(all(samp[, "class2_rate_log_sd"] > 0))

  # (d) Values evolve: not all identical (moves were accepted at least sometimes)
  # Use a relaxed test — at least one unique value beyond the starting value.
  expect_gt(length(unique(round(samp[, "class1_rate_log_sd"], 8))), 1L)
  expect_gt(length(unique(round(samp[, "class2_rate_log_sd"], 8))), 1L)
})


# ==============================================================================
# 2. unlink = "ratemultiplier"
# ==============================================================================

test_that("unlink=ratemultiplier chain runs, w columns present on simplex", {
  d <- .setup_2class_data()
  set.seed(43L)
  result <- RunMkPrime(
    data      = d$mkd,
    tree      = d$tree,
    mcmc      = .short_mcmc(),
    partition = d$part,
    unlink    = "ratemultiplier"
  )

  samp <- result$samples

  # (a) All LLs finite
  expect_true(all(is.finite(samp[, "log_likelihood"])))
  expect_gt(result$nSamples, 0L)

  # (b) Per-class columns present
  expect_true("w_1" %in% colnames(samp))
  expect_true("w_2" %in% colnames(samp))

  # (c) Simplex constraint: w_1 + w_2 ~= 1 every sample
  w_sum <- samp[, "w_1"] + samp[, "w_2"]
  expect_true(all(abs(w_sum - 1.0) < 1e-8))

  # (d) Weights are in (0, 1)
  expect_true(all(samp[, "w_1"] > 0 & samp[, "w_1"] < 1))
  expect_true(all(samp[, "w_2"] > 0 & samp[, "w_2"] < 1))

  # (e) Values evolve (at least some acceptance)
  expect_gt(length(unique(round(samp[, "w_1"], 8))), 1L)
})


# ==============================================================================
# 3. unlink = c("shape", "ratemultiplier")
# ==============================================================================

test_that("unlink=c(shape,ratemultiplier) chain runs, all per-class columns present", {
  d <- .setup_2class_data()
  set.seed(44L)
  result <- RunMkPrime(
    data      = d$mkd,
    tree      = d$tree,
    mcmc      = .short_mcmc(),
    partition = d$part,
    unlink    = c("shape", "ratemultiplier")
  )

  samp <- result$samples

  # (a) All LLs finite
  expect_true(all(is.finite(samp[, "log_likelihood"])))

  # (b) All per-class columns present
  for (col in c("class1_rate_log_sd", "class2_rate_log_sd", "w_1", "w_2")) {
    expect_true(col %in% colnames(samp), label = paste("column", col, "present"))
  }

  # (c) Simplex constraint holds
  w_sum <- samp[, "w_1"] + samp[, "w_2"]
  expect_true(all(abs(w_sum - 1.0) < 1e-8))

  # (d) Per-class shape values are positive
  expect_true(all(samp[, "class1_rate_log_sd"] > 0))
  expect_true(all(samp[, "class2_rate_log_sd"] > 0))
})
