# Tests for Beta-Geometric hyperparameter moves — (s, r) reparameterisation.
#
# Pre-2026-05-28: axis-aligned Bactrian + slice on raw (log α, log β) — could
# not traverse the BG posterior ridge; failed T9 SBC at AD p = 0 on both axes.
#
# Now: univariate slice on (s, r) = (log(α + β), logit(α/(α+β))) decorrelates
# the ridge across the parameter space. State persists as canonical (α, β);
# (s, r) is computed inside the move and the derived (α, β) is stored back.
skip_if_not_installed("phangorn")

# --- Helper -------------------------------------------------------------------
make_bg_fixture <- function() {
  tree <- ape::rtree(6)
  tree$edge.length <- tree$edge.length / sum(tree$edge.length) * 3
  nChar <- 5L
  mat <- matrix(sample(0:2, 6 * nChar, replace = TRUE), nrow = 6,
                dimnames = list(tree$tip.label, NULL))
  pd <- phangorn::phyDat(mat, type = "USER", levels = as.character(0:2))
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(kPrimePrior = "beta_geometric", treeLengthRate = 0.5)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  tree <- TreeTools::Preorder(tree)
  list(tree = tree, mkd = mkd, model = model)
}

# --- Slice sampler move construction ------------------------------------------

test_that(".BuildMoves registers slice_kprime_s and slice_kprime_r for BG", {
  f <- make_bg_fixture()
  nTrans <- sum(f$mkd$type == "transformational")
  nEdge <- nrow(f$tree$edge)
  mcmc <- suppressWarnings(MkPrimeMCMC(nIter = 100))

  moves <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE, mcmc = mcmc,
                                  kPrimePrior = "beta_geometric")
  moveNames <- vapply(moves, `[[`, character(1), "name")

  expect_true("slice_kprime_s" %in% moveNames)
  expect_true("slice_kprime_r" %in% moveNames)

  # The retired axis-aligned moves must NOT be registered.
  expect_false("kprime_alpha" %in% moveNames)
  expect_false("kprime_beta" %in% moveNames)
  expect_false("slice_kprime_alpha" %in% moveNames)
  expect_false("slice_kprime_beta" %in% moveNames)
})

test_that("slice_kprime_s and slice_kprime_r carry correct sliceParamIdx", {
  f <- make_bg_fixture()
  nTrans <- sum(f$mkd$type == "transformational")
  nEdge <- nrow(f$tree$edge)
  mcmc <- suppressWarnings(MkPrimeMCMC(nIter = 100))

  moves <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE, mcmc = mcmc,
                                  kPrimePrior = "beta_geometric")
  moveNames <- vapply(moves, `[[`, character(1), "name")

  sMove <- moves[[which(moveNames == "slice_kprime_s")]]
  rMove <- moves[[which(moveNames == "slice_kprime_r")]]

  # 0 → s = log(α+β); 1 → r = log(α/β) — must match C++ paramCode dispatch.
  expect_equal(sMove$sliceParamIdx, 0L)
  expect_equal(rMove$sliceParamIdx, 1L)

  # Both share the same C++ implementation (slice_kprime_hyper).
  expect_equal(sMove$type, "slice_kprime_hyper")
  expect_equal(rMove$type, "slice_kprime_hyper")
})

test_that(".BuildMoves does NOT include BG slice for geometric prior", {
  f <- make_bg_fixture()
  nTrans <- sum(f$mkd$type == "transformational")
  nEdge <- nrow(f$tree$edge)
  mcmc <- suppressWarnings(MkPrimeMCMC(nIter = 100))

  moves <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE, mcmc = mcmc,
                                  kPrimePrior = "geometric")
  moveNames <- vapply(moves, `[[`, character(1), "name")

  expect_false("slice_kprime_s" %in% moveNames)
  expect_false("slice_kprime_r" %in% moveNames)
  expect_true("p" %in% moveNames)
})

# --- Slice width adaptation ---------------------------------------------------

test_that(".AdaptSliceWidths handles slice_kprime_s/r widths", {
  moves <- list(
    list(name = "slice_kprime_s", type = "slice_kprime_hyper",
         target = "kprime_s", weight = 1, dim = 1L,
         sliceParamIdx = 0L),
    list(name = "slice_kprime_r", type = "slice_kprime_hyper",
         target = "kprime_r", weight = 1, dim = 1L,
         sliceParamIdx = 1L)
  )
  tuning <- list(slice_width_kprime_s = 1.0, slice_width_kprime_r = 1.0)
  proposeCount <- c(slice_kprime_s = 50L, slice_kprime_r = 50L)
  # High expansion count on both → widths too narrow → should increase.
  sliceExpCount <- c(slice_kprime_s = 300, slice_kprime_r = 300)

  result <- MkPrime:::.AdaptSliceWidths(
    tuning, proposeCount, sliceExpCount, moves
  )
  expect_gt(result$slice_width_kprime_s, tuning$slice_width_kprime_s)
  expect_gt(result$slice_width_kprime_r, tuning$slice_width_kprime_r)
})

test_that("MkPrimeMCMC defaults expose slice_width_kprime_s/r", {
  mcmc <- suppressWarnings(MkPrimeMCMC(nIter = 100))
  expect_equal(mcmc$tuning$slice_width_kprime_s, 0.5)
  expect_equal(mcmc$tuning$slice_width_kprime_r, 1.0)
  # The retired tuning keys must be gone.
  expect_null(mcmc$tuning$scale_kprime_alpha)
  expect_null(mcmc$tuning$scale_kprime_beta)
  expect_null(mcmc$tuning$slice_width_kprime_alpha)
  expect_null(mcmc$tuning$slice_width_kprime_beta)
})

# --- Integration: short BG run produces non-constant α, β --------------------

test_that("Short BG run moves kprime_alpha and kprime_beta", {
  set.seed(6147)
  f <- make_bg_fixture()
  mcmc <- suppressWarnings(MkPrimeMCMC(
    nIter = 500, thin = 1, nChains = 1L,
    maxWarmup = 100, minWarmup = 50
  ))

  post <- suppressWarnings(RunMkPrime(
    f$mkd, f$tree, model = f$model, mcmc = mcmc
  ))
  samp <- post$samples

  alpha_vals <- samp[, "kprime_alpha"]
  beta_vals <- samp[, "kprime_beta"]

  # Must have more than 1 unique value (not frozen).
  expect_gt(length(unique(alpha_vals)), 1,
            label = "kprime_alpha is not frozen at a single value")
  expect_gt(length(unique(beta_vals)), 1,
            label = "kprime_beta is not frozen at a single value")
})

# --- Integration: β drives down toward small-β posterior under informative data
#
# This is the regression detector for the (log α, log β) ridge that the pre-
# fix axis-aligned slice / Bactrian moves could not traverse.
#
# Construction: ~40 fully binary characters on a 10-tip tree. With Gibbs
# pulling every kPrime to kObs = 2, all u_i = 0, and the BG posterior on β
# has a likelihood factor ∝ (α / (α + β))^n that — at α ≈ 1, n ≈ 40 — is
# sharply peaked at β → 0 (modal β ≈ 1/n ≈ 0.025). Posterior mean of log β
# should land below −1 (β < 0.37) given an Exp(1) prior on β.
#
# Init at (α, β) = (1, 1) starts the chain at log β = 0, a factor of ~14×
# above the posterior mode. The axis-aligned sampler (cases 27/28 + axis-
# slice) needs O(κ²) iterations to traverse this ridge — at n = 40, κ on
# (log α, log β) is ~10–15, so 2000 iter is not enough. The (s, r)
# reparameterised sampler reaches the small-β regime within hundreds.
test_that("β drives down to small-β regime under informative binary data", {
  set.seed(6148)
  nTip  <- 10L
  nChar <- 40L
  tree <- ape::rtree(nTip)
  tree$edge.length <- tree$edge.length / sum(tree$edge.length) * 3
  # Binary characters with mild variation — every column has kObs = 2.
  mat <- matrix(sample(0:1, nTip * nChar, replace = TRUE), nrow = nTip,
                dimnames = list(tree$tip.label, NULL))
  # Force at least one of each state per column (guarantees kObs = 2).
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1, j] <- 1L - mat[1, j]
  }
  pd <- phangorn::phyDat(mat, type = "USER", levels = as.character(0:1))
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(kPrimePrior = "beta_geometric", treeLengthRate = 0.5)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  tree <- TreeTools::Preorder(tree)

  mcmc <- suppressWarnings(MkPrimeMCMC(
    nIter = 2000, thin = 1, nChains = 1L,
    maxWarmup = 500, minWarmup = 200,
    autoTune = FALSE
  ))

  post <- suppressWarnings(RunMkPrime(
    mkd, tree, model = model, mcmc = mcmc
  ))
  beta_vals <- post$samples[, "kprime_beta"]
  # Drop warmup-tail samples (first 25%) before computing the mean.
  keep <- beta_vals[round(length(beta_vals) * 0.25):length(beta_vals)]

  expect_lt(mean(log(keep)), -1.0,
            label = paste0("posterior mean of log β stuck at ",
                           sprintf("%.2f", mean(log(keep))),
                           " — sampler not traversing (log α, log β) ridge"))
})
