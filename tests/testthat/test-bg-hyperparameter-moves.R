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

# --- Integration: β can escape init at 1.0 toward small-β regime -------------
#
# Pre-fix axis-aligned regression: with init at (α, β) = (1, 1) and a posterior
# that favours small β, the chain could not traverse the (log α, log β) ridge
# and β remained stuck near 1. Post-fix, the s-slice moves the concentration
# down, the r-slice moves the shape away from α = β, and β can land
# significantly below its starting value.
test_that("β escapes init = 1.0 within a short BG run", {
  set.seed(6148)
  f <- make_bg_fixture()
  mcmc <- suppressWarnings(MkPrimeMCMC(
    nIter = 2000, thin = 1, nChains = 1L,
    maxWarmup = 500, minWarmup = 200,
    autoTune = FALSE
  ))

  post <- suppressWarnings(RunMkPrime(
    f$mkd, f$tree, model = f$model, mcmc = mcmc
  ))
  beta_vals <- post$samples[, "kprime_beta"]
  # log range must span more than 1 unit (factor of e ≈ 2.7) — the pre-fix
  # axis-aligned sampler typically gives < 0.3 here at the same init.
  expect_gt(diff(range(log(beta_vals))), 1.0,
            label = "β log-range too narrow — sampler stuck on ridge")
})
