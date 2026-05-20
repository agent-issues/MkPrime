# Tests for Beta-Geometric hyperparameter (α, β) moves — M-163
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

# --- Adaptation fix -----------------------------------------------------------

test_that(".AdaptTuning handles kprime_alpha/beta targets", {
  moves <- list(
    list(name = "tree_length", type = "scale", target = "tree_length",
         weight = 1, dim = 1L),
    list(name = "kprime_alpha", type = "kprime_alpha",
         target = "kprime_alpha", weight = 0.05, dim = 1L),
    list(name = "kprime_beta", type = "kprime_beta",
         target = "kprime_beta", weight = 0.05, dim = 1L)
  )

  tuning <- list(
    scale_tree_length = 0.5,
    scale_kprime_alpha = 0.05,
    scale_kprime_beta = 0.05
  )

  # 50% acceptance — above 35% target → scale should increase
  acceptCount <- c(tree_length = 25L, kprime_alpha = 25L, kprime_beta = 25L)
  proposeCount <- c(tree_length = 50L, kprime_alpha = 50L, kprime_beta = 50L)

  result <- MkPrime:::.AdaptTuning(tuning, acceptCount, proposeCount, moves)

  # Key: tuning must not be NA (was the bug before M-163)
  expect_false(is.na(result$scale_kprime_alpha))
  expect_false(is.na(result$scale_kprime_beta))
  expect_true(is.finite(result$scale_kprime_alpha))
  expect_true(is.finite(result$scale_kprime_beta))

  # With 50% acceptance > 35% target, scale should increase
  expect_gt(result$scale_kprime_alpha, tuning$scale_kprime_alpha)
  expect_gt(result$scale_kprime_beta, tuning$scale_kprime_beta)
})

test_that(".AdaptTuning shrinks kprime scales when acceptance too low", {
  moves <- list(
    list(name = "kprime_alpha", type = "kprime_alpha",
         target = "kprime_alpha", weight = 0.05, dim = 1L)
  )

  tuning <- list(scale_kprime_alpha = 0.5)
  # 5% acceptance — well below 35% target
  acceptCount <- c(kprime_alpha = 5L)
  proposeCount <- c(kprime_alpha = 100L)

  result <- MkPrime:::.AdaptTuning(tuning, acceptCount, proposeCount, moves)
  expect_lt(result$scale_kprime_alpha, tuning$scale_kprime_alpha)
})

test_that("MkPrimeMCMC default scale_kprime_alpha/beta matches tuned defaults", {
  mcmc <- suppressWarnings(MkPrimeMCMC(nIter = 100))
  # M-165: increased from 0.05 to better match posterior concentration
  expect_equal(mcmc$tuning$scale_kprime_alpha, 0.3)
  expect_equal(mcmc$tuning$scale_kprime_beta, 0.5)
})

# --- Slice sampler move construction ------------------------------------------

test_that(".BuildMoves includes slice_kprime_alpha/beta for BG prior", {
  f <- make_bg_fixture()
  nTrans <- sum(f$mkd$type == "transformational")
  nEdge <- nrow(f$tree$edge)
  mcmc <- suppressWarnings(MkPrimeMCMC(nIter = 100))

  moves <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE, mcmc = mcmc,
                                  kPrimePrior = "beta_geometric")
  moveNames <- vapply(moves, `[[`, character(1), "name")

  expect_true("slice_kprime_alpha" %in% moveNames)
  expect_true("slice_kprime_beta" %in% moveNames)
  # Scale proposals also present
  expect_true("kprime_alpha" %in% moveNames)
  expect_true("kprime_beta" %in% moveNames)
})

test_that("slice_kprime_alpha/beta have correct sliceParamIdx", {
  f <- make_bg_fixture()
  nTrans <- sum(f$mkd$type == "transformational")
  nEdge <- nrow(f$tree$edge)
  mcmc <- suppressWarnings(MkPrimeMCMC(nIter = 100))

  moves <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE, mcmc = mcmc,
                                  kPrimePrior = "beta_geometric")
  moveNames <- vapply(moves, `[[`, character(1), "name")

  alphaMove <- moves[[which(moveNames == "slice_kprime_alpha")]]
  betaMove <- moves[[which(moveNames == "slice_kprime_beta")]]

  expect_equal(alphaMove$sliceParamIdx, 0L)
  expect_equal(betaMove$sliceParamIdx, 1L)
})

test_that(".BuildMoves does NOT include BG slice for geometric prior", {
  f <- make_bg_fixture()
  nTrans <- sum(f$mkd$type == "transformational")
  nEdge <- nrow(f$tree$edge)
  mcmc <- suppressWarnings(MkPrimeMCMC(nIter = 100))

  moves <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE, mcmc = mcmc,
                                  kPrimePrior = "geometric")
  moveNames <- vapply(moves, `[[`, character(1), "name")

  expect_false("slice_kprime_alpha" %in% moveNames)
  expect_false("slice_kprime_beta" %in% moveNames)
  expect_false("kprime_alpha" %in% moveNames)
  expect_false("kprime_beta" %in% moveNames)
  expect_true("p" %in% moveNames)
})

# --- Slice width adaptation ---------------------------------------------------

test_that(".AdaptSliceWidths handles kprime hyper slice widths", {
  moves <- list(
    list(name = "slice_kprime_alpha", type = "slice_kprime_hyper",
         target = "kprime_alpha", weight = 1, dim = 1L,
         sliceParamIdx = 0L)
  )
  tuning <- list(slice_width_kprime_alpha = 1.0)
  proposeCount <- c(slice_kprime_alpha = 50L)
  # High expansion count → width too narrow → should increase
  sliceExpCount <- c(slice_kprime_alpha = 300)

  result <- MkPrime:::.AdaptSliceWidths(
    tuning, proposeCount, sliceExpCount, moves
  )
  expect_gt(result$slice_width_kprime_alpha, tuning$slice_width_kprime_alpha)
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

  # Must have more than 1 unique value (not frozen)
  expect_gt(length(unique(alpha_vals)), 1,
            label = "kprime_alpha is not frozen at a single value")
  expect_gt(length(unique(beta_vals)), 1,
            label = "kprime_beta is not frozen at a single value")
})
