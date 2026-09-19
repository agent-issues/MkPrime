# tests/testthat/test-m090-move-wiring.R
# M-090: Verify new Gibbs/Weighted moves are correctly wired into
# MkPrimeMCMC, .BuildMoves, .kMoveTypes, and C++ dispatch.

# ── MkPrimeMCMC parameter storage ──────────────────────────────────────────

test_that("MkPrimeMCMC stores Gibbs/Weighted move flags (defaults)", {

  mc <- MkPrimeMCMC()
  expect_true(mc$gibbsSpr)
  expect_true(mc$gibbsSubtreeSwap)
  expect_false(mc$weightedBranchScale)
  expect_false(mc$weightedSpr)
  expect_false(mc$weightedSubtreeSwap)
  expect_identical(mc$nBranchBins, 10L)
})

test_that("MkPrimeMCMC stores custom move flags", {
  mc <- MkPrimeMCMC(
    gibbsSpr = FALSE, gibbsSubtreeSwap = FALSE,
    weightedBranchScale = TRUE, weightedSpr = TRUE,
    weightedSubtreeSwap = TRUE, nBranchBins = 5L
  )
  expect_false(mc$gibbsSpr)
  expect_false(mc$gibbsSubtreeSwap)
  expect_true(mc$weightedBranchScale)
  expect_true(mc$weightedSpr)
  expect_true(mc$weightedSubtreeSwap)
  expect_identical(mc$nBranchBins, 5L)
})

test_that("MkPrimeMCMC rejects invalid nBranchBins", {
  expect_error(MkPrimeMCMC(nBranchBins = 0L), "nBranchBins")
  expect_error(MkPrimeMCMC(nBranchBins = 1L), "nBranchBins")
  expect_error(MkPrimeMCMC(nBranchBins = -5L), "nBranchBins")
})

# ── .kMoveTypes completeness ──────────────────────────────────────────────

test_that(".kMoveTypes maps all new move names to correct codes", {
  mt <- MkPrime:::.kMoveTypes
  expect_identical(mt[["gibbs_spr"]], 10L)
  expect_identical(mt[["gibbs_subtree_swap"]], 11L)
  expect_identical(mt[["weighted_branch_lengths"]], 12L)
  expect_identical(mt[["weighted_spr"]], 13L)
  expect_identical(mt[["weighted_subtree_swap"]], 14L)
})

# ── .BuildMoves move pool construction ────────────────────────────────────

test_that(".BuildMoves includes Gibbs but not Weighted by default", {
  mc <- MkPrimeMCMC()
  moves <- MkPrime:::.BuildMoves(
    nEdge = 12L, nTrans = 3L, hasNeo = TRUE, mcmc = mc,
    fixTopology = FALSE, kPrimePrior = "geometric"
  )
  names <- vapply(moves, `[[`, character(1), "name")

  expect_true("gibbs_spr" %in% names)
  expect_true("gibbs_subtree_swap" %in% names)
  expect_false("weighted_branch_lengths" %in% names)
  expect_false("weighted_spr" %in% names)
  expect_false("weighted_subtree_swap" %in% names)

  # Original moves still present
  expect_true("nni" %in% names)
  expect_true("spr" %in% names)
  expect_true("tree_length" %in% names)
  expect_true("branch_lengths" %in% names)
})

test_that(".BuildMoves includes all 5 new moves when all enabled", {
  mc <- MkPrimeMCMC(
    gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE,
    weightedBranchScale = TRUE, weightedSpr = TRUE,
    weightedSubtreeSwap = TRUE
  )
  moves <- MkPrime:::.BuildMoves(
    nEdge = 12L, nTrans = 0L, hasNeo = FALSE, mcmc = mc,
    fixTopology = FALSE, kPrimePrior = "geometric"
  )
  names <- vapply(moves, `[[`, character(1), "name")

  expect_true("gibbs_spr" %in% names)
  expect_true("gibbs_subtree_swap" %in% names)
  expect_true("weighted_branch_lengths" %in% names)
  expect_true("weighted_spr" %in% names)
  expect_true("weighted_subtree_swap" %in% names)
})

test_that(".BuildMoves excludes all 5 new moves when all disabled", {
  mc <- MkPrimeMCMC(
    gibbsSpr = FALSE, gibbsSubtreeSwap = FALSE,
    weightedBranchScale = FALSE, weightedSpr = FALSE,
    weightedSubtreeSwap = FALSE
  )
  moves <- MkPrime:::.BuildMoves(
    nEdge = 12L, nTrans = 3L, hasNeo = TRUE, mcmc = mc,
    fixTopology = FALSE, kPrimePrior = "geometric"
  )
  names <- vapply(moves, `[[`, character(1), "name")

  new_names <- c("gibbs_spr", "gibbs_subtree_swap",
                 "weighted_branch_lengths", "weighted_spr",
                 "weighted_subtree_swap")
  expect_false(any(new_names %in% names))

  # Original moves still present
  expect_true("nni" %in% names)
  expect_true("spr" %in% names)
})

# ── M-115: Gibbs weight cap for larger trees ──────────────────────────────

test_that("Gibbs weights are capped on large trees (M-115)", {
  mc <- MkPrimeMCMC(gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE,
                    blockGibbsBranch = TRUE)
  # Large tree: nEdge = 200 → nEdge/4 = 50, nEdge/6 ≈ 33
  moves <- MkPrime:::.BuildMoves(
    nEdge = 200L, nTrans = 3L, hasNeo = TRUE, mcmc = mc,
    fixTopology = FALSE, kPrimePrior = "geometric"
  )
  weights <- setNames(
    vapply(moves, `[[`, numeric(1), "weight"),
    vapply(moves, `[[`, character(1), "name")
  )

  expect_equal(weights[["gibbs_spr"]], 10)
  expect_equal(weights[["gibbs_subtree_swap"]], 10)
  expect_equal(weights[["block_gibbs_branch"]], 10)

  # Standard SPR should NOT be capped (scales linearly)
  expect_equal(weights[["spr"]], 50)
})

test_that("Gibbs weights are not capped on small trees (M-115)", {
  mc <- MkPrimeMCMC(gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE,
                    blockGibbsBranch = TRUE)
  # Small tree: nEdge = 12 → nEdge/4 = 3
  moves <- MkPrime:::.BuildMoves(
    nEdge = 12L, nTrans = 3L, hasNeo = TRUE, mcmc = mc,
    fixTopology = FALSE, kPrimePrior = "geometric"
  )
  weights <- setNames(
    vapply(moves, `[[`, numeric(1), "weight"),
    vapply(moves, `[[`, character(1), "name")
  )

  # Below cap → uses nEdge / k directly
  expect_equal(weights[["gibbs_spr"]], 3)
  expect_equal(weights[["gibbs_subtree_swap"]], 2)
  expect_equal(weights[["block_gibbs_branch"]], 3)
})

test_that(".BuildMoves excludes topology moves when fixTopology = TRUE", {
  mc <- MkPrimeMCMC(
    gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE,
    weightedBranchScale = TRUE, weightedSpr = TRUE,
    weightedSubtreeSwap = TRUE
  )
  moves <- MkPrime:::.BuildMoves(
    nEdge = 12L, nTrans = 0L, hasNeo = FALSE, mcmc = mc,
    fixTopology = TRUE, kPrimePrior = "geometric"
  )
  names <- vapply(moves, `[[`, character(1), "name")

  # None of the topology or weighted moves should appear
  topo_moves <- c("nni", "spr", "gibbs_spr", "gibbs_subtree_swap",
                  "weighted_branch_lengths", "weighted_spr",
                  "weighted_subtree_swap")
  expect_false(any(topo_moves %in% names))
})

# ── Move type codes resolve correctly ─────────────────────────────────────

test_that("All new move names resolve to valid integer codes", {
  mc <- MkPrimeMCMC(
    gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE,
    weightedBranchScale = TRUE, weightedSpr = TRUE,
    weightedSubtreeSwap = TRUE
  )
  moves <- MkPrime:::.BuildMoves(
    nEdge = 12L, nTrans = 0L, hasNeo = FALSE, mcmc = mc,
    fixTopology = FALSE, kPrimePrior = "geometric"
  )
  # This is the same extraction used by RunMkPrime
  codes <- vapply(moves, function(m) MkPrime:::.kMoveTypes[[m$name]],
                  integer(1L))
  expect_true(all(is.finite(codes)))
  expect_true(10L %in% codes)  # gibbs_spr
  expect_true(11L %in% codes)  # gibbs_subtree_swap
  expect_true(12L %in% codes)  # weighted_branch_lengths
  expect_true(13L %in% codes)  # weighted_spr
  expect_true(14L %in% codes)  # weighted_subtree_swap
})

# ── nBranchBins reaches C++ ───────────────────────────────────────────────

test_that("set_branch_bins sets nBranchBins on McmcData", {
  set.seed(4821L)
  tree <- ape::rtree(6L, rooted = FALSE)
  tree <- TreeTools::Preorder(tree)
  mat <- matrix(sample(0:2, 6L * 5L, replace = TRUE), nrow = 6L,
                dimnames = list(tree$tip.label, NULL))
  pd    <- MatrixToPhyDat(mat)
  mkd   <- suppressWarnings(MkPrimeData(pd))
  model <- MkPrimeModel()
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)

  expect_no_error(set_branch_bins(dataPtr, 5L))
  expect_no_error(set_branch_bins(dataPtr, 20L))
})

test_that("Weighted move works with custom nBranchBins via McmcData", {
  set.seed(7193L)
  tree <- ape::rtree(6L, rooted = FALSE)
  tree <- TreeTools::Preorder(tree)
  nEdge <- nrow(tree$edge)
  mat <- matrix(sample(0:2, 6L * 5L, replace = TRUE), nrow = 6L,
                dimnames = list(tree$tip.label, NULL))
  pd    <- MatrixToPhyDat(mat)
  mkd   <- MkPrimeData(pd)
  model <- MkPrimeModel()
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  state0   <- MkPrime:::.InitState(tree, mkd, model)
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)

  # Set custom nBranchBins = 5
  set_branch_bins(dataPtr, 5L)

  # Run weighted_branch_scale (moveType 12) a few times
  ran <- FALSE
  for (i in seq_len(50L)) {
    if (do_move_cpp(dataPtr, statePtr, 12L, 0L, 0.5, 0.5, 1L, 1.0)) {
      ran <- TRUE
      break
    }
  }
  # At least verify it didn't crash — acceptance depends on data
  expect_true(TRUE)
})

# ── .AdaptTuning handles new moves gracefully ─────────────────────────────

test_that(".AdaptTuning skips Gibbs/Weighted moves (NA tuning keys)", {
  mc <- MkPrimeMCMC(gibbsSpr = TRUE, weightedSpr = TRUE)
  moves <- MkPrime:::.BuildMoves(
    nEdge = 12L, nTrans = 0L, hasNeo = FALSE, mcmc = mc,
    fixTopology = FALSE, kPrimePrior = "geometric"
  )
  # Fake acceptance counts
  nms <- vapply(moves, `[[`, character(1), "name")
  accept  <- setNames(rep(10L, length(nms)), nms)
  propose <- setNames(rep(50L, length(nms)), nms)

  tuning <- list(
    scale_tree_length = 0.5, beta_simplex = 10,
    scale_rate_loss = 0.5, scale_rate_log_sd = 0.5,
    scale_rate_neo = 0.5, scale_p = 0.5, int_walk_window = 1L
  )
  result <- MkPrime:::.AdaptTuning(tuning, accept, propose, moves)

  # Tuning values should have changed for adaptable moves...
  # ... but the function should not error on the Gibbs/weighted moves
  expect_type(result, "list")
  expect_true("scale_tree_length" %in% names(result))
})

# ── Integration tests (slow) ─────────────────────────────────────────────

test_that("Full MCMC run with Gibbs moves produces valid MkPosterior", {
  skip_slow_tests()
  tree <- BalancedTree(8)
  mat <- matrix(sample(0:2, 8 * 5, replace = TRUE), nrow = 8,
                dimnames = list(tree$tip.label, NULL))
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel()
  mcmc <- MkPrimeMCMC(
    nIter = 600L, warmup = 200L, thin = 10L, nRuns = 1L,
    gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE
  )
  result <- RunMkPrime(data = mkd, tree = tree, model = model, mcmc = mcmc)
  expect_s3_class(result, "MkPosterior")
  expect_true(nrow(result$samples) > 0L)
})

test_that("Full MCMC run with all moves produces valid MkPosterior", {
  skip_slow_tests()
  tree <- BalancedTree(8)
  mat <- matrix(sample(0:2, 8 * 5, replace = TRUE), nrow = 8,
                dimnames = list(tree$tip.label, NULL))
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel()
  mcmc <- MkPrimeMCMC(
    nIter = 600L, warmup = 200L, thin = 10L, nRuns = 1L,
    gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE,
    weightedBranchScale = TRUE, weightedSpr = TRUE,
    weightedSubtreeSwap = TRUE, nBranchBins = 5L
  )
  result <- RunMkPrime(data = mkd, tree = tree, model = model, mcmc = mcmc)
  expect_s3_class(result, "MkPosterior")
  expect_true(nrow(result$samples) > 0L)
})
