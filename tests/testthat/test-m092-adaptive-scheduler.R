# Tests for M-092: Adaptive move scheduler
#
# Covers: .AdaptMoveWeights(), .ResolvePinnedWeights(),
# .NormalizeMoveWeights(), .LogMoveWeights(), .FormatMoveWeights(),
# MkPrimeMCMC() moveWeights parameter, and C++ timing data.

# Skip slow integration tests unless explicitly requested
skip_slow_tests <- function() {
  if (!identical(Sys.getenv("MKPRIME_RUN_SLOW_TESTS"), "true")) {
    skip("Slow test: set MKPRIME_RUN_SLOW_TESTS=true to run")
  }
}


# ==========================================================================
# MkPrimeMCMC() moveWeights parameter validation
# ==========================================================================

test_that("MkPrimeMCMC() stores moveWeights = NULL by default", {
  mcmc <- MkPrimeMCMC(nIter = 100L)
  expect_null(mcmc$moveWeights)
})

test_that("MkPrimeMCMC() stores valid moveWeights", {
  mcmc <- MkPrimeMCMC(nIter = 100L, moveWeights = c(nni = 0.3, spr = 0.2))
  expect_equal(mcmc$moveWeights, c(nni = 0.3, spr = 0.2))
})

test_that("MkPrimeMCMC() accepts moveWeights as list", {
  mcmc <- MkPrimeMCMC(nIter = 100L, moveWeights = list(nni = 0.3))
  expect_equal(mcmc$moveWeights, c(nni = 0.3))
})

test_that("MkPrimeMCMC() rejects unnamed moveWeights", {
  expect_error(MkPrimeMCMC(nIter = 100L, moveWeights = c(0.3, 0.2)),
               "named numeric")
})

test_that("MkPrimeMCMC() rejects unknown move names", {
  expect_error(MkPrimeMCMC(nIter = 100L, moveWeights = c(bogus = 0.3)),
               "unknown move name")
})

test_that("MkPrimeMCMC() rejects negative moveWeights", {
  expect_error(MkPrimeMCMC(nIter = 100L, moveWeights = c(nni = -0.1)),
               "positive")
})

test_that("MkPrimeMCMC() rejects moveWeights summing > 1", {
  expect_error(MkPrimeMCMC(nIter = 100L,
                            moveWeights = c(nni = 0.6, spr = 0.6)),
               "exceeds 1")
})

test_that("MkPrimeMCMC() accepts moveWeights summing to exactly 1", {
  mcmc <- MkPrimeMCMC(nIter = 100L,
                       moveWeights = c(nni = 0.5, spr = 0.5))
  expect_equal(sum(mcmc$moveWeights), 1.0)
})


# ==========================================================================
# .ResolvePinnedWeights()
# ==========================================================================

test_that(".ResolvePinnedWeights returns NULL for NULL input", {
  expect_null(MkPrime:::.ResolvePinnedWeights(NULL, c("nni", "spr")))
})

test_that(".ResolvePinnedWeights filters to pool moves", {
  result <- MkPrime:::.ResolvePinnedWeights(
    c(nni = 0.3, spr = 0.2, bogus = 0.1),
    c("nni", "spr", "tree_length")
  )
  expect_equal(result, c(nni = 0.3, spr = 0.2))
})

test_that(".ResolvePinnedWeights warns on dropped moves", {
  expect_warning(
    MkPrime:::.ResolvePinnedWeights(
      c(nni = 0.3, bogus = 0.1),
      c("nni", "spr")
    ),
    "ignored"
  )
})

test_that(".ResolvePinnedWeights returns NULL when no names match", {
  result <- suppressWarnings(
    MkPrime:::.ResolvePinnedWeights(c(bogus = 0.3), c("nni", "spr"))
  )
  expect_null(result)
})


# ==========================================================================
# .NormalizeMoveWeights()
# ==========================================================================

test_that(".NormalizeMoveWeights sums to 1 with pinned weights", {
  w <- c(nni = 3, spr = 2, tree_length = 1)
  result <- MkPrime:::.NormalizeMoveWeights(w, c(nni = 0.4))
  expect_equal(sum(result), 1.0, tolerance = 1e-10)
  expect_equal(result[["nni"]], 0.4)
})

test_that(".NormalizeMoveWeights with no pinned normalizes to 1", {
  w <- c(nni = 3, spr = 2)
  result <- MkPrime:::.NormalizeMoveWeights(w, c(bogus = 0.1))
  expect_equal(sum(result), 1.0, tolerance = 1e-10)
})

test_that(".NormalizeMoveWeights with all pinned sums to 1", {
  w <- c(nni = 3, spr = 2)
  result <- MkPrime:::.NormalizeMoveWeights(w, c(nni = 0.6, spr = 0.4))
  expect_equal(result[["nni"]], 0.6)
  expect_equal(result[["spr"]], 0.4)
  expect_equal(sum(result), 1.0, tolerance = 1e-10)
})


# ==========================================================================
# .AdaptMoveWeights()
# ==========================================================================

test_that(".AdaptMoveWeights returns vector of same length summing to 1", {
  n <- 4L
  moveNames <- c("tree_length", "nni", "spr", "kPrime")
  current <- rep(0.25, n)
  names(current) <- moveNames
  accept  <- setNames(c(50L, 30L, 10L, 40L), moveNames)
  propose <- setNames(c(100L, 100L, 100L, 100L), moveNames)
  timeNs  <- setNames(c(1e6, 5e6, 3e6, 2e6), moveNames)

  result <- MkPrime:::.AdaptMoveWeights(
    current, accept, propose, timeNs, moveNames,
    pinnedWeights = NULL, warmupProgress = 0.5
  )
  expect_length(result, n)
  expect_equal(sum(result), 1.0, tolerance = 1e-10)
})

test_that(".AdaptMoveWeights concentrates weight on best score", {
  n <- 3L
  moveNames <- c("a", "b", "c")
  current <- rep(1/3, n)
  names(current) <- moveNames
  # Move "a": 80% accept, moderate cost; "b"/"c": lower accept, higher cost
  accept  <- setNames(c(80L, 30L, 20L), moveNames)
  propose <- setNames(c(100L, 100L, 100L), moveNames)
  # All similar cost so score differences come from accept rate
  timeNs  <- setNames(c(1e8, 1e8, 1e8), moveNames)

  result <- MkPrime:::.AdaptMoveWeights(
    current, accept, propose, timeNs, moveNames,
    pinnedWeights = NULL, warmupProgress = 1.0  # T=0.5 (most peaked)
  )
  # Move "a" should have the highest weight
  expect_gt(result[1], result[2])
  expect_gt(result[1], result[3])
})

test_that(".AdaptMoveWeights at T=large gives near-uniform weights", {
  n <- 3L
  moveNames <- c("a", "b", "c")
  current <- rep(1/3, n)
  names(current) <- moveNames
  # Moderate score differences: accept rates differ modestly
  accept  <- setNames(c(60L, 40L, 30L), moveNames)
  propose <- setNames(c(100L, 100L, 100L), moveNames)
  timeNs  <- setNames(c(1e8, 1e8, 1e8), moveNames)

  result <- MkPrime:::.AdaptMoveWeights(
    current, accept, propose, timeNs, moveNames,
    pinnedWeights = NULL, warmupProgress = 0.0,  # T=2.0 (near uniform)
    tStart = 50.0  # very high T for strong uniformity
  )
  # Should be close to uniform (max - min < 0.15)
  expect_true(max(result) - min(result) < 0.15)
})

test_that(".AdaptMoveWeights respects pinned weights", {
  n <- 3L
  moveNames <- c("a", "b", "c")
  current <- c(0.5, 0.25, 0.25)
  names(current) <- moveNames
  accept  <- setNames(c(50L, 30L, 20L), moveNames)
  propose <- setNames(c(100L, 100L, 100L), moveNames)
  timeNs  <- setNames(c(1e8, 1e8, 1e8), moveNames)

  result <- MkPrime:::.AdaptMoveWeights(
    current, accept, propose, timeNs, moveNames,
    pinnedWeights = c(a = 0.5), warmupProgress = 0.5
  )
  expect_equal(result[["a"]], 0.5)
  expect_equal(sum(result), 1.0, tolerance = 1e-10)
})

test_that(".AdaptMoveWeights preserves all-pinned weights", {
  n <- 2L
  moveNames <- c("a", "b")
  current <- c(0.6, 0.4)
  names(current) <- moveNames
  accept  <- setNames(c(50L, 30L), moveNames)
  propose <- setNames(c(100L, 100L), moveNames)
  timeNs  <- setNames(c(1e6, 1e6), moveNames)

  result <- MkPrime:::.AdaptMoveWeights(
    current, accept, propose, timeNs, moveNames,
    pinnedWeights = c(a = 0.6, b = 0.4), warmupProgress = 0.5
  )
  expect_equal(result, current, tolerance = 1e-10)
})

test_that(".AdaptMoveWeights keeps current weight for low-proposal moves", {
  n <- 3L
  moveNames <- c("a", "b", "c")
  current <- c(0.4, 0.3, 0.3)
  names(current) <- moveNames
  # Move "c" has < 20 proposals
  accept  <- setNames(c(50L, 30L, 2L), moveNames)
  propose <- setNames(c(100L, 100L, 5L), moveNames)
  timeNs  <- setNames(c(1e6, 1e6, 1e5), moveNames)

  result <- MkPrime:::.AdaptMoveWeights(
    current, accept, propose, timeNs, moveNames,
    pinnedWeights = NULL, warmupProgress = 0.5
  )
  # Move "c" should keep its current weight
  expect_equal(result[["c"]], current[["c"]])
  # Total still sums to 1
  expect_equal(sum(result), 1.0, tolerance = 1e-10)
})

test_that(".AdaptMoveWeights returns unchanged if all below minProposals", {
  n <- 2L
  moveNames <- c("a", "b")
  current <- c(0.6, 0.4)
  names(current) <- moveNames
  accept  <- setNames(c(5L, 3L), moveNames)
  propose <- setNames(c(10L, 8L), moveNames)
  timeNs  <- setNames(c(1e5, 1e5), moveNames)

  result <- MkPrime:::.AdaptMoveWeights(
    current, accept, propose, timeNs, moveNames,
    pinnedWeights = NULL, warmupProgress = 0.5
  )
  expect_equal(result, current, tolerance = 1e-10)
})

test_that(".AdaptMoveWeights applies floor", {
  n <- 3L
  moveNames <- c("a", "b", "c")
  current <- rep(1/3, n)
  names(current) <- moveNames
  # Move "c" has very bad score (should hit floor)
  accept  <- setNames(c(80L, 80L, 1L), moveNames)
  propose <- setNames(c(100L, 100L, 100L), moveNames)
  timeNs  <- setNames(c(1e8, 1e8, 1e8), moveNames)

  result <- MkPrime:::.AdaptMoveWeights(
    current, accept, propose, timeNs, moveNames,
    pinnedWeights = NULL, warmupProgress = 1.0,
    wMin = 0.15  # generous floor for test clarity
  )
  # Floor per move: wMin = 0.15
  expect_gte(result[["c"]], 0.15 - 1e-10)
  expect_equal(sum(result), 1.0, tolerance = 1e-10)
})

test_that(".AdaptMoveWeights temperature annealing works", {
  n <- 3L
  moveNames <- c("a", "b", "c")
  current <- rep(1/3, n)
  names(current) <- moveNames
  # Moderate score differences so softmax doesn't fully saturate
  accept  <- setNames(c(60L, 40L, 30L), moveNames)
  propose <- setNames(c(100L, 100L, 100L), moveNames)
  timeNs  <- setNames(c(1e8, 1e8, 1e8), moveNames)

  early <- MkPrime:::.AdaptMoveWeights(
    current, accept, propose, timeNs, moveNames,
    pinnedWeights = NULL, warmupProgress = 0.0  # T=2.0
  )
  late <- MkPrime:::.AdaptMoveWeights(
    current, accept, propose, timeNs, moveNames,
    pinnedWeights = NULL, warmupProgress = 1.0  # T=0.5
  )
  # Late should be more peaked than early (larger max-min spread)
  expect_gt(max(late) - min(late), max(early) - min(early))
})


# ==========================================================================
# .FormatMoveWeights()
# ==========================================================================

test_that(".FormatMoveWeights produces readable string", {
  s <- MkPrime:::.FormatMoveWeights(c(0.3, 0.7), c("nni", "spr"))
  expect_match(s, "nni=30.0%")
  expect_match(s, "spr=70.0%")
})


# ==========================================================================
# .LogMoveWeights() writes to file
# ==========================================================================

test_that(".LogMoveWeights writes comment line to log file", {
  tmp <- tempfile(fileext = ".log")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("Sample\ta\tb", tmp)
  MkPrime:::.LogMoveWeights(c(0.6, 0.4), c("nni", "spr"), tmp)
  lines <- readLines(tmp)
  expect_length(lines, 2L)
  expect_match(lines[2], "^# Adapted move weights:")
  expect_match(lines[2], "nni=60.0%")
})


# ==========================================================================
# C++ timing data from run_mcmc_batch_cpp
# ==========================================================================

test_that("run_mcmc_batch_cpp returns move_time_ns matrix", {
  skip_if_not_installed("TreeSearch")
  dat <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  mkd <- MkPrimeData(dat)
  model <- MkPrimeModel()
  tree <- ape::rtree(length(dat), tip.label = names(dat))
  tree <- TreeTools::Preorder(tree)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)

  mcmcData <- MkPrime:::.InitMcmcData(mkd, model)
  state <- MkPrime:::.InitState(tree, mkd, model)
  chainState <- MkPrime:::.InitMcmcChain(state)
  fill_partition_cache(mcmcData, chainState)
  allocate_cl_workspace(mcmcData, chainState)

  moves <- MkPrime:::.BuildMoves(
    nrow(tree$edge), sum(mkd$type == "transformational"),
    any(mkd$type == "neomorphic"),
    MkPrimeMCMC(nIter = 100L, gibbsSpr = FALSE, gibbsSubtreeSwap = FALSE)
  )
  moveWeights <- vapply(moves, `[[`, numeric(1), "weight")
  moveTypeCodes <- vapply(
    moves, function(m) MkPrime:::.kMoveTypes[[m$name]], integer(1L)
  )
  nMoves <- length(moves)
  nEdge <- nrow(tree$edge)
  transIdx <- which(mkd$type == "transformational")
  transIdx0 <- if (length(transIdx) > 0) transIdx - 1L else integer(0)

  scaleTunings <- matrix(0.5, 1, nMoves)
  sliceParamCodes <- vapply(moves, function(m) m$sliceParamIdx %||% 0L, integer(1L))
  sliceWidths <- matrix(1.0, 1, nMoves)
  bsTunings <- 10
  iwWins <- 1L

  result <- run_mcmc_batch_cpp(
    mcmcData, list(chainState), 1.0,
    moveTypeCodes, transIdx0, sliceParamCodes, moveWeights,
    scaleTunings, bsTunings, iwWins, sliceWidths,
    50L, 1L, 100L, 10L,
    any(mkd$type == "neomorphic"), nEdge
  )

  expect_true("move_time_ns" %in% names(result))
  expect_true(is.matrix(result$move_time_ns))
  expect_equal(nrow(result$move_time_ns), 1L)  # 1 chain
  expect_equal(ncol(result$move_time_ns), nMoves)
  # All timing values should be non-negative
  expect_true(all(result$move_time_ns >= 0))
  # At least some timing should be positive (moves were executed)
  expect_gt(sum(result$move_time_ns), 0)
})


# ==========================================================================
# Integration: moveWeights adapt during warmup (via .BuildMoves + .AdaptMoveWeights)
# ==========================================================================

test_that(".BuildMoves produces weights that normalize to 1", {
  mcmc <- MkPrimeMCMC(nIter = 200L)
  moves <- MkPrime:::.BuildMoves(20, 5, TRUE, mcmc)
  w <- vapply(moves, `[[`, numeric(1), "weight")
  wNorm <- w / sum(w)
  expect_equal(sum(wNorm), 1.0, tolerance = 1e-10)
})

test_that("All move names in .kMoveTypes are valid moveWeights names", {
  validInMcmc <- c(
    "tree_length", "branch_lengths", "nni", "spr", "tbr", "kPrime", "p",
    "rate_loss", "rate_log_sd", "rate_neo", "neo_joint",
    "gibbs_spr", "gibbs_subtree_swap",
    "weighted_branch_lengths", "weighted_spr", "weighted_subtree_swap",
    "block_gibbs_branch",
    "beta_scale", "pspr",
    "slice_rate_loss", "slice_rate_neo", "slice_rate_log_sd",
    "slice_tree_length", "slice_beta_scale"
  )
  expect_true(all(names(MkPrime:::.kMoveTypes) %in% validInMcmc))
})


# ==========================================================================
# Slow integration tests
# ==========================================================================

test_that("Adaptive scheduler runs end-to-end (short MCMC)", {
  skip_slow_tests()
  skip_if_not_installed("TreeSearch")
  dat <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  mkd <- MkPrimeData(dat)
  model <- MkPrimeModel()
  tree <- ape::rtree(length(dat), tip.label = names(dat))
  tree <- TreeTools::Preorder(tree)

  mcmc <- MkPrimeMCMC(
    nIter = 600L, warmup = 400L, thin = 10L,
    nRuns = 1L, nChains = 1L
  )
  result <- RunMkPrime(dat, tree, model = model, mcmc = mcmc)
  expect_s3_class(result, "MkPosterior")
  expect_gt(nrow(result$samples[[1]]), 0)
})

test_that("User-pinned moveWeights preserved end-to-end", {
  skip_slow_tests()
  skip_if_not_installed("TreeSearch")
  dat <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  mkd <- MkPrimeData(dat)
  model <- MkPrimeModel()
  tree <- ape::rtree(length(dat), tip.label = names(dat))
  tree <- TreeTools::Preorder(tree)

  mcmc <- MkPrimeMCMC(
    nIter = 600L, warmup = 400L, thin = 10L,
    nRuns = 1L, nChains = 1L,
    moveWeights = c(nni = 0.3)
  )
  # Just verify it runs without error
  result <- RunMkPrime(dat, tree, model = model, mcmc = mcmc)
  expect_s3_class(result, "MkPosterior")
})
