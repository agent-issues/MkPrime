# Tests for M-092: Adaptive move scheduler
#
# Covers: .AdaptMoveWeights(), .ResolvePinnedWeights(),
# .NormalizeMoveWeights(), .LogMoveWeights(), .FormatMoveWeights(),
# MkPrimeMCMC() moveWeights parameter, and C++ timing data.

# ==========================================================================
# MkPrimeMCMC() moveWeights parameter validation
# ==========================================================================

test_that("MkPrimeMCMC() stores moveWeights = NULL by default", {
  mcmc <- MkPrimeMCMC(nIter = 100L, minWarmup = 50L)
  expect_null(mcmc$moveWeights)
})

test_that("MkPrimeMCMC() stores valid moveWeights", {
  mcmc <- MkPrimeMCMC(nIter = 100L, minWarmup = 50L, moveWeights = c(nni = 0.3, spr = 0.2))
  expect_equal(mcmc$moveWeights, c(nni = 0.3, spr = 0.2))
})

test_that("MkPrimeMCMC() accepts moveWeights as list", {
  mcmc <- MkPrimeMCMC(nIter = 100L, minWarmup = 50L, moveWeights = list(nni = 0.3))
  expect_equal(mcmc$moveWeights, c(nni = 0.3))
})

test_that("MkPrimeMCMC() rejects unnamed moveWeights", {
  expect_error(MkPrimeMCMC(nIter = 100L, minWarmup = 50L, moveWeights = c(0.3, 0.2)),
               "named numeric")
})

test_that("MkPrimeMCMC() rejects unknown move names", {
  expect_error(MkPrimeMCMC(nIter = 100L, minWarmup = 50L, moveWeights = c(bogus = 0.3)),
               "unknown move name")
})

# M-167: superseded by test-proposal-tuning.R, which derives the buildable
# move names from .BuildMoves()/.BuildMovesPartitioned() instead of restating
# them. A hand-copied list here had gone stale in exactly the way it was
# meant to prevent.

test_that("MkPrimeMCMC() rejects negative moveWeights", {
  expect_error(MkPrimeMCMC(nIter = 100L, minWarmup = 50L, moveWeights = c(nni = -0.1)),
               "positive")
})

test_that("MkPrimeMCMC() rejects moveWeights summing > 1", {
  expect_error(MkPrimeMCMC(nIter = 100L, minWarmup = 50L,
                            moveWeights = c(nni = 0.6, spr = 0.6)),
               "exceeds 1")
})

test_that("MkPrimeMCMC() accepts moveWeights summing to exactly 1", {
  mcmc <- MkPrimeMCMC(nIter = 100L, minWarmup = 50L,
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
  result <- suppressWarnings(MkPrime:::.ResolvePinnedWeights(
    c(nni = 0.3, spr = 0.2, bogus = 0.1),
    c("nni", "spr", "tree_length")
  ))
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
# SBC-WARMUP-002 regression: scalar-target moves must not be annealed to
# the global floor by the softmax `log(dim)` advantage of multi-parameter
# competitors.
# ==========================================================================

test_that(".AdaptMoveWeights honours wMinScalar for scalar-target moves", {
  # Reproduce the SBC-WARMUP-002 trace: a dim=1 scale move on a singleton
  # scalar parameter (`tree_length`) competes against high-dim moves
  # (`kPrime` with dim ~ nTrans) whose log(dim) term in the softmax score
  # crushes the singleton to `wMin` even when its acceptance rate is
  # healthy. Without the scalar floor the scheduler freezes that parameter.
  moveNames <- c("tree_length", "branch_lengths", "kPrime", "nni")
  current <- setNames(c(0.10, 0.30, 0.30, 0.30), moveNames)
  accept  <- setNames(c(35L, 86L, 94L, 23L) * 10L, moveNames)
  propose <- setNames(rep(1000L, 4L), moveNames)
  timeNs  <- setNames(c(1.4e10, 4.5e10, 1.4e10, 5e10), moveNames)
  moveDim <- setNames(c(1L, 1L, 30L, 1L), moveNames)

  # Without scalar floor: tree_length drops to global wMin = 0.01.
  res0 <- MkPrime:::.AdaptMoveWeights(
    current, accept, propose, timeNs, moveNames,
    moveDim = moveDim, pinnedWeights = NULL,
    warmupProgress = 1.0
  )
  expect_equal(res0[["tree_length"]], 0.01, tolerance = 1e-10)

  # With scalar floor: tree_length pinned at wMinScalar = 0.02.
  res1 <- MkPrime:::.AdaptMoveWeights(
    current, accept, propose, timeNs, moveNames,
    moveDim = moveDim, pinnedWeights = NULL,
    scalarFloorMoves = c("tree_length"),
    warmupProgress = 1.0
  )
  expect_gte(res1[["tree_length"]], 0.02 - 1e-10)
  expect_equal(sum(res1), 1.0, tolerance = 1e-10)
})

test_that("SBC-WARMUP-002-BRANCH: beta_simplex move stays ≥wMinScalar against high-dim competitor", {
  # Reproduce the branch_lengths symptom: dim=1 beta_simplex with 86%
  # acceptance is crushed to wMin=0.01 by kPrime (dim=30). After adding
  # "beta_simplex" to .kScalarFloorTypes, the scalarFloorMoves vector
  # must include "branch_lengths" and its floor must be wMinScalar=0.02.
  moveNames <- c("branch_lengths", "kPrime", "nni")
  current  <- setNames(c(0.30, 0.40, 0.30), moveNames)
  accept   <- setNames(c(860L, 940L, 230L), moveNames)
  propose  <- setNames(rep(1000L, 3L), moveNames)
  timeNs   <- setNames(c(4.5e10, 1.4e10, 5.0e10), moveNames)
  moveDim  <- setNames(c(1L, 30L, 1L), moveNames)

  # Without floor: branch_lengths drops to global wMin = 0.01.
  res0 <- MkPrime:::.AdaptMoveWeights(
    current, accept, propose, timeNs, moveNames,
    moveDim = moveDim, pinnedWeights = NULL,
    warmupProgress = 1.0
  )
  expect_equal(res0[["branch_lengths"]], 0.01, tolerance = 1e-10)

  # With floor: branch_lengths stays at ≥wMinScalar = 0.02.
  res1 <- MkPrime:::.AdaptMoveWeights(
    current, accept, propose, timeNs, moveNames,
    moveDim = moveDim, pinnedWeights = NULL,
    scalarFloorMoves = c("branch_lengths"),
    warmupProgress = 1.0
  )
  expect_gte(res1[["branch_lengths"]], 0.02 - 1e-10)
  expect_equal(sum(res1), 1.0, tolerance = 1e-10)
})

test_that("SBC-WARMUP-002: scalar moves keep ≥wMinScalar through warmup", {
  # End-to-end check that the scalar floor is wired through RunMkPrime
  # and stored in the final adapted weights. Uses a small simulated
  # dataset matching the SBC-WARMUP-002 reproducer dimensions. Kept fast
  # by limiting nIter; the key invariant is the weight floor, not the
  # trace.
  skip_on_cran()
  skip_if_not_installed("TreeTools")
  skip_if_not_installed("ape")

  set.seed(20262527L)
  nTip <- 8L
  tr <- ape::rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  tr$edge.length <- rep_len(0.1, nrow(tr$edge))
  tr <- Preorder(tr)
  mat <- matrix(sample.int(2L, nTip * 20L, replace = TRUE) - 1L,
                nTip, 20L, dimnames = list(tr$tip.label, NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(coding = "variable", kPrimePrior = "geometric",
                         expSteps = 50)
  mcmc <- MkPrimeMCMC(nIter = 3500L, thin = 50L,
                       minWarmup = 3000L, maxWarmup = 3000L,
                       autoTune = FALSE, nRuns = 1L, nChains = 1L)

  res <- allow_warning(
    RunMkPrime(mkd, tr, model = model, mcmc = mcmc,
               fixTopology = TRUE, overwrite = TRUE),
    "without stabilisation"
  )

  expect_true("tree_length" %in% names(res$moveWeights))
  expect_gte(res$moveWeights[["tree_length"]], 0.02 - 1e-6)
  expect_gte(res$moveWeights[["rate_log_sd"]], 0.02 - 1e-6)
  expect_gte(res$moveWeights[["branch_lengths"]], 0.02 - 1e-6)
})


# ==========================================================================
# .FormatMoveWeights()
# ==========================================================================

test_that(".FormatMoveWeights produces readable string", {
  s <- MkPrime:::.FormatMoveWeights(c(0.3, 0.7), c("nni", "spr"))
  expect_match(s, "nni:30.0%")
  expect_match(s, "spr:70.0%")
})


# ==========================================================================
# .LogMoveWeights() writes to file
# ==========================================================================

test_that(".LogMoveWeights writes categorized comment lines to log file", {
  tmp <- tempfile(fileext = ".log")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("Sample\ta\tb", tmp)
  MkPrime:::.LogMoveWeights(c(0.6, 0.4), c("nni", "spr"), tmp)
  lines <- readLines(tmp)
  # Header line + one category line (Topology)
  expect_length(lines, 2L)
  expect_match(lines[2], "^# Topology:")
  expect_match(lines[2], "nni:60.0%")
  expect_match(lines[2], "spr:40.0%")
})


# ==========================================================================
# C++ timing data from run_mcmc_batch_cpp
# ==========================================================================

test_that("run_mcmc_batch_cpp returns move_time_ns matrix", {
  skip_if_not_installed("TreeSearch")
  dat <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  mkd <- suppressWarnings(MkPrimeData(dat))
  model <- MkPrimeModel()
  tree <- ape::rtree(length(dat), tip.label = names(dat))
  tree <- Preorder(tree)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)

  mcmcData <- MkPrime:::.InitMcmcData(mkd, model)
  state <- MkPrime:::.InitState(tree, mkd, model)
  chainState <- MkPrime:::.InitMcmcChain(state)
  fill_partition_cache(mcmcData, chainState)
  allocate_cl_workspace(mcmcData, chainState)

  moves <- MkPrime:::.BuildMoves(
    nrow(tree$edge), sum(mkd$type == "transformational"),
    any(mkd$type == "neomorphic"),
    MkPrimeMCMC(nIter = 100L, minWarmup = 50L, gibbsSpr = FALSE, gibbsSubtreeSwap = FALSE)
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
  jointRhos <- matrix(0.0, 1, nMoves)
  bsTunings <- 10
  iwWins <- 1L
  moveIntPars <- integer(nMoves)

  result <- run_mcmc_batch_cpp(
    mcmcData, list(chainState), 1.0,
    moveTypeCodes, transIdx0, sliceParamCodes, moveWeights,
    scaleTunings, bsTunings, iwWins, moveIntPars, sliceWidths, jointRhos,
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

# Types that .BuildMoves' init-time scalar floor covers, and the schedule
# fields the tests below read back.
.moveWeights <- function(moves) {
  w <- vapply(moves, `[[`, numeric(1), "weight")
  names(w) <- vapply(moves, `[[`, character(1), "name")
  w
}
.scalarFloored <- function(moves) {
  scalarTypes <- c("scale", "int_walk", "gibbs_p", "scale_p", "logit_scale_p",
                   "slice", "beta_simplex", "gibbs_p_marginal")
  type <- vapply(moves, function(m) if (is.null(m$type)) m$name else m$type,
                 character(1))
  dim <- vapply(moves, `[[`, numeric(1), "dim")
  (dim == 1 & type %in% scalarTypes) | type == "joint_2d"
}

test_that(".BuildMoves floors scalar moves against a growing schedule", {
  # Normalising a vector and checking it sums to 1 constrains nothing.  What
  # the schedule must guarantee is that no scalar move is starved as the
  # topology and branch moves grow with nEdge.
  mcmc <- MkPrimeMCMC(nIter = 200L, minWarmup = 100L)
  small <- MkPrime:::.BuildMoves(20L, 5L, TRUE, mcmc)
  big   <- MkPrime:::.BuildMoves(200L, 5L, TRUE, mcmc)

  for (moves in list(small, big)) {
    w <- .moveWeights(moves)
    expect_true(all(is.finite(w)))
    expect_true(all(w > 0))
    # Without the floor the smallest scalar weight stays at its raw value
    # while the total grows with nEdge, so this share collapses.
    expect_gt(min(w[.scalarFloored(moves)]) / sum(w), 0.01)
  }

  # The floor is a share of the schedule, so it tracks the schedule's size.
  expect_gt(min(.moveWeights(big)[.scalarFloored(big)]),
            3 * min(.moveWeights(small)[.scalarFloored(small)]))
})

test_that(".BuildMoves builds a usable schedule for every k' arm", {
  # No test drove the empirical_geometric arm at all, so nothing asserted
  # that it schedules mh_logit_p -- the move that samples p there.
  mcmc <- MkPrimeMCMC(nIter = 200L, minWarmup = 100L)
  arms <- expand.grid(
    kPrimePrior = c("geometric", "empirical_geometric", "beta_geometric",
                    "logseries"),
    likelihoodMode = c("sampled_k", "marginal_k"),
    qHeterogeneity = c(FALSE, TRUE),
    fixTopology = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(arms))) {
    arm <- arms[i, ]
    label <- paste(unlist(arm), collapse = "/")
    moves <- suppressMessages(MkPrime:::.BuildMoves(
      20L, 5L, TRUE, mcmc,
      fixTopology = arm$fixTopology,
      kPrimePrior = arm$kPrimePrior,
      qHeterogeneity = arm$qHeterogeneity,
      likelihoodMode = arm$likelihoodMode
    ))
    w <- .moveWeights(moves)
    expect_true(all(is.finite(w) & w > 0), label = label)
    expect_false(anyDuplicated(names(w)) > 0L, label = label)
    expect_true(all(c("tree_length", "branch_lengths") %in% names(w)),
                label = label)
    expect_identical(any(c("nni", "spr") %in% names(w)), !arm$fixTopology,
                     label = label)
  }

  # p is sampled by mh_logit_p on both geometric arms; the beta_geometric arm
  # marginalises p away and must not schedule it.
  pMove <- function(prior) {
    "mh_logit_p" %in% names(.moveWeights(suppressMessages(
      MkPrime:::.BuildMoves(20L, 5L, TRUE, mcmc, kPrimePrior = prior))))
  }
  expect_true(pMove("empirical_geometric"))
  expect_true(pMove("geometric"))
  expect_false(pMove("beta_geometric"))
})

test_that("composed warmup adaptation keeps a positive, floored schedule", {
  # .AdaptMoveWeights, .DecayLowAcceptMoves and the gibbs cap were each
  # tested alone; the schedule the chain actually samples from is their
  # composition, applied batch after batch.
  mcmc <- MkPrimeMCMC(nIter = 200L, minWarmup = 100L)
  moves <- MkPrime:::.BuildMoves(20L, 5L, TRUE, mcmc)
  moveNames <- vapply(moves, `[[`, character(1), "name")
  moveDim <- as.integer(vapply(moves, `[[`, numeric(1), "dim"))
  floorMoves <- moveNames[.scalarFloored(moves)]
  gibbsKpIdx <- match("gibbs_kPrime", moveNames)

  initial <- .moveWeights(moves) / sum(.moveWeights(moves))
  # RunMkPrime auto-pins the always-accepting Gibbs and slice moves.
  pinned <- initial[c("gibbs_kPrime", "slice_rate_log_sd")]
  budget <- 1 - sum(pinned)

  starved <- "spr"
  propose <- initial
  propose[] <- 100
  accept <- propose * 0.4
  accept[[starved]] <- 0
  moveTimeNs <- propose
  moveTimeNs[] <- 1e6

  w <- initial
  for (batch in seq_len(6L)) {
    w <- MkPrime:::.AdaptMoveWeights(
      w, accept, propose, moveTimeNs, moveNames, moveDim,
      pinnedWeights = pinned, warmupProgress = batch / 6,
      scalarFloorMoves = floorMoves
    )
    w <- MkPrime:::.DecayLowAcceptMoves(
      w, accept, propose, initial, moveNames, pinnedWeights = pinned
    )
    expect_equal(sum(w), 1, tolerance = 1e-12)
    expect_true(all(w > 0))
    expect_equal(w[names(pinned)], pinned, tolerance = 1e-12)
    # wMinScalar (0.02) applies within the unpinned budget.
    expect_gte(min(w[floorMoves]), 0.02 * budget - 1e-9)

    capped <- MkPrime:::.WarmupGibbsCap(w, pinned, gibbsKpIdx)
    expect_equal(sum(capped), 1, tolerance = 1e-12)
    expect_true(all(capped > 0))
    expect_equal(capped[["gibbs_kPrime"]], pinned[["gibbs_kPrime"]] / 3,
                 tolerance = 1e-12)
    expect_equal(MkPrime:::.RestoreGibbsCap(capped, pinned, gibbsKpIdx), w,
                 tolerance = 1e-12)
  }

  # A move that never accepts must end below where it started, and must not
  # have been promoted along the way.
  expect_lt(w[[starved]], initial[[starved]])
})

# The real guard against move-name drift lives in test-proposal-tuning.R: it
# drives MkPrimeMCMC() itself rather than a hand-copied list of names.

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
  tree <- Preorder(tree)

  mcmc <- MkPrimeMCMC(
    nIter = 600L, maxWarmup = 400L, minWarmup = 400L, thin = 10L,
    autoTune = FALSE, nRuns = 1L, nChains = 1L
  )
  result <- RunMkPrime(dat, tree, model = model, mcmc = mcmc)
  expect_s3_class(result, "MkPosterior")
  expect_gt(nrow(result$samples), 0)
})

test_that("User-pinned moveWeights preserved end-to-end", {
  skip_slow_tests()
  skip_if_not_installed("TreeSearch")
  dat <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  mkd <- MkPrimeData(dat)
  model <- MkPrimeModel()
  tree <- ape::rtree(length(dat), tip.label = names(dat))
  tree <- Preorder(tree)

  mcmc <- MkPrimeMCMC(
    nIter = 600L, maxWarmup = 400L, minWarmup = 400L, thin = 10L,
    autoTune = FALSE, nRuns = 1L, nChains = 1L,
    moveWeights = c(nni = 0.3)
  )
  # Just verify it runs without error
  result <- RunMkPrime(dat, tree, model = model, mcmc = mcmc)
  expect_s3_class(result, "MkPosterior")
})


# --- #68: the move schedule a user reads back ---

test_that("the gibbs cap leaves a probability vector (#68)", {
  # With no free move to absorb the freed weight the vector is renormalized,
  # which is the schedule mcmc.cpp samples from either way -- so the cap must
  # still bind, and the total must still be 1.
  allPinned <- c(gibbs_kPrime = 0.6, nni = 0.4)          # no free move at all
  capped <- .WarmupGibbsCap(allPinned, allPinned, 1L)
  expect_equal(sum(capped), 1)
  expect_equal(capped[["gibbs_kPrime"]], 0.2 / 0.6)      # capped share, not 0.6

  noFreeWeight <- c(gibbs_kPrime = 1, nni = 0)           # free moves hold none
  onlyGibbs <- .WarmupGibbsCap(noFreeWeight, c(gibbs_kPrime = 1), 1L)
  expect_equal(sum(onlyGibbs), 1)

  tooLittleFree <- c(gibbs_kPrime = 0.2, nni = 1e-15)    # too little to reclaim
  restored <- .RestoreGibbsCap(tooLittleFree, c(gibbs_kPrime = 0.8), 1L)
  expect_equal(sum(restored), 1)
  expect_equal(restored[["gibbs_kPrime"]], 1, tolerance = 1e-12)
})


test_that("a feasible gibbs cap still caps, and restores (#68)", {
  weights <- c(gibbs_kPrime = 0.6, nni = 0.3, spr = 0.1)
  pinned  <- c(gibbs_kPrime = 0.6)

  capped <- .WarmupGibbsCap(weights, pinned, 1L, factor = 1 / 3)
  expect_equal(capped[["gibbs_kPrime"]], 0.2)
  expect_equal(sum(capped), 1)
  # Freed weight is shared in proportion to the free moves' existing weights.
  expect_equal(capped[["nni"]] / capped[["spr"]], 3)

  restored <- .RestoreGibbsCap(capped, pinned, 1L)
  expect_equal(restored[["gibbs_kPrime"]], 0.6)
  expect_equal(sum(restored), 1)
  expect_equal(unname(restored), unname(weights))
})


test_that("a fully pinned schedule survives the gibbs cap round trip", {
  pins <- c(nni = 0.3, gibbs_kPrime = 0.07895, spr = 0.42105, slice_p = 0.2)
  capped <- .WarmupGibbsCap(pins, pins, 2L)
  expect_lt(capped[["gibbs_kPrime"]], pins[["gibbs_kPrime"]])
  expect_equal(.RestoreGibbsCap(capped, pins, 2L), pins)

  # Pins exhaust the budget while a free move holds nothing.
  withFree <- c(pins, tbr = 0)
  capped <- .WarmupGibbsCap(withFree, pins, 2L)
  expect_equal(.RestoreGibbsCap(capped, pins, 2L), withFree)

  # Raw pins that do not sum to 1 restore to the normalized schedule.
  raw <- pins * 5
  capped <- replace(pins / (1 - 0.05), 2L, 0.02895 / 0.95)
  expect_equal(.RestoreGibbsCap(capped, raw, 2L), pins)
})


test_that(".BuildResult surfaces every run's frozen schedule (#68)", {
  # Runs adapt independently, so `$moveWeights` -- run 1's -- presents one
  # run's decisions as though they were the analysis's.
  paramNames <- c("log_posterior", "log_likelihood", "tree_length")
  Run <- function(nniWeight) {
    list(
      samples = matrix(0, 4L, length(paramNames),
                       dimnames = list(NULL, paramNames)),
      tree_samples = vector("list", 4L),
      saved_idx = 4L, tree_saved_idx = 4L,
      chain_accept = list(c(nni = 1)), chain_propose = list(c(nni = 2)),
      chain_tuning = list(list()), logPostHistory = numeric(0),
      moveWeights = c(nni = nniWeight, spr = 1 - nniWeight)
    )
  }
  mcmc <- list(nChains = 1L, nRuns = 2L, warmup = 0L, thin = 1L, treeThin = 1L)

  result <- .BuildResult(
    runs = list(Run(0.7), Run(0.4)), model = NULL, mkd = NULL, mcmc = mcmc,
    paramNames = paramNames, logFilePaths = NULL, actualIter = 4L,
    stopReason = "nIter"
  )

  expect_equal(result$moveWeights, c(nni = 0.7, spr = 0.3))
  expect_length(result$runMoveWeights, 2L)
  expect_equal(result$runMoveWeights[[1]], c(nni = 0.7, spr = 0.3))
  expect_equal(result$runMoveWeights[[2]], c(nni = 0.4, spr = 0.6))
})
