# T-017: OMP determinism + bit-identity tests for the parallel ecology
# partition pruner.
#
# These tests verify that cpp_log_likelihood_ecology returns BIT-identical
# results regardless of OMP_NUM_THREADS, within a single build.  The
# parallel implementation uses per-thread partial-sum slots followed by a
# fixed-order serial merge — the merge order is the work-item list build
# order, which is itself deterministic — so the floating-point sum is
# identical whether 1 or N threads ran the inner work.
#
# What we are NOT testing here:
#   - Bit-identity between a serial-only build (no -fopenmp) and an
#     OMP build.  Compile-with-OMP can change auto-vectorisation choices
#     at ULP level; that's documented in the design and out of scope.
#   - Bit-identity of MCMC sample trajectories — Phase 1 of T-017 does not
#     touch any RNG; existing per-move tests cover that.
#
# These are exercise tests for the orchestrator only.  Move-handler call
# sites of cpp_partition_log_likelihood_ecology are covered by the regular
# test-ecology-* tests, which we also run under both OMP=1 and OMP=4
# during CI to confirm safety end-to-end.

library("TreeTools")

.MakeOmpEcologyFixture <- function(seed = 20260521) {
  set.seed(seed)
  tips <- paste0("t", 1:6)
  mat <- matrix(c(
    0, 1, 2, 0, 1, 2,   # transformational, k = 3
    1, 0, 1, 2, 2, 0,   # transformational, k = 3
    2, 1, 0, 1, 0, 2,   # transformational, k = 3
    0, 1, 0, 1, 0, 1,   # neomorphic (binary)
    1, 0, 1, 0, 1, 0,   # neomorphic (binary)
    0, 0, 1, 1, 2, 2    # ecology (3 states)
  ), nrow = 6, ncol = 6, byrow = FALSE,
     dimnames = list(tips, NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = c(4L, 5L), ecology = 6L)
  tree <- Preorder(ape::rtree(6, tip.label = tips))
  list(mkd = mkd, tree = tree)
}

.MakeOmpEcoDataPtr <- function(mkd, model) {
  parts <- lapply(mkd$partitions, function(p) {
    list(type = p$type, k = p$k, kObs = p$kObs,
         char_indices = p$char_indices,
         tip_states = p$tip_states,
         unique_tip_states = p$unique_tip_states,
         pattern_index = p$pattern_index)
  })
  ecoTip <- as.integer(mkd$ecology)
  ecoTip[is.na(ecoTip)] <- -1L
  prepare_mcmc_data(
    parts, as.integer(mkd$kObs), mkd$type,
    any(mkd$type == "neomorphic"),
    model$nCat, model$coding, model$relabel,
    model$treeLengthShape, model$treeLengthRate,
    model$rateLossMeanlog, model$rateLossSdlog,
    model$rateLogSdShape, model$rateLogSdRate,
    model$rateNeoMeanlog, model$rateNeoSdlog,
    model$kprimeHyperA, model$kprimeHyperB,
    identical(model$kPrimePrior, "logseries"),
    model$kprimeLogseriesC %||% 0.7,
    FALSE, FALSE, 4L, 1.0, 1.0,
    FALSE, numeric(0), 1L, 0L, 0, -1e308,
    TRUE,
    ecoTip, as.integer(mkd$kEcology),
    model$magnitudeMode %||% "global",
    model$rho0Alpha %||% 75.0, model$rho0Beta %||% 25.0,
    model$sigmaPhi %||% 1.0, as.integer(model$gibbsZEvery %||% 50L),
    model$thetaAlpha %||% 1.0, model$thetaBeta %||% 1.0
  )
}

# Helper: evaluate cpp_log_likelihood_ecology for a fixture under a fixed
# argument set.  Returns the scalar logLik.  Caller is responsible for the
# OMP_NUM_THREADS env-var setting (already in effect when this is invoked).
.OmpEvalEco <- function(seedZ, mkd_tree, model_opts = list()) {
  f <- mkd_tree
  defaults <- list(
    ecologyAware = TRUE,
    expSteps     = 10,
    kPrimePrior  = "geometric",
    coding       = "none",
    nCat         = 4L          # exercise ACRV path
  )
  model_args <- modifyList(defaults, model_opts)
  model <- do.call(MkPrimeModel, model_args)
  dataPtr <- .MakeOmpEcoDataPtr(f$mkd, model)
  set.seed(seedZ)
  zMat <- matrix(as.integer(sample(0:2, f$mkd$nChar * (f$mkd$kEcology - 1L),
                                   replace = TRUE)),
                 nrow = f$mkd$nChar, ncol = f$mkd$kEcology - 1L)
  parent <- f$tree$edge[, 1]
  child  <- f$tree$edge[, 2]
  edgeLen <- f$tree$edge.length
  MkPrime:::.CppLogLikelihoodEcology(
    dataPtr, parent, child, edgeLen,
    kPrime    = as.integer(f$mkd$kObs),
    rateLoss  = 1.2,
    rateLogSd = 0.5,     # exercises useAcrv = TRUE path
    rateNeo   = 0.9,
    phi       = 2.0,
    zMatrix   = zMat,
    pi0       = 0.5
  )
}


# T-017 test 1: same seed, OMP_NUM_THREADS=1 vs =4 yields BIT-identical
# logLik (no tolerance — we want EXACT equality at IEEE-754 representation
# level).  The serial merge order is fixed regardless of thread count,
# so the post-merge sum should reproduce the serial reference exactly.
test_that("OMP_NUM_THREADS={1,4} produce bit-identical logLik (T-017)", {
  fix <- .MakeOmpEcologyFixture(seed = 1L)

  prev <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
  on.exit(if (is.na(prev)) Sys.unsetenv("OMP_NUM_THREADS")
          else Sys.setenv(OMP_NUM_THREADS = prev),
          add = TRUE)

  Sys.setenv(OMP_NUM_THREADS = "1")
  ll1 <- .OmpEvalEco(seedZ = 100L, mkd_tree = fix)
  Sys.setenv(OMP_NUM_THREADS = "4")
  ll4 <- .OmpEvalEco(seedZ = 100L, mkd_tree = fix)

  expect_identical(ll4, ll1)
})


# T-017 test 2: same seed, two independent OMP_NUM_THREADS=4 repeats
# yield BIT-identical logLik.  This catches any leaked per-thread mutable
# state that would make repeated parallel evaluations non-deterministic.
test_that("OMP_NUM_THREADS=4 is repeatable across calls (T-017)", {
  fix <- .MakeOmpEcologyFixture(seed = 2L)

  prev <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
  on.exit(if (is.na(prev)) Sys.unsetenv("OMP_NUM_THREADS")
          else Sys.setenv(OMP_NUM_THREADS = prev),
          add = TRUE)

  Sys.setenv(OMP_NUM_THREADS = "4")
  ll_a <- .OmpEvalEco(seedZ = 200L, mkd_tree = fix)
  ll_b <- .OmpEvalEco(seedZ = 200L, mkd_tree = fix)

  expect_identical(ll_a, ll_b)
})


# T-017 test 3: different z seeds produce different logLik (a sanity check
# that the fixture actually exercises non-trivial state — if the test above
# silently returned 0 from a degenerate fixture, this would also be 0 vs 0
# and would not catch the problem).
test_that("Different z seeds produce different logLik (sanity)", {
  fix <- .MakeOmpEcologyFixture(seed = 3L)
  ll_a <- .OmpEvalEco(seedZ = 1L,    mkd_tree = fix)
  ll_b <- .OmpEvalEco(seedZ = 1000L, mkd_tree = fix)
  expect_true(is.finite(ll_a))
  expect_true(is.finite(ll_b))
  expect_false(isTRUE(all.equal(ll_a, ll_b, tolerance = 0)))
})
