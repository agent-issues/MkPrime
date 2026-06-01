# Tests for Phase 9: Gibbs moves
#
# M-082: Gibbs update for p
# M-084: Gibbs update for k'_i (added here when implemented)

library("ape")
library("TreeTools")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.small_trans_tree <- function() {
  read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
}

.small_trans_pd <- function() {
  mat <- matrix(c(0, 1, 2, 0,
                  0, 1, 0, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  MatrixToPhyDat(mat)
}


# ---------------------------------------------------------------------------
# M-082: Gibbs update for p
# ---------------------------------------------------------------------------

test_that("geometric arm schedules mh_logit_p for p (truncated geometric is non-conjugate)", {
  # Stage 2 (MARGINAL-K-TRUNC-001): the geometric prior is truncated at K, whose
  # p-normaliser Z(p) breaks Beta conjugacy, so p is sampled via mh_logit_p
  # (case 30), NOT the legacy gibbs_p (case 9) -- mirroring marginal_k and
  # empirical_geometric. .BuildMoves defaults to the geometric, sampled_k arm.
  pd <- .small_trans_pd()
  mkd <- MkPrimeData(pd)
  nEdge <- 2 * length(.small_trans_tree()$tip.label) - 3L
  nTrans <- sum(mkd$type == "transformational")
  moves <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE,
                                 mcmc = MkPrimeMCMC(), fixTopology = FALSE)
  types   <- vapply(moves, `[[`, character(1L), "type")
  targets <- vapply(moves, function(m) if (is.null(m$target)) "" else m$target,
                    character(1L))
  p_moves <- moves[targets == "p"]
  expect_length(p_moves, 1L)
  expect_equal(p_moves[[1L]]$type, "logit_scale_p")
  expect_equal(p_moves[[1L]]$name, "mh_logit_p")
  expect_false("gibbs_p" %in% types)
})

test_that("gibbs_p move is absent when there are no transformational chars", {
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, knownStates = c("1" = 2L))
  nEdge <- 2 * 4L - 3L
  moves <- MkPrime:::.BuildMoves(nEdge, nTrans = 0L, hasNeo = FALSE,
                                 mcmc = MkPrimeMCMC(), fixTopology = FALSE)
  p_names <- vapply(moves, `[[`, character(1L), "name")
  expect_false("p" %in% p_names)
})

# The legacy conjugate-Beta-draw tests for gibbs_p (empirical Beta moments and
# acceptance == 1) validated a full conditional that holds ONLY for the
# UNtruncated geometric. Stage 2 (MARGINAL-K-TRUNC-001) truncates the geometric
# prior -- Z(p) is p-dependent -- so that conjugacy no longer holds for any
# active model, and the move is no longer scheduled. Those mechanism tests have
# been retired. What matters now is (a) the C++ engine REFUSES a misrouted
# gibbs_p so it cannot silently emit wrong samples, and (b) RunMkPrime samples p
# via mh_logit_p. Both are pinned below.

test_that("C++ engine rejects gibbs_p under the truncated geometric (non-conjugate guard)", {
  skip_if_not(requireNamespace("ape", quietly = TRUE))
  set.seed(2194)

  pd   <- .small_trans_pd()
  mkd  <- MkPrimeData(pd)
  model <- MkPrimeModel(kPrimePrior = "geometric")   # truncated at K (Stage 2)
  tree <- .small_trans_tree()

  model  <- MkPrime:::.FinalizeModel(model, tree, mkd)
  state0 <- MkPrime:::.InitState(tree, mkd, model)
  mcmcData <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(mcmcData, statePtr)
  allocate_cl_workspace(mcmcData, statePtr)

  s_before <- get_mcmc_state(statePtr)
  # Route a gibbs_p (case 9) call directly to the engine: it MUST be refused,
  # because the truncated geometric's p full conditional is not Beta.
  accepted <- do_move_cpp(mcmcData, statePtr,
    moveType = 9L, charIdx = 0L,
    scaleTuning = 0.5, betaSimplexTuning = 10.0,
    intWalkWindow = 1L, beta = 1.0)
  s_after <- get_mcmc_state(statePtr)

  expect_false(accepted)                            # non-conjugate -> rejected
  expect_equal(s_after$p, s_before$p)                # p untouched
  expect_equal(s_after$logPrior, s_before$logPrior)  # prior untouched
})

test_that("RunMkPrime samples p via mh_logit_p under the geometric prior (valid posterior)", {
  set.seed(7342)
  tree <- .small_trans_tree()
  pd   <- .small_trans_pd()

  result <- suppressWarnings(RunMkPrime(pd, tree,
    model = MkPrimeModel(kPrimePrior = "geometric"),
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 600L, thin = 5L,
                       maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE)))

  expect_s3_class(result, "MkPosterior")
  expect_true(all(is.finite(result$samples[, "log_posterior"])))
  p_vals <- result$samples[, "p"]
  expect_true(all(p_vals > 0 & p_vals < 1))

  # p is now driven by mh_logit_p (case 30), not gibbs_p. The acceptance table is
  # keyed by move name: mh_logit_p must be present and gibbs_p ("p") absent.
  expect_true("mh_logit_p" %in% names(result$acceptance))
  expect_false("p" %in% names(result$acceptance))
  expect_gt(stats::sd(p_vals), 0)                    # p actually explores
})
