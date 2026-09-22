# tests/testthat/test-gibbs-het.R
# M-114: Validate partial CL for Gibbs moves under Q-heterogeneity

# ---------------------------------------------------------------------------
# Helper: create a small dataset + model with qHeterogeneity
# ---------------------------------------------------------------------------
make_het_setup <- function(nTip = 8L, nChar = 10L, kMax = 3L,
                            seed = 5719L) {
  set.seed(seed)
  tree <- ape::rtree(nTip, rooted = FALSE)
  tree <- Preorder(tree)
  mat <- matrix(sample(0:(kMax - 1L), nTip * nChar, replace = TRUE),
                nrow = nTip,
                dimnames = list(tree$tip.label, paste0("c", seq_len(nChar))))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(qHeterogeneity = TRUE)
  list(tree = tree, mkd = mkd, model = model)
}

# C++ data + state pointers for the Q-het model, ready for do_move_cpp().
make_het_pts <- function(...) {
  setup <- make_het_setup(...)
  model <- MkPrime:::.FinalizeModel(setup$model, setup$tree, setup$mkd)
  state0   <- MkPrime:::.InitState(setup$tree, setup$mkd, model)
  dataPtr  <- MkPrime:::.InitMcmcData(setup$mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  list(dataPtr = dataPtr, statePtr = statePtr)
}

# ---------------------------------------------------------------------------
# Test: F81 transition unit test
# ---------------------------------------------------------------------------
test_that("f81_transition matches analytical formula", {
  # f81_transition (src/f81.h) evaluates (P(t) %*% cl) without ever forming
  # P(t), via P_ij(t) = pi_j + (delta_ij - pi_j) exp(-mu t).  Build P(t) here
  # and multiply it out, so the kernel is checked against the matrix it
  # claims to apply rather than against a rearrangement of its own formula.
  pi <- c(0.3, 0.5, 0.2)
  k <- 3L
  mu <- 1 / (1 - sum(pi^2))
  tVal <- 0.15

  Pmat <- outer(seq_len(k), seq_len(k),
                function(i, j) pi[j] + ((i == j) - pi[j]) * exp(-mu * tVal))
  expect_equal(rowSums(Pmat), rep(1, k), tolerance = 1e-14)
  expect_true(all(Pmat >= 0))

  # Column c of `cl` is one character's conditional likelihoods.
  cl <- matrix(c(0.4, 0.1, 0.5,
                 0.2, 0.7, 0.1), nrow = k)
  expect_equal(f81_transition_cpp(cl, pi, mu, tVal), Pmat %*% cl,
               tolerance = 1e-14)

  # t = 0 is the identity; t -> Inf leaves every state at pi %*% cl.
  expect_equal(f81_transition_cpp(cl, pi, mu, 0), cl, tolerance = 1e-14)
  expect_equal(f81_transition_cpp(cl, pi, mu, 1e6),
               matrix(rep(as.vector(pi %*% cl), each = k), nrow = k),
               tolerance = 1e-12)

  expect_error(f81_transition_cpp(cl, pi[1:2], mu, tVal), "one entry per row")
})

# ---------------------------------------------------------------------------
# Test: Q-het Gibbs partial CL matches full evaluation
# ---------------------------------------------------------------------------
test_that("Q-het Gibbs SPR partial CL matches full evaluation", {
  # The Gibbs moves score candidates from partial CLs propagated along the
  # regraft path, and write the winner's logLik into the state.  Re-running
  # the full pruning at the committed state must reproduce it exactly.
  set.seed(3847L)
  pts <- make_het_pts(nTip = 8L, nChar = 12L, seed = 3847L)
  expect_equal(get_state_log_lik(pts$statePtr),
               eval_full_loglik_cpp(pts$dataPtr, pts$statePtr),
               tolerance = 1e-10)

  nAcc <- 0L
  for (i in seq_len(100L)) {
    if (do_move_cpp(pts$dataPtr, pts$statePtr, 10L, 0L, 0.5, 0.5, 1L, 1.0))
      nAcc <- nAcc + 1L
    expect_equal(get_state_log_lik(pts$statePtr),
                 eval_full_loglik_cpp(pts$dataPtr, pts$statePtr),
                 tolerance = 1e-10)
  }
  expect_gt(nAcc, 0L)
})

test_that("Q-het Gibbs subtree swap partial CL matches full evaluation", {
  set.seed(6243L)
  pts <- make_het_pts(nTip = 8L, nChar = 12L, seed = 6243L)

  nAcc <- 0L
  for (i in seq_len(100L)) {
    if (do_move_cpp(pts$dataPtr, pts$statePtr, 11L, 0L, 0.5, 0.5, 1L, 1.0))
      nAcc <- nAcc + 1L
    expect_equal(get_state_log_lik(pts$statePtr),
                 eval_full_loglik_cpp(pts$dataPtr, pts$statePtr),
                 tolerance = 1e-10)
  }
  expect_gt(nAcc, 0L)
})

# ---------------------------------------------------------------------------
# Test: Q-het Gibbs SPR runs end to end
# ---------------------------------------------------------------------------
test_that("Q-het Gibbs SPR runs and produces reasonable results", {
  skip_under_memcheck()
  setup <- make_het_setup(nTip = 8L, nChar = 12L, seed = 3847L)
  mcmc <- MkPrimeMCMC(
    nIter = 500L, maxWarmup = 200L, minWarmup = 200L, thin = 5L,
    autoTune = FALSE, gibbsSpr = TRUE, gibbsSubtreeSwap = FALSE,
    nRuns = 1L
  )
  result <- allow_warning(
    RunMkPrime(data = setup$mkd, tree = setup$tree,
               model = setup$model, mcmc = mcmc),
    "without stabilisation"
  )
  expect_s3_class(result, "MkPosterior")
  expect_true(nrow(result$samples) > 0)
  expect_true(all(is.finite(result$samples[, "log_posterior"])))
})

# ---------------------------------------------------------------------------
# Test: Q-het Gibbs subtree swap runs and produces reasonable results
# ---------------------------------------------------------------------------
test_that("Q-het Gibbs subtree swap runs and produces reasonable results", {
  skip_under_memcheck()
  setup <- make_het_setup(nTip = 8L, nChar = 12L, seed = 6243L)
  mcmc <- MkPrimeMCMC(
    nIter = 500L, maxWarmup = 200L, minWarmup = 200L, thin = 5L,
    autoTune = FALSE, gibbsSpr = FALSE, gibbsSubtreeSwap = TRUE,
    nRuns = 1L
  )
  result <- allow_warning(
    RunMkPrime(data = setup$mkd, tree = setup$tree,
               model = setup$model, mcmc = mcmc),
    "without stabilisation"
  )
  expect_s3_class(result, "MkPosterior")
  expect_true(nrow(result$samples) > 0)
  expect_true(all(is.finite(result$samples[, "log_posterior"])))
})

# ---------------------------------------------------------------------------
# Test: Q-het Gibbs SPR + swap together
# ---------------------------------------------------------------------------
test_that("Q-het with both Gibbs moves produces valid MCMC", {
  skip_under_memcheck()
  setup <- make_het_setup(nTip = 8L, nChar = 12L, seed = 9102L)
  mcmc <- MkPrimeMCMC(
    nIter = 600L, maxWarmup = 200L, minWarmup = 200L, thin = 5L,
    autoTune = FALSE, gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE,
    nRuns = 1L
  )
  result <- allow_warning(
    RunMkPrime(data = setup$mkd, tree = setup$tree,
               model = setup$model, mcmc = mcmc),
    "without stabilisation"
  )
  expect_s3_class(result, "MkPosterior")

  # Both Gibbs moves must be scheduled, and both must have accepted at least
  # one proposal: `result$acceptance` reports accepted / proposed per move.
  rates <- result$acceptance
  expect_true(all(c("gibbs_spr", "gibbs_subtree_swap") %in% names(rates)))
  expect_true(all(rates >= 0 & rates <= 1))
  expect_gt(rates[["gibbs_spr"]], 0)
  expect_gt(rates[["gibbs_subtree_swap"]], 0)
})
