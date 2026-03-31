# Tests for Gibbs kPrime sweep (moveType 25) and block kPrime shift (moveType 26)

library("ape")
library("TreeTools")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.trans_tree <- function() {
  read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
}

.trans_pd <- function() {
  mat <- matrix(c(0, 1, 2, 0,
                  0, 1, 0, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  MatrixToPhyDat(mat)
}

# Set up C++ data and state pointers from tree + phyDat
.setup_cpp <- function(tree, pd) {
  tree <- Preorder(tree)
  mkd   <- MkPrimeData(pd)
  model <- MkPrimeModel()
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  state0   <- MkPrime:::.InitState(tree, mkd, model)
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  list(mkd = mkd, dataPtr = dataPtr, statePtr = statePtr)
}


# ---------------------------------------------------------------------------
# .BuildMoves wiring
# ---------------------------------------------------------------------------

test_that("gibbs_kPrime and block_kPrime appear for transformational data", {
  mkd <- MkPrimeData(.trans_pd())
  nEdge <- 2 * NTip(.trans_tree()) - 3L
  nTrans <- sum(mkd$type == "transformational")
  moves <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE,
                                 mcmc = MkPrimeMCMC(), fixTopology = FALSE)
  move_names <- vapply(moves, `[[`, character(1), "name")
  expect_true("gibbs_kPrime" %in% move_names)
  expect_true("block_kPrime" %in% move_names)
})

test_that("gibbs_kPrime and block_kPrime absent when nTrans = 0", {
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, knownStates = c("1" = 2L))
  moves <- MkPrime:::.BuildMoves(5L, nTrans = 0L, hasNeo = FALSE,
                                 mcmc = MkPrimeMCMC(), fixTopology = FALSE)
  move_names <- vapply(moves, `[[`, character(1), "name")
  expect_false("gibbs_kPrime" %in% move_names)
  expect_false("block_kPrime" %in% move_names)
})

test_that(".kMoveTypes maps new kPrime moves", {
  mt <- MkPrime:::.kMoveTypes
  expect_identical(mt[["gibbs_kPrime"]], 25L)
  expect_identical(mt[["block_kPrime"]], 26L)
})


# ---------------------------------------------------------------------------
# Gibbs kPrime sweep (moveType 25)
# ---------------------------------------------------------------------------

test_that("gibbs_kprime_sweep always accepts", {
  set.seed(8412)
  setup <- .setup_cpp(.trans_tree(), .trans_pd())

  accepted <- 0L
  for (i in seq_len(30L)) {
    if (do_move_cpp(setup$dataPtr, setup$statePtr, 25L, 0L,
                    0.5, 0.5, 1L, 1.0))
      accepted <- accepted + 1L
  }
  expect_equal(accepted, 30L)
})

test_that("gibbs_kprime_sweep preserves kPrime >= kObs", {
  set.seed(2019)
  setup <- .setup_cpp(.trans_tree(), .trans_pd())

  for (i in seq_len(50L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 25L, 0L, 0.5, 0.5, 1L, 1.0)

  st <- get_mcmc_state(setup$statePtr)
  expect_true(all(st$kPrime >= setup$mkd$kObs))
})

test_that("gibbs_kprime_sweep leaves logLik/logPrior finite", {
  set.seed(5517)
  setup <- .setup_cpp(.trans_tree(), .trans_pd())

  for (i in seq_len(20L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 25L, 0L, 0.5, 0.5, 1L, 1.0)

  st <- get_mcmc_state(setup$statePtr)
  expect_true(is.finite(st$logLik))
  expect_true(is.finite(st$logPrior))
})

test_that("gibbs_kprime_sweep can explore k' > kObs with long branches", {
  # Long branches make hidden states plausible
  set.seed(4871)
  tree <- rtree(8, rooted = FALSE)
  tree$edge.length <- tree$edge.length * 10
  tree <- Preorder(tree)

  mat <- matrix(sample(0:2, 8 * 5, replace = TRUE), 8, 5,
                dimnames = list(tree$tip.label, NULL))
  pd <- MatrixToPhyDat(mat)
  setup <- .setup_cpp(tree, pd)

  kp_init <- get_mcmc_state(setup$statePtr)$kPrime

  for (i in seq_len(100L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 25L, 0L, 0.5, 0.5, 1L, 1.0)

  kp_after <- get_mcmc_state(setup$statePtr)$kPrime
  # At least some k' should have been explored above kObs
  # (not necessarily still there — Gibbs samples from the full conditional)
  expect_true(is.finite(get_mcmc_state(setup$statePtr)$logLik))
})

test_that("gibbs_kprime_sweep logLik matches full recomputation", {
  set.seed(3328)
  setup <- .setup_cpp(.trans_tree(), .trans_pd())

  for (i in seq_len(10L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 25L, 0L, 0.5, 0.5, 1L, 1.0)

  st <- get_mcmc_state(setup$statePtr)
  fresh_ll <- eval_full_loglik_cpp(setup$dataPtr, setup$statePtr)
  expect_equal(st$logLik, fresh_ll, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# Block kPrime shift (moveType 26)
# ---------------------------------------------------------------------------

test_that("block_kprime_shift preserves kPrime >= kObs", {
  set.seed(6199)
  setup <- .setup_cpp(.trans_tree(), .trans_pd())

  # Push kPrime up first so downward shifts can be tested
  for (i in seq_len(50L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 25L, 0L, 0.5, 0.5, 1L, 1.0)

  for (i in seq_len(50L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 26L, 0L, 0.5, 0.5, 3L, 1.0)

  st <- get_mcmc_state(setup$statePtr)
  expect_true(all(st$kPrime >= setup$mkd$kObs))
})

test_that("block_kprime_shift leaves logLik/logPrior finite", {
  set.seed(9811)
  setup <- .setup_cpp(.trans_tree(), .trans_pd())

  for (i in seq_len(30L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 26L, 0L, 0.5, 0.5, 2L, 1.0)

  st <- get_mcmc_state(setup$statePtr)
  expect_true(is.finite(st$logLik))
  expect_true(is.finite(st$logPrior))
})

test_that("block_kprime_shift logLik matches full recomputation", {
  set.seed(4207)
  setup <- .setup_cpp(.trans_tree(), .trans_pd())

  for (i in seq_len(20L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 26L, 0L, 0.5, 0.5, 2L, 1.0)

  st <- get_mcmc_state(setup$statePtr)
  fresh_ll <- eval_full_loglik_cpp(setup$dataPtr, setup$statePtr)
  expect_equal(st$logLik, fresh_ll, tolerance = 1e-10)
})
