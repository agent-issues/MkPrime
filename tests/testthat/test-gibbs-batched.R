# Tests for M-155: batched Gibbs kPrime sweep.
# Verifies that the batched precomputation + sampling produces the same
# logLik as full recomputation, across various partition configurations.

library("ape")
library("TreeTools")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.setup_cpp <- function(tree, pd, model = MkPrimeModel()) {
  tree <- Preorder(tree)
  mkd   <- MkPrimeData(pd)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  state0   <- MkPrime:::.InitState(tree, mkd, model)
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  list(mkd = mkd, dataPtr = dataPtr, statePtr = statePtr)
}

# Multi-partition trans data: chars with kObs=2, 3, 4 → 3 trans partitions
.multi_kobs_pd <- function(nTip = 10) {
  labs <- paste0("t", seq_len(nTip))
  # 5 binary chars (kObs=2), 3 ternary (kObs=3), 2 quaternary (kObs=4)
  m <- cbind(
    matrix(sample(0:1, nTip * 5, TRUE), nTip, 5),
    matrix(sample(0:2, nTip * 3, TRUE), nTip, 3),
    matrix(sample(0:3, nTip * 2, TRUE), nTip, 2)
  )
  dimnames(m) <- list(labs, NULL)
  MatrixToPhyDat(m)
}


# ---------------------------------------------------------------------------
# Batched Gibbs sweep — logLik consistency
# ---------------------------------------------------------------------------

test_that("batched Gibbs sweep logLik matches recomputation (multi-kObs)", {
  set.seed(6394)
  tree <- rtree(10, rooted = FALSE)
  tree$edge.length <- tree$edge.length * 3
  tree <- Preorder(tree)
  pd <- .multi_kobs_pd(10)
  setup <- .setup_cpp(tree, pd)

  for (i in seq_len(30L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 25L, 0L, 0.5, 0.5, 1L, 1.0)

  st <- get_mcmc_state(setup$statePtr)
  fresh_ll <- eval_full_loglik_cpp(setup$dataPtr, setup$statePtr)
  expect_equal(st$logLik, fresh_ll, tolerance = 1e-10)
  expect_true(all(st$kPrime >= setup$mkd$kObs))
})

test_that("batched Gibbs sweep logLik matches after many sweeps", {
  set.seed(2781)
  tree <- rtree(12, rooted = FALSE)
  tree$edge.length <- tree$edge.length * 5
  tree <- Preorder(tree)
  pd <- .multi_kobs_pd(12)
  setup <- .setup_cpp(tree, pd)

  for (i in seq_len(100L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 25L, 0L, 0.5, 0.5, 1L, 1.0)

  st <- get_mcmc_state(setup$statePtr)
  fresh_ll <- eval_full_loglik_cpp(setup$dataPtr, setup$statePtr)
  expect_equal(st$logLik, fresh_ll, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# Interleaving Gibbs sweeps with other moves
# ---------------------------------------------------------------------------

test_that("Gibbs sweep + tree moves maintain consistency", {
  set.seed(7043)
  tree <- rtree(10, rooted = FALSE)
  tree$edge.length <- tree$edge.length * 2
  tree <- Preorder(tree)
  pd <- .multi_kobs_pd(10)
  setup <- .setup_cpp(tree, pd)

  # Interleave Gibbs (25), int_walk kPrime (7), NNI (1), scale treeLength (3)
  move_seq <- rep(c(25L, 7L, 1L, 3L, 25L), 10L)
  for (mt in move_seq)
    do_move_cpp(setup$dataPtr, setup$statePtr, mt, 0L, 0.5, 0.5, 1L, 1.0)

  st <- get_mcmc_state(setup$statePtr)
  fresh_ll <- eval_full_loglik_cpp(setup$dataPtr, setup$statePtr)
  expect_equal(st$logLik, fresh_ll, tolerance = 1e-10)
  expect_true(all(st$kPrime >= setup$mkd$kObs))
})


# ---------------------------------------------------------------------------
# With ACRV (nCat > 1)
# ---------------------------------------------------------------------------

test_that("batched Gibbs sweep correct with ACRV nCat=4", {
  set.seed(5019)
  tree <- rtree(8, rooted = FALSE)
  tree$edge.length <- tree$edge.length * 4
  tree <- Preorder(tree)

  m <- cbind(
    matrix(sample(0:1, 8 * 4, TRUE), 8, 4),
    matrix(sample(0:2, 8 * 3, TRUE), 8, 3)
  )
  dimnames(m) <- list(tree$tip.label, NULL)
  pd <- MatrixToPhyDat(m)

  model <- MkPrimeModel(nCat = 4L)
  setup <- .setup_cpp(tree, pd, model)

  for (i in seq_len(40L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 25L, 0L, 0.5, 0.5, 1L, 1.0)

  st <- get_mcmc_state(setup$statePtr)
  fresh_ll <- eval_full_loglik_cpp(setup$dataPtr, setup$statePtr)
  expect_equal(st$logLik, fresh_ll, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# Edge case: single transformational character
# ---------------------------------------------------------------------------

test_that("batched Gibbs sweep works with single trans char", {
  set.seed(1847)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  tree <- Preorder(tree)

  mat <- matrix(c(0, 1, 2, 0), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)
  setup <- .setup_cpp(tree, pd)

  for (i in seq_len(20L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 25L, 0L, 0.5, 0.5, 1L, 1.0)

  st <- get_mcmc_state(setup$statePtr)
  fresh_ll <- eval_full_loglik_cpp(setup$dataPtr, setup$statePtr)
  expect_equal(st$logLik, fresh_ll, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# M-172: pattern deduplication in Gibbs sweep
# ---------------------------------------------------------------------------

test_that("Gibbs sweep logLik matches recomputation with fully duplicated patterns", {
  # All transformational characters share one of two patterns — maximises dedup
  set.seed(3861)
  nTip <- 10
  labs <- paste0("t", seq_len(nTip))
  pat1 <- sample(0:1, nTip, replace = TRUE)
  pat2 <- sample(0:1, nTip, replace = TRUE)
  # 8 characters: alternating pat1 / pat2
  m <- matrix(rep(c(pat1, pat2), 4)[seq_len(nTip * 8)],
              nrow = nTip, ncol = 8)
  for (j in seq_len(8)) m[, j] <- if (j %% 2 == 1) pat1 else pat2
  dimnames(m) <- list(labs, NULL)
  pd <- MatrixToPhyDat(m)

  tree <- rtree(nTip, rooted = FALSE)
  tree$edge.length <- tree$edge.length * 3
  tree <- Preorder(tree)
  setup <- .setup_cpp(tree, pd)

  # Confirm partition has deduplication
  mkd <- setup$mkd
  part <- mkd$partitions[[1]]
  expect_lt(ncol(part$unique_tip_states), part$nChar)

  for (i in seq_len(40L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 25L, 0L, 0.5, 0.5, 1L, 1.0)

  st <- get_mcmc_state(setup$statePtr)
  fresh_ll <- eval_full_loglik_cpp(setup$dataPtr, setup$statePtr)
  expect_equal(st$logLik, fresh_ll, tolerance = 1e-10)
  expect_true(all(st$kPrime >= setup$mkd$kObs))
})

test_that("Gibbs sweep logLik correct with dedup + ACRV nCat=6", {
  set.seed(7752)
  nTip <- 8
  labs <- paste0("t", seq_len(nTip))
  base_pat <- sample(0:2, nTip, replace = TRUE)
  # 6 characters: first 3 identical, last 3 distinct
  m <- cbind(
    matrix(rep(base_pat, 3), nTip, 3),
    matrix(sample(0:2, nTip * 3, replace = TRUE), nTip, 3)
  )
  dimnames(m) <- list(labs, NULL)
  pd <- MatrixToPhyDat(m)

  tree <- rtree(nTip, rooted = FALSE)
  tree$edge.length <- tree$edge.length * 2
  tree <- Preorder(tree)
  model <- MkPrimeModel(nCat = 6L)
  setup <- .setup_cpp(tree, pd, model)

  for (i in seq_len(50L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 25L, 0L, 0.5, 0.5, 1L, 1.0)

  st <- get_mcmc_state(setup$statePtr)
  fresh_ll <- eval_full_loglik_cpp(setup$dataPtr, setup$statePtr)
  expect_equal(st$logLik, fresh_ll, tolerance = 1e-10)
})

test_that("Gibbs sweep with dedup + ascertainment correction is correct", {
  set.seed(6103)
  nTip <- 10
  labs <- paste0("t", seq_len(nTip))
  pat  <- sample(0:1, nTip, replace = TRUE)
  # 6 identical binary patterns + 2 unique
  m <- cbind(
    matrix(rep(pat, 6), nTip, 6),
    matrix(sample(0:1, nTip * 2, replace = TRUE), nTip, 2)
  )
  dimnames(m) <- list(labs, NULL)
  pd <- MatrixToPhyDat(m)

  tree <- rtree(nTip, rooted = FALSE)
  tree <- Preorder(tree)
  model <- MkPrimeModel(coding = "variable")
  setup <- .setup_cpp(tree, pd, model)

  for (i in seq_len(30L))
    do_move_cpp(setup$dataPtr, setup$statePtr, 25L, 0L, 0.5, 0.5, 1L, 1.0)

  st <- get_mcmc_state(setup$statePtr)
  fresh_ll <- eval_full_loglik_cpp(setup$dataPtr, setup$statePtr)
  expect_equal(st$logLik, fresh_ll, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Larger dataset (Sun2018 subset) if available
# ---------------------------------------------------------------------------

test_that("batched Gibbs sweep on Sun2018 subset matches recomputation", {
  nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
  skip_if(nexFile == "", message = "TreeSearch not available")

  set.seed(9273)
  pd <- TreeTools::ReadAsPhyDat(nexFile)
  # Subset to 20 taxa for speed
  keep <- sample(names(pd), 20)
  pd <- pd[keep]

  tree <- rtree(20, tip.label = keep, rooted = FALSE)
  tree$edge.length <- tree$edge.length * 3
  tree <- Preorder(tree)

  neo <- MkPrime::AutoDetectNeomorphic(pd)
  mkd <- MkPrimeData(pd, neomorphic = neo)
  nTrans <- sum(mkd$type == "transformational")
  skip_if(nTrans == 0, "No transformational characters in subset")

  model <- MkPrimeModel(nCat = 4L)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  state0   <- MkPrime:::.InitState(tree, mkd, model)
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)

  for (i in seq_len(20L))
    do_move_cpp(dataPtr, statePtr, 25L, 0L, 0.5, 0.5, 1L, 1.0)

  st <- get_mcmc_state(statePtr)
  fresh_ll <- eval_full_loglik_cpp(dataPtr, statePtr)
  expect_equal(st$logLik, fresh_ll, tolerance = 1e-10)
  expect_true(all(st$kPrime >= mkd$kObs))
})
