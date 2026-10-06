# Moves that cannot score the marginal_k target must refuse it in C++ as well
# as being left out by .BuildMoves: gibbs_p_marginal (case 35) under a heated
# chain (#269), the fixed-k' moves 11, 25 and 26 (#270), and int_walk, case 7,
# whose walk on k' nothing bounds once k' is summed out (#367).

.GuardTree <- function() {
  Preorder(read.tree(text = paste0(
    "(((t1:0.05,t2:0.07):0.04,(t3:0.06,t4:0.05):0.03):0.05,",
    "((t5:0.04,t6:0.06):0.05,(t7:0.05,t8:0.04):0.06):0.04);"
  )))
}

.GuardData <- function() {
  set.seed(4L)
  tips <- paste0("t", 1:8)
  mat <- matrix(sample(0:3, 8 * 20, replace = TRUE), 8, 20,
                dimnames = list(tips, NULL))
  MkPrimeData(MatrixToPhyDat(mat))
}

.GuardBuild <- function(p) {
  tree <- .GuardTree()
  mkd <- .GuardData()
  model <- MkPrime:::.FinalizeModel(
    MkPrimeModel(kPrimePrior = "geometric", likelihoodMode = "marginal_k"),
    tree, mkd
  )
  state0 <- MkPrime:::.InitState(tree, mkd, model)
  state0$p <- p
  state0$rate_log_sd <- 0
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)
  state0$log_prior <- eval_log_prior_cpp(
    dataPtr, MkPrime:::.InitMcmcChain(state0)
  )
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  list(dataPtr = dataPtr, statePtr = statePtr, kObs = mkd$kObs)
}

.GuardMove <- function(x, moveType, beta = 1) {
  do_move_cpp(x$dataPtr, x$statePtr, moveType = moveType, charIdx = 0L,
              scaleTuning = 1.5, betaSimplexTuning = 0, intWalkWindow = 2L,
              beta = beta)
}

test_that("case 35 leaves a heated chain on its tempered target (#269)", {
  beta <- 0.2
  # Exact tempered posterior of p on a grid, from cold marginals at each p;
  # the Beta(1, 1) hyperprior is flat.
  grid <- seq(0.0025, 0.9975, by = 0.005)
  logLik <- vapply(grid, function(p) {
    get_mcmc_state(.GuardBuild(p)$statePtr)$logLik
  }, 0)
  w <- exp(beta * (logLik - max(logLik)))
  w <- w / sum(w)
  exactMean <- sum(w * grid)
  exactSd <- sqrt(sum(w * grid ^ 2) - exactMean ^ 2)

  x <- .GuardBuild(0.3)
  set.seed(1L)
  nIter <- 5000L
  p <- numeric(nIter)
  accepted <- logical(nIter)
  for (i in seq_len(nIter)) {
    .GuardMove(x, 30L, beta)
    accepted[i] <- .GuardMove(x, 35L, beta)
    p[i] <- get_mcmc_state(x$statePtr)$p
  }
  expect_false(any(accepted))
  batchMeans <- colMeans(matrix(p, ncol = 50L))
  mcse <- stats::sd(batchMeans) / sqrt(50L)
  # Untempered, case 35 pulled the mean 0.044 low and halved the sd.
  expect_lt(abs(mean(p) - exactMean), 4 * mcse)
  expect_lt(abs(stats::sd(p) / exactSd - 1), 0.15)

  # The cold chain still uses it.
  cold <- .GuardBuild(0.3)
  set.seed(2L)
  expect_true(any(replicate(5L, .GuardMove(cold, 35L, beta = 1))))
})

test_that("k' moves refuse marginal_k on direct dispatch (#270, #367)", {
  for (moveType in c(7L, 11L, 25L, 26L)) {
    x <- .GuardBuild(0.3)
    before <- get_mcmc_state(x$statePtr)
    set.seed(3L)
    accepted <- vapply(seq_len(30L), function(i) .GuardMove(x, moveType),
                       logical(1))
    after <- get_mcmc_state(x$statePtr)
    expect_false(any(accepted))
    expect_equal(after$kPrime, x$kObs)
    expect_identical(after$edge, before$edge)
    expect_identical(after$relBrLengths, before$relBrLengths)
    expect_equal(after$logLik,
                 get_mcmc_state(.GuardBuild(after$p)$statePtr)$logLik,
                 tolerance = 1e-10)
  }
})
