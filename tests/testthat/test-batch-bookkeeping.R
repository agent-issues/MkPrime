# What run_mcmc_batch_cpp() carries between batches and per chain (#304,
# #306), and how it hands back an interrupt (#398, #317).

.BatchFixture <- function(nChains, moveTypes = c(26L, 0L)) {
  tree <- Preorder(read.tree(
    text = "(((t1:0.1,t2:0.2):0.1,(t3:0.1,t4:0.3):0.2):0.1,(t5:0.2,t6:0.1):0.1);"
  ))
  mat <- matrix(c(0, 1, 2, 0, 1, 2,
                  0, 1, 0, 1, 1, 0,
                  2, 0, 1, 1, 0, 2), 6,
                dimnames = list(paste0("t", 1:6), NULL))
  mkd <- MkPrimeData(MatrixToPhyDat(mat))
  model <- .FinalizeModel(MkPrimeModel(), tree, mkd)
  state0 <- .InitState(tree, mkd, model)
  dataPtr <- .InitMcmcData(mkd, model)
  states <- lapply(seq_len(nChains), function(ch) {
    statePtr <- .InitMcmcChain(state0)
    fill_partition_cache(dataPtr, statePtr)
    allocate_cl_workspace(dataPtr, statePtr)
    statePtr
  })
  nMoves <- length(moveTypes)
  Batch <- function(nBatch, startIter, betas = 0.8^(seq_len(nChains) - 1),
                    blockKpWins = integer(0), ptCarry = NULL, warmup = 0L,
                    thin = 1L, run = run_mcmc_batch_cpp) {
    run(
      dataPtr, states, betas,
      moveTypes, which(mkd$type == "transformational") - 1L,
      integer(nMoves), rep(1, nMoves),
      matrix(0.5, nChains, nMoves), rep(10, nChains), rep(1L, nChains),
      integer(nMoves), matrix(1, nChains, nMoves), matrix(0, nChains, nMoves),
      nBatch, startIter, warmup, thin, FALSE, nrow(tree$edge), 1,
      blockKpWins, ptCarry
    )
  }
  list(Batch = Batch, states = states,
       nTrans = sum(mkd$type == "transformational"), nEdge = nrow(tree$edge))
}

test_that("each chain shifts block_kPrime by its own window (#304)", {
  set.seed(1304)
  fx <- .BatchFixture(1L, moveTypes = 26L)
  shifts <- integer(0)
  for (i in seq_len(200)) {
    before <- get_mcmc_state(fx$states[[1]])$kPrime + 0L
    fx$Batch(1L, i, betas = 0, blockKpWins = 6L)
    shifts <- c(shifts, get_mcmc_state(fx$states[[1]])$kPrime[1] - before[1])
  }
  # The chain-level int_walk window is 1; only the block window allows more.
  expect_gt(max(abs(shifts)), 1L)
  expect_lte(max(abs(shifts)), 6L)
})

test_that("round trips and cold swaps carry across batch boundaries (#306)", {
  whole <- .BatchFixture(3L)
  split <- .BatchFixture(3L)
  set.seed(2306)
  one <- whole$Batch(400L, 1L, thin = 7L)
  set.seed(2306)
  first <- split$Batch(137L, 1L, thin = 7L)
  second <- split$Batch(263L, 138L, ptCarry = first$pt_carry, thin = 7L)

  expect_gt(one$round_trip_count, 0L)
  expect_identical(first$round_trip_count + second$round_trip_count,
                   one$round_trip_count)
  # swap_cold precedes topo_hash and the k' and branch columns.
  swapCold <- ncol(one$scalar_samples) - whole$nTrans - whole$nEdge - 1L
  expect_identical(
    c(first$scalar_samples[, swapCold], second$scalar_samples[, swapCold]),
    one$scalar_samples[, swapCold]
  )
})

test_that("an interrupted batch unwinds normally and advances the seed (#398, #317)", {
  skip_on_cran()
  skip_under_memcheck()
  fx <- .BatchFixture(1L)
  set.seed(398)
  seed <- .Random.seed
  setTimeLimit(elapsed = 0.5, transient = TRUE)
  on.exit(setTimeLimit(), add = TRUE)
  err <- tryCatch(fx$Batch(1e9L, 1L, thin = 1e6L, run = .RunBatch),
                  error = identity)
  setTimeLimit()
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "reached elapsed time limit")
  # Previously the longjmp skipped RNGScope, so the next R draw replayed the
  # uniforms the batch had consumed.
  expect_false(identical(.Random.seed, seed))
})

test_that("a heat outside the adaptation clamp warns (#306)", {
  expect_warning(MkPrimeMCMC(nChains = 4L, heat = 0.8),
                 "outside \\[0.01, 0.5\\]")
  expect_warning(MkPrimeMCMC(nChains = 4L, heat = 0.005), "0.01")
  expect_no_warning(MkPrimeMCMC(nChains = 4L, heat = 0.5))
  expect_no_warning(MkPrimeMCMC(nChains = 1L, heat = 0.8))
})

test_that("the ladder adapts on the swaps since it last moved (#306)", {
  skip_on_cran()
  skip_under_memcheck()
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  seen <- new.env()
  seen$propose <- numeric(0)
  seen$moved <- logical(0)
  adapt <- MkPrime:::.AdaptTemperatures
  local_mocked_bindings(
    .AdaptTemperatures = function(betas, swapAccept, swapPropose, ...) {
      out <- adapt(betas, swapAccept, swapPropose, ...)
      seen$propose <- c(seen$propose, sum(swapPropose))
      seen$moved <- c(seen$moved,
                      !identical(as.numeric(out), as.numeric(betas)))
      out
    }
  )
  set.seed(306)
  allow_warning(
    result <- RunMkPrime(MatrixToPhyDat(mat), tree, mcmc = MkPrimeMCMC(
      nRuns = 1L, nIter = 3500L, thin = 10L, maxWarmup = 3000L,
      minWarmup = 3000L, autoTune = FALSE, nChains = 3L, heat = 0.3,
      maxTime = 60
    )),
    "without stabilisation"
  )
  expect_length(seen$propose, 6L)
  expect_true(any(seen$moved[-6L]))
  # One swap is proposed per iteration, over 500-iteration batches: a ladder
  # that moved starts its count afresh, one that held keeps accumulating.
  expected <- numeric(6L)
  window <- 0
  for (k in seq_len(6L)) {
    window <- window + 500
    expected[k] <- window
    if (seen$moved[k]) window <- 0
  }
  expect_equal(seen$propose, expected)
})
