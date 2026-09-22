# Issue #66: the k' Gibbs sweep's work must track the prior it is sampling,
# not the compile-time candidate cap.

SweepFixture <- function() {
  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L,
                  0L, 1L, 1L, 0L),
                nrow = 4,
                dimnames = list(paste0("t", 1:4), NULL))
  tree <- Preorder(
    ape::read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  )
  list(tree = tree, mkd = MkPrimeData(MatrixToPhyDat(mat)))
}


test_that("the k' Gibbs sweep stops at the geometric truncation cap", {
  fx <- SweepFixture()

  Candidates <- function(truncK, qHet = FALSE) {
    model <- MkPrime:::.FinalizeModel(
      MkPrimeModel(kPrimePrior = "geometric", kprimeTruncK = truncK,
                   qHeterogeneity = qHet),
      fx$tree, fx$mkd)
    state <- MkPrime:::.InitState(fx$tree, fx$mkd, model)
    # A near-flat k' prior: nothing but the cap can bound the candidate range.
    state$p <- 1e-3
    mcmcData <- MkPrime:::.InitMcmcData(fx$mkd, model)
    statePtr <- MkPrime:::.InitMcmcChain(state)
    MkPrime:::fill_partition_cache(mcmcData, statePtr)
    MkPrime:::allocate_cl_workspace(mcmcData, statePtr)
    MkPrime:::kprime_sweep_candidates(mcmcData, statePtr, 0)
  }

  kObs <- fx$mkd$kObs[fx$mkd$type == "transformational"]
  expect_equal(Candidates(12L), 12L - kObs + 1L)
  expect_equal(Candidates(24L), 24L - kObs + 1L)
  # qHeterogeneity is the regime the cap exists for: each extra candidate
  # costs an O(k^2) pruning pass there.
  expect_equal(Candidates(12L, qHet = TRUE), 12L - kObs + 1L)
})


test_that("k' draws are unchanged by a cap above the prior's live range", {
  fx <- SweepFixture()

  Sweep <- function(truncK) {
    model <- MkPrime:::.FinalizeModel(
      MkPrimeModel(kPrimePrior = "geometric", kprimeTruncK = truncK),
      fx$tree, fx$mkd)
    state <- MkPrime:::.InitState(fx$tree, fx$mkd, model)
    mcmcData <- MkPrime:::.InitMcmcData(fx$mkd, model)
    statePtr <- MkPrime:::.InitMcmcChain(state)
    MkPrime:::fill_partition_cache(mcmcData, statePtr)
    MkPrime:::allocate_cl_workspace(mcmcData, statePtr)
    set.seed(508)
    for (i in seq_len(5)) {
      MkPrime:::do_move_cpp(mcmcData, statePtr, 25L, 0L, 0.5, 0.5, 1L, 1)
    }
    MkPrime:::get_mcmc_state(statePtr)$kPrime
  }

  # 256 == kMaxKprimeCand, so this contrasts an effectively uncapped
  # enumeration with the shipped default. At the initial p the weight cutoff
  # bites far below either cap, so the two must draw identically.
  expect_equal(Sweep(256L), Sweep(200L))
})
