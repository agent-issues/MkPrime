# Issue #65: each stone must be sampled with the schedule the caller asked for.

FourTipData <- function() {
  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L,
                  0L, 1L, 1L, 0L),
                nrow = 4,
                dimnames = list(paste0("t", 1:4), NULL))
  MatrixToPhyDat(mat)
}

FourTipTree <- function() {
  Preorder(
    ape::read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  )
}

# Names of the moves a stepping-stone run actually dispatches.
DispatchedMoves <- function(pd, tree, model = NULL, mcmc = NULL,
                            nIter = 200L) {
  realDoMove <- MkPrime:::.DoMove
  seen <- character(0)
  local_mocked_bindings(
    .DoMove = function(move, ...) {
      seen <<- c(seen, move$name)
      realDoMove(move, ...)
    },
    .package = "MkPrime"
  )
  set.seed(3607)
  mkp_stepping_stone(pd, tree, model = model, mcmc = mcmc, nStones = 2L,
                     nIter = nIter, warmup = 0L, verbose = FALSE)
  seen
}


test_that("stepping stone dispatches the moves `model` and `mcmc` specify", {
  pd <- FourTipData()
  tree <- FourTipTree()
  Dispatched <- function(...) unique(DispatchedMoves(pd, tree, ...))

  atDefaults <- Dispatched()
  expect_true(all(c("tbr", "pspr", "gibbs_spr", "gibbs_subtree_swap",
                    "dirichlet_branch", "local_dirichlet") %in% atDefaults))
  expect_false("beta_scale" %in% atDefaults)

  expect_true(all(c("beta_scale", "slice_beta_scale") %in%
                    Dispatched(model = MkPrimeModel(qHeterogeneity = TRUE))))

  expect_false("tbr" %in% Dispatched(mcmc = MkPrimeMCMC(tbr = FALSE)))
})


test_that("stepping stone follows the model's k' arm", {
  pd <- FourTipData()
  tree <- FourTipTree()
  Dispatched <- function(...) unique(DispatchedMoves(pd, tree, ...))

  bg <- Dispatched(model = MkPrimeModel(kPrimePrior = "beta_geometric"))
  expect_true(all(c("slice_kprime_s", "slice_kprime_r") %in% bg))
  expect_false("mh_logit_p" %in% bg)

  marginal <- Dispatched(model = MkPrimeModel(kPrimePrior = "geometric",
                                              likelihoodMode = "marginal_k"))
  expect_true("mh_logit_p" %in% marginal)
  expect_false(any(c("kPrime", "gibbs_kPrime", "block_kPrime") %in% marginal))
})


test_that("stepping stone declines the untempered marginal-k p move", {
  pd <- FourTipData()
  tree <- FourTipTree()
  model <- MkPrimeModel(kPrimePrior = "geometric",
                        likelihoodMode = "marginal_k")
  mcmc <- MkPrimeMCMC(gibbsPMarginal = TRUE)

  expect_warning(
    seen <- DispatchedMoves(pd, tree, model = model, mcmc = mcmc,
                            nIter = 40L),
    "gibbsPMarginal"
  )
  expect_false("gibbs_p_marginal" %in% seen)
})


test_that("stepping stone honours pinned move weights", {
  seen <- DispatchedMoves(FourTipData(), FourTipTree(),
                          mcmc = MkPrimeMCMC(moveWeights = c(nni = 0.9)))
  expect_gt(mean(seen == "nni"), 0.7)
})


test_that("stepping stone wires nBranchBins into the C++ data struct", {
  seen <- NULL
  realSetBins <- MkPrime:::set_branch_bins
  local_mocked_bindings(
    set_branch_bins = function(dataPtr, nBins) {
      seen <<- nBins
      realSetBins(dataPtr, nBins)
    },
    .package = "MkPrime"
  )
  mkp_stepping_stone(FourTipData(), FourTipTree(),
                     mcmc = MkPrimeMCMC(nBranchBins = 7L),
                     nStones = 2L, nIter = 20L, warmup = 0L, verbose = FALSE)
  expect_equal(seen, 7L)
})


test_that("a stone's slice moves sample their own parameter", {
  tree <- FourTipTree()
  mkd <- MkPrimeData(FourTipData())
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(), tree, mkd)
  state <- MkPrime:::.InitState(tree, mkd, model)
  mcmcData <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  MkPrime:::fill_partition_cache(mcmcData, statePtr)
  MkPrime:::allocate_cl_workspace(mcmcData, statePtr)

  move <- list(name = "slice_rate_log_sd", type = "slice",
               target = "rate_log_sd", weight = 1, dim = 1L,
               sliceParamIdx = 2L)
  before <- MkPrime:::get_mcmc_state(statePtr)
  set.seed(4410)
  accepted <- vapply(1:20, function(i) {
    MkPrime:::.DoMove(move, statePtr, tuning = MkPrimeMCMC()$tuning,
                      beta = 0.5, mcmcData = mcmcData)$accept
  }, logical(1))
  after <- MkPrime:::get_mcmc_state(statePtr)

  expect_true(any(accepted))
  expect_false(isTRUE(all.equal(before$rateLogSd, after$rateLogSd)))
  # paramIdx 2 is rate_log_sd; tree_length (paramIdx 0) must be untouched.
  expect_equal(before$treeLength, after$treeLength)
})
