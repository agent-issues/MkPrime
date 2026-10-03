# The C++ Metropolis-Hastings step scales the log-likelihood ratio by beta
# (power posterior, #331). Driven through the XPtr state, as production
# tempering and stepping stone do.

HeatedAcceptCount <- function(beta, nMoves = 400L) {
  tree <- ape::read.tree(
    text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  mkd <- MkPrimeData(MatrixToPhyDat(mat))
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(), tree, mkd)
  tree <- TreeTools::Preorder(tree)
  state <- MkPrime:::.InitState(tree, mkd, model)
  mcmcData <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  fill_partition_cache(mcmcData, statePtr)

  tuning <- MkPrimeMCMC()$tuning
  tuning$scale_tree_length <- 6
  move <- list(name = "tree_length", type = "scale",
               target = "tree_length", weight = 1)
  set.seed(7562)
  sum(vapply(seq_len(nMoves), function(i) {
    MkPrime:::.DoMove(move, statePtr, tuning = tuning, beta = beta,
                      mcmcData = mcmcData)$accept
  }, logical(1)))
}

test_that("C++ MH step: heated chain accepts more than cold (#331)", {
  expect_gt(HeatedAcceptCount(0.05), HeatedAcceptCount(1))
})

test_that("C++ MH step: beta = 0 ignores the likelihood (#331)", {
  expect_gt(HeatedAcceptCount(0), HeatedAcceptCount(1))
})
