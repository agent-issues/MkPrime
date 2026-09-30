# #255: the charLL warm cache keeps raw per-(char, k') LLs, which do not
# depend on p, but the candidate support they are summed over does. A cache
# filled at high p stops short of the tail a lower p favours, so chained
# case-30 (mh_logit_p) moves must re-derive the support at their own p and
# refill when the fill never evaluated it. The oracle is the forced-cold
# logLik that fill_partition_cache stores on a fresh pointer at the same p;
# eval_full_loglik_cpp would go through the warm path under test.

.SupportTree <- function() {
  Preorder(read.tree(text = paste0(
    "(((t1:0.05,t2:0.07):0.04,(t3:0.06,t4:0.05):0.03):0.05,",
    "((t5:0.04,t6:0.06):0.05,(t7:0.05,t8:0.04):0.06):0.04);"
  )))
}

.SupportData <- function() {
  set.seed(4L)
  tips <- paste0("t", 1:8)
  mat <- matrix(sample(0:3, 8 * 20, replace = TRUE), 8, 20,
                dimnames = list(tips, NULL))
  MkPrimeData(MatrixToPhyDat(mat))
}

# relabel = TRUE (the default) is what makes the corrected LL rise in k' and
# the truncation error O(1). kprimeHyperB only drives p down from the fill; the
# LL, and so the oracle, does not depend on it.
.SupportBuild <- function(p, hyperB = 1000) {
  tree <- .SupportTree()
  mkd <- .SupportData()
  model <- MkPrimeModel(kPrimePrior = "geometric",
                        likelihoodMode = "marginal_k",
                        kprimeHyperB = hyperB)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  state0 <- MkPrime:::.InitState(tree, mkd, model)
  state0$p <- p
  state0$rate_log_sd <- 0
  state0$tree_length <- sum(tree$edge.length)
  state0$rel_br_lengths <- tree$edge.length / state0$tree_length
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)
  state0$log_prior <- eval_log_prior_cpp(
    dataPtr, MkPrime:::.InitMcmcChain(state0)
  )
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  list(dataPtr = dataPtr, statePtr = statePtr)
}

.ColdLogLik <- function(p) {
  get_mcmc_state(.SupportBuild(p)$statePtr)$logLik
}

.PMove <- function(x, moveType = 30L) {
  do_move_cpp(x$dataPtr, x$statePtr, moveType = moveType, charIdx = 0L,
              scaleTuning = 0.5, betaSimplexTuning = 0, intWalkWindow = 1L,
              beta = 1)
}

# Fill at p = 0.95, then chain case-30 moves down towards p ~ 0.07.
.StaleChain <- function(nMove = 400L) {
  x <- .SupportBuild(0.95)
  set.seed(1L)
  errors <- numeric(0)
  for (i in seq_len(nMove)) {
    if (.PMove(x)) {
      st <- get_mcmc_state(x$statePtr)
      errors <- c(errors, st$logLik - .ColdLogLik(st$p))
    }
  }
  list(x = x, errors = errors)
}

test_that("chained case-30 moves score the cold marginal at every p (#255)", {
  chain <- .StaleChain()
  expect_gt(length(chain$errors), 20L)
  expect_lt(get_mcmc_state(chain$x$statePtr)$p, 0.3)
  expect_lt(max(abs(chain$errors)), 1e-8)

  # A repeat evaluation at the final p reproduces the stored value.
  expect_equal(eval_full_loglik_cpp(chain$x$dataPtr, chain$x$statePtr),
               get_mcmc_state(chain$x$statePtr)$logLik, tolerance = 1e-12)
})

test_that("case 35 imputes over the support at its own p (#255)", {
  stale <- .StaleChain()$x
  p <- get_mcmc_state(stale$statePtr)$p
  cold <- .SupportBuild(p)

  # Both pointers draw nTrans uniforms and one rbeta, so with the same seed
  # and the same support the proposed (and nearly always accepted) p* agree.
  set.seed(7L)
  .PMove(stale, 35L)
  set.seed(7L)
  .PMove(cold, 35L)
  pStar <- get_mcmc_state(cold$statePtr)$p
  expect_false(pStar == p)
  expect_equal(get_mcmc_state(stale$statePtr)$p, pStar, tolerance = 1e-12)
})
