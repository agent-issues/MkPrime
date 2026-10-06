# Issue #281: stepping out must split its budget between the two sides at
# random (Neal 2003, Fig. 3). A fixed cap per side is not reversible once it
# binds: with a width far below the slice's, a point near one end of the slice
# reaches that end while the other side stops at its cap, and no point near
# the other end generates the same interval. At beta = 0 the target is the
# prior, so starting from exact draws tests pi K^n = pi directly. The fixed
# cap gives KS p = 4e-5 and 3e-7 on the two samplers here.

SliceFixture <- function(kPrimePrior = "geometric") {
  mat <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L), nrow = 4,
                dimnames = list(paste0("t", 1:4), NULL))
  tree <- Preorder(
    ape::read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  )
  mkd <- MkPrimeData(MatrixToPhyDat(mat))
  list(tree = tree, mkd = mkd,
       model = MkPrime:::.FinalizeModel(
         MkPrimeModel(kPrimePrior = kPrimePrior), tree, mkd))
}

test_that("narrow slice_scalar holds the rate_log_sd prior at beta = 0", {
  skip_under_memcheck()
  fx <- SliceFixture()
  shape <- fx$model$rateLogSdShape
  rate <- fx$model$rateLogSdRate
  set.seed(2810)
  x <- vapply(seq_len(1000), function(r) {
    state <- MkPrime:::.InitState(fx$tree, fx$mkd, fx$model)
    state$rate_log_sd <- rgamma(1, shape, rate)
    state$log_prior <- MkPrime:::LogPrior(state, fx$model, fx$mkd)
    dataPtr <- MkPrime:::.InitMcmcData(fx$mkd, fx$model)
    statePtr <- MkPrime:::.InitMcmcChain(state)
    MkPrime:::fill_partition_cache(dataPtr, statePtr)
    MkPrime:::allocate_cl_workspace(dataPtr, statePtr)
    for (i in 1:50) {
      MkPrime:::do_move_cpp(dataPtr, statePtr, 19L, 2L, 0.05, 0.5, 1L, 0)
    }
    MkPrime:::get_mcmc_state(statePtr)$rateLogSd
  }, numeric(1))
  p <- ks.test(x, "pgamma", shape, rate)$p.value
  expect_gt(p, 1e-4, label = paste("KS p =", signif(p, 2)))
})

test_that("narrow slice_kprime_hyper holds the beta_geometric prior", {
  skip_under_memcheck()
  fx <- SliceFixture("beta_geometric")
  trans <- which(fx$mkd$type == "transformational")
  # Exact joint draws of (alpha, beta, u); the move updates (alpha, beta)
  # given u, so conditioning the draws on u keeps them exact.
  Draw <- function() {
    repeat {
      a <- rexp(1)
      b <- rexp(1)
      u <- floor(log(runif(length(trans))) /
                   log1p(-rbeta(length(trans), a, b)))
      if (all(u < 1e6)) return(list(a = a, b = b, u = u))
    }
  }
  set.seed(2811)
  s <- vapply(seq_len(3000), function(r) {
    d <- Draw()
    state <- MkPrime:::.InitState(fx$tree, fx$mkd, fx$model)
    state$kprime_alpha <- d$a
    state$kprime_beta <- d$b
    state$kPrime[trans] <- fx$mkd$kObs[trans] + as.integer(d$u)
    state$log_prior <- MkPrime:::LogPrior(state, fx$model, fx$mkd)
    dataPtr <- MkPrime:::.InitMcmcData(fx$mkd, fx$model)
    statePtr <- MkPrime:::.InitMcmcChain(state)
    # paramCode 0: s = log(alpha + beta), which the cap binds on most.
    for (i in 1:50) {
      MkPrime:::do_move_cpp(dataPtr, statePtr, 29L, 0L, 0.05, 0.5, 1L, 0)
    }
    st <- MkPrime:::get_mcmc_state(statePtr)
    log(st$kprimeAlpha + st$kprimeBeta)
  }, numeric(1))
  reference <- vapply(seq_len(30000), function(i) {
    d <- Draw()
    log(d$a + d$b)
  }, numeric(1))
  p <- ks.test(s, reference)$p.value
  expect_gt(p, 1e-4, label = paste("KS p =", signif(p, 2)))
})
