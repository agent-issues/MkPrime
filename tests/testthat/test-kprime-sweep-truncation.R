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


# Issue #278: the sweep enumerates at most kMaxKprimeCand = 256 candidates, but
# beta_geometric, logseries and empirical_geometric put mass beyond them. A
# k' outside the window must be left for int_walk and block_kprime_shift;
# drawing it back inside loses the tail (13% of beta_geometric mass).

BetaGeomSweepChain <- function(kPrime, alpha = 0.2, beta = 1) {
  nChar <- length(kPrime)
  mat <- matrix(rep(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L),
                    length.out = 4L * nChar),
                nrow = 4, dimnames = list(paste0("t", 1:4), NULL))
  tree <- Preorder(
    ape::read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  )
  mkd <- MkPrimeData(MatrixToPhyDat(mat))
  model <- MkPrime:::.FinalizeModel(
    MkPrimeModel(kPrimePrior = "beta_geometric", kprimeAlpha = alpha,
                 kprimeBeta = beta),
    tree, mkd)
  trans <- which(mkd$type == "transformational")
  state <- MkPrime:::.InitState(tree, mkd, model)
  state$kprime_alpha <- alpha
  state$kprime_beta <- beta
  state$kPrime[trans] <- mkd$kObs[trans] + as.integer(kPrime)
  state$log_prior <- MkPrime:::LogPrior(state, model, mkd)
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  MkPrime:::fill_partition_cache(dataPtr, statePtr)
  MkPrime:::allocate_cl_workspace(dataPtr, statePtr)
  list(data = dataPtr, state = statePtr, trans = trans,
       kObs = mkd$kObs[trans])
}

test_that("the k' sweep leaves a k' beyond its window in place", {
  ch <- BetaGeomSweepChain(c(0L, 300L, 5L))
  set.seed(278)
  MkPrime:::do_move_cpp(ch$data, ch$state, 25L, 0L, 0.5, 0.5, 1L, 0)
  u <- MkPrime:::get_mcmc_state(ch$state)$kPrime[ch$trans] - ch$kObs
  expect_equal(u[2], 300L)
})

test_that("sweep and int_walk hold the beta_geometric k' prior at beta = 0", {
  skip_under_memcheck()
  # Exact draws from P(u | alpha, beta) given u < uMax. Both moves leave
  # P(u >= 256) invariant from this start: int_walk's MH step sees the
  # truncation only within its window of uMax, far above 256.
  alpha <- 0.2
  uMax <- 1000L
  Draw <- function(n) {
    u <- integer(0)
    while (length(u) < n) {
      d <- floor(log(runif(n)) / log1p(-rbeta(n, alpha, 1)))
      u <- c(u, d[d < uMax])
    }
    u[seq_len(n)]
  }
  Tail <- function(n) exp(lbeta(alpha, 1 + n) - lbeta(alpha, 1))
  expected <- (Tail(256) - Tail(uMax)) / (1 - Tail(uMax))

  set.seed(2780)
  nChar <- 60L
  u <- unlist(lapply(seq_len(20), function(r) {
    ch <- BetaGeomSweepChain(Draw(nChar), alpha = alpha)
    for (round in 1:3) {
      MkPrime:::do_move_cpp(ch$data, ch$state, 25L, 0L, 0.5, 0.5, 1L, 0)
      for (i in ch$trans) {
        MkPrime:::do_move_cpp(ch$data, ch$state, 7L, i - 1L, 0.5, 0.5, 5L, 0)
      }
    }
    MkPrime:::get_mcmc_state(ch$state)$kPrime[ch$trans] - ch$kObs
  }))
  # Characters are independent at beta = 0, so the count is binomial.
  se <- sqrt(expected * (1 - expected) / length(u))
  expect_lt(abs(mean(u >= 256) - expected), 4 * se,
            label = paste("P(u >= 256) =", mean(u >= 256), "vs", signif(expected, 3)))
})
