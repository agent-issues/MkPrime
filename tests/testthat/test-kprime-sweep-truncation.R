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


# Issue #392: every k' prior is truncated at K and renormalised over its capped
# support, so the sweep enumerates the whole support (#278) and int_walk never
# wanders into a tail it cannot leave.

BinaryChain <- function(kPrimePrior, nChar, truncK, ...) {
  mat <- matrix(rep(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 0L),
                    length.out = 4L * nChar),
                nrow = 4, dimnames = list(paste0("t", 1:4), NULL))
  tree <- Preorder(
    ape::read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  )
  mkd <- MkPrimeData(MatrixToPhyDat(mat))
  model <- MkPrime:::.FinalizeModel(
    MkPrimeModel(kPrimePrior = kPrimePrior, kprimeTruncK = truncK, ...),
    tree, mkd)
  list(model = model, mkd = mkd,
       state = MkPrime:::.InitState(tree, mkd, model),
       trans = which(mkd$type == "transformational"))
}

StartChain <- function(ch, state) {
  state$log_prior <- MkPrime:::LogPrior(state, ch$model, ch$mkd)
  dataPtr <- MkPrime:::.InitMcmcData(ch$mkd, ch$model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  MkPrime:::fill_partition_cache(dataPtr, statePtr)
  MkPrime:::allocate_cl_workspace(dataPtr, statePtr)
  list(data = dataPtr, state = statePtr)
}

test_that("the k' sweep stops at K under every prior", {
  Candidates <- function(kPrimePrior, Flatten) {
    ch <- BinaryChain(kPrimePrior, 3L, 12L)
    # A near-flat k' prior: nothing but the cap can bound the candidate range.
    ptr <- StartChain(ch, Flatten(ch$state))
    MkPrime:::kprime_sweep_candidates(ptr$data, ptr$state, 0)
  }
  expect_equal(Candidates("beta_geometric", function(st) {
    st$kprime_alpha <- 1e-3
    st
  }), rep(11L, 3))
  expect_equal(Candidates("empirical_geometric", function(st) {
    st$p <- 1e-3
    st
  }), rep(11L, 3))
})

test_that("the k' sweep redraws a k' above K into the support", {
  # Only a checkpoint written before the cap holds such a state.
  ch <- BinaryChain("beta_geometric", 3L, 30L)
  state <- ch$state
  state$kPrime[ch$trans] <- c(2L, 300L, 7L)
  ptr <- StartChain(ch, state)
  set.seed(392)
  MkPrime:::do_move_cpp(ptr$data, ptr$state, 25L, 0L, 0.5, 0.5, 1L, 1)
  expect_lte(max(MkPrime:::get_mcmc_state(ptr$state)$kPrime), 30L)
})

# At beta = 0 the target is the capped prior, and starting from exact joint
# draws tests pi K^n = pi. The hyperparameter moves see the cap only through
# its normaliser, so their marginals detect a missing or wrong one.
test_that("k' and (alpha, beta) moves hold the capped beta_geometric prior", {
  skip_under_memcheck()
  truncK <- 30L
  ch <- BinaryChain("beta_geometric", 20L, truncK)
  n <- truncK - ch$mkd$kObs[ch$trans] + 1L
  # u | (alpha, beta) by rejection from the untruncated mixture, exactly.
  DrawU <- function(a, b) {
    u <- rep(NA_real_, length(n))
    while (anyNA(u)) {
      todo <- which(is.na(u))
      d <- floor(log(runif(length(todo))) /
                   log1p(-rbeta(length(todo), a, b)))
      u[todo] <- ifelse(d < n[todo], d, NA_real_)
    }
    u
  }
  set.seed(3920)
  hyper <- t(vapply(seq_len(400), function(r) {
    state <- ch$state
    state$kprime_alpha <- rexp(1)
    state$kprime_beta <- rexp(1)
    state$kPrime[ch$trans] <- ch$mkd$kObs[ch$trans] +
      as.integer(DrawU(state$kprime_alpha, state$kprime_beta))
    ptr <- StartChain(ch, state)
    for (round in 1:3) {
      MkPrime:::do_move_cpp(ptr$data, ptr$state, 25L, 0L, 0.5, 0.5, 1L, 0)
      for (i in ch$trans) {
        MkPrime:::do_move_cpp(ptr$data, ptr$state, 7L, i - 1L, 0.5, 0.5, 3L, 0)
      }
      for (param in 0:1) for (i in 1:5) {
        MkPrime:::do_move_cpp(ptr$data, ptr$state, 29L, param, 0.5, 0.5, 1L, 0)
      }
    }
    st <- MkPrime:::get_mcmc_state(ptr$state)
    c(st$kprimeAlpha, st$kprimeBeta)
  }, numeric(2)))
  pA <- ks.test(hyper[, 1], "pexp")$p.value
  pB <- ks.test(hyper[, 2], "pexp")$p.value
  expect_gt(pA, 1e-3, label = paste("alpha: KS p =", signif(pA, 2)))
  expect_gt(pB, 1e-3, label = paste("beta: KS p =", signif(pB, 2)))
})

test_that("k' and p moves hold the capped empirical_geometric prior", {
  skip_under_memcheck()
  truncK <- 10L
  ch <- BinaryChain("empirical_geometric", 30L, truncK)
  a <- ch$model$kprimeHyperA
  b <- ch$model$kprimeHyperB
  emp <- ch$model$empiricalNObs
  set.seed(3921)
  p <- vapply(seq_len(400), function(r) {
    state <- ch$state
    state$p <- rbeta(1, a, b)
    # Every kObs is 2, so the Model A prior on [2, K] is the conditional.
    pmf <- vapply(2:truncK, function(k) {
      MkPrime:::.LogPriorEmpiricalGeometric(k, emp, state$p, K = truncK,
                                            unconditional = TRUE)
    }, 0)
    state$kPrime[ch$trans] <- sample(2:truncK, length(ch$trans),
                                     replace = TRUE, prob = exp(pmf))
    ptr <- StartChain(ch, state)
    for (round in 1:3) {
      MkPrime:::do_move_cpp(ptr$data, ptr$state, 25L, 0L, 0.5, 0.5, 1L, 0)
      for (i in 1:10) {
        MkPrime:::do_move_cpp(ptr$data, ptr$state, 30L, 0L, 1, 0.5, 1L, 0)
      }
    }
    MkPrime:::get_mcmc_state(ptr$state)$p
  }, numeric(1))
  pKs <- ks.test(p, "pbeta", a, b)$p.value
  expect_gt(pKs, 1e-3, label = paste("p: KS p =", signif(pKs, 2)))
})
