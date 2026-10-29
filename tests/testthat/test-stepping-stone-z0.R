# #267: under sampled_k no move proposes k'_i < kObs_i, so the beta = 0 stone
# samples the prior restricted to k' >= kObs and the stones estimate
# log Z1 - log Z0. The oracle is the exact log marginal likelihood by
# enumeration over k' and quadrature over the k' prior's hyperparameters, at
# the fixed tree the stones never leave (only k' and hyperparameter moves are
# weighted).

.Z0Tree <- function() {
  Preorder(read.tree(text = "((t1:0.1,t2:0.2):0.05,(t3:0.15,t4:0.1):0.05);"))
}

.Z0Data <- function() {
  mat <- matrix(c(0, 1, 2, 2,  0, 1, 2, 0,  0, 1, 2, 3), 4, 3,
                dimnames = list(paste0("t", 1:4), NULL))
  MkPrimeData(MatrixToPhyDat(mat))
}

.Z0Model <- function(mode = "sampled_k", prior = "geometric") {
  MkPrimeModel(kPrimePrior = prior, kprimeTruncK = 8L, kprimeHyperA = 1,
               kprimeHyperB = 1, likelihoodMode = mode)
}

.LogSumExp <- function(x) max(x) + log(sum(exp(x - max(x))))

# log pi(k) for k = 2..K.
.Z0LogPrior <- function(prior, hyper, K, emp) {
  logW <- switch(prior,
    geometric = (2:K - 2) * log1p(-hyper$p),
    empirical_geometric = vapply(2:K, function(k) {
      MkPrime:::.LogPriorEmpiricalGeometric(k, emp, hyper$p)
    }, 0),
    beta_geometric = lbeta(hyper$kprime_alpha + 1, hyper$kprime_beta + 2:K - 2)
  )
  logW - .LogSumExp(logW)
}

# log of sum_{k'} exp(LL(k' | h)) prod_i P(k'_i | h), integrated against the
# prior on h: p ~ Beta(1, 1), or alpha, beta ~ Exp(1) for beta_geometric. h
# enters LL only through the joint ascertainment term, which k' does not.
.Z0ExactLogMl <- function(prior = "geometric") {
  tree <- .Z0Tree()
  mkd <- .Z0Data()
  model <- MkPrime:::.FinalizeModel(.Z0Model(prior = prior), tree, mkd)
  state <- MkPrime:::.InitState(tree, mkd, model)
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)
  LogLik <- function(kPrime, hyper = list()) {
    state$kPrime <- kPrime
    state[names(hyper)] <- hyper
    statePtr <- MkPrime:::.InitMcmcChain(state)
    MkPrime:::fill_partition_cache(dataPtr, statePtr)
    MkPrime:::eval_full_loglik_cpp(dataPtr, statePtr)
  }
  K <- 8L
  kObs <- mkd$kObs
  grid <- as.matrix(expand.grid(lapply(kObs, seq.int, to = K)))
  logLiks <- apply(grid, 1, LogLik)
  logLik0 <- max(logLiks)
  base <- LogLik(kObs)
  LogTerm <- function(hyper) {
    lp <- .Z0LogPrior(prior, hyper, K, model$empiricalNObs)
    .LogSumExp(logLiks - logLik0 + rowSums(matrix(lp[grid - 1L], nrow(grid)))) +
      LogLik(kObs, hyper) - base
  }
  if (prior == "beta_geometric") {
    # The trapezium rule in (log alpha, log beta), whose error falls
    # exponentially with the step for so smooth an integrand.
    step <- 0.3
    x <- seq(-14, 4, by = step)
    logTerms <- outer(x, x, Vectorize(function(xa, xb) {
      a <- exp(xa)
      b <- exp(xb)
      LogTerm(list(kprime_alpha = a, kprime_beta = b)) + xa + xb - a - b
    }))
    return(logLik0 + .LogSumExp(logTerms) + 2 * log(step))
  }
  Integrand <- function(p) {
    vapply(p, function(pp) exp(LogTerm(list(p = pp))), 0)
  }
  logLik0 + log(stats::integrate(Integrand, 0, 1, rel.tol = 1e-10)$value)
}

.Z0SteppingStone <- function(mode, prior = "geometric") {
  weights <- if (prior == "beta_geometric") {
    c(slice_kprime_s = 0.2, slice_kprime_r = 0.2, kPrime = 0.2,
      gibbs_kPrime = 0.3, block_kPrime = 0.1)
  } else if (mode == "sampled_k") {
    c(mh_logit_p = 0.4, kPrime = 0.2, gibbs_kPrime = 0.3, block_kPrime = 0.1)
  } else {
    c(mh_logit_p = 1)
  }
  # The oracle needs the tree and branch lengths frozen, which pins summing to
  # one would do but mkp_stepping_stone() refuses (#277).
  local_mocked_bindings(
    .SchedulePins = function(moveWeights, moveTypes, userPins) userPins,
    .package = "MkPrime"
  )
  set.seed(1L)
  mkp_stepping_stone(.Z0Data(), .Z0Tree(), model = .Z0Model(mode, prior),
                     mcmc = MkPrimeMCMC(moveWeights = weights),
                     nStones = 20L, nIter = 2000L, warmup = 200L,
                     fixTopology = TRUE, nRuns = 1L, verbose = FALSE)
}

test_that("stepping stone matches the exact log ML in both modes (#267)", {
  exact <- .Z0ExactLogMl()
  for (mode in c("sampled_k", "marginal_k")) {
    ss <- .Z0SteppingStone(mode)
    # Before the fix, sampled_k sat -log Z0 = 1.975 above the exact value.
    # The tolerance allows for seed spread.
    expect_lt(abs(ss$log_marginal - exact), 0.25)
  }
})

test_that("empirical_geometric stepping stone matches the exact log ML (#377)", {
  # The beta = 0 stone samples the prior capped at K, so log Z0 must integrate
  # the capped tail: integrating it to infinity was 2 nats out (#377).
  exact <- .Z0ExactLogMl("empirical_geometric")
  ss <- .Z0SteppingStone("sampled_k", "empirical_geometric")
  expect_lt(abs(ss$log_marginal - exact), 0.25)
})

test_that("beta_geometric stepping stone matches the exact log ML", {
  skip_under_memcheck()
  exact <- .Z0ExactLogMl("beta_geometric")
  ss <- .Z0SteppingStone("sampled_k", "beta_geometric")
  expect_lt(abs(ss$log_marginal - exact), 0.25)
})

# log Z0 = log integral of prod_i P(k'_i >= kObs_i | p) dBeta(p), from the
# complement of the pmf below kObs rather than the tail-sum forms .LogZ0 uses.
.Z0ByComplement <- function(LogPmf, kObs, a = 1, b = 1) {
  Integrand <- function(p) {
    vapply(p, function(pp) {
      tails <- vapply(kObs, function(k) {
        1 - sum(exp(LogPmf(seq.int(2L, k - 1L), pp)))
      }, 0)
      prod(tails) * stats::dbeta(pp, a, b)
    }, 0)
  }
  log(stats::integrate(Integrand, 0, 1, rel.tol = 1e-10)$value)
}

test_that(".LogZ0 integrates the tail mass of each k' prior (#267)", {
  mkd <- list(kObs = c(3L, 3L, 4L, 2L, 6L),
              type = c(rep("transformational", 4), "neomorphic"))
  kObs <- c(3L, 3L, 4L)

  geom <- .Z0Model()
  expect_equal(MkPrime:::.LogZ0(geom, mkd), -1.975117, tolerance = 1e-6)
  expect_equal(
    MkPrime:::.LogZ0(geom, mkd),
    .Z0ByComplement(function(k, p) {
      log(p) + (k - 2) * log1p(-p) - log1p(-(1 - p) ^ 7)
    }, kObs),
    tolerance = 1e-8
  )

  # The empirical_geometric pmf is renormalised on [2, K] (#377).
  eg <- .Z0Model(prior = "empirical_geometric")
  emp <- MkPrime:::.EmpiricalNObs()
  LogPmf <- function(k, p) {
    vapply(k, function(kk) {
      MkPrime:::.LogPriorEmpiricalGeometric(kk, emp, p)
    }, 0)
  }
  expect_equal(
    MkPrime:::.LogZ0(eg, mkd),
    .Z0ByComplement(function(k, p) {
      LogPmf(k, p) - log(sum(exp(LogPmf(2:8, p))))
    }, kObs),
    tolerance = 1e-8
  )

  logseries <- MkPrimeModel(kPrimePrior = "logseries")
  lsC <- logseries$kprimeLogseriesC
  lsK <- logseries$kprimeTruncK
  lsTail <- vapply(kObs, function(k) {
    1 - sum(lsC ^ (2:(k - 1)) / (2:(k - 1))) / sum(lsC ^ (2:lsK) / (2:lsK))
  }, 0)
  expect_equal(MkPrime:::.LogZ0(logseries, mkd), sum(log(lsTail)),
               tolerance = 1e-12)

  # Alpha and beta carry Exp(1) priors; a grid in their logs, with the
  # tails summed directly.
  bg <- .Z0Model(prior = "beta_geometric")
  step <- 0.05
  x <- seq(-30, 5, by = step)
  ab <- expand.grid(a = exp(x), b = exp(x))
  logW <- vapply(2:8, function(k) lbeta(ab$a + 1, ab$b + k - 2), ab$a)
  top <- apply(logW, 1, max)
  logTotal <- log(rowSums(exp(logW - top)))
  logF <- log(ab$a) + log(ab$b) - ab$a - ab$b
  for (k in kObs) {
    logF <- logF +
      log(rowSums(exp(logW[, k:8 - 1L, drop = FALSE] - top))) - logTotal
  }
  expect_equal(MkPrime:::.LogZ0(bg, mkd),
               .LogSumExp(logF) + 2 * log(step), tolerance = 1e-8)

  # Z0 = 1 under marginal_k, which carries the mass itself, or where no kObs
  # exceeds 2.
  expect_equal(MkPrime:::.LogZ0(.Z0Model("marginal_k"), mkd), 0)
  expect_equal(
    MkPrime:::.LogZ0(geom, list(kObs = c(2L, 2L),
                                type = rep("transformational", 2))),
    0
  )
})

test_that("default logseries subtracts log Z0 (#365)", {
  mkd <- list(kObs = c(3L, 3L, 4L, 2L, 6L),
              type = c(rep("transformational", 4), "neomorphic"))
  default <- MkPrimeModel(kPrimePrior = "logseries")
  expect_lt(MkPrime:::.LogZ0(default, mkd), 0)
  expect_equal(MkPrime:::.LogZ0(default, .Z0Data()), -2.579877,
               tolerance = 1e-6)
})

test_that("default logseries stepping stone matches the exact log ML (#365)", {
  # Exact value by enumerating k'_i in [kObs_i, 30] under the untruncated
  # logseries pmf (c = 0.7); the omitted tail is 3e-6 per character.
  exact <- -13.63322
  local_mocked_bindings(
    .SchedulePins = function(moveWeights, moveTypes, userPins) userPins,
    .package = "MkPrime"
  )
  set.seed(1L)
  ss <- mkp_stepping_stone(
    .Z0Data(), .Z0Tree(), model = MkPrimeModel(kPrimePrior = "logseries"),
    mcmc = MkPrimeMCMC(moveWeights = c(kPrime = 0.4, gibbs_kPrime = 0.6)),
    nStones = 20L, nIter = 2000L, warmup = 200L, fixTopology = TRUE,
    nRuns = 1L, verbose = FALSE
  )
  # Before the fix it sat -log Z0 = 2.58 above the exact value; seeds spread
  # by about 0.01.
  expect_lt(abs(ss$log_marginal - exact), 0.25)
})
