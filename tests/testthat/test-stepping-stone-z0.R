# #267: under sampled_k and the unconditional (Model A) k' prior, no move
# proposes k'_i < kObs_i, so the beta = 0 stone samples the prior restricted to
# k' >= kObs and the stones estimate log Z1 - log Z0. The oracle is the exact
# log marginal likelihood by enumeration over k' and quadrature over p, at the
# fixed tree the stones never leave (only p and k' moves are weighted).

.Z0Tree <- function() {
  Preorder(read.tree(text = "((t1:0.1,t2:0.2):0.05,(t3:0.15,t4:0.1):0.05);"))
}

.Z0Data <- function() {
  mat <- matrix(c(0, 1, 2, 2,  0, 1, 2, 0,  0, 1, 2, 3), 4, 3,
                dimnames = list(paste0("t", 1:4), NULL))
  MkPrimeData(MatrixToPhyDat(mat))
}

.Z0Model <- function(mode = "sampled_k", prior = "geometric") {
  MkPrimeModel(kPrimePrior = prior, kprimeTruncK = 8L,
               priorVariant = "unconditional", kprimeHyperA = 1,
               kprimeHyperB = 1, likelihoodMode = mode)
}

# log of sum_{k'} exp(LL(k')) prod_i P(k'_i | p), integrated against Beta(1, 1)
.Z0ExactLogMl <- function() {
  tree <- .Z0Tree()
  mkd <- .Z0Data()
  model <- MkPrime:::.FinalizeModel(.Z0Model(), tree, mkd)
  state <- MkPrime:::.InitState(tree, mkd, model)
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)
  LogLik <- function(kPrime) {
    state$kPrime <- kPrime
    statePtr <- MkPrime:::.InitMcmcChain(state)
    MkPrime:::fill_partition_cache(dataPtr, statePtr)
    MkPrime:::eval_full_loglik_cpp(dataPtr, statePtr)
  }
  K <- 8L
  kObs <- mkd$kObs
  grid <- as.matrix(expand.grid(lapply(kObs, seq.int, to = K)))
  logLiks <- apply(grid, 1, LogLik)
  logLik0 <- max(logLiks)
  LogPk <- function(k, p) {
    log(p) + (k - 2) * log1p(-p) - log1p(-(1 - p) ^ (K - 1))
  }
  Integrand <- function(p) {
    vapply(p, function(pp) {
      sum(exp(logLiks - logLik0 + rowSums(LogPk(grid, pp))))
    }, 0)
  }
  logLik0 + log(stats::integrate(Integrand, 0, 1, rel.tol = 1e-10)$value)
}

.Z0SteppingStone <- function(mode) {
  weights <- if (mode == "sampled_k") {
    c(mh_logit_p = 0.4, kPrime = 0.2, gibbs_kPrime = 0.3, block_kPrime = 0.1)
  } else {
    c(mh_logit_p = 1)
  }
  set.seed(1L)
  mkp_stepping_stone(.Z0Data(), .Z0Tree(), model = .Z0Model(mode),
                     mcmc = MkPrimeMCMC(moveWeights = weights),
                     nStones = 20L, nIter = 2000L, warmup = 200L,
                     fixTopology = TRUE, verbose = FALSE)
}

test_that("stepping stone matches the exact log ML in both modes (#267)", {
  exact <- .Z0ExactLogMl()
  for (mode in c("sampled_k", "marginal_k")) {
    ss <- .Z0SteppingStone(mode)
    # Before the fix, sampled_k sat -log Z0 = 1.975 above the exact value. The
    # delta-method se misses the stones' shared chain, so allow for seed spread.
    expect_lt(abs(ss$log_marginal - exact), 0.25)
  }
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

test_that(".LogZ0 integrates the Model A tail mass of each k' prior (#267)", {
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

  eg <- .Z0Model(prior = "empirical_geometric")
  emp <- MkPrime:::.EmpiricalNObs()
  expect_equal(
    MkPrime:::.LogZ0(eg, mkd),
    .Z0ByComplement(function(k, p) {
      vapply(k, function(kk) {
        MkPrime:::.LogPriorEmpiricalGeometric(kk, emp, p, unconditional = TRUE)
      }, 0)
    }, kObs),
    tolerance = 1e-8
  )

  logseries <- MkPrimeModel(kPrimePrior = "logseries",
                            priorVariant = "unconditional")
  lsC <- logseries$kprimeLogseriesC
  lsTail <- vapply(kObs, function(k) {
    1 - sum(lsC ^ (2:(k - 1)) / (2:(k - 1))) / (-log1p(-lsC) - lsC)
  }, 0)
  expect_equal(MkPrime:::.LogZ0(logseries, mkd), sum(log(lsTail)),
               tolerance = 1e-12)

  # Z0 = 1 where the prior already lives on k' >= kObs, or no kObs exceeds 2.
  expect_equal(MkPrime:::.LogZ0(.Z0Model("marginal_k"), mkd), 0)
  expect_equal(
    MkPrime:::.LogZ0(MkPrimeModel(kPrimePrior = "geometric",
                                  priorVariant = "conditional"), mkd),
    0
  )
  expect_equal(
    MkPrime:::.LogZ0(geom, list(kObs = c(2L, 2L),
                                type = rep("transformational", 2))),
    0
  )
})
