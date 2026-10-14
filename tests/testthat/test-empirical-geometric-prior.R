# Tests for the empirical_geometric kPrime prior (convolution of empirical
# pmf on N_obs with Geometric(p) on N_unobs).

test_that("MkPrimeEmpiricalPrior constructor normalises body + tail to 1", {
  emp <- MkPrimeEmpiricalPrior(body = c(0.5, 0.3, 0.2), tail_decay = 0)
  expect_s3_class(emp, "MkPrimeEmpiricalPrior")
  expect_equal(sum(emp$body), 1, tolerance = 1e-12)
  expect_equal(emp$tail_start_p, 0)

  emp2 <- MkPrimeEmpiricalPrior(body = c(0.6, 0.3), tail_decay = 0.4)
  totalMass <- sum(emp2$body) + emp2$tail_start_p / (1 - emp2$tail_decay)
  expect_equal(totalMass, 1, tolerance = 1e-12)
})


test_that("MkPrimeEmpiricalPrior rejects invalid inputs", {
  expect_error(MkPrimeEmpiricalPrior(body = c(-0.1, 0.5)))
  expect_error(MkPrimeEmpiricalPrior(body = c(0.5, 0.5), tail_decay = 1.5))
  expect_error(MkPrimeEmpiricalPrior(body = c(0.5, 0.3), tail_decay = 0))
})


test_that("MkPrimeEmpiricalPrior rejects body with zero total mass", {
  # body = c(0, 0) with tail_decay = 0 sums to 0, which is already caught by
  # the "must sum to 1" check -- but with tail_decay > 0 the old code divided
  # by a zero total (anchor = 0, tailMass = 0, total = 0) and returned an
  # all-NaN pmf silently, poisoning every later log-prior.
  expect_error(MkPrimeEmpiricalPrior(body = c(0, 0), tail_decay = 0.4),
               "positive total mass")
  expect_error(MkPrimeEmpiricalPrior(body = c(0, 0), tail_decay = 0),
               "positive total mass")
})


test_that(".LogPriorEmpiricalGeometric matches hand-computed convolution", {
  emp <- MkPrimeEmpiricalPrior(body = c(0.6, 0.3, 0.1), tail_decay = 0)
  p <- 0.7
  # P(k' = m) = sum_{j=2..m} P_emp(j) * p * (1-p)^(m - j)
  # m=2: 0.6 * 0.7
  # m=3: 0.6 * 0.21 + 0.3 * 0.7
  # m=4: 0.6 * 0.063 + 0.3 * 0.21 + 0.1 * 0.7
  expected2 <- log(0.6 * 0.7)
  expected3 <- log(0.6 * 0.21 + 0.3 * 0.7)
  expected4 <- log(0.6 * 0.063 + 0.3 * 0.21 + 0.1 * 0.7)

  expect_equal(MkPrime:::.LogPriorEmpiricalGeometric(2L, emp, p),
               expected2, tolerance = 1e-12)
  expect_equal(MkPrime:::.LogPriorEmpiricalGeometric(3L, emp, p),
               expected3, tolerance = 1e-12)
  expect_equal(MkPrime:::.LogPriorEmpiricalGeometric(4L, emp, p),
               expected4, tolerance = 1e-12)
  expect_equal(MkPrime:::.LogPriorEmpiricalGeometric(c(2L, 3L, 4L), emp, p),
               expected2 + expected3 + expected4, tolerance = 1e-12)
})


test_that(".LogPriorEmpiricalGeometric returns -Inf for invalid inputs", {
  emp <- MkPrimeEmpiricalPrior(body = c(0.6, 0.4), tail_decay = 0)
  expect_equal(MkPrime:::.LogPriorEmpiricalGeometric(2L, emp, 0), -Inf)
  expect_equal(MkPrime:::.LogPriorEmpiricalGeometric(2L, emp, 1), -Inf)
  expect_equal(MkPrime:::.LogPriorEmpiricalGeometric(1L, emp, 0.5), -Inf)
})


test_that("convolution places less prior mass at u = 0 than pure geometric", {
  # The motivation for this prior: counter the tendency of the simple
  # geometric prior to concentrate mass at u = 0 (no unobserved states).
  # `P(k' = 2) < p` alone reduces to `P_emp(2) < 1`, true of any non-degenerate
  # pmf, so pin the whole distribution instead: the convolution is the
  # geometric shifted by N_obs, which moves E[k'] up by E[N_obs] - 2.
  emp <- empiricalNObs  # package data
  p <- 0.7
  ks <- 2:400

  pEmp <- exp(MkPrime:::.LogPemp(max(ks), emp))
  expect_equal(sum(pEmp), 1, tolerance = 1e-9)

  pEG <- exp(vapply(ks, MkPrime:::.LogPriorEmpiricalGeometric,
                    numeric(1), emp = emp, p = p))
  expect_equal(sum(pEG), 1, tolerance = 1e-9)

  # At the minimum the convolution has the single term P_emp(2) * p.
  expect_equal(pEG[[1L]], pEmp[[1L]] * p, tolerance = 1e-12)
  expect_lt(pEG[[1L]], p)

  pGeo <- p * (1 - p) ^ (ks - 2)
  expect_equal(sum(ks * pEG) - sum(ks * pGeo),
               sum(ks * pEmp) - 2,
               tolerance = 1e-8)
})


test_that("MkPrimeEmpiricalPrior joins body to tail without shifting mass", {
  # The constructor's own normalisation is what `sum(body) + tail mass == 1`
  # re-derives, so it cannot see a mis-placed junction. Read the pmf back
  # through .LogPemp() instead: the tail must start one past `body` and
  # continue it at rate `tail_decay`.
  body <- c(0.6, 0.3)
  q <- 0.4
  emp <- MkPrimeEmpiricalPrior(body = body, tail_decay = q)
  expect_equal(emp$tail_start_k, length(body) + 2L)

  pmf <- exp(MkPrime:::.LogPemp(300L, emp))
  expect_equal(sum(pmf), 1, tolerance = 1e-12)
  # Body entries keep their relative sizes; the first tail entry continues
  # them, so every consecutive ratio from k = 3 on is `tail_decay`.
  expect_equal(pmf[[2L]] / pmf[[1L]], body[[2L]] / body[[1L]],
               tolerance = 1e-12)
  expect_equal(pmf[3:20] / pmf[2:19], rep(q, 18L), tolerance = 1e-12)

  # Passing the one admissible `tail_start_k` explicitly must give the same
  # object as letting it default; the rejections live in the test named
  # "MkPrimeEmpiricalPrior rejects a tail detached from the body".
  expect_identical(
    MkPrimeEmpiricalPrior(body = body, tail_decay = q,
                          tail_start_k = length(body) + 2L),
    emp
  )
})


test_that("MkPrimeModel(kPrimePrior='empirical_geometric') is the default", {
  model <- MkPrimeModel()
  expect_identical(model$kPrimePrior, "empirical_geometric")
})


test_that("LogPrior under empirical_geometric prior is finite for valid state", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(expSteps = 10, kPrimePrior = "empirical_geometric")

  state <- list(
    tree_length = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss = 1.5,
    rate_log_sd = 0.3,
    kPrime = 3L,
    p = 0.4
  )
  expect_true(is.finite(MkPrime:::LogPrior(state, model, mkd)))

  # boundary: k' < kObs returns -Inf
  state$kPrime <- 1L
  expect_equal(MkPrime:::LogPrior(state, model, mkd), -Inf)
})


test_that("user-supplied empiricalNObs overrides package default", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  customEmp <- MkPrimeEmpiricalPrior(body = c(0.99, 0.01), tail_decay = 0)
  model <- MkPrimeModel(expSteps = 10, kPrimePrior = "empirical_geometric",
                         empiricalNObs = customEmp)
  expect_identical(model$empiricalNObs, customEmp)

  state <- list(
    tree_length = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss = 1.0,
    rate_log_sd = 0.0,
    kPrime = 2L,
    p = 0.5
  )
  lp_custom <- MkPrime:::LogPrior(state, model, mkd)
  # With custom emp body c(0.99, 0.01) and p = 0.5:
  # P(k' = 2) = 0.99 * 0.5 = 0.495
  expected_kp <- log(0.99 * 0.5)
  # Strip non-k' terms by computing with same state under custom emp
  # Just check finite and that the k' contribution differs from default
  expect_true(is.finite(lp_custom))
  model_default <- MkPrimeModel(expSteps = 10, kPrimePrior = "empirical_geometric")
  lp_default <- MkPrime:::LogPrior(state, model_default, mkd)
  expect_false(isTRUE(all.equal(lp_custom, lp_default)))
})


test_that("convolution prior remains positive in the tail (no hard cutoff)", {
  # Construct a small empirical with a known geometric tail; verify
  # the prior remains > 0 (log > -Inf) for k' values far beyond the body.
  emp <- MkPrimeEmpiricalPrior(body = c(0.7, 0.2), tail_decay = 0.5)
  for (k in c(10L, 30L, 60L)) {
    expect_true(is.finite(MkPrime:::.LogPriorEmpiricalGeometric(k, emp, 0.5)))
  }
})


test_that("R and C++ log priors agree numerically under empirical_geometric", {
  set.seed(42)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  tree <- Preorder(tree)
  # Two transformational chars with different kObs (2 and 3).
  mat <- matrix(c(0, 1, 0, 1,
                  0, 1, 2, 0), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(expSteps = 10, kPrimePrior = "empirical_geometric")
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)

  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)

  for (p in c(0.2, 0.5, 0.8)) {
    # empiricalNObs has a 15-entry body, so its tail starts at k = 17: the
    # last two pairs carry the convolution past the end of empLogBody and into
    # the geometric tail, which the body-length cases never reach.
    for (kp in list(c(2L, 3L), c(3L, 5L), c(4L, 8L),
                    c(2L, 25L), c(20L, 30L))) {
      state <- list(
        tree = tree,
        tree_length = 0.5,
        rel_br_lengths = tree$edge.length / sum(tree$edge.length),
        rate_loss = 1.0, rate_log_sd = 0.2,
        rate_neo = 1.0,
        kPrime = kp, p = p,
        log_lik = 0.0, log_prior = 0.0
      )
      lp_r <- MkPrime:::LogPrior(state, model, mkd)

      statePtr <- MkPrime:::.InitMcmcChain(state)
      lp_cpp <- eval_log_prior_cpp(dataPtr, statePtr)

      expect_equal(lp_cpp, lp_r, tolerance = 1e-9,
                   info = sprintf("R vs C++ at p=%.2f, kp=%s",
                                  p, paste(kp, collapse = ",")))
    }
  }
})


test_that("empirical_geometric prior runs short MCMC end-to-end", {
  set.seed(11)
  tree <- rtree(5, tip.label = paste0("t", 1:5))
  tree$edge.length <- runif(nrow(tree$edge), 0.05, 0.25)
  mat <- matrix(sample(0:2, 5 * 12, replace = TRUE), 5, 12,
                dimnames = list(tree$tip.label, NULL))
  variable <- apply(mat, 2, function(x) length(unique(x)) > 1)
  mat <- mat[, variable, drop = FALSE]
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(kPrimePrior = "empirical_geometric")
  res <- allow_warning(
    RunMkPrime(
      mkd, tree,
      model = model,
      mcmc = MkPrimeMCMC(nIter = 200L, thin = 10L,
                         maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE)
    ),
    "without stabilisation"
  )
  expect_true("p" %in% colnames(res$samples))
  p_samples <- res$samples[, "p"]
  expect_true(all(p_samples > 0 & p_samples < 1))
  # Empirical_geometric must schedule `mh_logit_p` (the logit-scale MH on p)
  # rather than the broken multiplicative `mh_p` move.
  expect_true("mh_logit_p" %in% names(res$acceptance))
  expect_false("mh_p" %in% names(res$acceptance))
  # Acceptance should be non-trivial even on this tiny dataset.  We don't
  # assert a tight lower bound because the chain is only 200 iters, but a
  # totally stuck move would show ~0%.
  expect_gt(res$acceptance[["mh_logit_p"]], 0.01)
})


test_that("empirical_geometric keeps posterior u off zero when true k' > kObs", {
  skip_slow_tests()
  # Construct a scenario where the true number of states (k' = 5) exceeds
  # the typically observed count: with only 6 tips on a short tree, JC(5)
  # rarely realises all 5 states in a single character.  The empirical
  # prior should refuse to collapse u to 0.
  set.seed(2026)
  nTip <- 6
  nChar <- 80
  kTrue <- 5L
  true_tree <- rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  true_tree$edge.length <- runif(nrow(true_tree$edge), 0.05, 0.2)
  true_tree <- Preorder(true_tree)

  sim_mat <- matrix(NA_integer_, nTip, nChar,
                     dimnames = list(true_tree$tip.label, NULL))
  for (ch in seq_len(nChar)) {
    node_states <- integer(2 * nTip - 1)
    root_idx <- nTip + 1L
    node_states[root_idx] <- sample.int(kTrue, 1L) - 1L
    edges <- true_tree$edge
    el    <- true_tree$edge.length
    for (e in rev(seq_len(nrow(edges)))) {
      pa <- edges[e, 1]; ch2 <- edges[e, 2]; t <- el[e]
      pSame <- 1 / kTrue + (1 - 1 / kTrue) * exp(-kTrue * t)
      if (runif(1) < pSame) {
        node_states[ch2] <- node_states[pa]
      } else {
        node_states[ch2] <- sample(setdiff(seq.int(0, kTrue - 1L),
                                            node_states[pa]), 1L)
      }
    }
    sim_mat[, ch] <- node_states[seq_len(nTip)]
  }
  variable <- apply(sim_mat, 2, function(x) length(unique(x)) > 1)
  sim_mat <- sim_mat[, variable, drop = FALSE]

  pd <- MatrixToPhyDat(sim_mat)
  mkd <- MkPrimeData(pd)
  # Restrict to characters where kObs < kTrue — the interesting ones.
  partialObs <- which(mkd$kObs < kTrue)
  expect_true(length(partialObs) > 5L,
              info = "simulated dataset should contain partial-observation chars")

  run <- function(prior, seed) {
    set.seed(seed)
    allow_warning(
      RunMkPrime(
        mkd, true_tree,
        model = MkPrimeModel(kPrimePrior = prior,
                             expSteps = sum(true_tree$edge.length)),
        mcmc = MkPrimeMCMC(nIter = 3000L, thin = 10L,
                           maxWarmup = 1500L, minWarmup = 1500L,
                           autoTune = FALSE)
      ),
      "without stabilisation"
    )
  }

  # `kPrime_<i>` carries the GLOBAL character index, so look the column up by
  # name: positional indexing into the transformational-only column vector
  # agrees only while every character is transformational.
  uMedian <- function(res) {
    vapply(partialObs, function(i) {
      median(res$samples[, paste0("kPrime_", i)] - mkd$kObs[[i]])
    }, numeric(1))
  }

  resEmp <- run("empirical_geometric", 101L)
  resGeo <- run("geometric", 101L)
  uMedEmp <- uMedian(resEmp)
  uMedGeo <- uMedian(resGeo)

  # Every draw stays in support, and the prior's stated purpose holds: it does
  # not let u collapse to 0 on a character whose states are under-observed.
  for (i in partialObs) {
    expect_true(all(resEmp$samples[, paste0("kPrime_", i)] >= mkd$kObs[[i]]))
  }
  expect_gt(mean(uMedEmp > 0), 0.9)

  # The choice of prior must reach the sampler: the two arms disagree on some
  # character by at least a whole state.
  expect_gte(max(abs(uMedEmp - uMedGeo)), 1)
})


test_that("the capped empirical_geometric prior normalises over [2, K]", {
  # The prior is pre-data, so kObs never enters it: the likelihood enforces
  # k' >= kObs.
  emp <- MkPrimeEmpiricalPrior(body = c(0.6, 0.3, 0.1), tail_decay = 0)
  p <- 0.7
  eg <- function(k, K = Inf) MkPrime:::.LogPriorEmpiricalGeometric(k, emp, p, K)

  P <- exp(vapply(2:4, eg, 0))
  expect_equal(P, c(0.42, 0.336, 0.1708), tolerance = 1e-12)
  expect_equal(eg(2:4, 4L), sum(log(P / sum(P))), tolerance = 1e-12)
  for (K in c(3L, 5L, 30L)) {
    expect_equal(sum(exp(vapply(2:K, eg, 0, K = K))), 1, tolerance = 1e-12,
                 info = K)
  }
  expect_equal(eg(6L, 5L), -Inf)
})


test_that("LogPrior reports which character carries a missing k' or kObs", {
  # EG-004: a bare any() on a vector containing NA made `if()` raise
  # "missing value where TRUE/FALSE needed", naming the prior rather than the
  # character whose kObs failed to be ingested.
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 1, 1, 0), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  # Pinned, not defaulted: the frozen value below is sensitive to
  # `treeLengthShape` (3.0 nats) and `rateLogSdRate` (0.39 nats), twelve
  # orders outside its 1e-12 tolerance.
  model <- MkPrimeModel(
    expSteps        = 10,
    kPrimePrior     = "empirical_geometric",
    treeLengthShape = 2,
    rateLossMeanlog = 0,
    rateLossSdlog   = 2,
    rateLogSdShape  = 1,
    rateLogSdRate   = 1,
    kprimeTruncK    = 200L
  )

  state <- list(
    tree_length = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss = 1.5,
    rate_log_sd = 0.3,
    kPrime = c(3L, 3L),
    p = 0.4
  )
  # Value measured on origin/main before the guard existed: it must be a
  # strict no-op on NA-free input, not merely finite.
  expect_equal(MkPrime:::LogPrior(state, model, mkd), -2.32046728147438,
               tolerance = 1e-12)

  state$kPrime <- c(3L, NA_integer_)
  expect_error(MkPrime:::LogPrior(state, model, mkd), "character 2")

  state$kPrime <- c(3L, 3L)
  mkd$kObs[[1]] <- NA_integer_
  expect_error(MkPrime:::LogPrior(state, model, mkd), "character 1")

  state$kPrime <- c(3L, NA_integer_)
  expect_error(MkPrime:::LogPrior(state, model, mkd), "characters 1 and 2")
})


test_that("the packaged empiricalNObs satisfies the tail-join invariant", {
  expect_identical(as.integer(empiricalNObs$tail_start_k),
                   length(empiricalNObs$body) + 2L)
})


test_that("MkPrimeEmpiricalPrior rejects a tail detached from the body", {
  # EG-005: the tail's mass is anchored on body[nBody], so a tail_start_k
  # beyond nBody + 2 translates that mass outward instead of rescaling it.
  # The pmf still sums to 1 and no evaluation errors, so the reshaped prior
  # would reach the posterior silently.
  expect_error(
    MkPrimeEmpiricalPrior(body = c(0.6, 0.3), tail_decay = 0.4,
                          tail_start_k = 10L),
    "tail_start_k"
  )
  # The value one past the body is the only admissible one.
  expect_s3_class(
    MkPrimeEmpiricalPrior(body = c(0.6, 0.3), tail_decay = 0.4,
                          tail_start_k = 4L),
    "MkPrimeEmpiricalPrior"
  )
})


test_that("prepare_mcmc_data derives the empirical body length from the body", {
  # EG-006: empBodyLastK duplicated 2 + length(empLogBody) across the Rcpp
  # boundary but was never read; a second source of truth that could silently
  # disagree with the first.
  expect_false("empBodyLastK" %in%
                 names(formals(MkPrime:::prepare_mcmc_data)))
})


test_that("R and C++ share one body/tail precedence rule for P_emp", {
  tree <- Preorder(
    read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);"))
  mat <- matrix(c(0, 1, 0, 1,
                  0, 1, 2, 0), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  mkd <- MkPrimeData(MatrixToPhyDat(mat))

  wellFormed <- MkPrimeEmpiricalPrior(body = c(0.6, 0.3, 0.05, 0.05),
                                      tail_decay = 0.4)
  # The constructor refuses a tail that overlaps the body, so build one by
  # hand: both sides must still read the pmf the same way.
  overlapping <- wellFormed
  overlapping$tail_start_k <- 3L

  for (emp in list(wellFormed, overlapping)) {
    model <- MkPrimeModel(expSteps = 10, kPrimePrior = "empirical_geometric",
                          empiricalNObs = emp)
    model <- MkPrime:::.FinalizeModel(model, tree, mkd)
    dataPtr <- MkPrime:::.InitMcmcData(mkd, model)
    for (kp in list(c(3L, 4L), c(5L, 3L), c(2L, 9L))) {
      state <- list(
        tree = tree, tree_length = 0.5,
        rel_br_lengths = tree$edge.length / sum(tree$edge.length),
        rate_loss = 1.0, rate_log_sd = 0.2, rate_neo = 1.0,
        kPrime = kp, p = 0.5, log_lik = 0.0, log_prior = 0.0
      )
      expect_equal(
        eval_log_prior_cpp(dataPtr, MkPrime:::.InitMcmcChain(state)),
        MkPrime:::LogPrior(state, model, mkd),
        tolerance = 1e-10,
        info = sprintf("tail_start_k = %d, kp = %s", emp$tail_start_k,
                       paste(kp, collapse = ","))
      )
    }
  }
})


test_that("P_emp normalises over k >= 2", {
  for (emp in list(MkPrimeEmpiricalPrior(body = c(0.6, 0.3, 0.05, 0.05),
                                         tail_decay = 0.4),
                   MkPrimeEmpiricalPrior(body = rep(0.25, 4)),
                   MkPrime:::.EmpiricalNObs())) {
    expect_equal(sum(exp(MkPrime:::.LogPemp(2000L, emp))), 1,
                 tolerance = 1e-12)
  }
})


test_that("the empirical_geometric prior is exact far into its tail (#259)", {
  emp <- MkPrime::empiricalNObs
  K <- 200L
  logP <- function(k, p) MkPrime:::.LogPriorEmpiricalGeometric(k, emp, p)
  # The untruncated convolution pmf, summed over [2, K].
  ExactLogZ <- function(p) {
    lp <- vapply(2:K, logP, double(1), p = p)
    mx <- max(lp)
    mx + log(sum(exp(lp - mx)))
  }
  kPrime <- c(50L, 60L)

  tree <- Preorder(read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);"))
  mat <- matrix(c(0, 1, 0, 1, 0, 1, 2, 0), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  mkd <- MkPrimeData(MatrixToPhyDat(mat))
  model <- MkPrimeModel(expSteps = 10, kPrimePrior = "empirical_geometric",
                        kprimeTruncK = K)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)

  for (p in c(0.5, 0.9, 0.99)) {
    exact <- sum(vapply(kPrime, logP, double(1), p = p)) - 2 * ExactLogZ(p)
    lpR <- MkPrime:::.LogPriorEmpiricalGeometric(kPrime, emp, p, K)
    expect_equal(lpR, exact, tolerance = 1e-12, info = sprintf("R, p = %g", p))

    state <- list(
      tree = tree, tree_length = 0.5,
      rel_br_lengths = tree$edge.length / sum(tree$edge.length),
      rate_loss = 1, rate_log_sd = 0.2, rate_neo = 1,
      kPrime = kPrime, p = p, log_lik = 0, log_prior = 0
    )
    expect_equal(
      eval_log_prior_cpp(dataPtr, MkPrime:::.InitMcmcChain(state)),
      MkPrime:::LogPrior(state, model, mkd),
      tolerance = 1e-12, info = sprintf("C++ vs R, p = %g", p)
    )
  }
})
