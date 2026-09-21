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
  emp <- empiricalNObs  # package data
  p <- 0.7
  # P(k' = 2 | empirical_geometric) — the convolution at the minimum
  log_eg <- MkPrime:::.LogPriorEmpiricalGeometric(2L, emp, p)
  # P(u = 0 | geometric) = p
  log_geom <- log(p)
  expect_lt(log_eg, log_geom)
})


test_that("MkPrimeModel(kPrimePrior='empirical_geometric') is the default", {
  model <- MkPrimeModel()
  expect_identical(model$kPrimePrior, "empirical_geometric")
})


test_that("both geometric arms default to the unconditional (Model A) prior", {
  # A prior is pre-data: kObs is an observation and must not enter it.  Both
  # empirical_geometric and geometric therefore default to unconditional (Model A).
  # The geometric SBC harness draws k' ~ 2 + Geo(p), so the default aligns
  # inference with the validated SBC forward model. Explicit values are honoured.
  expect_identical(MkPrimeModel(kPrimePrior = "empirical_geometric")$priorVariant,
                   "unconditional")
  expect_identical(MkPrimeModel()$priorVariant, "unconditional")  # default arm
  expect_identical(MkPrimeModel(kPrimePrior = "geometric")$priorVariant,
                   "unconditional")
  expect_identical(
    MkPrimeModel(kPrimePrior = "empirical_geometric",
                 priorVariant = "conditional")$priorVariant,
    "conditional")
})


test_that("LogPrior under empirical_geometric prior is finite for valid state", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
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
  pd <- TreeTools::MatrixToPhyDat(mat)
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
  tree <- TreeTools::Preorder(tree)
  # Two transformational chars with different kObs (2 and 3).
  mat <- matrix(c(0, 1, 0, 1,
                  0, 1, 2, 0), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
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
  pd <- TreeTools::MatrixToPhyDat(mat)
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


test_that("empirical_geometric posterior on u beats geometric when true k' > kObs", {
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
  true_tree <- TreeTools::Preorder(true_tree)

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

  pd <- TreeTools::MatrixToPhyDat(sim_mat)
  mkd <- MkPrimeData(pd)
  # Restrict to characters where kObs < kTrue — the interesting ones.
  partialObs <- which(mkd$kObs < kTrue)
  expect_true(length(partialObs) > 5L,
              info = "simulated dataset should contain partial-observation chars")

  run <- function(prior, seed) {
    set.seed(seed)
    RunMkPrime(
      mkd, true_tree,
      model = MkPrimeModel(kPrimePrior = prior, expSteps = sum(true_tree$edge.length)),
      mcmc = MkPrimeMCMC(nIter = 3000L, thin = 10L,
                         maxWarmup = 1500L, minWarmup = 1500L,
                         autoTune = FALSE)
    )
  }
  res_emp <- run("empirical_geometric", 101)
  res_geo <- run("geometric", 101)

  kPrimeCols <- grep("^kPrime_", colnames(res_emp$samples), value = TRUE)
  expect_true(length(kPrimeCols) >= length(partialObs))

  # Per-character posterior median of u = kPrime - kObs
  uMedEmp <- vapply(seq_along(partialObs), function(i) {
    median(res_emp$samples[, kPrimeCols[partialObs[i]]] - mkd$kObs[partialObs[i]])
  }, numeric(1))
  uMedGeo <- vapply(seq_along(partialObs), function(i) {
    median(res_geo$samples[, kPrimeCols[partialObs[i]]] - mkd$kObs[partialObs[i]])
  }, numeric(1))

  # The empirical_geometric prior should keep posterior u above zero more
  # often than the simple geometric prior.  Sum of posterior medians is a
  # robust proxy: empirical should give a larger total inferred u.
  expect_gt(sum(uMedEmp), sum(uMedGeo))
})


test_that("EG-001: per-character truncation normaliser Z_i(p) is correct", {
  # Regression guard for EG-001 (proofs/kprime-priors.md s4.3). LogPrior
  # enforces k'_i >= kObs_i, so the empirical_geometric density must be
  # renormalised by Z_i(p) = sum_{k>=kObs_i} P(k|p) = 1 - sum_{m=2}^{kObs_i-1}
  # P(m|p). No-op at kObs = 2.
  emp <- MkPrimeEmpiricalPrior(body = c(0.6, 0.3, 0.1), tail_decay = 0)
  p <- 0.7
  eg <- function(k, kObs) MkPrime:::.LogPriorEmpiricalGeometric(k, emp, p, kObs)

  # Untruncated convolution pmf (kObs = 2 => Z_i = 1).
  P2 <- exp(eg(2L, 2L)); P3 <- exp(eg(3L, 2L)); P4 <- exp(eg(4L, 2L))
  expect_equal(P2, 0.42,   tolerance = 1e-12)   # 0.6 * 0.7
  expect_equal(P3, 0.336,  tolerance = 1e-12)
  expect_equal(P4, 0.1708, tolerance = 1e-12)

  # kObs = 2 is a strict no-op vs the (default-arg) untruncated form.
  for (k in 2:6) {
    expect_equal(eg(k, 2L), MkPrime:::.LogPriorEmpiricalGeometric(k, emp, p),
                 tolerance = 1e-14)
  }

  # Hand-verified Z_i at kObs = 3 and 4.
  Z3 <- 1 - P2            # 0.58
  Z4 <- 1 - P2 - P3       # 0.244
  expect_equal(eg(3L, 3L), log(P3) - log(Z3), tolerance = 1e-12)
  expect_equal(eg(4L, 3L), log(P4) - log(Z3), tolerance = 1e-12)
  expect_equal(eg(4L, 4L), log(P4) - log(Z4), tolerance = 1e-12)

  # Truncating raises the per-state density (mass redistributed; Z_i <= 1).
  expect_gte(eg(5L, 4L), eg(5L, 2L))

  # The renormalised prior integrates to 1 over its support k >= kObs.
  for (kObs in c(2L, 3L, 4L, 5L)) {
    mass <- sum(vapply(kObs:300L, function(k) exp(eg(k, kObs)), numeric(1)))
    expect_equal(mass, 1, tolerance = 1e-8,
                 info = sprintf("renormalised mass at kObs=%d", kObs))
  }

  # Vectorised per-character kObs.
  expect_equal(
    MkPrime:::.LogPriorEmpiricalGeometric(c(3L, 4L), emp, p, c(3L, 4L)),
    (log(P3) - log(Z3)) + (log(P4) - log(Z4)), tolerance = 1e-12)

  # Boundary: if mass below kObs is ~1 (Z_i -> 0), prior is -Inf, not NaN.
  empLow <- MkPrimeEmpiricalPrior(body = c(0.999, 0.001), tail_decay = 0)
  v <- MkPrime:::.LogPriorEmpiricalGeometric(8L, empLow, 0.999, 8L)
  expect_true(is.finite(v) || v == -Inf)   # never NaN
})


test_that("EG-001 Model A: priorVariant='unconditional' drops the Z_i correction", {
  # Model A (unconditional) places the empirical_geometric prior on the full
  # support k' >= 2 with no per-character Z_i(p) truncation correction; the
  # likelihood enforces the k' >= kObs floor. It must therefore equal the
  # untruncated convolution for every kObs, and differ from the conditional
  # (Model B) variant by exactly + sum_i log Z_i(p).
  emp <- MkPrimeEmpiricalPrior(body = c(0.6, 0.3, 0.1), tail_decay = 0)
  p <- 0.7
  egC <- function(k, kObs)
    MkPrime:::.LogPriorEmpiricalGeometric(k, emp, p, kObs)
  egA <- function(k, kObs)
    MkPrime:::.LogPriorEmpiricalGeometric(k, emp, p, kObs, unconditional = TRUE)

  # Model A == untruncated convolution (kObs = 2 default => Z_i = 1) for any kObs.
  for (kObs in 2:5) {
    for (k in kObs:6) {
      expect_equal(egA(k, kObs),
                   MkPrime:::.LogPriorEmpiricalGeometric(k, emp, p),
                   tolerance = 1e-14,
                   info = sprintf("Model A k=%d kObs=%d", k, kObs))
    }
  }

  # The two variants coincide at kObs = 2 (Z_i = 1 either way).
  for (k in 2:6) expect_equal(egA(k, 2L), egC(k, 2L), tolerance = 1e-14)

  # At kObs > 2 they differ by exactly + log Z_i (Model B subtracts log Z_i).
  P2 <- exp(egC(2L, 2L)); P3 <- exp(egC(3L, 2L))
  Z3 <- 1 - P2; Z4 <- 1 - P2 - P3
  expect_equal(egA(3L, 3L) - egC(3L, 3L), log(Z3), tolerance = 1e-12)
  expect_equal(egA(4L, 4L) - egC(4L, 4L), log(Z4), tolerance = 1e-12)

  # Vectorised: difference == sum of per-character log Z_i.
  expect_equal(
    MkPrime:::.LogPriorEmpiricalGeometric(c(3L, 4L), emp, p, c(3L, 4L),
                                          unconditional = TRUE) -
      MkPrime:::.LogPriorEmpiricalGeometric(c(3L, 4L), emp, p, c(3L, 4L)),
    log(Z3) + log(Z4), tolerance = 1e-12)
})


test_that("EG-001 Model A: R and C++ EG priors agree under priorVariant='unconditional'", {
  tree <- TreeTools::Preorder(
    read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);"))
  # Two transformational chars, kObs = 2 and 3 (the second exercises Z_i).
  mat <- matrix(c(0, 1, 0, 1,
                  0, 1, 2, 0), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- suppressMessages(
    MkPrimeModel(expSteps = 10, kPrimePrior = "empirical_geometric",
                 priorVariant = "unconditional"))
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)

  for (p in c(0.2, 0.5, 0.8)) {
    # empiricalNObs has a 15-entry body, so its tail starts at k = 17: the
    # last two pairs carry the convolution past the end of empLogBody and into
    # the geometric tail, which the body-length cases never reach.
    for (kp in list(c(2L, 3L), c(3L, 5L), c(4L, 8L),
                    c(2L, 25L), c(20L, 30L))) {
      state <- list(
        tree = tree, tree_length = 0.5,
        rel_br_lengths = tree$edge.length / sum(tree$edge.length),
        rate_loss = 1.0, rate_log_sd = 0.2, rate_neo = 1.0,
        kPrime = kp, p = p, log_lik = 0.0, log_prior = 0.0
      )
      lp_r   <- MkPrime:::LogPrior(state, model, mkd)
      lp_cpp <- eval_log_prior_cpp(dataPtr, MkPrime:::.InitMcmcChain(state))
      expect_equal(lp_cpp, lp_r, tolerance = 1e-9,
                   info = sprintf("Model A R vs C++ at p=%.2f, kp=%s",
                                  p, paste(kp, collapse = ",")))
    }
  }
})


test_that("LogPrior reports which character carries a missing k' or kObs", {
  # EG-004: a bare any() on a vector containing NA made `if()` raise
  # "missing value where TRUE/FALSE needed", naming the prior rather than the
  # character whose kObs failed to be ingested.
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 1, 1, 0), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  # Pinned, not defaulted: the frozen value below is sensitive to
  # `treeLengthShape` (3.0 nats) and `rateLogSdRate` (0.39 nats), twelve
  # orders outside its 1e-12 tolerance.
  model <- MkPrimeModel(
    expSteps        = 10,
    kPrimePrior     = "empirical_geometric",
    priorVariant    = "unconditional",
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
