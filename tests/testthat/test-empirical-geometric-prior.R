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


test_that("LogPrior under empirical_geometric prior is finite for valid state", {
  library("ape")
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
  library("ape")
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
  library("ape")
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
    for (kp in list(c(2L, 3L), c(3L, 5L), c(4L, 8L))) {
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


test_that("block_kPrime move is omitted under empirical_geometric", {
  # The block_kPrime move proposes a uniform integer shift in k', which
  # has ~0 % acceptance under the empirical_geometric prior (the prior is
  # heterogeneous in k', so a uniform shift lands every character in a
  # wildly different prior region).  mh_logit_p + joint_p_kprime carry
  # the (p, k') coordination instead.  Verify the move is not scheduled
  # for this prior, but is still scheduled for the plain geometric prior.
  pd <- TreeTools::MatrixToPhyDat(matrix(
    c(0, 1, 2, 0,
      0, 1, 0, 1), 4, 2,
    dimnames = list(paste0("t", 1:4), NULL)
  ))
  mkd <- MkPrimeData(pd)
  nTrans <- sum(mkd$type == "transformational")

  moves_eg <- MkPrime:::.BuildMoves(
    nEdge = 5L, nTrans = nTrans, hasNeo = FALSE,
    mcmc = MkPrimeMCMC(), fixTopology = FALSE,
    kPrimePrior = "empirical_geometric"
  )
  expect_false("block_kPrime" %in% vapply(moves_eg, `[[`, character(1L), "name"))

  moves_g <- MkPrime:::.BuildMoves(
    nEdge = 5L, nTrans = nTrans, hasNeo = FALSE,
    mcmc = MkPrimeMCMC(), fixTopology = FALSE,
    kPrimePrior = "geometric"
  )
  expect_true("block_kPrime" %in% vapply(moves_g, `[[`, character(1L), "name"))
})


test_that("empirical_geometric prior runs short MCMC end-to-end", {
  library("ape")
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
  res <- RunMkPrime(
    mkd, tree,
    model = model,
    mcmc = MkPrimeMCMC(nIter = 200L, thin = 10L,
                       maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE)
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
  library("ape")
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


test_that("empirical_geometric posterior on Sigma u is two-sided bounded near truth", {
  # Two-sided guard against the over-recovery regression diagnosed in
  # dev/plans/2026-05-13-1430-empirical-geometric-prior-overrecovery-diagnosis.md
  # (Σu posterior was ~345% of truth before the Beta(15, 1) hyperprior).
  # Catches both under- and over-shrinkage on a synthetic dataset whose
  # ground-truth Σu is known.
  skip_slow_tests()
  library("ape")
  set.seed(2026)
  nTip <- 8L
  nCharTotal <- 80L
  true_tree <- rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  true_tree$edge.length <- runif(nrow(true_tree$edge), 0.05, 0.25)
  true_tree <- TreeTools::Preorder(true_tree)

  kTrueAll <- sample(c(2L, 3L, 4L, 5L), nCharTotal, replace = TRUE,
                      prob = c(0.40, 0.30, 0.20, 0.10))
  sim_mat <- matrix(NA_integer_, nTip, nCharTotal,
                     dimnames = list(true_tree$tip.label, NULL))
  for (ch in seq_len(nCharTotal)) {
    kTrue <- kTrueAll[ch]
    node_states <- integer(2 * nTip - 1)
    rootIdx <- nTip + 1L
    node_states[rootIdx] <- sample.int(kTrue, 1L) - 1L
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
  kTrueAll <- kTrueAll[variable]

  pd <- TreeTools::MatrixToPhyDat(sim_mat)
  mkd <- MkPrimeData(pd)
  unseen <- sum(kTrueAll - mkd$kObs)
  start_tree <- TreeTools::NJTree(pd, edgeLengths = TRUE)

  set.seed(101)
  suppressWarnings(suppressMessages({
    res <- RunMkPrime(
      mkd, start_tree,
      model = MkPrimeModel(coding = "variable",
                           kPrimePrior = "empirical_geometric",
                           expSteps = sum(true_tree$edge.length)),
      mcmc = MkPrimeMCMC(nIter = 3000L, thin = 20L,
                         maxWarmup = 1500L, minWarmup = 1500L,
                         autoTune = FALSE,
                         nRuns = 1L, nChains = 1L)
    )
  }))
  kpCols <- grep("^kPrime_", colnames(res$samples), value = TRUE)
  kPostMean <- colMeans(res$samples[, kpCols, drop = FALSE])
  kPrimeIdx <- as.integer(sub("kPrime_", "", kpCols))
  ord <- order(kPrimeIdx)
  sigmaUPost <- sum((kPostMean - mkd$kObs)[ord])

  # Tolerance: 50% — wide enough to ride out chain noise on a short run,
  # tight enough to flag a 2x over-recovery (the pre-fix symptom).
  tol <- 0.5
  expect_lt(abs(sigmaUPost - unseen), tol * unseen,
            label = sprintf("posterior Sigma u = %.1f; truth = %d",
                            sigmaUPost, unseen))
})


test_that("empirical_geometric posterior Sigma u is bounded under truth = 0", {
  # Regression guard for the over-recovery pathology on a dataset whose
  # ground-truth Σu is 0 (every character truly binary, kObs = kTrue).
  #
  # Structural note: Fix B (Beta(15, 1) hyperprior) cannot drive Σu to
  # zero on this dataset.  At the prior asymptote p -> 1 the per-character
  # k' = 2 -> 3 prior log-ratio is log(P_emp(3)/P_emp(2)) ≈ -1.139, while
  # the relabel correction lgamma(k'+1) - lgamma(k'-kObs+1) shifts by
  # log(3) ≈ +1.099.  The two nearly cancel and the residual posterior
  # mass on k' > kObs is set by the (weak) likelihood signal on short
  # branches.  With the current Beta(15, 1) default the chain converges
  # at p ≈ 0.95 and Σu posterior settles around 30 on this 60-character
  # dataset; even Beta(200, 1) only pushes it to ~28.  Driving Σu below
  # ~5 would require Fix A on top of Fix B for characters with kObs > 2,
  # plus reconsideration of the relabel correction.  See
  # dev/plans/2026-05-13-1430-empirical-geometric-prior-overrecovery-diagnosis.md.
  #
  # What this test guards against: a revert of Fix B (Σu jumps back to
  # ~42 under Beta(1, 1) on this seed) or any future change that breaks
  # the partial counteraction Fix B provides.
  skip_slow_tests()
  library("ape")
  set.seed(7)
  nTip <- 7L
  nChar <- 60L
  tree <- rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  tree$edge.length <- runif(nrow(tree$edge), 0.05, 0.15)
  tree <- TreeTools::Preorder(tree)

  # Binary characters: simulate JC(k=2) on short branches and discard any
  # that fail to realise both states (we want kObs = 2 = kTrue for every
  # char, so Σu_true = 0).
  kTrue <- 2L
  buildBinary <- function() {
    node_states <- integer(2 * nTip - 1)
    rootIdx <- nTip + 1L
    node_states[rootIdx] <- sample.int(kTrue, 1L) - 1L
    edges <- tree$edge
    el    <- tree$edge.length
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
    node_states[seq_len(nTip)]
  }
  sim_mat <- matrix(NA_integer_, nTip, 0L,
                     dimnames = list(tree$tip.label, NULL))
  attempts <- 0L
  while (ncol(sim_mat) < nChar && attempts < 10L * nChar) {
    attempts <- attempts + 1L
    col <- buildBinary()
    if (length(unique(col)) == kTrue) {
      sim_mat <- cbind(sim_mat, col)
    }
  }
  expect_gte(ncol(sim_mat), nChar)
  sim_mat <- sim_mat[, seq_len(nChar)]
  colnames(sim_mat) <- NULL

  pd <- TreeTools::MatrixToPhyDat(sim_mat)
  mkd <- MkPrimeData(pd)
  expect_true(all(mkd$kObs[mkd$type == "transformational"] == 2L))

  set.seed(202)
  suppressWarnings(suppressMessages({
    res <- RunMkPrime(
      mkd, tree,
      model = MkPrimeModel(coding = "variable",
                           kPrimePrior = "empirical_geometric",
                           expSteps = sum(tree$edge.length)),
      mcmc = MkPrimeMCMC(nIter = 8000L, thin = 20L,
                         maxWarmup = 4000L, minWarmup = 4000L,
                         autoTune = FALSE,
                         nRuns = 1L, nChains = 1L)
    )
  }))
  kpCols <- grep("^kPrime_", colnames(res$samples), value = TRUE)
  kPostMean <- colMeans(res$samples[, kpCols, drop = FALSE])
  kPrimeIdx <- as.integer(sub("kPrime_", "", kpCols))
  ord <- order(kPrimeIdx)
  sigmaUPost <- sum((kPostMean - mkd$kObs)[ord])
  pPost <- mean(res$samples[, "p"])

  # Threshold = 40 catches a revert to Beta(1, 1) (≈42 on this seed) with
  # comfortable margin above the post-fix value (≈32).
  expect_lt(sigmaUPost, 40,
            label = sprintf("posterior Sigma u = %.2f (truth = 0)",
                            sigmaUPost))
  # p posterior should sit near the prior mean (≈0.94) — chain converging
  # to prior on weak-likelihood data is itself a useful invariant.
  expect_gt(pPost, 0.85,
            label = sprintf("posterior p mean = %.3f (prior E[p] = 0.9375)",
                            pPost))
})


test_that("empirical_geometric prior at k' = kObs matches analytic convolution (R == C++)", {
  # Direct parity check at the boundary case k' = kObs.  At p = 1 the
  # geometric on N_unobs collapses to a point mass at 0, so
  #   lim_{p -> 1} π(k' = kObs | p) = P_emp(kObs)
  # For p < 1 the convolution sum has multiple non-zero terms; we test
  # against an explicit closed-form expansion for kObs in {2, 3, 4} and
  # check that:
  #   (a) the R helper .LogPriorEmpiricalGeometric matches the closed form,
  #   (b) the full R LogPrior matches the C++ XPtr engine, and
  #   (c) the convolution converges to P_emp(kObs) as p -> 1.
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  tree <- TreeTools::Preorder(tree)
  # Three transformational chars with kObs in {2, 3, 4}.
  mat <- matrix(c(0, 1, 0, 1,
                  0, 1, 2, 0,
                  0, 1, 2, 3), 4, 3,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(expSteps = 10, kPrimePrior = "empirical_geometric")
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  emp <- model$empiricalNObs
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)

  pEmp <- as.numeric(emp$body[seq_len(3L)])  # P_emp(2), P_emp(3), P_emp(4)

  # Explicit convolution sums for kObs in {2, 3, 4}, summing j from 2 to kObs:
  closedForm <- function(p) {
    q <- 1 - p
    log(pEmp[1] * p) +                                           # kObs = 2
      log(pEmp[1] * p * q + pEmp[2] * p) +                       # kObs = 3
      log(pEmp[1] * p * q^2 + pEmp[2] * p * q + pEmp[3] * p)     # kObs = 4
  }

  for (p in c(0.05, 0.5, 0.9, 0.99)) {
    state <- list(
      tree = tree,
      tree_length = 0.5,
      rel_br_lengths = tree$edge.length / sum(tree$edge.length),
      rate_loss = 1.0, rate_log_sd = 0.2, rate_neo = 1.0,
      kPrime = as.integer(mkd$kObs),  # k' = kObs (boundary)
      p = p,
      log_lik = 0.0, log_prior = 0.0
    )
    lp_r <- MkPrime:::LogPrior(state, model, mkd)
    statePtr <- MkPrime:::.InitMcmcChain(state)
    lp_cpp <- eval_log_prior_cpp(dataPtr, statePtr)
    expect_equal(lp_cpp, lp_r, tolerance = 1e-10,
                 info = sprintf("R vs C++ at p = %.3f, k' = kObs", p))

    lp_helper <- MkPrime:::.LogPriorEmpiricalGeometric(state$kPrime, emp, p)
    expect_equal(lp_helper, closedForm(p), tolerance = 1e-12,
                 info = sprintf(".LogPriorEmpiricalGeometric at p = %.3f", p))
  }

  # Asymptotic check: as p -> 1, the convolution reduces to log(P_emp(kObs))
  # summed over characters (the geometric weight on j = kObs goes to 1, all
  # other j -> 0).
  ptest <- 1 - 1e-12
  lp_high <- MkPrime:::.LogPriorEmpiricalGeometric(
    as.integer(mkd$kObs), emp, ptest)
  expect_equal(lp_high, sum(log(pEmp)), tolerance = 1e-8)
})
