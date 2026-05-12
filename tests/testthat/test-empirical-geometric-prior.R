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
