# The ascertainment correction conditions each character on the event that
# selected it: variation among its OBSERVED tips. Its ?/- tips must be
# marginalised, not pinned to the constant state (#213).

.kMissTips <- c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE)

.MissTree <- function(seed = 1L, nTip = 6L) {
  set.seed(seed)
  tree <- Preorder(ape::rtree(nTip, rooted = FALSE))
  tree$edge.length <- runif(nrow(tree$edge), 0.05, 0.8)
  tree
}

# P(constant [or parsimony-uninformative] among the observed tips), by
# enumerating every completion of the missing tips: independent of the
# masked kernels.
.BruteAscProb <- function(tree, missing, k, neo, rateLoss, rates,
                          informative) {
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  nTip <- NTip(tree)
  obs <- which(!missing)
  mis <- which(missing)
  SiteLik <- function(pattern) {
    m <- matrix(as.integer(pattern), nTip, 1)
    exp(if (neo) {
      MkPrime:::pruning_mkn_acrv(parent, child, tree$edge.length, m, rateLoss,
                                 MkPrime:::mkn_stationary_freqs(rateLoss),
                                 rates)
    } else {
      MkPrime:::pruning_jc_acrv(parent, child, tree$edge.length, m, k,
                                rep(1 / k, k), rates)
    })
  }
  # Uninformative: fewer than two states each occur twice or more.
  obsPatterns <- as.matrix(expand.grid(rep(list(seq_len(k) - 1L),
                                           length(obs))))
  keep <- apply(obsPatterns, 1, function(x) {
    counts <- tabulate(x + 1L, k)
    if (informative) sum(counts >= 2L) < 2L else sum(counts > 0L) == 1L
  })
  patterns <- lapply(which(keep), function(i) {
    p <- integer(nTip)
    p[obs] <- obsPatterns[i, ]
    p
  })
  fills <- as.matrix(expand.grid(rep(list(seq_len(k) - 1L), length(mis))))
  # Return:
  sum(vapply(patterns, function(p) {
    sum(apply(fills, 1, function(f) {
      p[mis] <- f
      SiteLik(p)
    }))
  }, double(1)))
}

test_that("masked ascertainment probability marginalises the missing tips", {
  tree <- .MissTree()
  rates <- c(0.3, 0.8, 1.9)
  for (cfg in list(list(k = 3L, neo = FALSE, r = 1),
                   list(k = 5L, neo = FALSE, r = 1),
                   list(k = 2L, neo = TRUE, r = 0.3),
                   list(k = 2L, neo = TRUE, r = 4))) {
    for (informative in c(FALSE, TRUE)) {
      expect_equal(
        MkPrime:::asc_site_prob_missing(
          tree$edge[, 1], tree$edge[, 2], tree$edge.length, NTip(tree),
          cfg$k, cfg$neo, cfg$r, rates, .kMissTips, informative),
        .BruteAscProb(tree, .kMissTips, cfg$k, cfg$neo, cfg$r, rates,
                      informative),
        tolerance = 1e-12, info = paste(cfg, informative)
      )
    }
  }
})

test_that("masked ascertainment reduces to the all-observed probability", {
  tree <- .MissTree(2L)
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  el <- tree$edge.length
  rates <- c(0.5, 1.5)
  none <- logical(NTip(tree))
  expect_identical(
    MkPrime:::asc_site_prob_missing(parent, child, el, 6L, 4L, FALSE, 1,
                                    rates, none, FALSE),
    MkPrime:::constant_site_prob_jc(parent, child, el, 6L, 4L, rep(0.25, 4),
                                    rates))
  rf <- MkPrime:::mkn_stationary_freqs(0.3)
  expect_identical(
    MkPrime:::asc_site_prob_missing(parent, child, el, 6L, 2L, TRUE, 0.3,
                                    rates, none, TRUE),
    MkPrime:::constant_site_prob_mkn(parent, child, el, 6L, 0.3, rf, rates) +
      MkPrime:::singleton_site_prob_mkn(parent, child, el, 6L, 0.3, rf, rates))
})

# Mixed data with missing cells in every partition type: 4 neomorphic,
# 8 transformational and 3 known-state characters (else 11 transformational).
.MissingData <- function(seed = 4L, nTip = 9L, pMissing = 0.25,
                         known = TRUE) {
  set.seed(seed)
  tree <- Preorder(ape::rtree(nTip, rooted = FALSE,
                              br = function(n) runif(n, 0.05, 0.6)))
  nChar <- 15L
  repeat {
    mat <- cbind(matrix(sample(0:1, nTip * 4L, TRUE), nTip),
                 matrix(sample(0:2, nTip * 11L, TRUE), nTip))
    mat[matrix(runif(nTip * nChar) < pMissing, nTip)] <- NA
    kObs <- apply(mat, 2, function(x) length(unique(x[!is.na(x)])))
    if (all(kObs[1:4] == 2L) && all(kObs >= 2L)) break
  }
  tokens <- matrix(as.character(mat), nTip,
                   dimnames = list(tree$tip.label, NULL))
  tokens[is.na(mat)] <- sample(c("?", "-"), sum(is.na(mat)), TRUE)
  knownStates <- if (known) c("13" = 4L, "14" = 4L, "15" = 3L)
  mkd <- MkPrimeData(MatrixToPhyDat(tokens), neomorphic = 1:4,
                     knownStates = knownStates)
  list(tree = tree, mkd = mkd)
}

test_that("fixture has missing cells in every partition", {
  d <- .MissingData()
  expect_true(all(vapply(d$mkd$partitions,
                         function(p) anyNA(p$tip_states), logical(1))))
  expect_setequal(unique(d$mkd$type),
                  c("neomorphic", "transformational", "known"))
})

.DirectLogLik <- function(d, coding, rateLoss, rateLogSd, rateNeo,
                          kPrime = d$mkd$kObs) {
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(coding = coding), d$tree,
                                    d$mkd)
  dataPtr <- MkPrime:::.InitMcmcData(d$mkd, model)
  # Return:
  cpp_log_likelihood_xptr(dataPtr, d$tree$edge[, 1], d$tree$edge[, 2],
                          d$tree$edge.length, as.integer(kPrime),
                          rateLoss = rateLoss, rateLogSd = rateLogSd,
                          rateNeo = rateNeo)
}

test_that("each character is conditioned on its own observed tips", {
  d <- .MissingData()
  mkd <- d$mkd
  tree <- d$tree
  nCat <- MkPrimeModel()$nCat
  kPrime <- mkd$kObs
  kPrime[mkd$type == "transformational"][1:3] <- 4L
  isNeo <- mkd$type == "neomorphic"
  for (par in list(c(rl = 0.4, sd = 0, rn = 1.7),
                   c(rl = 2.5, sd = 0.6, rn = 0.5))) {
    rates <- DiscreteLognormalRates(par[["sd"]], nCat)
    nTotal <- length(isNeo)
    neoScale <- par[["rn"]] / (1 + par[["rn"]]) * nTotal / sum(isNeo)
    transScale <- 1 / (1 + par[["rn"]]) * nTotal / sum(!isNeo)
    for (coding in c("variable", "informative")) {
      correction <- sum(vapply(seq_along(isNeo), function(i) {
        part <- mkd$partitions[[which(vapply(
          mkd$partitions, function(p) i %in% p$char_indices, logical(1)))]]
        col <- part$tip_states[, match(i, part$char_indices)]
        scaled <- tree
        scaled$edge.length <- tree$edge.length *
          if (isNeo[i]) neoScale else transScale
        k <- if (mkd$type[i] == "known") part$k else kPrime[i]
        log1p(-.BruteAscProb(scaled, is.na(col), k, isNeo[i], par[["rl"]],
                             rates, coding == "informative"))
      }, double(1)))
      raw <- .DirectLogLik(d, "none", par[["rl"]], par[["sd"]], par[["rn"]],
                           kPrime)
      expect_equal(
        .DirectLogLik(d, coding, par[["rl"]], par[["sd"]], par[["rn"]],
                      kPrime),
        raw - correction, tolerance = 1e-9, info = coding)
      expect_equal(
        MkpLogLikelihood(tree, mkd, kPrime = kPrime, rate_loss = par[["rl"]],
                         rate_log_sd = par[["sd"]], rate_neo = par[["rn"]],
                         coding = coding),
        raw - correction, tolerance = 1e-9, info = coding)
    }
  }
})

# Cached, partial-CL and Gibbs paths each rebuild the correction; every one
# must agree with a fresh full evaluation after every move.
.MissingDrift <- function(d, model, moveNames, nMove = 300L) {
  model <- MkPrime:::.FinalizeModel(model, d$tree, d$mkd)
  dataPtr <- MkPrime:::.InitMcmcData(d$mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(
    MkPrime:::.InitState(d$tree, d$mkd, model))
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  moves <- MkPrime:::.kMoveTypes[moveNames]
  drift <- 0
  for (move in sample(moves, nMove, replace = TRUE)) {
    do_move_cpp(dataPtr, statePtr, move, 0L, 0.5, 10, 1L, 1.0)
    drift <- max(drift, abs(get_state_log_lik(statePtr) -
                              eval_full_loglik_cpp(dataPtr, statePtr)))
  }
  # Return:
  drift
}

test_that("cached and Gibbs paths carry the per-mask correction", {
  d <- .MissingData()
  moveNames <- c("nni", "branch_lengths", "dirichlet_branch", "gibbs_spr",
                 "gibbs_subtree_swap", "gibbs_kPrime")
  set.seed(8812)
  for (coding in c("variable", "informative")) {
    expect_lt(.MissingDrift(d, MkPrimeModel(coding = coding), moveNames),
              1e-6, label = coding)
  }
  expect_lt(.MissingDrift(d, MkPrimeModel(qHeterogeneity = TRUE),
                          c("gibbs_spr", "gibbs_subtree_swap",
                            "gibbs_kPrime", "branch_lengths")),
            1e-6, label = "qHeterogeneity")
})

# The marginal-k evaluator scores k' through the Gibbs sweep's per-(k, mask)
# correction cache; sampled-k scores it through the full evaluator.
test_that("the k' sweep conditions each character on its observed tips", {
  d <- .MissingData(known = FALSE)
  p <- 0.7
  Ptrs <- function(mode) {
    model <- MkPrime:::.FinalizeModel(
      MkPrimeModel(kPrimePrior = "geometric", likelihoodMode = mode,
                   priorVariant = "conditional", coding = "variable",
                   relabel = FALSE),
      d$tree, d$mkd)
    state <- MkPrime:::.InitState(d$tree, d$mkd, model)
    state$p <- p
    state$rate_log_sd <- 0
    state$tree_length <- sum(d$tree$edge.length)
    state$rel_br_lengths <- d$tree$edge.length / state$tree_length
    dataPtr <- MkPrime:::.InitMcmcData(d$mkd, model)
    statePtr <- MkPrime:::.InitMcmcChain(state)
    fill_partition_cache(dataPtr, statePtr)
    list(dataPtr = dataPtr, statePtr = statePtr, state = state)
  }
  marg <- Ptrs("marginal_k")
  samp <- Ptrs("sampled_k")
  kObs <- d$mkd$kObs
  LogLik <- function(kPrime) {
    cpp_log_likelihood_xptr(samp$dataPtr, d$tree$edge[, 1], d$tree$edge[, 2],
                            d$tree$edge.length, as.integer(kPrime),
                            rateLoss = samp$state$rate_loss, rateLogSd = 0,
                            rateNeo = samp$state$rate_neo, betaScale = 1)
  }
  base <- LogLik(kObs)
  # Weights p (1 - p)^u, u = k' - kObs; the tail beyond u = 14 is negligible.
  u <- 0:14
  perChar <- vapply(which(d$mkd$type == "transformational"), function(i) {
    terms <- vapply(u, function(ui) {
      kPrime <- kObs
      kPrime[i] <- kObs[i] + ui
      LogLik(kPrime) - base
    }, double(1)) + log(p) + u * log1p(-p)
    max(terms) + log(sum(exp(terms - max(terms))))
  }, double(1))
  expect_equal(eval_full_loglik_cpp(marg$dataPtr, marg$statePtr),
               base + sum(perChar), tolerance = 1e-7)
})
