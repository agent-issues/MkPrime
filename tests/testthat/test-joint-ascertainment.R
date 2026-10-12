# The joint ascertainment normaliser D = sum_k pi(k) Pr(V | k) needs
# Pr(V | k) at every k in the prior's support, from a lumped kernel that
# serves all missing-data masks at once and interpolates the k tail.

.JointTree <- function(seed, nTip, rooted = FALSE) {
  set.seed(seed)
  Preorder(ape::rtree(nTip, rooted = rooted,
                      br = function(n) runif(n, 0.05, 0.8)))
}

# Mask 1 has no missing tip; masks 3 and 4 differ at one tip; mask 5 leaves
# one tip observed.
.JointMasks <- function(nTip, nMask, seed) {
  set.seed(seed)
  masks <- matrix(0L, nTip, nMask)
  for (j in seq_len(nMask)[-1]) {
    masks[, j] <- as.integer(runif(nTip) < runif(1, 0.1, 0.7))
  }
  masks[, 4] <- masks[, 3]
  masks[1, 4] <- 1L - masks[1, 3]
  masks[, 5] <- 1L
  masks[2, 5] <- 0L
  masks
}

.PerKConstProb <- function(tree, masks, k, rates) {
  vapply(seq_len(ncol(masks)), function(m) {
    MkPrime:::asc_site_prob_missing(
      tree$edge[, 1], tree$edge[, 2], tree$edge.length, NTip(tree), k, FALSE,
      1, rates, as.logical(masks[, m]), FALSE)
  }, double(1))
}

test_that("lumped JC kernel matches the per-k constant-site probability", {
  for (rooted in c(FALSE, TRUE)) {
    for (rates in list(1, c(0.3, 0.9, 1.8))) {
      tree <- .JointTree(1, 8L, rooted)
      masks <- .JointMasks(8L, 6L, 2)
      got <- MkPrime:::lumped_const_probs(
        tree$edge[, 1], tree$edge[, 2], tree$edge.length, 8L, rates, masks,
        as.numeric(2:7))
      expected <- vapply(2:7, function(k) .PerKConstProb(tree, masks, k, rates),
                         double(ncol(masks)))
      expect_equal(got, expected, tolerance = 1e-12)
      expect_equal(
        got[1, ],
        vapply(2:7, function(k) {
          MkPrime:::constant_site_prob_jc(
            tree$edge[, 1], tree$edge[, 2], tree$edge.length, 8L, k,
            rep(1 / k, k), rates)
        }, double(1)),
        tolerance = 1e-12)
    }
  }
})

test_that("masks beyond the dense class table share partials by hashing", {
  tree <- .JointTree(11, 40L)
  masks <- .JointMasks(40L, 300L, 99)
  got <- MkPrime:::lumped_const_probs(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length, 40L, c(0.5, 1.5), masks,
    c(2, 5))
  expected <- vapply(c(2, 5), function(k) {
    .PerKConstProb(tree, masks, k, c(0.5, 1.5))
  }, double(ncol(masks)))
  expect_equal(got, expected, tolerance = 1e-12)
})

# Brute-force F81 pruning over every constant pattern, Beta bins as
# compute_het_bins: independent of the lumped kernel.
.HetBins <- function(alpha, k, nBins) {
  b <- (k - 1) * alpha
  vapply(seq_len(nBins) - 1, function(i) {
    lo <- qbeta(i / nBins, alpha, b)
    hi <- qbeta((i + 1) / nBins, alpha, b)
    alpha / (alpha + b) * (pbeta(hi, alpha + 1, b) - pbeta(lo, alpha + 1, b)) /
      (pbeta(hi, alpha, b) - pbeta(lo, alpha, b))
  }, double(1))
}

.BruteHetConst <- function(tree, missing, k, alpha, nBins, rates) {
  nTip <- NTip(tree)
  total <- 0
  for (beta in .HetBins(alpha, k, nBins)) for (rate in rates) {
    pi <- c(beta, rep((1 - beta) / (k - 1), k - 1))
    mu <- 1 / (1 - sum(pi^2))
    for (s in seq_len(k)) {
      cl <- matrix(1, max(tree$edge), k)
      cl[which(!missing), ] <- 0
      cl[which(!missing), s] <- 1
      for (e in rev(seq_len(nrow(tree$edge)))) {
        d <- exp(-mu * tree$edge.length[e] * rate)
        P <- (1 - d) * matrix(pi, k, k, byrow = TRUE) + d * diag(k)
        cl[tree$edge[e, 1], ] <- cl[tree$edge[e, 1], ] *
          as.vector(P %*% cl[tree$edge[e, 2], ])
      }
      total <- total + sum(pi * cl[nTip + 1, ])
    }
  }
  # Return:
  total / (nBins * length(rates))
}

test_that("lumped het kernel matches F81 pruning over constant patterns", {
  tree <- .JointTree(3, 6L)
  masks <- .JointMasks(6L, 5L, 5)
  rates <- c(0.4, 1.6)
  for (alpha in c(0.3, 2.5)) {
    got <- MkPrime:::lumped_const_probs(
      tree$edge[, 1], tree$edge[, 2], tree$edge.length, 6L, rates, masks,
      as.numeric(2:4), betaScale = alpha, nBetaCat = 4L)
    expected <- vapply(2:4, function(k) vapply(seq_len(ncol(masks)), function(m) {
      .BruteHetConst(tree, as.logical(masks[, m]), k, alpha, 4L, rates)
    }, double(1)), double(ncol(masks)))
    expect_equal(got, expected, tolerance = 1e-12)
    expect_equal(
      got[1, ],
      vapply(2:4, function(k) {
        MkPrime:::het_const_site_prob(
          tree$edge[, 1], tree$edge[, 2], tree$edge.length, 6L, k, 1,
          .HetBins(alpha, k, 4L), rates)
      }, double(1)),
      tolerance = 1e-12)
  }
})

.JointData <- function() {
  set.seed(4)
  nTip <- 9L
  tree <- Preorder(ape::rtree(nTip, rooted = FALSE,
                              br = function(n) runif(n, 0.05, 0.6)))
  repeat {
    mat <- matrix(sample(0:3, nTip * 24L, TRUE, prob = 4:1), nTip)
    mat[matrix(runif(nTip * 24L) < 0.25, nTip)] <- NA
    if (all(apply(mat, 2, function(x) length(unique(x[!is.na(x)]))) >= 2L)) {
      break
    }
  }
  tokens <- matrix(as.character(mat), nTip,
                   dimnames = list(tree$tip.label, NULL))
  tokens[is.na(mat)] <- "?"
  list(tree = tree, mkd = MkPrimeData(MatrixToPhyDat(tokens)))
}

.JointDataPtr <- function(d, ...) {
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(...), d$tree, d$mkd)
  MkPrime:::.InitMcmcData(d$mkd, model)
}

.LogPriorTable <- function(lw) lw - log(sum(exp(lw - max(lw)))) - max(lw)

test_that("k' prior table matches each arm on its support", {
  d <- .JointData()
  K <- 30L
  p <- 0.3
  for (kObs in 2:4) {
    lo <- c(A = 2L, B = kObs)
    k <- list(A = 2:K, B = kObs:K)
    for (variant in c("unconditional", "conditional")) {
      v <- if (variant == "conditional") "B" else "A"
      ptr <- .JointDataPtr(d, kPrimePrior = "geometric",
                           priorVariant = variant, kprimeTruncK = K)
      got <- MkPrime:::kprime_log_prior_table_cpp(ptr, kObs, p, 1, 1)
      expect_equal(attr(got, "lo"), lo[[v]])
      expect_equal(as.numeric(got),
                   .LogPriorTable((k[[v]] - lo[[v]]) * log1p(-p)))
    }
    ptr <- .JointDataPtr(d, kPrimePrior = "beta_geometric", kprimeTruncK = K)
    expect_equal(
      as.numeric(MkPrime:::kprime_log_prior_table_cpp(ptr, kObs, p, 0.7, 2)),
      .LogPriorTable(lbeta(1.7, 2 + k$B - kObs)))
    ptr <- .JointDataPtr(d, kPrimePrior = "logseries", kprimeTruncK = K)
    expect_equal(
      as.numeric(MkPrime:::kprime_log_prior_table_cpp(ptr, kObs, p, 1, 1)),
      .LogPriorTable(k$A * log(0.7) - log(k$A)))
    model <- MkPrime:::.FinalizeModel(
      MkPrimeModel(kPrimePrior = "empirical_geometric",
                   priorVariant = "conditional", kprimeTruncK = K),
      d$tree, d$mkd)
    ptr <- MkPrime:::.InitMcmcData(d$mkd, model)
    expect_equal(
      as.numeric(MkPrime:::kprime_log_prior_table_cpp(ptr, kObs, p, 1, 1)),
      .LogPriorTable(
        MkPrime:::.LogEmpGeomCapped(model$empiricalNObs, p, K)$logPk[k$B]))
  }
})

test_that("joint normaliser sums Pr(V | k) over the k' prior", {
  d <- .JointData()
  K <- 20L
  ptr <- .JointDataPtr(d, kPrimePrior = "geometric",
                       priorVariant = "conditional", kprimeTruncK = K)
  rates <- MkPrime:::DiscreteLognormalRates(0.8, 4L)
  for (partIdx in seq_along(d$mkd$partitions)) {
    part <- d$mkd$partitions[[partIdx]]
    if (part$type != "transformational") next
    missing <- is.na(part$tip_states)
    key <- apply(missing, 2, paste, collapse = "")
    k <- part$kObs:K
    prior <- exp(.LogPriorTable((k - part$kObs) * log1p(-0.3)))
    expected <- sum(vapply(unique(key), function(mask) {
      masks <- matrix(as.integer(missing[, match(mask, key)]))
      pV <- 1 - vapply(k, function(kk) {
        .PerKConstProb(d$tree, masks, kk, rates)
      }, double(1))
      sum(key == mask) * log(sum(prior * pV))
    }, double(1)))
    expect_equal(
      as.numeric(MkPrime:::joint_asc_log_cpp(
        ptr, partIdx, d$tree$edge[, 1], d$tree$edge[, 2], d$tree$edge.length,
        rates, 0, 0.3, 1, 1, direct = TRUE)),
      expected, tolerance = 1e-12)
  }
})

test_that("interpolated k tail matches the direct sum", {
  d <- .JointData()
  for (het in c(FALSE, TRUE)) {
    for (arm in c("geometric", "beta_geometric")) {
      ptr <- .JointDataPtr(d, kPrimePrior = arm, kprimeTruncK = 200L,
                           qHeterogeneity = het)
      for (scale in c(0.05, 1, 6)) {
        args <- list(ptr, 2L, d$tree$edge[, 1], d$tree$edge[, 2],
                     d$tree$edge.length * scale, c(0.5, 1.5),
                     if (het) 0.7 else 0, 0.02, 0.5, 2)
        quad <- do.call(MkPrime:::joint_asc_log_cpp, args)
        expect_gt(attr(quad, "nNode"), 0L)
        expect_equal(
          as.numeric(quad),
          as.numeric(do.call(MkPrime:::joint_asc_log_cpp,
                             c(args, direct = TRUE))),
          tolerance = 1e-12, info = paste(het, arm, scale))
      }
    }
  }
})
