# coding = "informative" must remove exactly the parsimony-uninformative
# patterns: those in which fewer than two states each occur twice or more.
# For k = 2 that is constant + singleton; for k >= 3 it also holds patterns
# such as (0, 0, 0, 1, 2), which the singleton term alone misses (#154).

# Exhaustive oracle: Felsenstein pruning under JC(k) over all k^nObs
# patterns of the observed tips, with missing tips marginalised.
.JcPatternProb <- function(parent, child, el, nTip, k, tipState, rates) {
  nNode <- max(parent)
  mean(vapply(rates, function(rate) {
    cl <- matrix(1, nNode, k)
    for (tip in seq_len(nTip)) {
      if (!is.na(tipState[tip])) {
        cl[tip, ] <- 0
        cl[tip, tipState[tip] + 1L] <- 1
      }
    }
    for (e in rev(seq_along(parent))) {
      expTerm <- exp(-k * el[e] * rate / (k - 1))
      pSame <- 1 / k + (1 - 1 / k) * expTerm
      pDiff <- (1 - expTerm) / k
      pMat <- matrix(pDiff, k, k)
      diag(pMat) <- pSame
      cl[parent[e], ] <- cl[parent[e], ] * (pMat %*% cl[child[e], ])
    }
    sum(cl[nTip + 1L, ]) / k
  }, double(1)))
}

.UninformativeMass <- function(tree, k, rates = 1, missing = NULL) {
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  nTip <- length(tree$tip.label)
  observed <- setdiff(seq_len(nTip), missing)
  patterns <- as.matrix(expand.grid(rep(list(seq_len(k) - 1L),
                                        length(observed))))
  total <- c(const = 0, uninf = 0, all = 0)
  for (i in seq_len(nrow(patterns))) {
    tipState <- rep(NA_integer_, nTip)
    tipState[observed] <- patterns[i, ]
    p <- .JcPatternProb(parent, child, tree$edge.length, nTip, k, tipState,
                        rates)
    counts <- tabulate(patterns[i, ] + 1L, k)
    total["all"] <- total["all"] + p
    if (sum(counts > 0) == 1L) total["const"] <- total["const"] + p
    if (sum(counts >= 2L) < 2L) total["uninf"] <- total["uninf"] + p
  }
  total
}

test_that("uninf_nonconst_prob_jc matches exhaustive enumeration", {
  set.seed(1540)
  cases <- list(
    list(nTip = 5L, k = 2L), list(nTip = 5L, k = 3L),
    list(nTip = 5L, k = 4L), list(nTip = 6L, k = 3L),
    list(nTip = 4L, k = 5L), list(nTip = 3L, k = 3L)
  )
  for (case in cases) {
    for (rooted in c(TRUE, FALSE)) {
      tree <- Preorder(ape::rtree(case$nTip, rooted = rooted,
                                  br = function(n) runif(n, 0.05, 0.6)))
      for (rates in list(1, c(0.3, 0.9, 1.8))) {
        oracle <- .UninformativeMass(tree, case$k, rates)
        expect_equal(oracle[["all"]], 1, tolerance = 1e-12)
        parent <- tree$edge[, 1]
        child <- tree$edge[, 2]
        el <- tree$edge.length
        pConst <- constant_site_prob_jc(parent, child, el, case$nTip, case$k,
                                        rep(1 / case$k, case$k), rates)
        pRest <- uninf_nonconst_prob_jc(parent, child, el, case$nTip,
                                        case$k, rates)
        info <- sprintf("nTip = %d, k = %d, rooted = %s, nCat = %d",
                        case$nTip, case$k, rooted, length(rates))
        expect_equal(pConst, oracle[["const"]], tolerance = 1e-12, info = info)
        expect_equal(pConst + pRest, oracle[["uninf"]], tolerance = 1e-12,
                     info = info)
      }
    }
  }
})

test_that("informative ascertainment marginalises missing tips exactly", {
  set.seed(1541)
  tree <- Preorder(ape::rtree(7L, rooted = FALSE,
                              br = function(n) runif(n, 0.05, 0.6)))
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  for (k in 3:4) {
    for (missing in list(c(2L, 5L), 1:4, 1:6)) {
      oracle <- .UninformativeMass(tree, k, missing = missing)
      miss <- seq_len(7L) %in% missing
      p <- asc_site_prob_missing(parent, child, tree$edge.length, 7L, k,
                                 FALSE, 1, 1, miss, TRUE)
      expect_equal(p, oracle[["uninf"]], tolerance = 1e-12,
                   info = sprintf("k = %d, missing = %s", k,
                                  paste(missing, collapse = ",")))
    }
  }
})

test_that("coding = 'informative' likelihood conditions on the exact set", {
  set.seed(1542)
  tree <- Preorder(ape::rtree(5L, rooted = FALSE,
                              br = function(n) runif(n, 0.05, 0.6)))
  mat <- matrix(c(0, 0, 1, 1, 2,
                  0, 1, 1, 2, 2,
                  1, 1, 0, 0, 0), nrow = 5L,
                dimnames = list(tree$tip.label, paste0("c", 1:3)))
  mkd <- MkPrimeData(MatrixToPhyDat(mat), knownStates = c(`1` = 3L, `2` = 3L, `3` = 3L))
  llVar <- MkpLogLikelihood(tree, mkd, coding = "variable")
  llInf <- MkpLogLikelihood(tree, mkd, coding = "informative")
  oracle <- .UninformativeMass(tree, 3L)
  pVar <- 1 - oracle[["const"]]
  pInf <- 1 - oracle[["uninf"]]
  expect_equal(llInf - llVar, -3 * (log(pInf) - log(pVar)), tolerance = 1e-10)
})

test_that("k = 2 informative mass marginalises missing tips exactly", {
  set.seed(1543)
  tree <- Preorder(ape::rtree(7L, rooted = FALSE,
                              br = function(n) runif(n, 0.05, 0.6)))
  rates <- c(0.4, 1.6)
  for (missing in list(integer(0), c(2L, 5L), 1:4)) {
    oracle <- .UninformativeMass(tree, 2L, rates, missing = missing)
    miss <- seq_len(7L) %in% missing
    expect_equal(
      asc_site_prob_missing(tree$edge[, 1], tree$edge[, 2], tree$edge.length,
                            7L, 2L, FALSE, 1, rates, miss, TRUE),
      oracle[["uninf"]], tolerance = 1e-12,
      info = paste(missing, collapse = ","))
  }
})

# #252 replaced the singleton kernels' one-pseudo-character-per-tip pruning
# with a two-partial pass shared across masks. Reference values are from the
# per-tip kernels.
.ParityTree <- function(nTip) {
  tree <- Preorder(BalancedTree(nTip))
  nEdge <- nrow(tree$edge)
  tree$edge.length <- 0.02 + 0.4 * ((seq_len(nEdge) * 0.618034) %% 1)
  tree
}

test_that("masked informative mass matches the per-tip kernels", {
  tree <- .ParityTree(24L)
  rates <- c(0.35, 0.9, 1.75)
  masks <- list(seq_len(24L) > 3L, (seq_len(24L) * 7L) %% 5L > 1L,
                seq_len(24L) > 11L, logical(24L))
  cfgs <- list(list(k = 2L, neo = FALSE, rl = 1),
               list(k = 3L, neo = FALSE, rl = 1),
               list(k = 4L, neo = FALSE, rl = 1),
               list(k = 2L, neo = TRUE, rl = 0.3),
               list(k = 2L, neo = TRUE, rl = 4))
  got <- t(vapply(masks, function(miss) {
    vapply(cfgs, function(cfg) {
      asc_site_prob_missing(tree$edge[, 1], tree$edge[, 2], tree$edge.length,
                            24L, cfg$k, cfg$neo, cfg$rl, rates, miss, TRUE)
    }, double(1))
  }, double(length(cfgs))))
  expected <- matrix(c(
    0.99999999999999989, 1.0000000000000018, 1, 1, 1.0000000000000007,
    0.16524578611100146, 0.18279127818222329, 0.20878860777916805,
    0.41299734949258726, 0.48267023856705127,
    0.18594945557368181, 0.19180540998096099, 0.20044918455134883,
    0.36917821700091774, 0.4259048656366955,
    0.034519821646537677, 0.040302537569187655, 0.044469679306904111,
    0.1010723288314237, 0.12706618949063464), 4L, byrow = TRUE)
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("informative Sun2018 likelihood matches the per-tip kernels", {
  nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
  skip_if(nexFile == "", message = "TreeSearch not available")
  sun <- ReadAsPhyDat(nexFile)
  tree <- .ParityTree(length(sun))
  tree$tip.label <- names(sun)
  binary <- c(1L, 3L, 5L, 6L, 7L, 9L, 10L, 11L, 12L, 13L, 14L, 15L)
  ll <- vapply(list(integer(0), binary), function(neo) {
    mkd <- suppressWarnings(MkPrimeData(sun, neomorphic = neo))
    model <- MkPrime:::.FinalizeModel(MkPrimeModel(coding = "informative"),
                                      tree, mkd)
    cpp_log_likelihood_xptr(MkPrime:::.InitMcmcData(mkd, model),
                            tree$edge[, 1], tree$edge[, 2], tree$edge.length,
                            as.integer(mkd$kObs), rateLoss = 0.6,
                            rateLogSd = 0.4, rateNeo = 1.3)
  }, double(1))
  expect_lt(max(abs(ll - c(-2829.4843325754368, -2826.5877050786489))),
            1e-10)
})
