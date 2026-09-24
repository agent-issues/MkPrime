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
