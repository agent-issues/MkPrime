test_that("Constant site prob is in [0, 1]", {
  tree <- read.tree(text = "((t1:0.1,t2:0.1):0.1,t3:0.2);")
  tree <- Preorder(tree)

  for (k in c(2, 3, 5)) {
    pconst <- MkPrime:::constant_site_prob_jc(
      tree$edge[, 1], tree$edge[, 2], tree$edge.length,
      3L, k, rep(1 / k, k), 1.0
    )
    expect_true(pconst > 0 && pconst < 1,
                label = sprintf("k=%d, pconst=%g", k, pconst))
  }
})


test_that("Constant site prob → 1 with zero-length branches", {
  # Near-zero branch lengths: tips are nearly identical to root
  tree <- read.tree(text = "((t1:1e-10,t2:1e-10):1e-10,t3:1e-10);")
  tree <- Preorder(tree)

  pconst <- MkPrime:::constant_site_prob_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    3L, 2L, c(0.5, 0.5), 1.0
  )

  # With near-zero branches, almost all sites are constant

  expect_gt(pconst, 0.99)
})


test_that("Constant site prob → 1/k with very long branches", {
  tree <- read.tree(text = "((t1:100,t2:100):100,t3:100);")
  tree <- Preorder(tree)

  for (k in c(2, 3, 4)) {
    pconst <- MkPrime:::constant_site_prob_jc(
      tree$edge[, 1], tree$edge[, 2], tree$edge.length,
      3L, k, rep(1 / k, k), 1.0
    )
    # With very long branches, states are independent at each tip
    # P(all same) = sum_s (1/k)^nTip = k * (1/k)^3 = 1/k^2 for 3 tips
    expected <- k * (1 / k)^3
    expect_equal(pconst, expected, tolerance = 0.01,
                 label = sprintf("k=%d", k))
  }
})


test_that("Ascertainment correction makes likelihood more negative", {
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  tree <- Preorder(tree)
  root_freqs <- rep(1 / 3, 3)
  tips <- matrix(c(0L, 1L, 2L, 0L), ncol = 1)

  ll_uncorrected <- MkPrime:::pruning_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    tips, 3L, root_freqs
  )

  pconst <- MkPrime:::constant_site_prob_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 3L, root_freqs, 1.0
  )

  nChar <- 1L
  ll_corrected <- ll_uncorrected - nChar * log(1 - pconst)

  # Ascertainment correction makes likelihood LESS negative:
  # L_corr = L / P(variable), and P(variable) < 1, so L_corr > L
  expect_gt(ll_corrected, ll_uncorrected)
})


test_that("MkN constant site prob with rate_loss=1 matches JC(2)", {
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  tree <- Preorder(tree)
  root_freqs <- c(0.5, 0.5)
  rates <- MkPrime:::DiscreteLognormalRates(0.5, 6L)

  pconst_jc <- MkPrime:::constant_site_prob_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 2L, root_freqs, rates
  )

  pconst_mkn <- MkPrime:::constant_site_prob_mkn(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 1.0, root_freqs, rates
  )

  expect_equal(pconst_mkn, pconst_jc, tolerance = 1e-12)
})


# --- Singleton (autapomorphy) probability tests ---

test_that("Singleton site prob JC is in (0, 1)", {
  tree <- read.tree(text = "((t1:0.1,t2:0.1):0.1,t3:0.2);")
  tree <- Preorder(tree)

  for (k in c(2, 3, 5)) {
    psingle <- MkPrime:::singleton_site_prob_jc(
      tree$edge[, 1], tree$edge[, 2], tree$edge.length,
      3L, k, rep(1 / k, k), 1.0
    )
    expect_true(psingle > 0 && psingle < 1,
                label = sprintf("k=%d, psingle=%g", k, psingle))
  }
})


test_that("Singleton prob → 0 with zero-length branches", {
  tree <- read.tree(text = "((t1:1e-10,t2:1e-10):1e-10,t3:1e-10);")
  tree <- Preorder(tree)

  psingle <- MkPrime:::singleton_site_prob_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    3L, 2L, c(0.5, 0.5), 1.0
  )

  # With near-zero branches, singletons are extremely rare
  expect_lt(psingle, 0.01)
})


test_that("P(constant) + P(singleton) < 1", {
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  tree <- Preorder(tree)

  for (k in c(2, 3, 5)) {
    rf <- rep(1 / k, k)
    pconst <- MkPrime:::constant_site_prob_jc(
      tree$edge[, 1], tree$edge[, 2], tree$edge.length,
      4L, k, rf, 1.0
    )
    psingle <- MkPrime:::singleton_site_prob_jc(
      tree$edge[, 1], tree$edge[, 2], tree$edge.length,
      4L, k, rf, 1.0
    )
    expect_lt(pconst + psingle, 1.0,
              label = sprintf("k=%d, P(uninf)=%g", k, pconst + psingle))
    expect_gt(pconst + psingle, 0.0)
  }
})


test_that("MkN singleton prob with rate_loss=1 matches JC(2)", {
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  tree <- Preorder(tree)
  root_freqs <- c(0.5, 0.5)

  psingle_jc <- MkPrime:::singleton_site_prob_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 2L, root_freqs, 1.0
  )

  psingle_mkn <- MkPrime:::singleton_site_prob_mkn(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 1.0, root_freqs, 1.0
  )

  expect_equal(psingle_mkn, psingle_jc, tolerance = 1e-12)
})


test_that("MkN singleton prob varies with rate_loss", {
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  tree <- Preorder(tree)

  # Use rate_loss values that aren't reciprocals (avoids state-swap symmetry)
  ps1 <- MkPrime:::singleton_site_prob_mkn(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 1.0, c(0.5, 0.5), 1.0
  )

  ps2 <- MkPrime:::singleton_site_prob_mkn(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 3.0, c(3 / 4, 1 / 4), 1.0
  )

  expect_false(isTRUE(all.equal(ps1, ps2)))
})


test_that("Singleton prob works with ACRV", {
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  tree <- Preorder(tree)
  rf <- rep(1 / 3, 3)

  psingle_no_acrv <- MkPrime:::singleton_site_prob_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 3L, rf, 1.0
  )

  rates <- MkPrime:::DiscreteLognormalRates(1.0, 6L)
  psingle_acrv <- MkPrime:::singleton_site_prob_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 3L, rf, rates
  )

  expect_false(isTRUE(all.equal(psingle_no_acrv, psingle_acrv)))
  expect_true(psingle_acrv > 0 && psingle_acrv < 1)
})


test_that("P(uninf) → 1 with very short branches", {
  tree <- read.tree(text = "((t1:0.001,t2:0.001):0.001,(t3:0.001,t4:0.001):0.001);")
  tree <- Preorder(tree)

  for (k in c(2, 3)) {
    rf <- rep(1 / k, k)
    pconst <- MkPrime:::constant_site_prob_jc(
      tree$edge[, 1], tree$edge[, 2], tree$edge.length,
      4L, k, rf, 1.0
    )
    psingle <- MkPrime:::singleton_site_prob_jc(
      tree$edge[, 1], tree$edge[, 2], tree$edge.length,
      4L, k, rf, 1.0
    )
    expect_gt(pconst + psingle, 0.99,
              label = sprintf("k=%d, P(uninf)=%g", k, pconst + psingle))
  }
})


# --- coding = "informative" through MkpLogLikelihood ---

test_that("MkpLogLikelihood accepts coding='informative'", {

  set.seed(4281)
  tree <- rtree(6, br = runif)
  tree <- unroot(tree)

  mat <- matrix(c(0L, 0L, 1L, 1L, 0L, 1L,
                  0L, 0L, 0L, 1L, 1L, 1L,
                  0L, 0L, 0L, 0L, 0L, 1L),
                nrow = 6,
                dimnames = list(tree$tip.label, NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  ll_none <- MkpLogLikelihood(tree, mkd, coding = "none", relabel = FALSE)
  ll_var  <- MkpLogLikelihood(tree, mkd, coding = "variable", relabel = FALSE)
  ll_inf  <- MkpLogLikelihood(tree, mkd, coding = "informative", relabel = FALSE)

  # Correction increases log-likelihood (conditions on observable patterns)
  expect_gt(ll_var, ll_none)
  expect_gt(ll_inf, ll_var)
})


test_that("Informative correction matches manual calculation", {

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")

  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L),
                nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  ll_none <- MkpLogLikelihood(tree, mkd, coding = "none", relabel = FALSE)
  ll_inf  <- MkpLogLikelihood(tree, mkd, coding = "informative", relabel = FALSE)

  # Manual calculation
  tree2 <- Preorder(tree)
  p <- tree2$edge[, 1]; ch <- tree2$edge[, 2]; el <- tree2$edge.length
  pconst <- MkPrime:::constant_site_prob_jc(p, ch, el, 4L, 2L, c(0.5, 0.5), 1.0)
  psingle <- MkPrime:::singleton_site_prob_jc(p, ch, el, 4L, 2L, c(0.5, 0.5), 1.0)
  expected <- ll_none - mkd$nChar * log(1 - pconst - psingle)

  expect_equal(ll_inf, expected, tolerance = 1e-10)
})


test_that("Informative correction works with neomorphic chars", {

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")

  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L),
                nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = 1:2)

  ll_var <- MkpLogLikelihood(tree, mkd, rate_loss = 1.5,
                              coding = "variable", relabel = FALSE)
  ll_inf <- MkpLogLikelihood(tree, mkd, rate_loss = 1.5,
                              coding = "informative", relabel = FALSE)

  expect_gt(ll_inf, ll_var)
})


test_that("MkPrimeModel accepts coding='informative'", {
  mod <- MkPrimeModel(coding = "informative")
  expect_equal(mod$coding, "informative")
})


test_that("Constant site prob works with ACRV", {
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  tree <- Preorder(tree)
  root_freqs <- rep(1 / 3, 3)

  pconst_no_acrv <- MkPrime:::constant_site_prob_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 3L, root_freqs, 1.0
  )

  rates <- MkPrime:::DiscreteLognormalRates(1.0, 6L)
  pconst_acrv <- MkPrime:::constant_site_prob_jc(
    tree$edge[, 1], tree$edge[, 2], tree$edge.length,
    4L, 3L, root_freqs, rates
  )

  # ACRV should change the constant site probability
  expect_false(isTRUE(all.equal(pconst_no_acrv, pconst_acrv)))
  # But it should still be in [0, 1]
  expect_true(pconst_acrv > 0 && pconst_acrv < 1)
})


test_that("het constant-site probability matches a brute-force reference", {
  tree <- Preorder(
    ape::read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);"))
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  edgeLength <- tree$edge.length
  nTip <- 4L
  baseRL <- 1
  rates <- c(0.4, 1.0, 1.9)

  # Every spike rotation and every constant pattern, evaluated separately.
  Reference <- function(k, betaBins) {
    nRot <- if (k == 2L) 1L else k
    maxNode <- 2L * nTip - 1L
    total <- 0
    for (rate in rates) for (b in betaBins) for (rot in seq_len(nRot)) {
      stationary <- if (k == 2L) {
        gain <- 2 * b / (1 + baseRL)
        loss <- 2 * (1 - b) * baseRL / (1 + baseRL)
        c(1 - gain / (gain + loss), gain / (gain + loss))
      } else {
        f <- rep((1 - b) / (k - 1), k)
        f[rot] <- b
        f
      }
      mu <- 1 / (1 - sum(stationary^2))
      for (s in seq_len(k)) {
        cl <- matrix(0, maxNode, k)
        cl[seq_len(nTip), s] <- 1
        seeded <- c(rep(TRUE, nTip), rep(FALSE, maxNode - nTip))
        for (e in rev(seq_along(parent))) {
          par <- parent[e]
          d <- exp(-mu * edgeLength[e] * rate)
          v <- (1 - d) * sum(stationary * cl[child[e], ]) + d * cl[child[e], ]
          cl[par, ] <- if (seeded[par]) cl[par, ] * v else v
          seeded[par] <- TRUE
        }
        total <- total + sum(stationary * cl[nTip + 1L, ])
      }
    }
    total / (length(rates) * length(betaBins) * nRot)
  }

  for (k in c(2L, 3L, 5L, 9L)) {
    betaBins <- seq(1 / k, 0.9, length.out = 4)
    expect_equal(
      MkPrime:::het_const_site_prob(parent, child, edgeLength, nTip, k,
                                    baseRL, betaBins, rates),
      Reference(k, betaBins)
    )
  }
})
