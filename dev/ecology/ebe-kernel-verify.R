# EBE kernel verification (T1 graceful degradation, T2 sanity).
# Compiles the package in THIS subprocess (pkgbuild + load_all), then calls the
# R entry points directly with hand-built small inputs.  Not a testthat file —
# a standalone smoke driver owned by the kernel agent.

suppressMessages({
  pkgbuild::compile_dll(".")
  devtools::load_all(".", quiet = TRUE)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Plain-R baseline MkN pruner (NO ecology), spec §2 P-formula with lambda = 2.
# This is the T1 reference (independent of the C++ kernel).
# ---------------------------------------------------------------------------
baseline_mkn_loglik <- function(parent, child, edgeLen, tipStates,
                                 rateLoss, rateMultipliers) {
  nTip  <- nrow(tipStates)
  nChar <- ncol(tipStates)
  nEdge <- length(parent)
  maxNode <- 2 * nTip - 1
  root <- nTip + 1
  pi1 <- 1 / (1 + rateLoss)
  pi0 <- rateLoss / (1 + rateLoss)
  nCat <- length(rateMultipliers)
  siteLik <- numeric(nChar)
  for (cat in seq_len(nCat)) {
    rate <- rateMultipliers[cat]
    cl <- matrix(0, nrow = maxNode, ncol = nChar * 2)  # node x (char*2)
    init <- logical(maxNode)
    # tips
    for (tip in seq_len(nTip)) {
      for (c in seq_len(nChar)) {
        st <- tipStates[tip, c]
        off <- (c - 1) * 2
        if (st < 0) { cl[tip, off + 1] <- 1; cl[tip, off + 2] <- 1 }
        else cl[tip, off + 1 + st] <- 1
      }
      init[tip] <- TRUE
    }
    for (e in seq(nEdge, 1)) {
      par <- parent[e]; ch <- child[e]
      t <- edgeLen[e] * rate
      ex <- exp(-2 * t)
      P00 <- pi0 + pi1 * ex; P01 <- pi1 - pi1 * ex
      P10 <- pi0 - pi0 * ex; P11 <- pi1 + pi0 * ex
      for (c in seq_len(nChar)) {
        off <- (c - 1) * 2
        cl0 <- cl[ch, off + 1]; cl1 <- cl[ch, off + 2]
        v0 <- P00 * cl0 + P01 * cl1
        v1 <- P10 * cl0 + P11 * cl1
        if (!init[par]) { cl[par, off + 1] <- v0; cl[par, off + 2] <- v1 }
        else { cl[par, off + 1] <- cl[par, off + 1] * v0
               cl[par, off + 2] <- cl[par, off + 2] * v1 }
      }
      init[par] <- TRUE
    }
    for (c in seq_len(nChar)) {
      off <- (c - 1) * 2
      siteLik[c] <- siteLik[c] + pi0 * cl[root, off + 1] + pi1 * cl[root, off + 2]
    }
  }
  sum(log(siteLik / nCat))
}

# Plain-R baseline constant-site correction (coding = variable), for T1 coding 1.
baseline_const_corr <- function(parent, child, edgeLen, nTip, nChar,
                                rateLoss, rateMultipliers) {
  corr <- 0
  for (c in seq_len(nChar)) {
    s0 <- matrix(0L, nrow = nTip, ncol = 1)
    s1 <- matrix(1L, nrow = nTip, ncol = 1)
    p0 <- exp(baseline_mkn_loglik(parent, child, edgeLen, s0, rateLoss, rateMultipliers))
    p1 <- exp(baseline_mkn_loglik(parent, child, edgeLen, s1, rateLoss, rateMultipliers))
    corr <- corr - log(1 - (p0 + p1))
  }
  corr
}

# ---------------------------------------------------------------------------
# Build a tiny ecology setup (4 tips, kEco = 2).
# ---------------------------------------------------------------------------
make_setup <- function(kEco = 2L, rateLoss = 1.3,
                       ecoTip = c(0L, 1L, 0L, 1L)) {
  # Simple 4-tip rooted tree, preorder edges.
  # nodes: 1..4 tips, 5 = root, 6 = internal
  parent <- c(5L, 6L, 6L, 5L)
  child  <- c(6L, 1L, 2L, 3L)
  # That is only 3 tips reachable; build a proper 4-tip balanced tree instead.
  # ((1,2),(3,4)) : root=5, c6=(1,2), c7=(3,4)
  parent <- c(5L, 6L, 6L, 5L, 7L, 7L)
  child  <- c(6L, 1L, 2L, 7L, 3L, 4L)
  edgeLen <- c(0.3, 0.2, 0.25, 0.35, 0.15, 0.4)
  marg <- MkPrime:::.EcologyNodeMarginals(parent, child, edgeLen,
                                          as.integer(ecoTip), kEco)
  wEdge <- MkPrime:::.EcologyEdgeWeights(marg, parent, child)
  pi1 <- 1 / (1 + rateLoss); pi0 <- rateLoss / (1 + rateLoss)
  list(parent = parent, child = child, edgeLen = edgeLen, wEdge = wEdge,
       kEco = kEco, rateLoss = rateLoss, rootFreqs = c(pi0, pi1))
}

cat("=== EBE kernel verification ===\n")

# ---- T1: all z = 0, phi != 1  =>  EBE == baseline MkN (coding 0 & 1) ----
s <- make_setup(kEco = 2L, rateLoss = 1.3)
nTip <- 4L
tipStates <- matrix(c(0L, 1L, 0L, 1L,
                      1L, 0L, 1L, 0L,
                      0L, 0L, 1L, 1L), nrow = 4, byrow = FALSE)
nChar <- ncol(tipStates)
zZero <- matrix(0L, nrow = nChar, ncol = s$kEco - 1L)
rateMults <- c(0.5, 1.0, 1.7)

phi_test <- 2.5  # phi != 1 — the key EBE guarantee

ll_ebe_c0 <- MkPrime:::.PruningMknEcology(
  s$parent, s$child, s$edgeLen, tipStates,
  s$rateLoss, s$rootFreqs, rateMults,
  s$wEdge, zZero, phi = phi_test, mode = 0L, refEcology = 0L)
ll_base_c0 <- baseline_mkn_loglik(s$parent, s$child, s$edgeLen, tipStates,
                                  s$rateLoss, rateMults)

cat(sprintf("T1 coding0:  EBE = %.15f\n", ll_ebe_c0))
cat(sprintf("T1 coding0:  base= %.15f\n", ll_base_c0))
cat(sprintf("T1 coding0:  |diff| = %.3e\n", abs(ll_ebe_c0 - ll_base_c0)))

# coding = 1: add constant-site correction via the R wrapper (same raw pruner).
corr_ebe <- 0
for (c in seq_len(nChar)) {
  pc <- MkPrime:::.ConstSiteProbMknEcology(
    s$parent, s$child, s$edgeLen, nTip,
    s$rateLoss, rateMults, s$wEdge,
    zZero[c, ], phi = phi_test, mode = 0L, refEcology = 0L)
  corr_ebe <- corr_ebe - log(1 - pc)
}
ll_ebe_c1 <- ll_ebe_c0 + corr_ebe
corr_base <- baseline_const_corr(s$parent, s$child, s$edgeLen, nTip, nChar,
                                 s$rateLoss, rateMults)
ll_base_c1 <- ll_base_c0 + corr_base

cat(sprintf("T1 coding1:  EBE = %.15f\n", ll_ebe_c1))
cat(sprintf("T1 coding1:  base= %.15f\n", ll_base_c1))
cat(sprintf("T1 coding1:  |diff| = %.3e\n", abs(ll_ebe_c1 - ll_base_c1)))

t1_pass <- abs(ll_ebe_c0 - ll_base_c0) < 1e-10 &&
           abs(ll_ebe_c1 - ll_base_c1) < 1e-10
cat(sprintf("T1 %s (tol 1e-10)\n", if (t1_pass) "PASS" else "FAIL"))

# ---- T1 per_ecology mode (mode 1): also must degrade with phi != 1 ----
s3 <- make_setup(kEco = 3L, rateLoss = 0.7, ecoTip = c(0L, 1L, 2L, 0L))
tip3 <- matrix(c(0L, 1L, 0L, 1L,
                 1L, 0L, 1L, 0L), nrow = 4, byrow = FALSE)
z3 <- matrix(0L, nrow = ncol(tip3), ncol = s3$kEco - 1L)
phi3 <- c(1.2, 2.5, 0.6)
ll_ebe3 <- MkPrime:::.PruningMknEcology(
  s3$parent, s3$child, s3$edgeLen, tip3,
  s3$rateLoss, s3$rootFreqs, rateMults,
  s3$wEdge, z3, phi = phi3, mode = 1L, refEcology = 0L)
ll_base3 <- baseline_mkn_loglik(s3$parent, s3$child, s3$edgeLen, tip3,
                                s3$rateLoss, rateMults)
cat(sprintf("T1 mode1:    |diff| = %.3e  %s\n",
            abs(ll_ebe3 - ll_base3),
            if (abs(ll_ebe3 - ll_base3) < 1e-10) "PASS" else "FAIL"))

# ---- T2 sanity: mixed z, phi != 1  =>  finite AND different from baseline ----
zMix <- matrix(c(1L, 2L, 0L), nrow = nChar, ncol = 1L)
ll_mix <- MkPrime:::.PruningMknEcology(
  s$parent, s$child, s$edgeLen, tipStates,
  s$rateLoss, s$rootFreqs, rateMults,
  s$wEdge, zMix, phi = phi_test, mode = 0L, refEcology = 0L)
cat(sprintf("T2 mixed z:  ll = %.15f  finite=%s  differs=%s\n",
            ll_mix, is.finite(ll_mix),
            abs(ll_mix - ll_base_c0) > 1e-6))
t2_pass <- is.finite(ll_mix) && abs(ll_mix - ll_base_c0) > 1e-6
cat(sprintf("T2 %s\n", if (t2_pass) "PASS" else "FAIL"))

# ---- T3 sanity: tilt direction.  z=1 (toward present) should RAISE the
# likelihood of an all-present character vs z=0; z=2 should lower it. ----
tipPres <- matrix(rep(1L, 4), nrow = 4, ncol = 1L)  # all present
zc0 <- matrix(0L, 1, 1); zc1 <- matrix(1L, 1, 1); zc2 <- matrix(2L, 1, 1)
ll_p0 <- MkPrime:::.PruningMknEcology(s$parent, s$child, s$edgeLen, tipPres,
  s$rateLoss, s$rootFreqs, 1, s$wEdge, zc0, phi = phi_test, mode = 0L, refEcology = 0L)
ll_p1 <- MkPrime:::.PruningMknEcology(s$parent, s$child, s$edgeLen, tipPres,
  s$rateLoss, s$rootFreqs, 1, s$wEdge, zc1, phi = phi_test, mode = 0L, refEcology = 0L)
ll_p2 <- MkPrime:::.PruningMknEcology(s$parent, s$child, s$edgeLen, tipPres,
  s$rateLoss, s$rootFreqs, 1, s$wEdge, zc2, phi = phi_test, mode = 0L, refEcology = 0L)
cat(sprintf("T-tilt all-present: z0=%.6f z1(present)=%.6f z2(absent)=%.6f\n",
            ll_p0, ll_p1, ll_p2))
cat(sprintf("T-tilt direction %s (expect z1 > z0 > z2)\n",
            if (ll_p1 > ll_p0 && ll_p0 > ll_p2) "PASS" else "CHECK"))

stopifnot(t1_pass, t2_pass)
cat("=== ALL GATING CHECKS PASS ===\n")
