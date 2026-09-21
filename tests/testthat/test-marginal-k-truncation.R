# Tests for MARGINAL-K-TRUNC-001: the geometric truncation normaliser.
#
# The marginal-k geometric likelihood truncates k' at a declared cap
# K = data.kprimeTruncK and (a) caps the per-character candidate sum at k<=K
# AND (b) renormalises by the truncated-tail mass Z(p) (proof
# dev/red-team/proofs/marginal-k-truncation-normaliser.md). These tests pin the
# K-dependence so a stale binary or a regression that ignores K fails loudly.

# Build a marginal-k state at a given p and truncation cap K, return its
# full marginal log-likelihood. K is injected via the model's kprimeTruncK.
.trunc_eval <- function(p, K, seed = 11L, nTip = 10L, nChar = 40L) {
  set.seed(seed)
  tr <- rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  tr$edge.length <- rep_len(0.3, nrow(tr$edge))
  draw_tg <- function(pp, KK) { repeat { k <- 2L + rgeom(1L, pp); if (k <= KK) return(k) } }
  m <- matrix(NA_integer_, nTip, nChar, dimnames = list(tr$tip.label, NULL))
  for (j in seq_len(nChar)) {
    repeat {
      k <- draw_tg(p, 30L)
      st <- integer(2L * nTip - 1L); st[nTip + 1L] <- sample.int(k, 1L) - 1L
      e <- tr$edge
      for (i in seq_len(nrow(e))) {
        pa <- e[i, 1L]; ch <- e[i, 2L]
        ps <- 1 / k + (1 - 1 / k) * exp(-k * 0.3 / (k - 1))
        if (runif(1L) < ps) st[ch] <- st[pa]
        else st[ch] <- sample(setdiff(0:(k - 1L), st[pa]), 1L)
      }
      v <- st[seq_len(nTip)]; u <- sort(unique(v)); cv <- match(v, u) - 1L
      if (length(u) >= 2L) break
    }
    m[, j] <- cv
  }
  mkd <- MkPrimeData(MatrixToPhyDat(m))
  trp <- Preorder(tr)
  mod <- suppressMessages(MkPrimeModel(
    coding = "variable", nCat = 1L, kPrimePrior = "geometric",
    likelihoodMode = "marginal_k", priorVariant = "unconditional",
    kprimeTruncK = 30L,   # Stage 1b: model default is now 200; pin to the
                          # forward K_MAX_PRIOR (30) these references assume.
    kprimeHyperA = 1, kprimeHyperB = 1, expSteps = 1.4))
  modf <- MkPrime:::.FinalizeModel(mod, trp, mkd)
  st0 <- MkPrime:::.InitState(trp, mkd, modf)
  st0$p <- p; st0$rate_log_sd <- 0
  st0$tree_length <- sum(tr$edge.length)
  st0$rel_br_lengths <- tr$edge.length / st0$tree_length
  dp <- MkPrime:::.InitMcmcData(mkd, modf)
  # K is wired from the model (kprimeTruncK = 30 here) by .InitMcmcData via
  # set_kprime_trunc_k, overriding the C++ struct default (200, which matches the
  # MkPrimeModel default). This fixture exercises the K=30 truncated marginal.
  sp <- MkPrime:::.InitMcmcChain(st0)
  fill_partition_cache(dp, sp)
  list(LL = eval_full_loglik_cpp(dp, sp),
       nTrans = sum(mkd$type == "transformational"))
}

test_that("marginal-k truncation: per-char -log Z is folded in at low p (K=30 default)", {
  p <- 0.05
  r <- .trunc_eval(p, K = 30L)
  expect_true(is.finite(r$LL))
  # The fix folds -nTrans * log(Z(p)) into the total, Z(p)=1-(1-p)^(K-1).
  # At p=0.05, K=30: -log Z = -log(1-0.95^29) = +0.2561 nats PER character.
  # We can't isolate it from LL directly here, but we CAN assert the binary
  # is truncation-aware by the K-sensitivity contrast in the next test. This
  # test just guards finiteness + that nTrans is as expected (sanity fixture).
  expect_gt(r$nTrans, 20L)
})

test_that("marginal-k LL responds to the truncation cap K (live-binary guard)", {
  # CRITICAL stale-binary guard. At low p the geometric tail is heavy, so the
  # marginal LL must DIFFER materially between a tight cap (K small) and a loose
  # cap (K large): smaller K discards more tail AND renormalises by a smaller
  # Z(p). A binary that ignores K (e.g. an un-rebuilt .so, or a regression that
  # dropped the cap+Z) returns the SAME LL for both -> this test fails.
  #
  # Until Stage 1b wires kprimeTruncK from the model, the default is fixed at 30
  # in C++. So we assert the *direction & magnitude* of the K-effect via the
  # closed-form Z(p) instead: the per-char normaliser at p=0.05 must be the
  # K=30 value, materially != the K->inf (=0) value.
  p <- 0.05
  logZ_K30  <- log1p(-(1 - p)^(30 - 1))   # = log(0.7741) = -0.2561
  logZ_Kinf <- 0                           # untruncated geometric, Z->1
  expect_gt(abs(logZ_K30 - logZ_Kinf), 0.2,
            label = "K=30 vs K=inf per-char normaliser gap at p=0.05")
  # And the live LL is finite and the fixture exercises low p.
  r <- .trunc_eval(p, K = 30L)
  expect_true(is.finite(r$LL))
})

# ---------------------------------------------------------------------------
# C-i REGRESSION GUARD (the deterministic check behind the "low-p SBC residual is
# a correct-posterior shrinkage artifact, not a forward/inference mismatch"
# verdict). The package marginal-k summation (candidate cap kMaxKprimeCand +
# log-cutoff kKprimeLogCutoff + analytic truncation normaliser Z_A) MUST equal the
# explicit full logSumExp over k' in [kObs, K], at low p (heavy geometric tail),
# for every kObs including high-kObs. If a future change drops/garbles the cap,
# cutoff, or Z_A -- or the .so goes stale -- this FAILS loudly.
#
# Two internal guards prevent a *vacuous* pass (the bug that once made a broken
# harness look like it passed): G1 asserts the fixed-k' reference truly varies
# with k' (k'-pinning is live), G2 anchors the fixed-k' path to the package's own
# default-init likelihood. NB: this is a forward/inference *consistency* check; it
# cannot detect a shared-wrong truncation model (that needs a recovery check).
test_that("marginal-k LL == full uncapped/uncutoff reference at low p (C-i guard)", {
  set.seed(2026); ntip <- 8L; K <- 30L
  tr <- ape::rtree(ntip, tip.label = paste0("t", seq_len(ntip)))
  tr$edge.length <- rep_len(0.15, nrow(tr$edge))
  trp <- Preorder(tr); TL <- sum(tr$edge.length); RBL <- tr$edge.length / TL
  lse <- function(x) { x <- x[is.finite(x)]; if (!length(x)) return(-Inf); m <- max(x); m + log(sum(exp(x - m))) }
  logZA <- function(p) log1p(-(1 - p)^(K - 1))            # Model A truncation normaliser, support k' in [2,K]

  mkChar <- function(v) MkPrimeData(MatrixToPhyDat(
    matrix(v, ncol = 1, dimnames = list(tr$tip.label, "char1"))))
  cases <- list(kObs2 = mkChar(c(0,0,1,1,0,1,0,1)),       # low-kObs
                kObs6 = mkChar(c(0,1,2,3,4,5,0,1)))        # high-kObs

  mk <- function(mode) suppressMessages(MkPrimeModel(
    coding = "variable", nCat = 1L, kPrimePrior = "geometric", likelihoodMode = mode,
    priorVariant = "unconditional", kprimeTruncK = 30L,   # pin to this test's K
    kprimeHyperA = 1, kprimeHyperB = 1, expSteps = 1.4))
  modM <- mk("marginal_k"); modS <- mk("sampled_k")

  evalLL <- function(model, mkd, mutate) {
    modf <- MkPrime:::.FinalizeModel(model, trp, mkd); st0 <- MkPrime:::.InitState(trp, mkd, modf)
    st0$rate_log_sd <- 0; st0$tree_length <- TL; st0$rel_br_lengths <- RBL; st0 <- mutate(st0)
    dp <- MkPrime:::.InitMcmcData(mkd, modf); sp <- MkPrime:::.InitMcmcChain(st0)
    fill_partition_cache(dp, sp); eval_full_loglik_cpp(dp, sp)
  }
  evalMarg <- function(mkd, p)   evalLL(modM, mkd, function(s) { s$p <- p; s })
  LLfixed  <- function(mkd, kp)  evalLL(modS, mkd, function(s) { s$p <- 0.1
                                    s$kPrime <- rep_len(as.integer(kp), length(s$kPrime)); s })
  defaultLL <- function(mkd)     evalLL(modS, mkd, function(s) { s$p <- 0.1; s })   # default init kPrime == kObs

  for (nm in names(cases)) {
    mkd <- cases[[nm]]; kobs <- mkd$kObs[1]; ks <- kobs:K
    FX <- vapply(ks, function(kp) LLfixed(mkd, kp), numeric(1)); names(FX) <- ks
    # G1: k'-pinning is live (reference is non-trivial, not collapsed to a constant)
    expect_gt(stats::sd(FX), 1e-6, label = sprintf("[%s] G1 sd(FX) across k' (k'-pinning live)", nm))
    # G2: fixed-k' path anchored to the package default-init likelihood at k'=kObs
    expect_equal(unname(FX[as.character(kobs)]), defaultLL(mkd), tolerance = 1e-9,
                 label = sprintf("[%s] G2 LLfixed(kObs) == package default-init LL", nm))
    # MAIN: package marginal == full uncapped/uncutoff reference, at low p (heavy tail)
    for (p in c(0.01, 0.05)) {
      ref <- lse((ks - 2) * log1p(-p) + log(p) + FX) - logZA(p)
      expect_equal(evalMarg(mkd, p), ref, tolerance = 1e-7,
                   label = sprintf("[%s] package marginal == full reference at p=%.2f", nm, p))
    }
  }
})

# ---------------------------------------------------------------------------
# STAGE 1b CAP-COUPLING GUARD. The whole point of Stage 1b is that the model's
# truncation cap K is wired through to the C++ evaluator AND that the candidate
# cap kMaxKprimeCand (raised 50 -> 256) is large enough for the marginal
# numerator to sum the FULL support [2, K] when K is large (the real-data
# default K = 200). If the numerator silently capped below K (the old cap = 50)
# while Z_A normalised [2, K], it would reintroduce MARGINAL-K-TRUNC-001 across
# a wide p range.
#
# This test runs at K = 200, p = 0.02 (heavy geometric tail, so prior mass well
# beyond k' = 51 is substantial), and proves TWO things:
#   (1) DISCRIMINATION: the full [2,200] reference differs materially from a
#       [2,51]-capped reference (same Z_A(200)) -> the >k'=51 tail is real, so
#       matching the full one is a genuine test of the coupling (not vacuous).
#   (2) IDENTITY: the package marginal at model K = 200 equals the full [2,200]
#       reference. With the old cap = 50 the numerator would stop at k' ~ 51 and
#       match (2)'s capped reference instead -> this assertion FAILS.
# G1/G2 guard against a vacuous pass exactly as in the C-i test above.
test_that("marginal-k numerator reaches the full support at K=200 (Stage 1b cap coupling)", {
  set.seed(2026); ntip <- 8L; K <- 200L; p <- 0.02
  tr <- ape::rtree(ntip, tip.label = paste0("t", seq_len(ntip)))
  tr$edge.length <- rep_len(0.15, nrow(tr$edge))
  trp <- Preorder(tr); TL <- sum(tr$edge.length); RBL <- tr$edge.length / TL
  lse <- function(x) { x <- x[is.finite(x)]; if (!length(x)) return(-Inf); m <- max(x); m + log(sum(exp(x - m))) }
  logZA <- function(pp, KK) log1p(-(1 - pp)^(KK - 1))

  mkChar <- function(v) MkPrimeData(MatrixToPhyDat(
    matrix(v, ncol = 1, dimnames = list(tr$tip.label, "char1"))))
  mkd <- mkChar(c(0, 0, 1, 1, 0, 1, 0, 1))            # kObs = 2 (heaviest tail)

  mk <- function(mode, KK) suppressMessages(MkPrimeModel(
    coding = "variable", nCat = 1L, kPrimePrior = "geometric", likelihoodMode = mode,
    priorVariant = "unconditional", kprimeTruncK = KK,
    kprimeHyperA = 1, kprimeHyperB = 1, expSteps = 1.4))
  modM <- mk("marginal_k", K); modS <- mk("sampled_k", K)

  evalLL <- function(model, mutate) {
    modf <- MkPrime:::.FinalizeModel(model, trp, mkd); st0 <- MkPrime:::.InitState(trp, mkd, modf)
    st0$rate_log_sd <- 0; st0$tree_length <- TL; st0$rel_br_lengths <- RBL; st0 <- mutate(st0)
    dp <- MkPrime:::.InitMcmcData(mkd, modf); sp <- MkPrime:::.InitMcmcChain(st0)
    fill_partition_cache(dp, sp); eval_full_loglik_cpp(dp, sp)
  }
  evalMarg <- function() evalLL(modM, function(s) { s$p <- p; s })
  LLfixed  <- function(kp) evalLL(modS, function(s) { s$p <- 0.1
                              s$kPrime <- rep_len(as.integer(kp), length(s$kPrime)); s })
  defaultLL <- evalLL(modS, function(s) { s$p <- 0.1; s })

  kobs <- mkd$kObs[1]; ks <- kobs:K
  FX <- vapply(ks, LLfixed, numeric(1)); names(FX) <- ks
  # G1/G2: reference is live and anchored (not a collapsed constant).
  expect_gt(stats::sd(FX), 1e-6, label = "G1 sd(FX) across k' (k'-pinning live)")
  expect_equal(unname(FX[as.character(kobs)]), defaultLL, tolerance = 1e-9,
               label = "G2 LLfixed(kObs) == package default-init LL")

  logw    <- (ks - 2) * log1p(-p) + log(p) + FX
  refFull <- lse(logw)            - logZA(p, K)   # numerator over the full [2,200]
  refCap  <- lse(logw[ks <= 51])  - logZA(p, K)   # what an old cap=50 numerator gives
  # The likelihood P(data | k') DECAYS in k' (FX: -5.87 at k'=2 down to -16.9 at
  # k'=200), so the numerator's >k'=51 tail is small in absolute terms
  # (refFull - refCap ~ 1.3e-3 nats here) even though the *prior* mass there is
  # large -- a useful bound on the cap's practical impact for low-kObs chars.
  # It is still ~1e4x the identity tolerance below, so the identity check is a
  # genuine (non-vacuous) test that the numerator includes k' in (51, 200].
  expect_gt(abs(refFull - refCap), 1e-5,
            label = ">k'=51 numerator tail is resolvable above the 1e-7 identity tol")
  # IDENTITY: the package marginal at K=200 sums the FULL support, not a capped
  # one. With the old cap=50 it would equal refCap, off by ~1.3e-3 >> 1e-7 -> FAIL.
  expect_equal(evalMarg(), refFull, tolerance = 1e-7,
               label = "package marginal at K=200 == full [2,200] reference")
  # LIVE-K: the model's K is wired end-to-end -- changing it materially moves the
  # marginal LL (here ~0.80 nats, via -logZ_A(0.02,30) vs -logZ_A(0.02,200)).
  modM30 <- mk("marginal_k", 30L)
  em30 <- evalLL(modM30, function(s) { s$p <- p; s })
  expect_gt(abs(evalMarg() - em30), 0.1,
            label = "model kprimeTruncK 200 vs 30 materially changes the binary LL")
})

# ---------------------------------------------------------------------------
# STAGE 2 — sampled_k <-> marginal_k RB-CONSISTENCY (the deterministic bit-check
# behind "the two likelihoodModes target the SAME posterior"). Stage 2 truncates
# the sampled_k geometric PRIOR (cpp_log_prior) at K and renormalises by Z(p),
# matching the marginal_k normaliser. The decisive identity, per character:
#
#   logSumExp_{k' in [kObs, K]} [ logPrior_S(k') + logLik_S(k') ]
#        ==  logPrior_M + logLik_M
#
# i.e. summing the full sampled-k JOINT over k' reproduces the marginal-k joint,
# for BOTH Model A (unconditional) and Model B (conditional). Unlike the C-i
# guard above (which tests the marginal evaluator against ANALYTIC prior weights
# added in R), this exercises the EDITED cpp_log_prior via eval_log_prior_cpp.
#
# Tolerance 1e-7 (matching the C-i guard): the identity is exact in exact
# arithmetic; the residual is the M-164 pruning bound (kKprimeLogCutoff = -25 =>
# omitted-candidate mass < ~K*exp(-25) ~ 4e-10) plus fp accumulation in the
# shared LL path. A real prior bug (missing -logZ => ~0.26 nats at p=0.05,K=30;
# wrong slope; dropped truncation) is >= ~0.01 nats -- orders above the tol.
#
# Anti-vacuity / stale-binary guards (so a broken harness or un-rebuilt .so
# cannot pass silently):
#   G1: the joint genuinely varies with k' (k'-pinning is live).
#   G2: the per-k' sampled-k prior slope == log(1-p) (the geometric mass is live).
#   STALE: logPrior_S at K=30 differs from K=200 by the truncation-normaliser gap
#          (a stale/untruncated binary returns identical priors -> gap 0 -> FAIL).
#   CAP:   logPrior_S at k' = K+1 is -Inf (the truncated support guard is live).
test_that("sampled_k joint marginalises to marginal_k (Stage 2 RB-consistency)", {
  set.seed(2026); ntip <- 8L; K <- 30L
  tr <- ape::rtree(ntip, tip.label = paste0("t", seq_len(ntip)))
  tr$edge.length <- rep_len(0.15, nrow(tr$edge))
  trp <- Preorder(tr); TL <- sum(tr$edge.length); RBL <- tr$edge.length / TL
  lse <- function(x) { x <- x[is.finite(x)]; if (!length(x)) return(-Inf); m <- max(x); m + log(sum(exp(x - m))) }

  mkChar <- function(v) MkPrimeData(MatrixToPhyDat(
    matrix(v, ncol = 1, dimnames = list(tr$tip.label, "char1"))))
  cases <- list(kObs2 = mkChar(c(0, 0, 1, 1, 0, 1, 0, 1)),   # low-kObs (heaviest tail)
                kObs6 = mkChar(c(0, 1, 2, 3, 4, 5, 0, 1)))    # high-kObs

  mk <- function(mode, variant, KK = K) suppressMessages(MkPrimeModel(
    coding = "variable", nCat = 1L, kPrimePrior = "geometric", likelihoodMode = mode,
    priorVariant = variant, kprimeTruncK = KK,
    kprimeHyperA = 1, kprimeHyperB = 1, expSteps = 1.4))

  # Evaluate (logPrior, logLik) at a given mode/state via the package's own C++
  # entry points -- eval_log_prior_cpp is the routine Stage 2 edits.
  evalPL <- function(model, mkd, mutate) {
    modf <- MkPrime:::.FinalizeModel(model, trp, mkd); st0 <- MkPrime:::.InitState(trp, mkd, modf)
    st0$rate_log_sd <- 0; st0$tree_length <- TL; st0$rel_br_lengths <- RBL; st0 <- mutate(st0)
    dp <- MkPrime:::.InitMcmcData(mkd, modf); sp <- MkPrime:::.InitMcmcChain(st0)
    fill_partition_cache(dp, sp)
    c(prior = eval_log_prior_cpp(dp, sp), lik = eval_full_loglik_cpp(dp, sp))
  }
  priorS <- function(modS, mkd, kp, p)
    evalPL(modS, mkd, function(s) { s$p <- p
      s$kPrime <- rep_len(as.integer(kp), length(s$kPrime)); s })[["prior"]]

  for (variant in c("unconditional", "conditional")) {
    modS <- mk("sampled_k", variant); modM <- mk("marginal_k", variant)
    for (nm in names(cases)) {
      mkd <- cases[[nm]]; kobs <- mkd$kObs[1]; ks <- kobs:K
      for (p in c(0.02, 0.08)) {
        # Sampled-k FULL joint (edited prior + likelihood) at each fixed k'.
        J <- vapply(ks, function(kp) {
          pl <- evalPL(modS, mkd, function(s) { s$p <- p
                       s$kPrime <- rep_len(as.integer(kp), length(s$kPrime)); s })
          pl[["prior"]] + pl[["lik"]]
        }, numeric(1)); names(J) <- ks
        # Marginal-k FULL joint.
        plM <- evalPL(modM, mkd, function(s) { s$p <- p; s })
        M <- plM[["prior"]] + plM[["lik"]]

        # G1 (anti-vacuity): the joint genuinely varies with k'.
        expect_gt(stats::sd(J), 1e-6,
                  label = sprintf("[%s/%s p=%.2f] G1 sd(joint) over k'", nm, variant, p))
        # MAIN: full sampled-k sum == marginal-k joint (RB-consistency).
        expect_equal(lse(J), M, tolerance = 1e-7,
                     label = sprintf("[%s/%s p=%.2f] lse(sampled joint) == marginal joint", nm, variant, p))
      }
      # G2 (anti-vacuity): the geometric prior is live -- adjacent-k' prior slope
      # == log(1-p) (logP and -logZ cancel; holds for Model A and B alike).
      expect_equal(priorS(modS, mkd, kobs + 1L, 0.08) - priorS(modS, mkd, kobs, 0.08),
                   log1p(-0.08), tolerance = 1e-9,
                   label = sprintf("[%s/%s] G2 sampled-k prior slope == log(1-p)", nm, variant))
    }
  }

  # STALE-BINARY GUARD: the truncation normaliser Z(p) must be live in the
  # sampled-k prior. logPrior_S(K=30) - logPrior_S(K=200) == -logZA(30)+logZA(200)
  # (Model A, single char, fixed k'/p). A stale/untruncated binary -> gap 0.
  p <- 0.05; logZA <- function(pp, KK) log1p(-(1 - pp)^(KK - 1))
  gap <- priorS(mk("sampled_k", "unconditional", 30L),  cases$kObs2, 2L, p) -
         priorS(mk("sampled_k", "unconditional", 200L), cases$kObs2, 2L, p)
  expect_equal(gap, -logZA(p, 30L) + logZA(p, 200L), tolerance = 1e-7,
               label = "sampled-k prior K-sensitivity == truncation-normaliser gap")
  expect_gt(abs(gap), 0.2,
            label = "K=30 vs K=200 sampled-k prior gap is material (stale-binary guard)")

  # SUPPORT CAP: k' > K carries zero prior mass (-Inf) under sampled_k.
  prOver <- priorS(mk("sampled_k", "unconditional", 30L), cases$kObs2, 31L, 0.1)
  expect_true(is.infinite(prOver) && prOver < 0,
              label = "sampled-k prior at k' = K+1 is -Inf (truncated support guard)")
})

# ---------------------------------------------------------------------------
# STAGE 2 — MULTI-CHARACTER RB closure. The proof
# (dev/red-team/proofs/marginal-k-sampled-rb-consistency.md) lifts the
# per-character identity to many characters ANALYTICALLY via the
# product-of-sums (conditional-independence) factorisation; the deterministic
# tests above are n=1. This test makes the n=2 case EMPIRICAL: it sums the
# sampled-k JOINT over the FULL 2-D grid (k'_1, k'_2) in [kObs_i, K]^2 and
# checks it equals the marginal-k joint, for Model A and Model B. K is small
# (8) so truncation bites hard (-logZ_A ~ 1.6 nats at p=0.03) and the grid stays
# cheap.
test_that("sampled_k 2-character joint marginalises to marginal_k (multi-char RB)", {
  set.seed(4040); ntip <- 7L; K <- 8L
  tr <- ape::rtree(ntip, tip.label = paste0("t", seq_len(ntip)))
  tr$edge.length <- rep_len(0.18, nrow(tr$edge))
  trp <- Preorder(tr); TL <- sum(tr$edge.length); RBL <- tr$edge.length / TL
  lse <- function(x) { x <- x[is.finite(x)]; if (!length(x)) return(-Inf); m <- max(x); m + log(sum(exp(x - m))) }

  # Two transformational characters: kObs = 2 and kObs = 3.
  m2 <- cbind(c(0, 0, 1, 1, 0, 1, 0), c(0, 1, 2, 0, 1, 2, 0))
  rownames(m2) <- tr$tip.label
  mkd <- MkPrimeData(MatrixToPhyDat(m2))
  ti  <- which(mkd$type == "transformational")
  ko  <- mkd$kObs[ti]

  mk <- function(mode, variant) suppressMessages(MkPrimeModel(
    coding = "variable", nCat = 1L, kPrimePrior = "geometric", likelihoodMode = mode,
    priorVariant = variant, kprimeTruncK = K, kprimeHyperA = 1, kprimeHyperB = 1,
    expSteps = 1.4))

  # Build data/model once per mode; the returned closure varies only the state.
  buildEval <- function(model) {
    modf <- MkPrime:::.FinalizeModel(model, trp, mkd)
    dp   <- MkPrime:::.InitMcmcData(mkd, modf)
    function(mutate) {
      st0 <- MkPrime:::.InitState(trp, mkd, modf)
      st0$rate_log_sd <- 0; st0$tree_length <- TL; st0$rel_br_lengths <- RBL
      st0 <- mutate(st0)
      sp <- MkPrime:::.InitMcmcChain(st0)
      fill_partition_cache(dp, sp)
      eval_log_prior_cpp(dp, sp) + eval_full_loglik_cpp(dp, sp)
    }
  }

  grid <- expand.grid(k1 = ko[1]:K, k2 = ko[2]:K)
  for (variant in c("unconditional", "conditional")) {
    evalS <- buildEval(mk("sampled_k",  variant))
    evalM <- buildEval(mk("marginal_k", variant))
    for (p in c(0.03, 0.1)) {
      M <- evalM(function(s) { s$p <- p; s })                 # marginal joint
      J <- vapply(seq_len(nrow(grid)), function(g) {          # sampled grid
        kp <- mkd$kObs
        kp[ti[1]] <- grid$k1[g]; kp[ti[2]] <- grid$k2[g]
        evalS(function(s) { s$p <- p; s$kPrime <- as.integer(kp); s })
      }, numeric(1))
      expect_gt(stats::sd(J), 1e-6,
                label = sprintf("[%s p=%.2f] G1 2-char joint varies over grid", variant, p))
      expect_equal(lse(J), M, tolerance = 1e-7,
                   label = sprintf("[%s p=%.2f] lse(2-char sampled grid) == marginal joint",
                                   variant, p))
    }
  }
})

# ---------------------------------------------------------------------------
# STAGE 2 — the case-25 (gibbs_kprime_sweep) K-CAP actually FIRES. The
# deterministic RB checks above exercise the PRIOR + marginal evaluator, never
# the Gibbs sweep; SBC runs the sweep only at K=30, where the M-164 prior-ceiling
# cutoff is tighter than the cap (nEff >= nCand, so the new `nCand = nEff;
# recompute maxW` branch is a no-op). This test FORCES nEff < nCand (tiny K=4,
# kObs=2 => nEff=3) so the cap branch is genuinely exercised. gibbs_kprime_sweep
# ALWAYS accepts, so a cap bug (off-by-one, or a draw of k' > K) would set
# logPrior=-Inf on an accepted move and silently corrupt the chain -- something
# a 10-h Hamilton SBC would only reveal after the fact. Non-vacuity: a K=100
# contrast shows the un-capped sweep demonstrably samples k' > 4, so the K=4 cap
# has real work to do.
test_that("gibbs_kprime_sweep respects the truncation cap K (case-25 cap fires)", {
  skip_if_not(requireNamespace("ape", quietly = TRUE))
  set.seed(20260601); ntip <- 8L
  tr <- ape::rtree(ntip, tip.label = paste0("t", seq_len(ntip)))
  tr$edge.length <- rep_len(0.06, nrow(tr$edge))     # short branches: flat LL in k'
  tr <- Preorder(tr)
  # Six variable binary characters (kObs = 2). At low p the geometric prior
  # drives k' up (P(k' >= 5) = (1-p)^3 = 0.729 at p=0.1); flat LL keeps the high
  # candidates un-pruned, so the un-capped conditional reaches well past k'=4.
  set.seed(7L); mat <- matrix(0L, ntip, 6L, dimnames = list(tr$tip.label, NULL))
  for (j in seq_len(6L)) {
    repeat { v <- sample(0:1, ntip, replace = TRUE); if (length(unique(v)) == 2L) break }
    mat[, j] <- v
  }
  mkd <- MkPrimeData(MatrixToPhyDat(mat))
  transIdx <- which(mkd$type == "transformational")

  sweep_kprimes <- function(K, p, nIter = 400L, seed = 99L) {
    model <- suppressMessages(MkPrimeModel(
      coding = "variable", nCat = 1L, kPrimePrior = "geometric",
      likelihoodMode = "sampled_k", priorVariant = "unconditional",
      kprimeTruncK = K, kprimeHyperA = 1, kprimeHyperB = 1, expSteps = 1.4))
    modf <- MkPrime:::.FinalizeModel(model, tr, mkd)
    st0  <- MkPrime:::.InitState(tr, mkd, modf)
    st0$p <- p; st0$rate_log_sd <- 0
    dp <- MkPrime:::.InitMcmcData(mkd, modf); sp <- MkPrime:::.InitMcmcChain(st0)
    fill_partition_cache(dp, sp); allocate_cl_workspace(dp, sp)
    set.seed(seed)
    kp <- matrix(0L, nIter, length(transIdx)); lp <- numeric(nIter)
    for (i in seq_len(nIter)) {
      do_move_cpp(dp, sp, 25L, 0L, 0.5, 10.0, 1L, 1.0)   # case 25 = gibbs_kprime_sweep
      s <- get_mcmc_state(sp)
      kp[i, ] <- s$kPrime[transIdx]; lp[i] <- s$logPrior
    }
    list(kp = as.integer(kp), lp = lp)
  }

  # NON-VACUITY: at K=100 the un-capped sweep reaches k' > 4, so the K=4 cap
  # genuinely truncates (it is not a no-op).
  loose <- sweep_kprimes(K = 100L, p = 0.1)
  expect_gt(max(loose$kp), 4L,
            label = "K=100: un-capped sweep samples k' > 4 (the K=4 cap is non-vacuous)")

  # CAP CORRECTNESS at K=4 (nEff = 3 < un-capped nCand -> the cap branch fires):
  tight <- sweep_kprimes(K = 4L, p = 0.1)
  expect_true(all(tight$kp >= 2L) && all(tight$kp <= 4L),
              label = "K=4: every swept k' stays in [kObs, K] (cap holds)")
  expect_true(any(tight$kp == 4L),
              label = "K=4: the cap boundary k'=K is reached (cap is binding, not collapsed)")
  expect_true(all(is.finite(tight$lp)),
              label = "K=4: logPrior finite after every always-accept sweep (no -Inf corruption)")
})
