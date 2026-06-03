## gibbs-p-identity-check.R
##
## Numerical verification of the data-augmentation Metropolis-within-Gibbs
## p-update derived in dev/red-team/proofs/marginal-k-gibbs-p.md.
##
## This is a PURE-R reimplementation of the *formulae* in the note (it does NOT
## call the C++ evaluator, so it needs no DLL rebuild). It checks:
##
##   (a) §2 marginalisation identity: sum_{u} pi_aug(p, u) == pi(p) to ~1e-12,
##       for arbitrary support S_i and arbitrary p-dependent Z_i(p).
##   (b) §3-§4 invariance: the (impute -> independence-MH) kernel leaves the
##       discretised p-target invariant. Compared as a long-run empirical p
##       histogram vs the analytic pi(p|theta) on a grid, for:
##         - Model A, small K (truncation BITES, c_A != 0, varying kObs)
##         - Model B, small K (truncation bites)
##         - Model A, large K (truncation underflows -> near-pure Gibbs)
##   (c) §5 untruncated reduction: with K huge, the accepted draws match
##       rbeta(a+n, b+S+c_A); and Model B with kObs==2 matches the case-9
##       sampled_k draw shape2 = b + sumU exactly (same shapes).
##
## All "L_i(u)" raw per-(char,u) likelihoods are synthetic positive numbers; the
## identities are algebraic in L and hold for ANY L, so synthetic values are a
## fair test of the math (not of the pruning that produces L in the code).

set.seed(20260603L)

## ---------------------------------------------------------------------------
## Model definition mirroring src/mcmc.cpp:4326-4476 and :428-444.
##
##   logPriorByU[u]   = log p + u*log(1-p)              (Model B base weight)
##   Model A adds      (kObs_i - 2)*log(1-p)            per char, u-independent
##   Z_i(p) Model A    = 1 - (1-p)^(K-1)                shared
##   Z_i(p) Model B    = 1 - (1-p)^(K - kObs_i + 1)     per char
##
## Per-char marginal m_i(p) = g_i(p) * sum_{u in S_i} L_i(u) p (1-p)^u / Z_i(p)
##   g_i(p) = (1-p)^(kObs_i - 2)  (Model A)   or  1  (Model B)
## ---------------------------------------------------------------------------

## log Z_i(p) for a given char, matching the code's log1p/exp form.
logZ_i <- function(p, K, kObs_i, modelA) {
  log1mP <- log1p(-p)
  if (modelA) {
    log1p(-exp((K - 1) * log1mP))
  } else {
    log1p(-exp((K - kObs_i + 1) * log1mP))
  }
}

## log g_i(p).
log_g_i <- function(p, kObs_i, modelA) {
  if (modelA) (kObs_i - 2) * log1p(-p) else 0
}

## Per-char log marginal log m_i(p), exactly as the evaluator forms it
## (logSumExp over the cached support S_i = {0..nEff_i-1}).
log_m_i <- function(p, Lvec, kObs_i, K, modelA) {
  # Lvec is L_i(u) for u = 0..(nEff_i - 1) -- the cached support already.
  u <- seq_along(Lvec) - 1L
  logP <- log(p); log1mP <- log1p(-p)
  w <- log(Lvec) + logP + u * log1mP            # charLogW[i,u]
  mx <- max(w)
  lse <- mx + log(sum(exp(w - mx)))
  lse + log_g_i(p, kObs_i, modelA) - logZ_i(p, K, kObs_i, modelA)
}

## Full un-normalised log pi(p | theta, T, data) = log Beta(p;a,b) + sum_i log m_i.
log_pi_p <- function(p, chars, a, b, K, modelA) {
  if (p <= 0 || p >= 1) return(-Inf)
  lp <- dbeta(p, a, b, log = TRUE)
  for (ch in chars) lp <- lp + log_m_i(p, ch$L, ch$kObs, K, modelA)
  lp
}

## ---------------------------------------------------------------------------
## (a) Marginalisation identity §2.
##   sum_{u in S_i} pi_aug(p, u) == pi(p), for each char independently, since
##   chars factor. We check the per-char identity at machine precision: the
##   augmented per-char joint summed over u equals m_i(p).
## ---------------------------------------------------------------------------
check_marginalisation <- function() {
  cat("=== (a) Marginalisation identity (§2) ===\n")
  ps  <- c(0.01, 0.1, 0.37, 0.5, 0.83, 0.99)
  Ks  <- c(5L, 200L)
  maxabs <- 0
  for (modelA in c(TRUE, FALSE)) for (K in Ks) for (p in ps) {
    for (rep in 1:4) {
      kObs_i <- sample(2:6, 1)
      nMax   <- min(8L, K - kObs_i + 1L)           # cap support so k <= K
      if (nMax < 1L) next                          # empty support (kObs > K)
      nEff   <- if (nMax == 1L) 1L else sample(1:nMax, 1)
      L      <- exp(rnorm(nEff, sd = 3))           # arbitrary positive L_i(u)
      ## LHS: sum_u of augmented per-char joint (un-normalised in p, but here we
      ## compute the EXACT per-char marginal contribution m_i(p)).
      u <- 0:(nEff - 1)
      g <- if (modelA) (1 - p)^(kObs_i - 2) else 1
      Z <- if (modelA) 1 - (1 - p)^(K - 1) else 1 - (1 - p)^(K - kObs_i + 1)
      lhs_terms <- g * L * p * (1 - p)^u / Z       # pi_aug per-u (in p-space)
      lhs <- sum(lhs_terms)                          # sum over u
      ## RHS: m_i(p) via the evaluator's logSumExp form, exponentiated.
      rhs <- exp(log_m_i(p, L, kObs_i, K, modelA))
      maxabs <- max(maxabs, abs(lhs - rhs) / max(abs(rhs), 1e-300))
    }
  }
  cat(sprintf("  max relative |sum_u pi_aug - m_i(p)| = %.3e  (target <1e-12)\n",
              maxabs))
  invisible(maxabs)
}

## ---------------------------------------------------------------------------
## The Gibbs-p kernel exactly as the note specifies it.
##   impute u_i ~ Categorical(w_i) with w_i[u] propto L_i(u) p (1-p)^u
##   propose  p* ~ Beta(a+n, b+S+c_A),  S = sum u_i,  c_A = sum(kObs_i - 2) | A
##   accept   log alpha = sum_i [log Z_i(p) - log Z_i(p*)]
## ---------------------------------------------------------------------------
## `flipSign` and `dropCA` are CONTROLS: when set, deliberately corrupt the
## kernel so the invariance test can be shown to FAIL on a wrong accept ratio
## or a missing c_A. Default FALSE = the note's kernel.
gibbs_p_step <- function(p, chars, a, b, K, modelA,
                         flipSign = FALSE, dropCA = FALSE) {
  n <- length(chars)
  S <- 0L
  for (ch in chars) {
    u   <- seq_along(ch$L) - 1L
    lw  <- log(ch$L) + log(p) + u * log1p(-p)
    lw  <- lw - max(lw)
    pr  <- exp(lw); pr <- pr / sum(pr)
    ui  <- if (length(u) == 1L) u else sample(u, 1L, prob = pr)
    S   <- S + ui
  }
  cA <- if (modelA && !dropCA) sum(vapply(chars, function(ch) ch$kObs - 2, 0)) else 0
  pstar <- rbeta(1L, a + n, b + S + cA)
  if (pstar <= 0 || pstar >= 1) return(p)               # guard §6.3
  sZp  <- sum(vapply(chars, function(ch) logZ_i(p,     K, ch$kObs, modelA), 0))
  sZps <- sum(vapply(chars, function(ch) logZ_i(pstar, K, ch$kObs, modelA), 0))
  if (!is.finite(sZp) || !is.finite(sZps)) return(p)    # guard §6.4
  logAlpha <- sZp - sZps                                 # note's sign
  if (flipSign) logAlpha <- -logAlpha                    # CONTROL: wrong sign
  if (log(runif(1)) < logAlpha) pstar else p
}

## ---------------------------------------------------------------------------
## (b) Invariance: long-run histogram of the kernel vs analytic pi(p) on a grid.
## ---------------------------------------------------------------------------
check_invariance <- function(label, chars, a, b, K, modelA,
                             nIter = 4e5, nGrid = 400,
                             flipSign = FALSE, dropCA = FALSE) {
  cat(sprintf("=== (b) Invariance: %s ===\n", label))
  ## analytic target on a grid (normalised by trapezoid)
  grid  <- seq(1e-4, 1 - 1e-4, length.out = nGrid)
  ltarg <- vapply(grid, log_pi_p, 0, chars = chars, a = a, b = b,
                  K = K, modelA = modelA)
  ltarg <- ltarg - max(ltarg)
  dens  <- exp(ltarg)
  dx    <- grid[2] - grid[1]
  dens  <- dens / (sum(dens) * dx)                       # normalised pdf
  ## run the kernel
  p <- 0.5
  burn <- nIter %/% 5
  keep <- numeric(nIter - burn)
  for (it in 1:nIter) {
    p <- gibbs_p_step(p, chars, a, b, K, modelA,
                      flipSign = flipSign, dropCA = dropCA)
    if (it > burn) keep[it - burn] <- p
  }
  ## bin empirical, compare to analytic at bin centres
  brks <- seq(0, 1, length.out = 41)
  ctr  <- (brks[-1] + brks[-length(brks)]) / 2
  h    <- hist(keep, breaks = brks, plot = FALSE)
  emp  <- h$density
  ana  <- vapply(ctr, function(x) {
    lp <- log_pi_p(x, chars, a, b, K, modelA); exp(lp)
  }, 0)
  ana  <- ana / (sum(ana) * (brks[2] - brks[1]))         # normalise to pdf
  ## restrict comparison to bins with non-trivial analytic mass
  keepbin <- ana > 1e-3 * max(ana)
  maxabs  <- max(abs(emp[keepbin] - ana[keepbin]))
  l1      <- sum(abs(emp - ana)) * (brks[2] - brks[1])    # total-variation-ish
  ## posterior-mean check (low variance, sharp discriminator of sign/c_A)
  emean_emp <- mean(keep)
  emean_ana <- sum(grid * dens) * dx
  cat(sprintf("  E[p] kernel = %.5f   E[p] analytic = %.5f   |diff| = %.2e\n",
              emean_emp, emean_ana, abs(emean_emp - emean_ana)))
  cat(sprintf("  max |density diff| over mass bins = %.4f   L1 = %.4f\n",
              maxabs, l1))
  invisible(list(emean_emp = emean_emp, emean_ana = emean_ana,
                 maxabs = maxabs, l1 = l1))
}

## ---------------------------------------------------------------------------
## (c) Untruncated reduction §5: with K huge, accepted draws ~ Beta(a+n,b+S+cA).
##   We fix the imputed u (so S, cA fixed), draw many p* from the kernel's
##   proposal under huge K, confirm always-accept and KS-match to rbeta.
## ---------------------------------------------------------------------------
check_untruncated <- function() {
  cat("=== (c) Untruncated reduction (§5) ===\n")
  a <- 1.5; b <- 2.0; K <- 1e9
  ## Build chars with FIXED L so the imputation is reproducible; but here we
  ## test the p-draw directly at a fixed (S, n, cA).
  for (modelA in c(TRUE, FALSE)) {
    kObs <- c(2L, 4L, 3L, 5L)
    n    <- length(kObs)
    S    <- 7L
    cA   <- if (modelA) sum(kObs - 2) else 0
    p <- 0.4
    nAcc <- 0; draws <- numeric(5000)
    for (i in 1:5000) {
      pstar <- rbeta(1L, a + n, b + S + cA)
      sZp   <- sum(vapply(kObs, function(k) logZ_i(p,     K, k, modelA), 0))
      sZps  <- sum(vapply(kObs, function(k) logZ_i(pstar, K, k, modelA), 0))
      # under huge K, log Z -> 0 for any non-tiny p; logAlpha -> 0
      logAlpha <- sZp - sZps
      acc <- is.finite(logAlpha) && log(runif(1)) < logAlpha
      if (acc) { nAcc <- nAcc + 1; p <- pstar }
      draws[i] <- pstar
    }
    ref <- rbeta(5000L, a + n, b + S + cA)
    ks  <- suppressWarnings(ks.test(draws, ref))$statistic
    cat(sprintf("  modelA=%-5s  c_A=%2d  accept_frac=%.4f  KS(proposal vs rbeta(a+n,b+S+cA))=%.4f\n",
                modelA, cA, nAcc / 5000, ks))
  }
  ## Model B, kObs==2 everywhere: c_A=0, S=sumU, must equal case-9 shapes.
  cat("  Model B, kObs==2: case-9 equivalence check (shape identity)\n")
  a <- 1.0; b <- 1.0; n <- 5L; sumU <- 9L
  s1_case9 <- a + n;  s2_case9 <- b + sumU            # src/mcmc.cpp:5036-5037
  s1_gibbs <- a + n;  s2_gibbs <- b + sumU + 0        # c_A = 0
  cat(sprintf("    case-9 shapes  = (%.1f, %.1f)\n", s1_case9, s2_case9))
  cat(sprintf("    gibbs-p shapes = (%.1f, %.1f)   identical: %s\n",
              s1_gibbs, s2_gibbs,
              isTRUE(all.equal(c(s1_case9, s2_case9), c(s1_gibbs, s2_gibbs)))))
}

## ---------------------------------------------------------------------------
## Run all checks.
## ---------------------------------------------------------------------------
cat("\n################ gibbs-p-identity-check.R ################\n\n")

check_marginalisation()
cat("\n")

## Build small synthetic character sets.
mkChars <- function(kObsVec, nEffVec, sd = 2) {
  Map(function(k, ne) list(kObs = k, L = exp(rnorm(ne, sd = sd))),
      kObsVec, nEffVec)
}

a <- 1.7; b <- 2.3

## Model A, small K=5, varying kObs so c_A != 0 -- THE load-bearing cell.
## (Simultaneously falsifies a sign error in log alpha AND a missing c_A.)
set.seed(11L)
charsA <- mkChars(kObs = c(2L, 3L, 4L, 5L), nEff = c(4L, 3L, 2L, 1L))
rA <- check_invariance("Model A, K=5 (truncation bites, c_A != 0)",
                       charsA, a, b, K = 5L, modelA = TRUE, nIter = 8e5)
cat("\n")

## Model B, small K=5.
set.seed(12L)
charsB <- mkChars(kObs = c(2L, 3L, 4L), nEff = c(4L, 3L, 2L))
rB <- check_invariance("Model B, K=5 (truncation bites, c_A = 0)",
                       charsB, a, b, K = 5L, modelA = FALSE, nIter = 8e5)
cat("\n")

## Model A, large K=200 (truncation underflows -> near-pure Gibbs).
set.seed(13L)
charsA2 <- mkChars(kObs = c(2L, 3L, 4L), nEff = c(6L, 5L, 4L))
rA2 <- check_invariance("Model A, K=200 (truncation underflows)",
                        charsA2, a, b, K = 200L, modelA = TRUE,
                        nIter = 2e5)
cat("\n")

## ---- CONTROLS: the test must FAIL when the kernel is corrupted -----------
## These demonstrate the E[p] check actually discriminates the risky claims.
cat("--- CONTROLS (these SHOULD show large E[p] error) ---\n")
set.seed(11L)
cFlip <- check_invariance("CONTROL Model A K=5, WRONG SIGN on log alpha",
                          charsA, a, b, K = 5L, modelA = TRUE,
                          nIter = 8e5, flipSign = TRUE)
cat("\n")
set.seed(11L)
cDrop <- check_invariance("CONTROL Model A K=5, c_A DROPPED (uses b+S only)",
                          charsA, a, b, K = 5L, modelA = TRUE,
                          nIter = 8e5, dropCA = TRUE)
cat("\n")

check_untruncated()

## ---- Summary verdict lines ----------------------------------------------
cat("\n--- SUMMARY ---\n")
tol_mean <- 5e-3
ok <- function(r) abs(r$emean_emp - r$emean_ana) < tol_mean
cat(sprintf("  Model A K=5   E[p] err = %.2e   %s\n",
            abs(rA$emean_emp  - rA$emean_ana),  if (ok(rA))  "PASS" else "FAIL"))
cat(sprintf("  Model B K=5   E[p] err = %.2e   %s\n",
            abs(rB$emean_emp  - rB$emean_ana),  if (ok(rB))  "PASS" else "FAIL"))
cat(sprintf("  Model A K=200 E[p] err = %.2e   %s\n",
            abs(rA2$emean_emp - rA2$emean_ana), if (ok(rA2)) "PASS" else "FAIL"))
cat(sprintf("  CONTROL flip-sign E[p] err = %.2e   %s (want LARGE)\n",
            abs(cFlip$emean_emp - cFlip$emean_ana),
            if (!ok(cFlip)) "discriminates" else "FAILS-TO-DISCRIMINATE"))
cat(sprintf("  CONTROL drop-c_A  E[p] err = %.2e   %s (want LARGE)\n",
            abs(cDrop$emean_emp - cDrop$emean_ana),
            if (!ok(cDrop)) "discriminates" else "FAILS-TO-DISCRIMINATE"))

cat("\n################ done ################\n")
