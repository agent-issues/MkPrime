# Red-team: shape of per-character k' posterior under Mk' priors
# ----------------------------------------------------------------
# Diagnostic for the Mk'-vs-mk_k40 tree-recovery gap. We test, on real
# binary characters from tree_01/rep_01:
#   (1) Prior PDFs over k' under geometric / EG / Beta(1,20) hyperpriors.
#   (2) Marginal likelihood profile L(k') with a fixed (true) tree &
#       reasonable branch lengths -- is the likelihood flat for k' > kObs?
#   (3) Implied posterior = prior x likelihood at a few branch-length
#       regimes.  Where does the mode sit?  Long right tail?
#
# Output: redteam_kprime_*.png + a tabular RDS into the analysis dir.

suppressPackageStartupMessages({
  library(MkPrime)
  library(TreeTools)
  library(ape)
})

OUT_DIR <- "dev/pilots/2026-05-12-prior-validation/analysis"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

K_GRID <- 2:50  # candidate k' values (Gibbs sweep is K_MAX_CAND=50)

# ---------------------------------------------------------------- prior PDFs --
# Geometric Beta hyperprior:   k' = kObs + u,  u ~ Geom(p),  p ~ Beta(a,b)
#   P(u | a, b) = Beta(a+1, b+u) / Beta(a, b)        (marginalised over p)
# This is *exactly* the beta_geometric prior in MkPrimeModel.R:498-504,
# and the "geometric" prior with the Beta hyperprior matches when we
# marginalise.
log_prior_beta_geom <- function(k, kObs, a, b) {
  u <- k - kObs
  u[u < 0] <- NA_real_
  out <- lbeta(a + 1, b + u) - lbeta(a, b)
  out[is.na(u)] <- -Inf
  out
}

# Empirical_geometric: k' = N_obs + N_unobs, with N_obs ~ empirical pmf,
# N_unobs ~ Geom(p). Marginal P(k'=K | kObs, a, b) is the convolution
# truncated to N_obs >= kObs (since k' >= kObs).
log_prior_emp_geom <- function(K, kObs, emp, a, b) {
  # Build a P(N_obs = n) lookup for n >= kObs.
  emp_body <- emp$body  # names "2"..."16"
  n_max_body <- as.integer(names(emp_body)[length(emp_body)])
  # Tail: P(N_obs = n) = tail_start_p * decay^(n - tail_start_k) for n >= tail_start_k
  # (linear interpretation; geometric tail).
  pn <- function(n) {
    if (n < 2L) return(0)
    if (n <= n_max_body) {
      v <- emp_body[as.character(n)]
      if (is.na(v)) return(0) else return(unname(v))
    }
    if (n >= emp$tail_start_k) {
      return(unname(emp$tail_start_p) * emp$tail_decay ^ (n - emp$tail_start_k))
    }
    0
  }
  # P(N_unobs = u) marginalised over p ~ Beta(a,b) = lbeta(a+1, b+u) / lbeta(a,b).
  log_pu <- function(u) lbeta(a + 1, b + u) - lbeta(a, b)
  # For each K, P(k' = K | kObs) ∝ sum_{n >= kObs, n <= K} P(N_obs = n) * P(N_unobs = K - n).
  # (truncation: kObs <= n <= K.)
  out <- vapply(K, function(K_i) {
    if (K_i < kObs) return(-Inf)
    ns <- seq.int(kObs, K_i)
    log_terms <- vapply(ns, function(n) {
      pp <- pn(n)
      if (pp <= 0) return(-Inf)
      log(pp) + log_pu(K_i - n)
    }, numeric(1))
    matrixStats <- NULL
    # logsumexp
    m <- max(log_terms)
    if (!is.finite(m)) return(-Inf)
    m + log(sum(exp(log_terms - m)))
  }, numeric(1))
  # Normalise over the K grid so we can read it as a PMF on the displayed range.
  m <- max(out)
  norm <- m + log(sum(exp(out - m)))
  out - norm
}

# ----------------------------------------------------------------- load data --
data(empiricalNObs, package = "MkPrime")
emp <- empiricalNObs

sim_dir <- file.path("C:/Users/pjjg18/GitHub/mkprime/tree-inference/tree_01/rep_01/sim")
fs <- list.files(sim_dir, pattern = "[.]nex$", full.names = TRUE)
cat("Reading", length(fs), "nexus files...\n")
# Each ReadAsPhyDat returns a 1-column phyDat; convert via PhyDatToMatrix
# (preserves tip labels), cbind into a wide matrix, MatrixToPhyDat back.
mats <- lapply(fs, function(f) {
  pd_i <- TreeTools::ReadAsPhyDat(f)
  TreeTools::PhyDatToMatrix(pd_i)
})
# Each is a single-column matrix [nTip x 1] with the same rownames.
charMat <- do.call(cbind, mats)
pd <- TreeTools::MatrixToPhyDat(charMat)
mkd <- MkPrimeData(pd)
cat("MkPrimeData: nTip =", mkd$nTip, " nChar =", mkd$nChar, "\n")
cat("kObs distribution:\n"); print(table(mkd$kObs))

tr <- read.tree("C:/Users/pjjg18/GitHub/mkprime/tree-inference/tree_01/tree.nwk")
tr$edge.length <- pmax(tr$edge.length, 1e-8)
cat("Tree: ", Ntip(tr), " tips, total length =",
    format(sum(tr$edge.length), digits = 3), "\n")

# Make a copy with TL = 1 (typical posterior median TL is roughly in this range)
tr_med <- tr; tr_med$edge.length <- tr$edge.length / sum(tr$edge.length) * 1.0
tr_short <- tr; tr_short$edge.length <- tr$edge.length / sum(tr$edge.length) * 0.3
tr_long  <- tr; tr_long$edge.length  <- tr$edge.length / sum(tr$edge.length) * 3.0

# --------------------------------------------------------- likelihood profile --
# For each binary character we evaluate L(k') on tr_med, tr_short, tr_long.
# We construct a single-character MkPrimeData (so other chars don't dilute the
# vector) and call MkpLogLikelihood with a fixed kPrime vector of length 1.

single_char_mkd <- function(mkd, char_i) {
  # mkd$matrix is the post-invariant-drop integer matrix [nTip x nChar].
  # Recode to character symbols ("0", "1", ..., "?") and rebuild via
  # MatrixToPhyDat -- preserves the actual data of the requested column.
  v_int <- mkd$matrix[, char_i]
  syms <- as.character(v_int)
  syms[is.na(v_int)] <- "?"
  m1 <- matrix(syms, ncol = 1L,
               dimnames = list(rownames(mkd$matrix), NULL))
  pd1 <- TreeTools::MatrixToPhyDat(m1)
  MkPrimeData(pd1)
}

profile_one_char <- function(mkd1, tr_local, K_grid, relabel = TRUE) {
  vapply(K_grid, function(K) {
    ll <- tryCatch(
      MkpLogLikelihood(tr_local, mkd1, kPrime = K,
                       coding = "variable", relabel = relabel),
      error = function(e) NA_real_
    )
    ll
  }, numeric(1))
}

# Limit to 12 binary characters for clarity (mix of patterns).
set.seed(1)
bin_idx_all <- which(mkd$kObs == 2L)
bin_idx <- sort(sample(bin_idx_all, min(12, length(bin_idx_all))))
cat("Profiling", length(bin_idx), "binary characters at 3 tree-length regimes\n")

prof    <- list()  # full Mk' likelihood (pruning+ascertainment+relabel)
prof_nr <- list()  # no relabel: pruning+ascertainment only (= "mk_kN-style")
for (regime in c("short", "med", "long")) {
  tr_use <- switch(regime, short = tr_short, med = tr_med, long = tr_long)
  mat   <- matrix(NA_real_, nrow = length(K_GRID), ncol = length(bin_idx),
                  dimnames = list(as.character(K_GRID), as.character(bin_idx)))
  matnr <- matrix(NA_real_, nrow = length(K_GRID), ncol = length(bin_idx),
                  dimnames = list(as.character(K_GRID), as.character(bin_idx)))
  for (j in seq_along(bin_idx)) {
    mkd1 <- single_char_mkd(mkd, bin_idx[j])
    mat[, j]   <- profile_one_char(mkd1, tr_use, K_GRID, relabel = TRUE)
    matnr[, j] <- profile_one_char(mkd1, tr_use, K_GRID, relabel = FALSE)
  }
  prof[[regime]]    <- mat
  prof_nr[[regime]] <- matnr
}

# Decomposition table at medium TL
relabel_term <- vapply(K_GRID, function(K)
  MkPrime:::mk_prime_relabel_log(as.integer(K), 2L), numeric(1))
cat("\n=== Likelihood decomposition at medium TL ===\n")
cat("Per-char slope from k=2 to k=10 -- WITH relabel (Mk') vs NO relabel (mk_kN-style):\n")
dec <- data.frame(
  char     = bin_idx,
  slope_with    = round(prof$med   ["10",] - prof$med   ["2",], 2),
  slope_pruning = round(prof_nr$med["10",] - prof_nr$med["2",], 2),
  relabel_delta = round(relabel_term[K_GRID == 10] - relabel_term[K_GRID == 2], 2)
)
print(dec)

# ------------------------------------------------------ priors over the grid --
# Geometric Beta(1,1) (mkp_geo default)
lp_geo  <- log_prior_beta_geom(K_GRID, kObs = 2, a = 1, b = 1)
# Geometric Beta(1,20) (mkp_highk: shrinks p -> 0, encourages large u)
lp_high <- log_prior_beta_geom(K_GRID, kObs = 2, a = 1, b = 20)
# Empirical_geometric Beta(1,1) (mkp_eg default)
lp_eg   <- log_prior_emp_geom(K_GRID, kObs = 2, emp = emp, a = 1, b = 1)

# Normalise each over K_GRID for plotting as a PMF.
norm_logpmf <- function(lp) {
  m <- max(lp); exp(lp - m - log(sum(exp(lp - m))))
}
pmf_geo  <- norm_logpmf(lp_geo)
pmf_high <- norm_logpmf(lp_high)
pmf_eg   <- norm_logpmf(lp_eg)

# ------------------------------------------------------------------- plot 1 --
png(file.path(OUT_DIR, "redteam_kprime_priors.png"),
    width = 1100, height = 500, res = 110)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
matplot(K_GRID, cbind(pmf_geo, pmf_high, pmf_eg), type = "l", lty = 1, lwd = 2,
        col = c("tomato", "darkgreen", "steelblue"),
        xlab = "k'", ylab = "P(k' | kObs=2)",
        main = "Prior PMFs (kObs = 2), linear")
legend("topright", c("geom Beta(1,1)  [mkp_geo]",
                     "geom Beta(1,20) [mkp_highk]",
                     "emp+geom Beta(1,1) [mkp_eg]"),
       col = c("tomato","darkgreen","steelblue"), lty = 1, lwd = 2, bty = "n")
matplot(K_GRID, log(cbind(pmf_geo, pmf_high, pmf_eg)), type = "l",
        lty = 1, lwd = 2,
        col = c("tomato", "darkgreen", "steelblue"),
        xlab = "k'", ylab = "log P(k' | kObs=2)",
        main = "Prior PMFs (kObs = 2), log-scale")
par(op); dev.off()
cat("Wrote redteam_kprime_priors.png\n")

# --------------------------- prior * likelihood (per-character posteriors) --
# For each binary character: posterior(k') propto prior(k') * L(k')
# Plot raw likelihood profiles + implied posteriors at the medium TL regime.
make_posterior <- function(loglik, logprior) {
  lp <- loglik + logprior
  m <- max(lp, na.rm = TRUE)
  out <- exp(lp - m); out / sum(out, na.rm = TRUE)
}

regimes <- c("short", "med", "long")
results <- list()
for (regime in regimes) {
  M <- prof[[regime]]
  post_geo  <- apply(M, 2, function(ll) make_posterior(ll, lp_geo))
  post_high <- apply(M, 2, function(ll) make_posterior(ll, lp_high))
  post_eg   <- apply(M, 2, function(ll) make_posterior(ll, lp_eg))
  results[[regime]] <- list(loglik = M,
                            post_geo = post_geo,
                            post_high = post_high,
                            post_eg = post_eg)
}

# Plot: likelihood profile (med) + posteriors under each prior, for each char.
png(file.path(OUT_DIR, "redteam_kprime_loglik.png"),
    width = 1100, height = 700, res = 110)
op <- par(mfrow = c(3, 4), mar = c(3.2, 3.2, 2, 0.5), mgp = c(2, 0.6, 0))
M <- results$med$loglik
for (j in seq_along(bin_idx)) {
  ll <- M[, j]
  plot(K_GRID, ll - max(ll, na.rm = TRUE), type = "l", lwd = 2, col = "grey20",
       xlab = "k'", ylab = "log L - max",
       main = sprintf("char %d (kObs=2, med TL)", bin_idx[j]),
       ylim = c(-5, 0.2))
  abline(h = 0, col = "grey80", lty = 3)
  abline(v = 2, lty = 2, col = "tomato")
}
par(op); dev.off()
cat("Wrote redteam_kprime_loglik.png\n")

png(file.path(OUT_DIR, "redteam_kprime_posteriors.png"),
    width = 1100, height = 700, res = 110)
op <- par(mfrow = c(3, 4), mar = c(3.2, 3.2, 2, 0.5), mgp = c(2, 0.6, 0))
for (j in seq_along(bin_idx)) {
  matplot(K_GRID,
          cbind(results$med$post_geo[, j],
                results$med$post_high[, j],
                results$med$post_eg[, j]),
          type = "l", lty = 1, lwd = 2,
          col = c("tomato", "darkgreen", "steelblue"),
          xlab = "k'", ylab = "P(k' | x, ...)",
          main = sprintf("char %d (kObs=2, med TL)", bin_idx[j]))
  if (j == 1L) legend("topright",
                      c("mkp_geo","mkp_highk","mkp_eg"),
                      col = c("tomato","darkgreen","steelblue"),
                      lty = 1, lwd = 2, bty = "n", cex = 0.7)
}
par(op); dev.off()
cat("Wrote redteam_kprime_posteriors.png\n")

# Tabulate posterior summaries (mean, median, P(k' > kObs+1), P(k' > kObs+5))
post_summary <- function(P, lbl) {
  do.call(rbind, lapply(seq_len(ncol(P)), function(j) {
    p <- P[, j]
    mean_k <- sum(K_GRID * p)
    csum <- cumsum(p)
    med_k <- K_GRID[which(csum >= 0.5)[1]]
    data.frame(prior = lbl, char = bin_idx[j],
               mean_k = mean_k, median_k = med_k,
               p_gt_3 = sum(p[K_GRID > 3]),
               p_gt_7 = sum(p[K_GRID > 7]),
               p_gt_22 = sum(p[K_GRID > 22]))
  }))
}

summary_df <- do.call(rbind, list(
  post_summary(results$med$post_geo,  "mkp_geo"),
  post_summary(results$med$post_high, "mkp_highk"),
  post_summary(results$med$post_eg,   "mkp_eg")
))
cat("\n=== Posterior summary on binary chars (medium TL) ===\n")
print(summary_df)

# Also: prior-only mean k' (no data) for context
prior_mean <- function(pmf) sum(K_GRID * pmf)
cat("\n=== Prior-only mean k' (kObs=2) ===\n")
cat(sprintf("  mkp_geo  Beta(1,1):  E[k'] = %.2f\n", prior_mean(pmf_geo)))
cat(sprintf("  mkp_highk Beta(1,20): E[k'] = %.2f\n", prior_mean(pmf_high)))
cat(sprintf("  mkp_eg   emp+Beta(1,1): E[k'] = %.2f\n", prior_mean(pmf_eg)))

saveRDS(list(K_GRID = K_GRID, bin_idx = bin_idx,
             prior_pmf = list(geo = pmf_geo, high = pmf_high, eg = pmf_eg),
             profiles = prof,
             results  = results,
             summary  = summary_df),
        file.path(OUT_DIR, "redteam_kprime_shape.rds"))
cat("\nSaved redteam_kprime_shape.rds\n")
