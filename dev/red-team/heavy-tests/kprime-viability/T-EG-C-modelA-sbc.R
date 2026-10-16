#!/usr/bin/env Rscript
# T-EG-C — Honest SBC for empirical_geometric under Model A (unconditional).
#
# Companion to T0-T4-honest-sbc.R (geometric arm).  Validates the EG arm's
# sampler+prior+likelihood pipeline under Model A, which has a coherent
# forward-generative joint (no per-character Z_i(p) subtraction).  Model A
# is the only k' prior since 186f5a9 retired `priorVariant`.
#
# Forward (per simulation):
#   1. p_true       ~ Beta(a, b)
#   2. tree         ~ rtree(N_TIP) with edge lengths from Gamma·Dirichlet
#   3. kTrue_i      = N_obs_i + N_unobs_i, N_obs_i ~ empiricalNObs,
#                     N_unobs_i ~ Geo(p_true), redrawn until
#                     kTrue_i <= K_MAX_PRIOR (truncation, not a clamp; #57)
#   4. y_i          ~ JC(kTrue_i) on the tree; an invariant (kTrue_i, y_i)
#                     pair is redrawn whole (joint Pr(K, y | V), 57628cf)
#   5. canonicalize y_i → kObs_i (= number of distinct realized states)
#
# Inference: MkPrime with kPrimePrior="empirical_geometric",
#            kprimeTruncK=K_MAX_PRIOR,
#            expSteps fixed.  Forward helpers: T-EG-C-forward.R.
#
# Rank statistics: kTrue_i (per-char, randomised ties), p (continuous),
# tree_length (continuous), rate_log_sd (continuous).
#
# Pass: AD p > 0.4 on p and tree_length at N_SIM = 200.  rate_log_sd is
# prior-only (nCat = 1, so it never enters the likelihood): its rank checks
# only that its moves sample the Gamma prior.  The kPrime_pooled rank suffers the same Talts boundary artefact
# as the geometric arm (kObs is at the support floor for many chars);
# diagnostic only.
#
# Modes:
#   --quick  : N_SIM=30, N_TIP=8, N_CHAR=30, short MCMC.  Local smoke test.
#   default  : N_SIM=200, N_TIP=8, N_CHAR=30, full MCMC.  Hamilton.
#
# Outputs: dev/red-team/heavy-tests/kprime-viability/T-EG-C/sims.rds,
#          ranks.rds, rank-histograms.png, verdict.txt.

suppressPackageStartupMessages({
  source("dev/red-team/heavy-tests/load-mkprime.R")
  LoadMkPrime(".")
})
source("dev/red-team/heavy-tests/kprime-viability/T-EG-C-forward.R")

OUT_DIR <- "dev/red-team/heavy-tests/kprime-viability"
CONFIG_DIR <- file.path(OUT_DIR, "T-EG-C")
dir.create(CONFIG_DIR, recursive = TRUE, showWarnings = FALSE)
# An aborted run must not leave an earlier verdict beside its sims.rds.
unlink(file.path(CONFIG_DIR, "verdict.txt"))
gitSha <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) "unknown", warning = function(w) "unknown"
)
runId <- Sys.getenv("SLURM_JOB_ID", "local")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (any(args %in% c("--quick", "-q"))) "quick" else "full"
cat(sprintf("[T-EG-C] mode=%s\n", mode))

if (mode == "quick") {
  N_SIM     <- 30L
  N_TIP     <- 8L
  N_CHAR    <- 30L
  N_ITER    <- 1000L
  N_THIN    <- 10L
  N_WARM    <- 500L
  maxTimeSim <- 120
} else {
  # Matches geometric Model A worktree run (STATUS.md lines 77-80): N_ITER
  # = 12000 was the configuration at which tree_length AD p reached 0.47.
  N_SIM     <- 200L
  N_TIP     <- 8L
  N_CHAR    <- 30L
  N_ITER    <- 12000L
  N_THIN    <- 60L
  N_WARM    <- 4000L
  maxTimeSim <- 600
}

EXPSTEPS_FIXED <- 50
TREE_SHAPE     <- 2
K_MAX_PRIOR    <- 30L          # truncation K, forward and inference
A_PRIOR        <- 1
B_PRIOR        <- 1
seedBase       <- 20260528L

EMP <- MkPrime::empiricalNObs

# -------- Talts randomised rank for discrete params --------------------
.rankRandomised <- function(true_val, post_samples) {
  if (!is.finite(true_val) || length(post_samples) == 0L) return(NA_real_)
  below <- sum(post_samples < true_val)
  ties  <- sum(post_samples == true_val)
  below + sample.int(ties + 1L, 1L) - 1L
}
.rankPlain <- function(true_val, post_samples) {
  if (!is.finite(true_val) || length(post_samples) == 0L) return(NA_real_)
  sum(post_samples < true_val)
}

# -------- One simulation ----------------------------------------------
.runOneSim <- function(sim_id, seed) {
  set.seed(seed)
  p_true         <- stats::rbeta(1, A_PRIOR, B_PRIOR)
  rateLogSd_true <- stats::rgamma(1, shape = 1, rate = 1)
  tr             <- .SimTree(N_TIP, TREE_SHAPE, EXPSTEPS_FIXED)
  tl_true        <- sum(tr$edge.length)

  sim_mat <- matrix(NA_integer_, N_TIP, N_CHAR,
                    dimnames = list(tr$tip.label, NULL))
  kTrue <- kObs <- integer(N_CHAR)
  for (j in seq_len(N_CHAR)) {
    ch <- .SimVariableChar(tr, p_true, EMP, K_MAX_PRIOR)
    if (is.null(ch)) return(list(skipped = TRUE, reason = "redraw_cap"))
    kTrue[j] <- ch$k
    sim_mat[, j] <- ch$states
    kObs[j] <- attr(ch$states, "kObs")
  }
  pd  <- TreeTools::MatrixToPhyDat(sim_mat)
  mkd <- MkPrimeData(pd)

  start_tree <- tr
  start_tree$edge.length <- rep_len(0.1, nrow(tr$edge))

  model <- suppressMessages(MkPrimeModel(
    coding       = "variable",
    nCat         = 1L,
    kPrimePrior  = "empirical_geometric",
    kprimeHyperA = A_PRIOR,
    kprimeHyperB = B_PRIOR,
    kprimeTruncK = K_MAX_PRIOR,
    empiricalNObs = EMP,
    expSteps     = EXPSTEPS_FIXED
  ))
  mcmc <- MkPrimeMCMC(
    nIter = N_ITER, thin = N_THIN,
    minWarmup = N_WARM, maxWarmup = N_WARM,
    autoTune = FALSE, nRuns = 1L, nChains = 1L, maxTime = maxTimeSim
  )
  t0 <- Sys.time()
  res <- tryCatch(suppressMessages(suppressWarnings(
    RunMkPrime(mkd, start_tree, model = model, mcmc = mcmc,
               fixTopology = TRUE, overwrite = TRUE)
  )), error = function(e) list(error = conditionMessage(e)))
  if (!is.null(res$error)) {
    return(list(skipped = TRUE, reason = paste0("mcmc:", res$error)))
  }
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  samples <- res$samples
  if (is.null(samples) || nrow(samples) < 10L) {
    return(list(skipped = TRUE, reason = "no_samples"))
  }

  kp_cols <- grep("^kPrime_", colnames(samples))
  rk <- list()
  rk$tree_length <- .rankPlain(tl_true,        samples[, "tree_length"])
  rk$rate_log_sd <- .rankPlain(rateLogSd_true, samples[, "rate_log_sd"])
  rk$p           <- .rankPlain(p_true,         samples[, "p"])
  rk$kPrime_per_char <- vapply(seq_len(N_CHAR), function(j) {
    char_post <- samples[, kp_cols[j]]
    .rankRandomised(kTrue[j], char_post)
  }, numeric(1))

  list(skipped = FALSE, L = nrow(samples), wall = dt,
       p_true = p_true, tl_true = tl_true, rateLogSd_true = rateLogSd_true,
       kTrue = kTrue, kObs = kObs,
       ranks = rk)
}

# -------- Driver -------------------------------------------------------
sims <- vector("list", N_SIM)
for (i in seq_len(N_SIM)) {
  seed <- seedBase + i
  cat(sprintf("[T-EG-C] sim %3d/%d ... ", i, N_SIM))
  s <- tryCatch(.runOneSim(i, seed),
                error = function(e) list(skipped = TRUE,
                                          reason = paste0("trycatch:",
                                                          conditionMessage(e))))
  if (isTRUE(s$skipped)) {
    cat(sprintf("SKIP (%s)\n", s$reason))
  } else {
    cat(sprintf("L=%d wall=%.1fs\n", s$L, s$wall))
  }
  sims[[i]] <- s
  saveRDS(sims[seq_len(i)], file.path(CONFIG_DIR, "sims.rds"))
}

good <- sapply(sims, function(s) !isTRUE(s$skipped))
cat(sprintf("\n[T-EG-C] Good sims: %d / %d\n", sum(good), N_SIM))

if (sum(good) < 20L) {
  cat("[T-EG-C] Too few good sims; aborting AD analysis.\n")
  quit(status = 1L)
}
goods <- sims[good]
chainLen <- vapply(goods, `[[`, integer(1), "L")
fullLen <- max(chainLen)
# A chain cut short by maxTime has fewer draws; slow sims are high-k', low-p
# ones, so these exclusions are not random and are reported in the verdict.
truncatedP <- vapply(goods[chainLen < fullLen], `[[`, numeric(1), "p_true")
goods <- goods[chainLen == fullLen]
L_used <- fullLen
if (length(goods) < 20L) {
  cat("[T-EG-C] Too few full-length chains; aborting AD analysis.\n")
  quit(status = 1L)
}

tl_ranks  <- sapply(goods, function(s) s$ranks$tree_length)
rls_ranks <- sapply(goods, function(s) s$ranks$rate_log_sd)
p_ranks   <- sapply(goods, function(s) s$ranks$p)
kp_ranks  <- unlist(lapply(goods, function(s) s$ranks$kPrime_per_char))

.normRanks <- function(r, L) (r + 0.5) / (L + 1)
.adP <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4L) return(NA_real_)
  if (requireNamespace("goftest", quietly = TRUE)) {
    suppressWarnings(goftest::ad.test(x, null = "punif")$p.value)
  } else {
    suppressWarnings(ks.test(x, "punif")$p.value)
  }
}
ad <- list(
  tree_length   = .adP(.normRanks(tl_ranks,  L_used)),
  rate_log_sd   = .adP(.normRanks(rls_ranks, L_used)),
  p             = .adP(.normRanks(p_ranks,   L_used)),
  kPrime_pooled = .adP(.normRanks(kp_ranks,  L_used))
)
cat(sprintf("\n[T-EG-C] AD p-values:\n"))
for (nm in names(ad)) {
  cat(sprintf("  %-18s: %.4f  %s\n", nm, ad[[nm]],
              if (is.na(ad[[nm]])) "NA"
              else if (ad[[nm]] > 0.4) "PASS"
              else if (ad[[nm]] > 0.01) "MARGINAL"
              else "FAIL"))
}

saveRDS(list(ranks = list(tree_length = tl_ranks, rate_log_sd = rls_ranks,
                          p = p_ranks, kPrime_pooled = kp_ranks),
             L = L_used, ad = ad),
        file.path(CONFIG_DIR, "ranks.rds"))

png(file.path(CONFIG_DIR, "rank-histograms.png"),
    width = 1200, height = 800, res = 110)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (nm in c("tree_length", "rate_log_sd", "p", "kPrime_pooled")) {
  r <- switch(nm,
              tree_length = tl_ranks, rate_log_sd = rls_ranks,
              p = p_ranks, kPrime_pooled = kp_ranks)
  hist(r, breaks = 30, freq = FALSE,
       main = sprintf("%s\nAD p=%.4f", nm, ad[[nm]]),
       xlab = "rank", col = "lightgrey", border = "white")
  abline(h = 1 / L_used, col = "red", lwd = 2)
}
par(op)
dev.off()

sink(file.path(CONFIG_DIR, "verdict.txt"))
cat(sprintf("Config:        T-EG-C — empirical_geometric Model A SBC\n"))
cat(sprintf("Commit:        %s  run_id: %s\n", gitSha, runId))
cat(sprintf("Mode:          %s\n", mode))
cat(sprintf("N_SIM:         %d (good = %d)\n", N_SIM, sum(good)))
cat(sprintf(
  "Truncated:     %d chain(s) cut by maxTime, excluded; p_true: %s\n",
  length(truncatedP), paste(signif(truncatedP, 3), collapse = ", ")
))
cat(sprintf("L_used:        %d\n", L_used))
cat(sprintf("N_TIP:         %d  N_CHAR (target): %d\n", N_TIP, N_CHAR))
cat(sprintf("Prior:         empirical_geometric (Model A)\n"))
cat(sprintf("Hyperprior:    p ~ Beta(%g, %g)\n", A_PRIOR, B_PRIOR))
cat(sprintf("Forward draw:  N_obs ~ empiricalNObs (corpus)\n"))
cat(sprintf("               N_unobs ~ Geo(p_true)\n"))
cat(sprintf(
  "               kTrue = N_obs + N_unobs, truncated to <= %d\n", K_MAX_PRIOR
))
cat("               invariant (kTrue, y) pairs redrawn whole\n")
cat(sprintf("Inference K:   kprimeTruncK = %d\n", K_MAX_PRIOR))
cat("\nAD p-values vs Uniform(0, 1) [pass gate: > 0.4]:\n")
for (nm in names(ad)) {
  cat(sprintf("  %-18s: %.4f  %s\n", nm, ad[[nm]],
              if (is.na(ad[[nm]])) "NA"
              else if (ad[[nm]] > 0.4) "PASS"
              else if (ad[[nm]] > 0.01) "MARGINAL"
              else "FAIL"))
}
cat("\nNote: kPrime_pooled is expected to show a Talts boundary artefact when\n")
cat("kObs == kTrue (the truth lies at the support floor) — diagnostic only;\n")
cat("the headline pass criterion is on tree_length and p.  rate_log_sd is\n")
cat("prior-only under nCat = 1: its rank checks only its moves.\n")
sink()

cat(sprintf("\n[T-EG-C] artefacts in %s\n", CONFIG_DIR))
