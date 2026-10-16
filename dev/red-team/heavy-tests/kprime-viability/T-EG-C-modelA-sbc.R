#!/usr/bin/env Rscript
# T-EG-C — Honest SBC for empirical_geometric under Model A (unconditional).
#
# Companion to T0-T4-honest-sbc.R (geometric arm).  Validates the EG arm's
# sampler+prior+likelihood pipeline by running the chain through the new
# `priorVariant = "unconditional"` flag (which drops the per-character
# Z_i(p) subtraction so the prior has a coherent forward-generative joint).
#
# Forward (per simulation):
#   1. p_true       ~ Beta(a, b)
#   2. tree         ~ rtree(N_TIP) with edge lengths from Gamma·Dirichlet
#   3. N_obs_i      ~ empiricalNObs           (i = 1..N_CHAR)
#   4. N_unobs_i    ~ Geo(p_true)
#   5. kTrue_i      = N_obs_i + N_unobs_i
#   6. y_i          ~ JC(kTrue_i) on the tree
#   7. canonicalise y_i → kObs_i (= number of distinct realised states)
#
# Inference: MkPrime with kPrimePrior="empirical_geometric",
#            priorVariant="unconditional", expSteps fixed.
#
# Rank statistics: kTrue_i (per-char, randomised ties), p (continuous),
# tree_length (continuous), rate_log_sd (continuous).
#
# Pass: AD p > 0.4 on p, tree_length, rate_log_sd at N_SIM = 200.  The
# kPrime_pooled rank suffers the same Talts boundary artefact as the
# geometric arm (kObs is at the support floor for many chars); diagnostic
# only.
#
# Modes:
#   --quick  : N_SIM=30, N_TIP=8, N_CHAR=30, short MCMC.  Local smoke test.
#   default  : N_SIM=200, N_TIP=8, N_CHAR=30, full MCMC.  Hamilton.
#
# Outputs: dev/red-team/heavy-tests/kprime-viability/T-EG-C/sims.rds,
#          ranks.rds, rank-histograms.png, verdict.txt.

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(ape)
  library(TreeTools)
})

OUT_DIR <- "dev/red-team/heavy-tests/kprime-viability"
CONFIG_DIR <- file.path(OUT_DIR, "T-EG-C")
dir.create(CONFIG_DIR, recursive = TRUE, showWarnings = FALSE)

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
} else {
  # Matches geometric Model A worktree run (STATUS.md lines 77-80): N_ITER
  # = 12000 was the configuration at which tree_length AD p reached 0.47.
  N_SIM     <- 200L
  N_TIP     <- 8L
  N_CHAR    <- 30L
  N_ITER    <- 12000L
  N_THIN    <- 60L
  N_WARM    <- 4000L
}

EXPSTEPS_FIXED <- 50
TREE_SHAPE     <- 2
K_MAX_PRIOR    <- 30L          # cap on kTrue (matches T0/T4)
A_PRIOR        <- 1
B_PRIOR        <- 1
seedBase       <- 20260528L

EMP <- MkPrime::empiricalNObs

# -------- Forward sim helpers ------------------------------------------
.simTree <- function(nTip) {
  tr <- ape::rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  rate_ <- TREE_SHAPE / EXPSTEPS_FIXED
  tl    <- stats::rgamma(1, shape = TREE_SHAPE, rate = rate_)
  g <- stats::rgamma(nrow(tr$edge), shape = 1)
  rels <- g / sum(g)
  tr$edge.length <- tl * rels
  tr
}
.simJCchar <- function(tree, kTrue) {
  nTip <- length(tree$tip.label)
  states <- integer(2L * nTip - 1L)
  rootIdx <- nTip + 1L
  states[rootIdx] <- sample.int(kTrue, 1L) - 1L
  edges <- tree$edge; el <- tree$edge.length
  for (e in seq_len(nrow(edges))) {
    pa <- edges[e, 1L]; ch <- edges[e, 2L]; t <- el[e]
    pSame <- 1 / kTrue + (1 - 1 / kTrue) * exp(-kTrue * t / (kTrue - 1))
    if (runif(1L) < pSame) {
      states[ch] <- states[pa]
    } else {
      states[ch] <- sample(setdiff(seq.int(0L, kTrue - 1L), states[pa]), 1L)
    }
  }
  states[seq_len(nTip)]
}
.canon <- function(v) {
  uvals <- sort(unique(v))
  out <- match(v, uvals) - 1L
  attr(out, "kObs") <- length(uvals)
  out
}

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
  tr             <- .simTree(N_TIP)
  tl_true        <- sum(tr$edge.length)

  # Joint Model A forward:
  #   N_obs_i  ~ empiricalNObs
  #   N_unobs  ~ Geo(p_true)
  #   kTrue    = N_obs + N_unobs   (capped at K_MAX_PRIOR for tractability)
  n_obs   <- MkPrime:::.SampleNObsEmpirical(N_CHAR, EMP)
  n_unobs <- stats::rgeom(N_CHAR, p_true)
  kTrue   <- pmin(n_obs + n_unobs, K_MAX_PRIOR)

  sim_mat <- matrix(NA_integer_, N_TIP, N_CHAR,
                    dimnames = list(tr$tip.label, NULL))
  kObs <- integer(N_CHAR)
  for (j in seq_len(N_CHAR)) {
    raw <- .simJCchar(tr, kTrue[j])
    cv  <- .canon(raw)
    sim_mat[, j] <- cv
    kObs[j] <- attr(cv, "kObs")
  }
  keep <- kObs >= 2L
  if (sum(keep) < 3L) return(list(skipped = TRUE, reason = "few_var_chars"))
  sim_mat <- sim_mat[, keep, drop = FALSE]
  kTrue   <- kTrue[keep]
  kObs    <- kObs[keep]
  n_char  <- sum(keep)
  pd  <- TreeTools::MatrixToPhyDat(sim_mat)
  mkd <- MkPrimeData(pd)

  start_tree <- tr
  start_tree$edge.length <- rep_len(0.1, nrow(tr$edge))

  model <- suppressMessages(MkPrimeModel(
    coding       = "variable",
    nCat         = 1L,
    kPrimePrior  = "empirical_geometric",
    priorVariant = "unconditional",
    kprimeHyperA = A_PRIOR,
    kprimeHyperB = B_PRIOR,
    expSteps     = EXPSTEPS_FIXED
  ))
  mcmc <- MkPrimeMCMC(
    nIter = N_ITER, thin = N_THIN,
    minWarmup = N_WARM, maxWarmup = N_WARM,
    autoTune = FALSE, nRuns = 1L, nChains = 1L
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
  rk$kPrime_per_char <- vapply(seq_len(n_char), function(j) {
    char_post <- samples[, kp_cols[j]]
    .rankRandomised(kTrue[j], char_post)
  }, numeric(1))

  list(skipped = FALSE, L = nrow(samples), wall = dt,
       p_true = p_true, tl_true = tl_true, rateLogSd_true = rateLogSd_true,
       kTrue = kTrue, kObs = kObs, n_obs_true = n_obs, n_unobs_true = n_unobs,
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
}

good <- sapply(sims, function(s) !isTRUE(s$skipped))
cat(sprintf("\n[T-EG-C] Good sims: %d / %d\n", sum(good), N_SIM))
saveRDS(sims, file.path(CONFIG_DIR, "sims.rds"))

if (sum(good) < 20L) {
  cat("[T-EG-C] Too few good sims; aborting AD analysis.\n")
  quit(status = 1L)
}
goods <- sims[good]
L_used <- max(sapply(goods, `[[`, "L"))

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
cat(sprintf("Mode:          %s\n", mode))
cat(sprintf("N_SIM:         %d (good = %d)\n", N_SIM, sum(good)))
cat(sprintf("L_used:        %d\n", L_used))
cat(sprintf("N_TIP:         %d  N_CHAR (target): %d\n", N_TIP, N_CHAR))
cat(sprintf("Prior:         empirical_geometric, priorVariant=unconditional\n"))
cat(sprintf("Hyperprior:    p ~ Beta(%g, %g)\n", A_PRIOR, B_PRIOR))
cat(sprintf("Forward draw:  N_obs ~ empiricalNObs (corpus)\n"))
cat(sprintf("               N_unobs ~ Geo(p_true)\n"))
cat(sprintf("               kTrue = pmin(N_obs + N_unobs, %d)\n", K_MAX_PRIOR))
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
cat("the headline pass criterion is on tree_length, rate_log_sd, and p.\n")
sink()

cat(sprintf("\n[T-EG-C] artefacts in %s\n", CONFIG_DIR))
