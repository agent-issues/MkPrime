#!/usr/bin/env Rscript
# T-SBC-sampled-geometric — Honest SBC for geometric arm under sampled-k.
#
# Plan reference: dev/notes/2026-05-28-marginal-k-plan.md §7.3 (Rung 3 SBC).
# Proof reference: dev/red-team/proofs/marginal-k-sampled-rb-consistency.md
#   (sampled_k <-> marginal_k RB-consistency under the truncated geometric) and
#   dev/red-team/proofs/marginal-k-truncation-normaliser.md (the Z(p) normaliser).
# This driver empirically validates the chain-level calibration of the sampled_k
# geometric arm (Stage 2: the truncated sampled_k prior matches marginal_k).
#
# Mirrors the marginal-k SBC (T-SBC-marginal-geometric.R) structurally; the
# only inference diff is likelihoodMode. The diffs vs the kprime-viability
# precedent are:
#   1. Forward: u_i ~ Geo(p_true), kTrue_i = 2 + u_i (pure geometric Model A,
#      no empiricalNObs convolution).
#   2. Inference: MkPrimeModel(kPrimePrior = "geometric",
#      likelihoodMode = "sampled_k").
#   3. kPrime IS a sampled state under sampled-k, but k'_pooled SBC is a Talts
#      boundary artefact (k' floored at kObs), so it is EXCLUDED from the pass
#      bar (project memory project_sbc_kprime_structural). Pass bar is on
#      (tree_length, rate_log_sd, p), as in the marginal-k driver.
#   4. Quick-mode trigger via env var MARGINAL_K_SBC_QUICK=1 (per PR-C spec).
#   5. Output dir: dev/red-team/heavy-tests/marginal-k/sbc-results-sampled/
#
# DO NOT MODIFY T-EG-C-modelA-sbc.R when iterating here — that harness is
# the validated sampled-k precedent and a moving target invalidates the
# comparison.
#
# Forward (per simulation):
#   1. p_true       ~ Beta(a, b)
#   2. rateLogSd    ~ Gamma(1, 1)
#   3. tree         ~ rtree(N_TIP) with edge lengths from Gamma · Dirichlet
#   4. u_i          ~ Geo(p_true)              (i = 1..N_CHAR)
#   5. kTrue_i      = pmin(2 + u_i, K_MAX_PRIOR)
#   6. y_i          ~ JC(kTrue_i) on the tree
#   7. canonicalise y_i → kObs_i (= number of distinct realised states)
#
# Inference:
#   MkPrimeModel(kPrimePrior = "geometric",
#                likelihoodMode = "marginal_k",
#                priorVariant = "unconditional",   # Model A — matches forward
#                expSteps = EXPSTEPS_FIXED,
#                kprimeHyperA = 1, kprimeHyperB = 1)
#
# Rank statistics (continuous parameters only):
#   tree_length, rate_log_sd, p.
#
# Pass bar (plan §7.3, full mode): AD p > 0.4 on each of tree_length,
# rate_log_sd, p at N_SIM = 200. MARGINAL between 0.01 and 0.4. FAIL below.
#
# Quick mode (env var): N_SIM=20, N_ITER=1000. Insufficient AD power —
# diagnostic only; PASS in quick mode is NOT evidence the chain is correct
# (see SBC-HARNESS-001..006 in dev/red-team/findings.md for the precedent
# pattern of harnesses that "pass quick" and fail full).
#
# Outputs: dev/red-team/heavy-tests/marginal-k/sbc-results-sampled/sims.rds,
#          ranks.rds, rank-histograms.png, verdict.txt.

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(ape)
  library(TreeTools)
})

OUT_DIR <- "dev/red-team/heavy-tests/marginal-k/sbc-results-sampled"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

quick_env <- Sys.getenv("MARGINAL_K_SBC_QUICK", unset = "")
mode <- if (identical(quick_env, "1")) "quick" else "full"
cat(sprintf("[T-SBC-sampled-geometric] mode=%s\n", mode))

if (mode == "quick") {
  # Local smoke-test config — N_SIM = 20 is too small for AD power.
  N_SIM     <- 20L
  N_TIP     <- 8L
  N_CHAR    <- 30L
  N_ITER    <- 1000L
  N_THIN    <- 10L
  N_WARM    <- 500L
} else {
  # Hamilton full config.
  #
  # 2026-05-29: regime moved off the saturated 8-tip / expSteps=50 toy
  # (mean tree length ~50 over 8 tips ⇒ mean edge ~3.6 ⇒ every character
  # randomised to noise). At saturation the chain cannot mix in
  # tree_length/p and SBC ranks spike at the extremes — a non-convergence
  # artefact, NOT a calibration defect (priors provably match: tree-length
  # Gamma(2, 2/expSteps) and the JC rate convention exp(−K/(K−1)·rate·t)).
  # Realistic low-homoplasy morphology (mean tree length ~1.4, 16 tips)
  # is both more relevant AND a more valid SBC test because the chain
  # actually mixes. N_CHAR raised to 100 to offset constant-character
  # (kObs<2) filtering at short branch lengths.
  N_SIM     <- 200L
  N_TIP     <- 16L
  N_CHAR    <- 100L
  N_ITER    <- 12000L
  N_THIN    <- 60L     # → ~133 retained / chain; well under 30k thin budget
  N_WARM    <- 4000L
}

EXPSTEPS_FIXED <- 1.4          # mean tree length (forward + inference prior)
TREE_SHAPE     <- 2
K_MAX_PRIOR    <- 30L          # cap on kTrue
A_PRIOR        <- 1
B_PRIOR        <- 1
# seedBase override (MARGINAL_K_SBC_SEEDBASE) lets an independent second batch
# run on disjoint seeds for the pre-registered 2-batch confirmation; default
# unchanged so the canonical run is reproducible.
seedBase       <- {
  sb_env <- Sys.getenv("MARGINAL_K_SBC_SEEDBASE", unset = "")
  if (nzchar(sb_env)) as.integer(sb_env) else 20260528L
}

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
  # Iterate edges in seq_len() order (NOT rev()): SBC-TL-MIX-001 fix.
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

  # Pure geometric Model A forward (per proof §2 + §7.4):
  #   u_i   ~ Geo(p_true)         (location-free; full support on u >= 0)
  #   kTrue = pmin(2 + u, K_MAX_PRIOR)
  #
  # 2026-05-29 LEWIS-MKV FORWARD (Model I-a). Mk' implements the standard
  # Lewis-Mkv ascertainment: k' ~ w_k(p) UNCONDITIONALLY (intrinsic property),
  # and each character's likelihood is conditioned on being variable GIVEN its
  # k, i.e. L(y|k)/a_k with a_k = P(variable|k,tree). Marginalising k gives
  # Σ_k w_k(p) L(y|k)/a_k (sum-of-ratios), which is exactly what the C++
  # marginal path computes — so the inference is CORRECT for this model (no
  # code change). Confirmed: MkPrimeData drops invariant characters
  # unconditionally, so the data path always conditions on variability.
  #
  # The MATCHING forward (what earlier runs got wrong by jointly redrawing k
  # and y, a different "filtered-pool" model giving ratio-of-sums): draw
  # kTrue_j ONCE from the unconditional prior, then redraw the DATA ONLY
  # (k held fixed) until the character is variable. This yields k ~ w_k and
  # y ~ P(y | k, variable) — Model I-a. All characters are variable, so
  # MkPrimeData drops nothing and the coding="variable" /a_k correction
  # applies cleanly. SBC must pass if the marginal-k core + Mkv are correct.
  # MARGINAL-K-TRUNC-001: draw k' from a TRUNCATED geometric on [2, K_MAX_PRIOR]
  # by rejection-redraw — NOT pmin(2+Geom, K) which piles a point mass at K.
  # The inference (cpp_log_likelihood_marginal) caps its marginal sum at k<=K
  # AND renormalises by Z(p)=1-(1-p)^(K-1); the forward must match for SBC
  # validity (proof dev/red-team/proofs/marginal-k-truncation-normaliser.md).
  kTrue <- integer(N_CHAR)
  for (jj in seq_len(N_CHAR)) {
    repeat { kk <- 2L + stats::rgeom(1L, p_true); if (kk <= K_MAX_PRIOR) break }
    kTrue[jj] <- kk
  }
  u_true <- kTrue - 2L                             # offset (k held FIXED below)
  MAX_REDRAW <- 100000L
  sim_mat <- matrix(NA_integer_, N_TIP, N_CHAR,
                    dimnames = list(tr$tip.label, NULL))
  kObs <- integer(N_CHAR)
  for (j in seq_len(N_CHAR)) {
    attempt <- 0L
    repeat {
      attempt <- attempt + 1L
      if (attempt > MAX_REDRAW) {
        return(list(skipped = TRUE, reason = "redraw_cap"))
      }
      cv <- .canon(.simJCchar(tr, kTrue[j]))        # redraw DATA only; k fixed
      if (attr(cv, "kObs") >= 2L) break
    }
    sim_mat[, j] <- cv
    kObs[j] <- attr(cv, "kObs")
  }
  n_char <- N_CHAR
  pd  <- TreeTools::MatrixToPhyDat(sim_mat)
  mkd <- MkPrimeData(pd)

  start_tree <- tr
  start_tree$edge.length <- rep_len(0.1, nrow(tr$edge))

  model <- suppressMessages(MkPrimeModel(
    # Lewis-Mkv (Model I-a): per-char likelihood conditioned on variable given
    # k. Matches the redraw-data-only forward; this is the model the C++
    # marginal path (sum-of-ratios Σ_k w_k L/a_k) actually implements.
    coding         = "variable",
    nCat           = 1L,
    kPrimePrior    = "geometric",
    likelihoodMode = "sampled_k",
    # Forward draws kTrue_i = 2 + u_i (Model A, unconditional on kObs_i), so
    # inference uses the unconditional marginal weights p (1-p)^(k - 2).
    priorVariant   = "unconditional",
    # Stage 1b: the model default truncation cap is now 200; the inference K
    # MUST equal the forward's K_MAX_PRIOR for calibration, so pin it here.
    kprimeTruncK   = K_MAX_PRIOR,
    kprimeHyperA   = A_PRIOR,
    kprimeHyperB   = B_PRIOR,
    expSteps       = EXPSTEPS_FIXED
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

  # Per plan §2 ("state vector schema does NOT shrink under marginal mode —
  # state->kPrime remains as a no-op field initialised to kObs"), kPrime_*
  # trace columns are still emitted under marginal-k but should sit
  # exactly at kObs. The plan §7.4 "no kPrime_* in trace" smoke test is
  # an aspirational follow-up (schema shrink), not a current correctness
  # check. We record presence + a kPrime==kObs spot check for the first
  # sim only.
  kp_cols <- grep("^kPrime_", colnames(samples), value = TRUE)
  if (sim_id == 1L && length(kp_cols) > 0L) {
    kp_first <- samples[1L, kp_cols]
    kobs_align <- as.integer(kp_first) == kObs[seq_along(kp_first)]
    cat(sprintf("[T-SBC-sampled-geometric] sim1 has %d kPrime_* cols (expected per plan §2); first-row matches kObs: %d/%d\n",
                length(kp_cols), sum(kobs_align), length(kobs_align)))
  }

  rk <- list()
  rk$tree_length <- .rankPlain(tl_true,        samples[, "tree_length"])
  rk$rate_log_sd <- .rankPlain(rateLogSd_true, samples[, "rate_log_sd"])
  rk$p           <- .rankPlain(p_true,         samples[, "p"])

  list(skipped = FALSE, L = nrow(samples), wall = dt,
       p_true = p_true, tl_true = tl_true, rateLogSd_true = rateLogSd_true,
       kTrue = kTrue, kObs = kObs, u_true = u_true,
       n_kPrime_cols = length(kp_cols),
       ranks = rk)
}

# -------- Sharding / aggregation control (backfill-friendly array job) --
# The full monolith needs ~N_SIM * per-sim-wall seconds; at ~3 min/sim that is
# ~10 h for N_SIM=200 — OVER the 8 h single-job walltime, and the verdict only
# computes after the whole loop, so a monolith was doomed to time out with NO
# verdict. Splitting into NSHARD array tasks (each a contiguous block of sims
# written to shards/) + one AGGREGATE pass also makes each task short enough to
# backfill under fair-share depression.
#
# Three execution paths, selected by env vars (all default to the monolith):
#   * AGGREGATE (MARGINAL_K_SBC_AGGREGATE=1): skip sims, merge shard files into
#     the full `sims` list, then run the identical AD/verdict code below.
#   * SHARD     (MARGINAL_K_SBC_NSHARD>1):    run this task's contiguous block
#     of sims -> shards/sims-shard-<id>.rds, then exit BEFORE aggregation.
#   * MONOLITH  (default):                    original behaviour, unchanged.
# Seeds are seedBase + i regardless of path, so the sharded run reproduces the
# monolith's sims EXACTLY; L_used and .normRanks are computed over all good
# sims in both, so the verdict is identical-by-construction.
nshard         <- max(1L, as.integer(Sys.getenv("MARGINAL_K_SBC_NSHARD", "1")))
shard_id       <- as.integer(Sys.getenv("MARGINAL_K_SBC_SHARD", "0"))   # 0-indexed
aggregate_mode <- identical(Sys.getenv("MARGINAL_K_SBC_AGGREGATE", ""), "1")
SHARD_DIR      <- file.path(OUT_DIR, "shards")
.shardFile     <- function(k) file.path(SHARD_DIR, sprintf("sims-shard-%03d.rds", k))

if (aggregate_mode) {
  # ---- AGGREGATE: reconstruct the full sims list from shard files --------
  sims <- vector("list", N_SIM)
  for (i in seq_len(N_SIM)) sims[[i]] <- list(skipped = TRUE, reason = "missing_shard")
  filled <- 0L
  for (k in 0:(nshard - 1L)) {
    f <- .shardFile(k)
    if (!file.exists(f)) { cat(sprintf("[aggregate] shard %03d MISSING\n", k)); next }
    sh <- readRDS(f)
    if (length(sh$idx)) { sims[sh$idx] <- sh$results; filled <- filled + length(sh$idx) }
  }
  cat(sprintf("[aggregate] merged %d shards; %d / %d sim slots filled\n", nshard, filled, N_SIM))
} else if (nshard > 1L) {
  # ---- SHARD: run my contiguous block, checkpoint, exit pre-aggregation --
  dir.create(SHARD_DIR, recursive = TRUE, showWarnings = FALSE)
  chunk   <- as.integer(ceiling(N_SIM / nshard))
  my_from <- shard_id * chunk + 1L
  my_to   <- min((shard_id + 1L) * chunk, N_SIM)
  my_idx  <- if (my_from <= N_SIM) seq.int(my_from, my_to) else integer(0)
  cat(sprintf("[shard %03d/%d] sims %s\n", shard_id, nshard,
              if (length(my_idx)) sprintf("%d..%d", my_from, my_to) else "(none)"))
  results <- vector("list", length(my_idx))
  for (j in seq_along(my_idx)) {
    i <- my_idx[j]; seed <- seedBase + i
    cat(sprintf("[shard %03d] sim %3d/%d ... ", shard_id, i, N_SIM))
    s <- tryCatch(.runOneSim(i, seed),
                  error = function(e) list(skipped = TRUE,
                            reason = paste0("trycatch:", conditionMessage(e))))
    if (isTRUE(s$skipped)) cat(sprintf("SKIP (%s)\n", s$reason)) else cat(sprintf("L=%d wall=%.1fs\n", s$L, s$wall))
    results[[j]] <- s
    # Checkpoint every sim so a per-task timeout keeps completed sims.
    saveRDS(list(idx = my_idx[seq_len(j)], results = results[seq_len(j)]), .shardFile(shard_id))
  }
  saveRDS(list(idx = my_idx, results = results), .shardFile(shard_id))
  cat(sprintf("[shard %03d/%d] done (%d sims)\n", shard_id, nshard, length(my_idx)))
  quit(save = "no", status = 0L)
} else {
  # ---- MONOLITH: original path (behaviour unchanged) ---------------------
  sims <- vector("list", N_SIM)
  for (i in seq_len(N_SIM)) {
    seed <- seedBase + i
    cat(sprintf("[T-SBC-sampled-geometric] sim %3d/%d ... ", i, N_SIM))
    s <- tryCatch(.runOneSim(i, seed),
                  error = function(e) list(skipped = TRUE,
                            reason = paste0("trycatch:", conditionMessage(e))))
    if (isTRUE(s$skipped)) {
      cat(sprintf("SKIP (%s)\n", s$reason))
    } else {
      cat(sprintf("L=%d wall=%.1fs\n", s$L, s$wall))
    }
    sims[[i]] <- s
    # Defensive: save partial sims rds every 10 sims so a Hamilton timeout
    # doesn't lose everything.
    if (i %% 10L == 0L) saveRDS(sims, file.path(OUT_DIR, "sims.rds"))
  }
}

good <- sapply(sims, function(s) !isTRUE(s$skipped))
cat(sprintf("\n[T-SBC-sampled-geometric] Good sims: %d / %d\n", sum(good), N_SIM))
saveRDS(sims, file.path(OUT_DIR, "sims.rds"))

min_good <- if (mode == "quick") 5L else 20L
if (sum(good) < min_good) {
  cat(sprintf("[T-SBC-sampled-geometric] Too few good sims (need %d); aborting AD analysis.\n",
              min_good))
  writeLines("INSUFFICIENT_SIMS",
             file.path(OUT_DIR, "verdict.txt"))
  quit(status = 1L)
}
goods <- sims[good]
L_used <- max(sapply(goods, `[[`, "L"))

tl_ranks  <- sapply(goods, function(s) s$ranks$tree_length)
rls_ranks <- sapply(goods, function(s) s$ranks$rate_log_sd)
p_ranks   <- sapply(goods, function(s) s$ranks$p)

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
  p             = .adP(.normRanks(p_ranks,   L_used))
)
.classify <- function(p) {
  if (is.na(p))       "NA"
  else if (p > 0.4)   "PASS"
  else if (p > 0.01)  "MARGINAL"
  else                "FAIL"
}
cat(sprintf("\n[T-SBC-sampled-geometric] AD p-values:\n"))
for (nm in names(ad)) {
  cat(sprintf("  %-14s: %.4f  %s\n", nm, ad[[nm]], .classify(ad[[nm]])))
}

saveRDS(list(ranks = list(tree_length = tl_ranks,
                          rate_log_sd = rls_ranks,
                          p           = p_ranks),
             L = L_used, ad = ad, mode = mode),
        file.path(OUT_DIR, "ranks.rds"))

png(file.path(OUT_DIR, "rank-histograms.png"),
    width = 1200, height = 400, res = 110)
op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
for (nm in c("tree_length", "rate_log_sd", "p")) {
  r <- switch(nm,
              tree_length = tl_ranks, rate_log_sd = rls_ranks, p = p_ranks)
  hist(r, breaks = 30, freq = FALSE,
       main = sprintf("%s\nAD p=%.4f", nm, ad[[nm]]),
       xlab = "rank", col = "lightgrey", border = "white")
  abline(h = 1 / L_used, col = "red", lwd = 2)
}
par(op)
dev.off()

# Headline verdict — PASS only if every continuous param's AD > 0.01 AND
# at least one is > 0.4. Quick-mode results are diagnostic only.
ad_finite <- vapply(ad, function(x) is.finite(x), logical(1L))
all_finite <- all(ad_finite)
all_marg_or_better <- all_finite && all(unlist(ad) > 0.01)
any_full_pass      <- all_finite && any(unlist(ad) > 0.4)
# STRICT gate = the real acceptance bar (plan §7.3 full mode): EVERY continuous
# parameter must individually clear AD > 0.4. This is STRICTER than the lenient
# "headline verdict" below, which prints PASS on any-one > 0.4 (e.g. AD =
# {tl 0.5, rls 0.02, p 0.02} prints "PASS" but is a strict FAIL). A cold-resume
# poll reads verdict.txt, so emit the strict result explicitly and unambiguously.
strict_pass <- all_finite && all(unlist(ad) > 0.4)
strict_gate <- if (mode == "quick") "QUICK_SANITY_ONLY" else if (!all_finite) "ERROR" else if (strict_pass) "PASS" else "FAIL"
verdict <- if (!all_finite) {
  "ERROR"
} else if (mode == "quick") {
  "QUICK_SANITY_ONLY"
} else if (all_marg_or_better && any_full_pass) {
  "PASS"
} else if (all_marg_or_better) {
  "MARGINAL"
} else {
  "FAIL"
}

sink(file.path(OUT_DIR, "verdict.txt"))
cat(sprintf("Driver:        T-SBC-sampled-geometric (PR-C / plan §7.3)\n"))
cat(sprintf("Mode:          %s\n", mode))
cat(sprintf("N_SIM:         %d (good = %d)\n", N_SIM, sum(good)))
cat(sprintf("L_used:        %d\n", L_used))
cat(sprintf("N_TIP:         %d  N_CHAR (target): %d\n", N_TIP, N_CHAR))
cat(sprintf("N_ITER:        %d  N_WARM: %d  N_THIN: %d\n", N_ITER, N_WARM, N_THIN))
cat(sprintf("Prior:         geometric (k'_i = kObs_i + Geo(p))\n"))
cat(sprintf("Mode flag:     likelihoodMode = 'marginal_k'\n"))
cat(sprintf("Prior variant: unconditional (Model A: k' = 2 + Geo(p))\n"))
cat(sprintf("Hyperprior:    p ~ Beta(%g, %g)\n", A_PRIOR, B_PRIOR))
cat(sprintf("Forward draw:  kTrue_i ~ TRUNCATED Geo(p_true) on [2, %d] via rejection\n",
            K_MAX_PRIOR))
cat("               redraw (NOT pmin-cap; MARGINAL-K-TRUNC-001), drawn ONCE\n")
cat("               (uncond.); then redraw DATA only (k fixed) until variable\n")
cat("               (Lewis-Mkv, Model I-a). Inference coding='variable',\n")
cat(sprintf("               marginal sum capped at k<=%d + renormalised by Z(p).\n", K_MAX_PRIOR))
cat("\nAD p-values vs Uniform(0, 1) [full-mode pass gate: > 0.4]:\n")
for (nm in names(ad)) {
  cat(sprintf("  %-14s: %.4f  %s\n", nm, ad[[nm]], .classify(ad[[nm]])))
}
cat(sprintf("\nHeadline verdict: %s  (lenient: all AD>0.01 AND any AD>0.4)\n", verdict))
cat(sprintf("STRICT gate:      %s  [the real bar: AD>0.4 on tree_length AND rate_log_sd AND p]\n", strict_gate))
if (all_finite && !strict_pass && mode != "quick") {
  cat(sprintf("  params below 0.4: %s\n", paste(names(ad)[unlist(ad) <= 0.4], collapse = ", ")))
}
if (mode == "quick") {
  cat("\nNote: QUICK mode (N_SIM=20) is INSUFFICIENT for AD power.\n")
  cat("Quick-mode AD numbers are noise; do not interpret as PASS/FAIL.\n")
  cat("Submit the full Hamilton run (sbatch submit-marginal-k-sbc.sh)\n")
  cat("for the validated §7.3 pass-bar check.\n")
}
cat("\nNote: kPrime_* rank tests are not reported under marginal-k —\n")
cat("kPrime is no longer a sampled state variable (proof §2). The pass\n")
cat("bar applies to tree_length, rate_log_sd, p only.\n")
sink()

writeLines(verdict, file.path(OUT_DIR, "verdict-headline.txt"))
writeLines(strict_gate, file.path(OUT_DIR, "verdict-strict.txt"))

cat(sprintf("\n[T-SBC-sampled-geometric] verdict: %s\n", verdict))
cat(sprintf("[T-SBC-sampled-geometric] artefacts in %s\n", OUT_DIR))
