#!/usr/bin/env Rscript
# T9 — Model A beta_geometric SBC (post s/r-reparameterisation fix).
#
# Forward (Model A joint generative):
#   alpha_true ~ Exp(1),  beta_true ~ Exp(1)
#   for each char i:
#     p_i ~ Beta(alpha_true, beta_true)
#     u_i ~ Geo(p_i)
#     kTrue_i = u_i + 2   (unconditional floor at 2)
#   tree ~ Gamma(2, 2/50) * Dirichlet
#   y_i ~ JC(kTrue_i)
#
# Inference uses Model A beta_geometric prior (u = k' - 2, same branch).
# The (α, β) hyperparameter sampler uses the s/r reparameterisation
# (fix/bg-sr-reparameterisation, c16be46) to traverse the BG posterior ridge.
#
# Rank tests: tree_length, rate_log_sd, kprime_alpha, kprime_beta, kPrime_pooled.
# Pass criterion: AD p > 0.01 for all five.
#
# This supersedes the pre-fix T9 run (tree_length 0.018, α/β/kPrime AD≈0).

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(ape); library(TreeTools)
})

OUT_DIR <- "dev/red-team/heavy-tests/kprime-viability"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
mode <- if (any(args %in% c("--quick", "-q"))) "quick" else "full"
cat(sprintf("[T9 Model A beta-geo] mode=%s\n", mode))

if (mode == "quick") {
  N_SIM  <- 30L
  N_TIP  <- 8L
  N_CHAR <- 30L
  N_ITER <- 1000L
  N_THIN <- 10L
  N_WARM <- 500L
} else {
  N_SIM  <- 200L
  N_TIP  <- 8L
  N_CHAR <- 30L
  N_ITER <- 12000L
  N_THIN <- 120L
  N_WARM <- 4000L
}
L_SAMPLES <- N_ITER %/% N_THIN

EXPSTEPS_FIXED <- 50
TREE_SHAPE     <- 2
K_MAX_PRIOR    <- 30L
seedBase       <- 20260528L + 200000L

.simTree <- function(nTip) {
  tr <- ape::rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  g <- rgamma(nrow(tr$edge), 1); rels <- g / sum(g)
  tl <- rgamma(1, TREE_SHAPE, TREE_SHAPE / EXPSTEPS_FIXED)
  tr$edge.length <- tl * rels; tr
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
    if (runif(1L) < pSame) states[ch] <- states[pa]
    else states[ch] <- sample(setdiff(seq.int(0L, kTrue - 1L), states[pa]), 1L)
  }
  states[seq_len(nTip)]
}
.canon <- function(v) {
  u <- sort(unique(v)); out <- match(v, u) - 1L
  attr(out, "kObs") <- length(u); out
}
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

.runOneSim <- function(sim_id, seed) {
  set.seed(seed)
  # Joint forward draw under Model A beta_geometric
  alpha_true     <- rexp(1, rate = 1)
  beta_true      <- rexp(1, rate = 1)
  rateLogSd_true <- rgamma(1, shape = 1, rate = 1)
  tr <- .simTree(N_TIP); tl_true <- sum(tr$edge.length)
  # Per-character kTrue: p_i ~ Beta(alpha, beta); u_i ~ Geo(p_i); k = u + 2
  pv    <- rbeta(N_CHAR, alpha_true, beta_true)
  pv    <- pmax(pmin(pv, 1 - 1e-9), 1e-9)
  u_true <- rgeom(N_CHAR, pv)
  kTrue  <- pmin(u_true + 2L, K_MAX_PRIOR)
  sim_mat <- matrix(NA_integer_, N_TIP, N_CHAR,
                    dimnames = list(tr$tip.label, NULL))
  kObs <- integer(N_CHAR)
  for (j in seq_len(N_CHAR)) {
    raw <- .simJCchar(tr, kTrue[j]); cv <- .canon(raw)
    sim_mat[, j] <- cv; kObs[j] <- attr(cv, "kObs")
  }
  keep <- kObs >= 2L
  if (sum(keep) < 5L) return(list(skipped = TRUE, reason = "few_var_chars"))
  sim_mat <- sim_mat[, keep, drop = FALSE]
  kTrue <- kTrue[keep]; kObs <- kObs[keep]; n_char <- sum(keep)
  pd <- TreeTools::MatrixToPhyDat(sim_mat)
  mkd <- MkPrimeData(pd)
  start_tree <- tr
  start_tree$edge.length <- rep_len(0.1, nrow(tr$edge))
  model <- MkPrimeModel(
    coding = "variable", nCat = 1L,
    kPrimePrior = "beta_geometric",
    kprimeAlpha = 1.0,
    kprimeBeta  = 1.0,
    expSteps = EXPSTEPS_FIXED
  )
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
  if (!is.null(res$error))
    return(list(skipped = TRUE, reason = paste0("mcmc:", res$error)))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  samples <- res$samples
  if (is.null(samples) || nrow(samples) < 10L)
    return(list(skipped = TRUE, reason = "no_samples"))
  kp_cols <- grep("^kPrime_", colnames(samples))
  rk <- list()
  rk$tree_length  <- .rankPlain(tl_true,         samples[, "tree_length"])
  rk$rate_log_sd  <- .rankPlain(rateLogSd_true,  samples[, "rate_log_sd"])
  rk$kprime_alpha <- .rankPlain(alpha_true,      samples[, "kprime_alpha"])
  rk$kprime_beta  <- .rankPlain(beta_true,       samples[, "kprime_beta"])
  rk$kPrime_per_char <- vapply(seq_len(n_char), function(j) {
    .rankRandomised(kTrue[j], samples[, kp_cols[j]])
  }, numeric(1))
  list(skipped = FALSE, L = nrow(samples), wall = dt,
       alpha_true = alpha_true, beta_true = beta_true,
       tl_true = tl_true, rateLogSd_true = rateLogSd_true,
       kTrue = kTrue, kObs = kObs, ranks = rk)
}

cat(sprintf("\n=== T9: Model A beta_geometric SBC (post s/r fix) ===\n"))
cat(sprintf("    N_SIM=%d, N_ITER=%d, N_WARM=%d, N_THIN=%d, L=%d\n",
            N_SIM, N_ITER, N_WARM, N_THIN, L_SAMPLES))
configDir <- file.path(OUT_DIR, "T9-betaGeo")
dir.create(configDir, recursive = TRUE, showWarnings = FALSE)
sims <- vector("list", N_SIM)
for (i in seq_len(N_SIM)) {
  seed <- seedBase + i
  cat(sprintf("  sim %3d/%d ... ", i, N_SIM))
  s <- tryCatch(.runOneSim(i, seed),
                error = function(e) list(skipped = TRUE,
                                          reason = paste0("trycatch:", conditionMessage(e))))
  if (isTRUE(s$skipped)) cat(sprintf("SKIP (%s)\n", s$reason))
  else                    cat(sprintf("L=%d wall=%.1fs\n", s$L, s$wall))
  sims[[i]] <- s
}
good <- sapply(sims, function(s) !isTRUE(s$skipped))
cat(sprintf("\nGood sims: %d / %d\n", sum(good), N_SIM))
if (sum(good) < 20L) {
  cat("Too few good sims; saving and exiting.\n")
  saveRDS(sims, file.path(configDir, "sims.rds"))
  quit(status = 1)
}
goods <- sims[good]
L_used <- max(sapply(goods, `[[`, "L"))
tl_ranks  <- sapply(goods, function(s) s$ranks$tree_length)
rls_ranks <- sapply(goods, function(s) s$ranks$rate_log_sd)
a_ranks   <- sapply(goods, function(s) s$ranks$kprime_alpha)
b_ranks   <- sapply(goods, function(s) s$ranks$kprime_beta)
kp_ranks  <- unlist(lapply(goods, function(s) s$ranks$kPrime_per_char))
.normRanks <- function(r, L) (r + 0.5) / (L + 1)
.adP <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4L) return(NA_real_)
  if (requireNamespace("goftest", quietly = TRUE))
    suppressWarnings(goftest::ad.test(x, null = "punif")$p.value)
  else suppressWarnings(ks.test(x, "punif")$p.value)
}
ad <- list(
  tree_length   = .adP(.normRanks(tl_ranks,  L_used)),
  rate_log_sd   = .adP(.normRanks(rls_ranks, L_used)),
  kprime_alpha  = .adP(.normRanks(a_ranks,   L_used)),
  kprime_beta   = .adP(.normRanks(b_ranks,   L_used)),
  kPrime_pooled = .adP(.normRanks(kp_ranks,  L_used))
)
cat(sprintf("\nT9 AD p-values (beta_geometric, Model A, post s/r fix):\n"))
for (nm in names(ad)) {
  cat(sprintf("  %-15s : %.4f  %s\n", nm, ad[[nm]],
              if (is.na(ad[[nm]])) "NA"
              else if (ad[[nm]] > 0.01) "PASS" else "FAIL"))
}
saveRDS(sims, file.path(configDir, "sims.rds"))
saveRDS(list(ranks = list(tree_length = tl_ranks, rate_log_sd = rls_ranks,
                          kprime_alpha = a_ranks, kprime_beta = b_ranks,
                          kPrime_pooled = kp_ranks),
             L = L_used, ad = ad),
        file.path(configDir, "ranks.rds"))
# Rank histograms
png(file.path(configDir, "rank-histograms.png"), width = 1500, height = 600,
    res = 110)
op <- par(mfrow = c(1, 5), mar = c(4, 4, 3, 1))
for (nm in c("tree_length", "rate_log_sd", "kprime_alpha", "kprime_beta",
             "kPrime_pooled")) {
  r <- switch(nm,
    tree_length   = tl_ranks,  rate_log_sd  = rls_ranks,
    kprime_alpha  = a_ranks,   kprime_beta  = b_ranks,
    kPrime_pooled = kp_ranks)
  hist(r, breaks = 30, freq = FALSE, col = "lightgrey", border = "white",
       main = sprintf("%s\nAD p=%.4f", nm, ad[[nm]]), xlab = "rank")
  abline(h = 1 / L_used, col = "red", lwd = 2)
}
par(op)
dev.off()
# Verdict
sink(file.path(OUT_DIR, "T9-summary.md"))
cat("# T9 — Model A beta_geometric SBC (post s/r-reparameterisation fix)\n\n")
cat(sprintf("Mode: %s.  N_SIM=%d (good=%d), N_TIP=%d, N_CHAR target=%d, L=%d.\n\n",
            mode, N_SIM, sum(good), N_TIP, N_CHAR, L_SAMPLES))
cat("Inference: Model A beta_geometric prior (u = k' − 2, unconditional).\n")
cat("Sampler: (s, r) reparameterised slice for (α, β) — fix/bg-sr-reparameterisation.\n")
cat("Forward: alpha_true, beta_true ~ Exp(1); p_i ~ Beta; u_i ~ Geo(p_i); kTrue = u+2.\n\n")
cat("## Results\n\n")
for (nm in names(ad)) {
  cat(sprintf("- %s : AD p = %.4f  %s\n", nm, ad[[nm]],
              if (is.na(ad[[nm]])) "NA"
              else if (ad[[nm]] > 0.01) "**PASS**" else "**FAIL**"))
}
cat(sprintf("\nOverall: %s\n",
            if (all(sapply(ad, function(x) is.finite(x) && x > 0.01)))
              "**PASS**" else "**FAIL**"))
sink()
cat(sprintf("\n[T9] output: %s/T9-summary.md  %s/T9-betaGeo/\n", OUT_DIR, OUT_DIR))
