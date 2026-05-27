#!/usr/bin/env Rscript
# SBC-TL-MIX-001 diagnostic.
#
# Sim 15 of v5 Mkp_geometric (seed 20262541) — 8 tips, 30 chars binary
# (28 informative kept), true topology fixed. Posterior on tree_length
# stationarises at ~1.4 vs tl_true ~ 4.13. This script discriminates
# between (1) sampler stuck and (2) likelihood/prior bug by:
#
#   (a) profiling logL + logPrior on a TL grid with rel_br fixed at truth
#   (b) instrumenting per-move accept rates and the adapted scale_tree_length
#       tuning on a longer chain
#
# Output:
#   dev/red-team/heavy-tests/sbc-tl-mix-001-diagnose.rds
#   dev/red-team/heavy-tests/sbc-tl-mix-001-diagnose.png
#   stdout summary

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(ape)
  library(TreeTools)
})

# ---- Setup: exact sim 15 reproduction (mirrors sbc-warmup-trace.R) ----
EXPSTEPS_FIXED <- 50
TREE_SHAPE     <- 2
N_TIP          <- 8L
N_CHAR         <- 30L
simIdx         <- 15L
seed           <- 20260526L + 1000L * 2L + simIdx   # = 20262541

cat(sprintf("[diag] seed=%d\n", seed))

.simTree <- function(nTip) {
  tr <- ape::rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  nEdge <- nrow(tr$edge)
  tl <- stats::rgamma(1, shape = TREE_SHAPE, rate = TREE_SHAPE / EXPSTEPS_FIXED)
  w <- stats::rexp(nEdge, rate = 1)
  tr$edge.length <- tl * w / sum(w)
  TreeTools::Preorder(tr)
}
.simJCchar <- function(tree, kTrue) {
  nTip <- length(tree$tip.label)
  states <- integer(2L * nTip - 1L)
  rootIdx <- nTip + 1L
  states[rootIdx] <- sample.int(kTrue, 1L) - 1L
  edges <- tree$edge
  el    <- tree$edge.length
  for (e in rev(seq_len(nrow(edges)))) {
    pa <- edges[e, 1L]; ch <- edges[e, 2L]; t <- el[e]
    pSame <- 1 / kTrue + (1 - 1 / kTrue) * exp(-kTrue * t / (kTrue - 1))
    if (runif(1L) < pSame) states[ch] <- states[pa]
    else states[ch] <- sample(setdiff(seq.int(0L, kTrue - 1L), states[pa]), 1L)
  }
  states[seq_len(nTip)]
}
.canonicaliseLabels <- function(vec) {
  uvals <- sort(unique(vec))
  out <- match(vec, uvals) - 1L
  attr(out, "kObs") <- length(uvals)
  out
}

set.seed(seed)
p_true         <- stats::rbeta(1, 1, 1)
rateLogSd_true <- stats::rgamma(1, shape = 1, rate = 1)
true_tree      <- .simTree(N_TIP)
tl_true        <- sum(true_tree$edge.length)
cat(sprintf("[diag] tl_true = %.4f  nEdge = %d\n",
            tl_true, nrow(true_tree$edge)))

kSim <- rep(2L, N_CHAR)
sim_mat <- matrix(NA_integer_, N_TIP, N_CHAR,
                  dimnames = list(true_tree$tip.label, NULL))
kObs <- integer(N_CHAR)
for (j in seq_len(N_CHAR)) {
  raw   <- .simJCchar(true_tree, kSim[j])
  canon <- .canonicaliseLabels(raw)
  sim_mat[, j] <- canon
  kObs[j] <- attr(canon, "kObs")
}
keep    <- kObs >= 2L
sim_mat <- sim_mat[, keep, drop = FALSE]
kObs    <- kObs[keep]
n_char  <- sum(keep)
cat(sprintf("[diag] n_char kept = %d/%d\n", n_char, N_CHAR))

pd  <- TreeTools::MatrixToPhyDat(sim_mat)
mkd <- MkPrimeData(pd)

# Parsimony score on truth (info-theoretic floor on TL signal)
fitchScore <- tryCatch(
  TreeTools::CharacterLength(true_tree, pd, compress = FALSE) |> sum(),
  error = function(e) NA_real_)
cat(sprintf("[diag] truth parsimony score = %s (over %d chars)\n",
            as.character(fitchScore), n_char))

# ============================================================
# PART A — likelihood profile on TL grid
# ============================================================
cat("\n[diag] ==== PART A: logL profile ====\n")

relBrTruth <- true_tree$edge.length / tl_true   # truth's normalised weights
stopifnot(abs(sum(relBrTruth) - 1) < 1e-9)

TL_grid <- c(0.5, 1.0, 1.2, 1.4, 1.5, 2.0, 2.5, 3.0, 4.0, tl_true,
             5.0, 6.0, 8.0, 12.0, 20.0, 50.0)
TL_grid <- sort(unique(TL_grid))

profile_at_TL <- function(TL, relabel = TRUE) {
  tr <- true_tree
  tr$edge.length <- TL * relBrTruth
  tr <- TreeTools::Preorder(tr)
  logL <- MkpLogLikelihood(
    tr, mkd,
    kPrime      = rep(2L, n_char),
    rate_loss   = 1.0,
    rate_log_sd = 0,
    nCat        = 1L,
    coding      = "variable",
    relabel     = relabel
  )
  logPrior_TL <- dgamma(TL, shape = TREE_SHAPE,
                        rate = TREE_SHAPE / EXPSTEPS_FIXED, log = TRUE)
  c(logL = logL, logPrior = logPrior_TL, logPost = logL + logPrior_TL)
}

profileMat <- t(vapply(TL_grid, profile_at_TL, numeric(3), relabel = TRUE))
rownames(profileMat) <- sprintf("%.3f", TL_grid)
print(round(profileMat, 3))

mleIdx  <- which.max(profileMat[, "logL"])
mapIdx  <- which.max(profileMat[, "logPost"])
cat(sprintf("\n[diag] MLE  TL (likelihood-only) = %.3f  (logL = %.3f)\n",
            TL_grid[mleIdx], profileMat[mleIdx, "logL"]))
cat(sprintf("[diag] MAP TL (logL + logPrior)  = %.3f  (logPost = %.3f)\n",
            TL_grid[mapIdx], profileMat[mapIdx, "logPost"]))
cat(sprintf("[diag] logL(truth=%.3f) = %.3f   logL(chain mode≈1.5) = %.3f   delta = %.3f\n",
            tl_true,
            profileMat[which.min(abs(TL_grid - tl_true)), "logL"],
            profileMat[which.min(abs(TL_grid - 1.5)),  "logL"],
            profileMat[which.min(abs(TL_grid - tl_true)), "logL"] -
              profileMat[which.min(abs(TL_grid - 1.5)),  "logL"]))

# Probe: same profile with relabel=FALSE (plain Mk) — discriminates
# whether the Mk' relabel correction is responsible.
profileMatPlainMk <- t(vapply(TL_grid, profile_at_TL, numeric(3),
                              relabel = FALSE))
rownames(profileMatPlainMk) <- sprintf("%.3f", TL_grid)
mleIdxPM <- which.max(profileMatPlainMk[, "logL"])
cat(sprintf("[diag] (relabel=FALSE / plain Mk) MLE TL = %.3f\n",
            TL_grid[mleIdxPM]))

# Refine the MLE with optimize() around the grid maximum
brentInt <- c(max(0.1, TL_grid[max(1, mleIdx - 1)]),
              TL_grid[min(length(TL_grid), mleIdx + 1)])
opt <- optimize(
  function(TL) profile_at_TL(TL, relabel = TRUE)["logL"],
  interval = brentInt, maximum = TRUE, tol = 1e-4
)
cat(sprintf("[diag] Brent-refined MLE TL ∈ [%.2f, %.2f] = %.4f (logL = %.3f)\n",
            brentInt[1], brentInt[2], opt$maximum, opt$objective))

# ============================================================
# PART B — instrumented MCMC: per-move accept rates + tuning
# ============================================================
cat("\n[diag] ==== PART B: accept rates + adapted scale ====\n")

# Long enough warmup that adaptive scheduler can move the proposal scale
N_WARM_INST <- 5000L
N_ITER_INST <- 20000L
N_THIN_INST <- 100L

st <- true_tree
st <- TreeTools::Preorder(st)

model <- MkPrimeModel(
  coding       = "variable",
  nCat         = 1L,
  kPrimePrior  = "geometric",
  expSteps     = EXPSTEPS_FIXED,
  kprimeHyperA = 1, kprimeHyperB = 1,
  kprimeAlpha  = 2, kprimeBeta  = 2
)
mcmc <- MkPrimeMCMC(
  nIter     = N_ITER_INST,
  thin      = N_THIN_INST,
  minWarmup = N_WARM_INST,
  maxWarmup = N_WARM_INST,
  autoTune  = TRUE,           # exercise the scale adaptation
  nRuns     = 1L, nChains = 1L
)

cat(sprintf("[diag] running RunMkPrime nIter=%d thin=%d warm=%d autoTune=TRUE\n",
            N_ITER_INST, N_THIN_INST, N_WARM_INST))
t0 <- Sys.time()
res <- suppressMessages(suppressWarnings(
  RunMkPrime(mkd, st, model = model, mcmc = mcmc,
             fixTopology = TRUE, overwrite = TRUE)
))
cat(sprintf("[diag] mcmc wall = %.1fs   nSamples = %d\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            nrow(res$samples)))

cat("\n[diag] --- per-move acceptance ---\n")
print(round(res$acceptance, 4))

cat("\n[diag] --- final adapted tuning (selected) ---\n")
tn <- res$tuning
sel <- intersect(c("scale_tree_length", "scale_rate_log_sd", "scale_rate_loss",
                   "beta_simplex"),
                 names(tn))
for (k in sel) {
  v <- tn[[k]]
  cat(sprintf("  %-22s = %s\n", k, format(v, digits = 5)))
}

cat("\n[diag] --- final adapted move weights (if available) ---\n")
mw <- res$mcmc$moveWeights %||% res$moveWeights
if (!is.null(mw)) print(round(mw, 4))

tl_samples <- res$samples[, "tree_length"]
cat(sprintf("\n[diag] tree_length samples: mean=%.3f sd=%.3f q025=%.3f q975=%.3f  unique=%d/%d\n",
            mean(tl_samples), sd(tl_samples),
            quantile(tl_samples, 0.025), quantile(tl_samples, 0.975),
            length(unique(tl_samples)), length(tl_samples)))
cat(sprintf("[diag] tl_true = %.3f  → posterior rank in kept = %d/%d\n",
            tl_true, sum(tl_samples < tl_true), length(tl_samples)))

# ============================================================
# Persist
# ============================================================
out <- list(
  seed          = seed,
  tl_true       = tl_true,
  n_char        = n_char,
  fitchScore    = fitchScore,
  TL_grid       = TL_grid,
  profile       = profileMat,
  profilePlainMk = profileMatPlainMk,
  mleIdx        = mleIdx,
  mleTL_Brent   = opt$maximum,
  mleLogL_Brent = opt$objective,
  acceptance    = res$acceptance,
  tuning        = res$tuning,
  moveWeights   = mw,
  tl_samples    = tl_samples
)
saveRDS(out, "dev/red-team/heavy-tests/sbc-tl-mix-001-diagnose.rds")

# Plot
pngPath <- "dev/red-team/heavy-tests/sbc-tl-mix-001-diagnose.png"
grDevices::png(pngPath, width = 1000, height = 700)
op <- graphics::par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot(TL_grid, profileMat[, "logL"], type = "b", pch = 19,
     log = "x",
     xlab = "tree_length (log)", ylab = "logL",
     main = sprintf("logL profile (rel_br fixed at truth, n_char=%d)", n_char))
graphics::abline(v = tl_true, col = "red", lty = 2, lwd = 2)
graphics::abline(v = opt$maximum, col = "blue", lty = 3, lwd = 2)
graphics::legend("bottomright", c("tl_true", "Brent MLE"),
                 col = c("red", "blue"), lty = c(2, 3), lwd = 2, bty = "n")

graphics::hist(tl_samples, breaks = 50, col = "steelblue",
               xlab = "tree_length", main = "posterior samples (autoTune=TRUE, fixTopology=TRUE)")
graphics::abline(v = tl_true, col = "red", lty = 2, lwd = 2)
graphics::abline(v = opt$maximum, col = "blue", lty = 3, lwd = 2)
graphics::par(op)
grDevices::dev.off()
cat(sprintf("[diag] saved %s and .rds\n", pngPath))
