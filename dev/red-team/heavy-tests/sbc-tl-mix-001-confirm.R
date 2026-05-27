#!/usr/bin/env Rscript
# SBC-TL-MIX-001 confirmation: forward simulator's `rev(seq_len(...))` is
# postorder iteration, which reads `states[pa]` before the parent state has
# been drawn. Standard tree-simulation idiom is preorder (forward iteration
# on a Preorder edge list). This script re-simulates sim 15 with the FIXED
# direction, refits the MLE profile, and reports whether the MLE now
# matches tl_true.

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(ape); library(TreeTools)
})

EXPSTEPS_FIXED <- 50
TREE_SHAPE     <- 2
N_TIP          <- 8L
N_CHAR         <- 30L
simIdx         <- 15L
seed           <- 20260526L + 1000L * 2L + simIdx

.simTree <- function(nTip) {
  tr <- ape::rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  nEdge <- nrow(tr$edge)
  tl <- stats::rgamma(1, shape = TREE_SHAPE, rate = TREE_SHAPE / EXPSTEPS_FIXED)
  w <- stats::rexp(nEdge, rate = 1)
  tr$edge.length <- tl * w / sum(w)
  TreeTools::Preorder(tr)
}

# BUGGY version (matches sbc.R:139)
.simJCchar_buggy <- function(tree, kTrue) {
  nTip <- length(tree$tip.label)
  states <- integer(2L * nTip - 1L)
  rootIdx <- nTip + 1L
  states[rootIdx] <- sample.int(kTrue, 1L) - 1L
  edges <- tree$edge; el <- tree$edge.length
  for (e in rev(seq_len(nrow(edges)))) {
    pa <- edges[e, 1L]; ch <- edges[e, 2L]; t <- el[e]
    pSame <- 1 / kTrue + (1 - 1 / kTrue) * exp(-kTrue * t / (kTrue - 1))
    if (runif(1L) < pSame) states[ch] <- states[pa]
    else states[ch] <- sample(setdiff(seq.int(0L, kTrue - 1L), states[pa]), 1L)
  }
  states[seq_len(nTip)]
}
# FIXED version: forward iteration in preorder.
.simJCchar_fixed <- function(tree, kTrue) {
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
.canonicaliseLabels <- function(vec) {
  uvals <- sort(unique(vec))
  out <- match(vec, uvals) - 1L
  attr(out, "kObs") <- length(uvals); out
}

# Reproduce sim 15 setup
set.seed(seed)
p_true <- stats::rbeta(1, 1, 1)
rateLogSd_true <- stats::rgamma(1, shape = 1, rate = 1)
true_tree <- .simTree(N_TIP)
tl_true <- sum(true_tree$edge.length)

# Sanity: print first few edges to confirm Preorder ordering (parents introduced
# before children — edges[1, 1] should be the root index = nTip + 1).
cat(sprintf("[confirm] first edge (root, child) = (%d, %d); root expected = %d\n",
            true_tree$edge[1L, 1L], true_tree$edge[1L, 2L], N_TIP + 1L))
cat(sprintf("[confirm] tl_true = %.4f\n", tl_true))

profile_TL <- function(mkd, n_char, TL, true_tree, tl_true) {
  relBrTruth <- true_tree$edge.length / tl_true
  tr <- true_tree; tr$edge.length <- TL * relBrTruth
  tr <- TreeTools::Preorder(tr)
  MkpLogLikelihood(tr, mkd, kPrime = rep(2L, n_char), rate_loss = 1.0,
                   rate_log_sd = 0, nCat = 1L, coding = "variable",
                   relabel = TRUE)
}

simAndProfile <- function(sim_fn, label, useScratchSeed = TRUE) {
  # Use a *separate* RNG stream for the per-character draws so we can compare
  # the two simulators on the same set of attempted characters.
  if (useScratchSeed) set.seed(seed + 999L)
  kSim <- rep(2L, N_CHAR)
  sim_mat <- matrix(NA_integer_, N_TIP, N_CHAR,
                    dimnames = list(true_tree$tip.label, NULL))
  kObs <- integer(N_CHAR)
  for (j in seq_len(N_CHAR)) {
    raw <- sim_fn(true_tree, kSim[j])
    canon <- .canonicaliseLabels(raw)
    sim_mat[, j] <- canon
    kObs[j] <- attr(canon, "kObs")
  }
  keep <- kObs >= 2L
  n_char <- sum(keep)
  sim_mat <- sim_mat[, keep, drop = FALSE]
  pd <- TreeTools::MatrixToPhyDat(sim_mat)
  mkd <- MkPrimeData(pd)

  TL_grid <- c(0.5, 1, 1.2, 1.5, 2, 3, 4, tl_true, 5, 8, 12, 20, 50)
  logL <- vapply(TL_grid, function(TL) profile_TL(mkd, n_char, TL, true_tree, tl_true),
                 numeric(1))
  mle_idx <- which.max(logL)
  opt <- optimize(function(TL) profile_TL(mkd, n_char, TL, true_tree, tl_true),
                  interval = c(0.3, 30), maximum = TRUE, tol = 1e-4)
  cat(sprintf("[%s] n_char kept = %d/%d  MLE_TL = %.4f  logL(MLE)=%.3f  logL(tl_true=%.3f)=%.3f\n",
              label, n_char, N_CHAR, opt$maximum, opt$objective, tl_true,
              profile_TL(mkd, n_char, tl_true, true_tree, tl_true)))
  invisible(list(mle = opt$maximum, n_char = n_char, logL = logL, grid = TL_grid))
}

cat("\n[confirm] === profile on data from BUGGY simulator (current sbc.R) ===\n")
buggy <- simAndProfile(.simJCchar_buggy, "buggy")

cat("\n[confirm] === profile on data from FIXED simulator (forward iter) ===\n")
fixed <- simAndProfile(.simJCchar_fixed, "fixed")

# Repeat across multiple seeds to check the bias is systematic
cat("\n[confirm] === systematic check: 30 distinct seeds, both simulators ===\n")
mleBuggy <- numeric(30); mleFixed <- numeric(30); tlTrues <- numeric(30)
for (s in seq_len(30)) {
  set.seed(20260526L + 1000L * 2L + s + 100L)
  p_true2         <- stats::rbeta(1, 1, 1)
  rateLogSd_true2 <- stats::rgamma(1, shape = 1, rate = 1)
  tr2 <- .simTree(N_TIP); tl2 <- sum(tr2$edge.length); tlTrues[s] <- tl2

  simOne <- function(simFn) {
    sm <- matrix(NA_integer_, N_TIP, N_CHAR, dimnames = list(tr2$tip.label, NULL))
    kO <- integer(N_CHAR)
    for (j in seq_len(N_CHAR)) {
      raw <- simFn(tr2, 2L)
      ca <- .canonicaliseLabels(raw); sm[, j] <- ca; kO[j] <- attr(ca, "kObs")
    }
    keep <- kO >= 2L
    pd <- TreeTools::MatrixToPhyDat(sm[, keep, drop = FALSE])
    md <- MkPrimeData(pd)
    opt <- optimize(function(TL) profile_TL(md, sum(keep), TL, tr2, tl2),
                    interval = c(0.1, 200), maximum = TRUE, tol = 1e-3)
    opt$maximum
  }
  scratch <- 20260526L + 1000L * 2L + s + 100L + 99999L
  set.seed(scratch); mleBuggy[s] <- simOne(.simJCchar_buggy)
  set.seed(scratch); mleFixed[s] <- simOne(.simJCchar_fixed)
}
cat(sprintf("[confirm] buggy: median MLE/tl_true = %.3f  (mean = %.3f)\n",
            median(mleBuggy / tlTrues), mean(mleBuggy / tlTrues)))
cat(sprintf("[confirm] fixed: median MLE/tl_true = %.3f  (mean = %.3f)\n",
            median(mleFixed / tlTrues), mean(mleFixed / tlTrues)))
cat(sprintf("[confirm] buggy: %d/30 sims have MLE < tl_true\n", sum(mleBuggy < tlTrues)))
cat(sprintf("[confirm] fixed: %d/30 sims have MLE < tl_true\n", sum(mleFixed < tlTrues)))

# Wilcoxon: is buggy MLE systematically below truth (one-sided)?
cat(sprintf("[confirm] Wilcoxon signed-rank (buggy MLE vs tl_true, alt=less): p = %.2e\n",
            tryCatch(wilcox.test(mleBuggy, tlTrues, paired = TRUE, alternative = "less")$p.value,
                     error = function(e) NA_real_)))
cat(sprintf("[confirm] Wilcoxon signed-rank (fixed MLE vs tl_true, two-sided): p = %.3f\n",
            tryCatch(wilcox.test(mleFixed, tlTrues, paired = TRUE)$p.value,
                     error = function(e) NA_real_)))

saveRDS(list(buggy = buggy, fixed = fixed,
             mleBuggy = mleBuggy, mleFixed = mleFixed, tlTrues = tlTrues),
        "dev/red-team/heavy-tests/sbc-tl-mix-001-confirm.rds")
cat("[confirm] saved sbc-tl-mix-001-confirm.rds\n")
