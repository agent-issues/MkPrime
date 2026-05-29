#!/usr/bin/env Rscript
# marginal-k-slice-freeze-repro.R
#
# Discriminating local test for MARGINAL-K-SLICE-001 (the marginal_k
# tree_length freeze). Runs ONE simulated dataset (16 tips, 100 chars,
# tl~1.4, Lewis-Mkv / Model I-a forward, fixTopology) through BOTH
# likelihoodMode = "sampled_k" and "marginal_k" with otherwise-identical
# config, and reports the tree_length / p / rate_log_sd posterior mean + sd.
#
# PRE-FIX  : marginal_k tree_length sd ~ 0.000 (frozen, 1 unique value).
# POST-FIX : marginal_k tree_length sd ~ 0.08, comparable to sampled_k, and
#            tree_length posterior mean near truth (not pinned too long).
#
# This mirrors dev/red-team/heavy-tests/marginal-k/T-SBC-marginal-geometric.R
# per-sim forward exactly (so a pass here predicts the Hamilton SBC).

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(ape)
  library(TreeTools)
})

EXPSTEPS_FIXED <- 1.4
TREE_SHAPE     <- 2
K_MAX_PRIOR    <- 30L
A_PRIOR        <- 1
B_PRIOR        <- 1
N_TIP          <- 16L
N_CHAR         <- 100L
# Modest local budget — enough to expose freeze-vs-mixing in minutes.
N_ITER         <- 6000L
N_THIN         <- 30L
N_WARM         <- 2000L

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

simulate_one <- function(seed) {
  set.seed(seed)
  p_true         <- stats::rbeta(1, A_PRIOR, B_PRIOR)
  rateLogSd_true <- stats::rgamma(1, shape = 1, rate = 1)
  tr             <- .simTree(N_TIP)
  tl_true        <- sum(tr$edge.length)
  u_true <- stats::rgeom(N_CHAR, p_true)
  kTrue  <- pmin(2L + u_true, K_MAX_PRIOR)
  sim_mat <- matrix(NA_integer_, N_TIP, N_CHAR,
                    dimnames = list(tr$tip.label, NULL))
  kObs <- integer(N_CHAR)
  for (j in seq_len(N_CHAR)) {
    repeat {
      cv <- .canon(.simJCchar(tr, kTrue[j]))
      if (attr(cv, "kObs") >= 2L) break
    }
    sim_mat[, j] <- cv
    kObs[j] <- attr(cv, "kObs")
  }
  list(tr = tr, tl_true = tl_true, p_true = p_true,
       rateLogSd_true = rateLogSd_true,
       mkd = MkPrimeData(TreeTools::MatrixToPhyDat(sim_mat)))
}

run_mode <- function(sim, mode, seed) {
  model <- suppressMessages(MkPrimeModel(
    coding         = "variable",
    nCat           = 1L,
    kPrimePrior    = "geometric",
    likelihoodMode = mode,
    priorVariant   = "unconditional",
    kprimeHyperA   = A_PRIOR,
    kprimeHyperB   = B_PRIOR,
    expSteps       = EXPSTEPS_FIXED
  ))
  mcmc <- MkPrimeMCMC(
    nIter = N_ITER, thin = N_THIN,
    minWarmup = N_WARM, maxWarmup = N_WARM,
    autoTune = FALSE, nRuns = 1L, nChains = 1L
  )
  start_tree <- sim$tr
  start_tree$edge.length <- rep_len(0.1, nrow(sim$tr$edge))
  set.seed(seed + 7L)
  res <- suppressMessages(suppressWarnings(
    RunMkPrime(sim$mkd, start_tree, model = model, mcmc = mcmc,
               fixTopology = TRUE, overwrite = TRUE)
  ))
  res$samples
}

summarise <- function(s, par) {
  if (is.null(s) || !(par %in% colnames(s))) return(c(mean = NA, sd = NA, nuniq = NA))
  v <- s[, par]
  c(mean = mean(v), sd = stats::sd(v), nuniq = length(unique(v)))
}

SEEDS <- c(20260528L, 20260529L, 20260530L)
cat(sprintf("%-10s %-10s | %-28s | %-22s | %-22s\n",
            "seed", "mode", "tree_length (mean/sd/nuniq)",
            "p (mean/sd)", "rate_log_sd (mean/sd)"))
cat(strrep("-", 104), "\n")
for (sd_ in SEEDS) {
  sim <- simulate_one(sd_)
  for (mode in c("sampled_k", "marginal_k")) {
    s <- run_mode(sim, mode, sd_)
    tl <- summarise(s, "tree_length")
    pp <- summarise(s, "p")
    rl <- summarise(s, "rate_log_sd")
    cat(sprintf("%-10d %-10s | %7.4f / %6.4f / %4d      | %6.4f / %6.4f      | %6.4f / %6.4f\n",
                sd_, mode, tl["mean"], tl["sd"], tl["nuniq"],
                pp["mean"], pp["sd"], rl["mean"], rl["sd"]))
  }
  cat(sprintf("  (truth: tl=%.4f  p=%.4f  rateLogSd=%.4f)\n",
              sim$tl_true, sim$p_true, sim$rateLogSd_true))
}
