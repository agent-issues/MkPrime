# T-OVL: Posterior overlap, sampled-k vs marginal-k geometric arm
# ===============================================================
#
# Plan reference: dev/notes/2026-05-28-marginal-k-plan.md §7.2 (Rung 2).
#
# Validates that the v1 marginal-k geometric mode delivers the same
# marginal posterior on (tree_length, rate_log_sd, p) as the legacy
# sampled-k mode. By Rao-Blackwell, the two posteriors must agree in the
# limit; this driver runs short matched-seed chains across a small
# (n_tip × n_char) grid and applies KS tests on the marginal posteriors.
#
# Pass bar: KS p > 0.01 on every parameter in every cell.
#
# DO NOT RUN INLINE — this is a heavy test. Submit to Hamilton via a
# companion submit-marginal-k-ovl.sh (TODO; not yet written). See
# `dev/red-team/heavy-tests/submit-sbc.sh` for the canonical pre-build
# pattern and the `feedback_pkgload_prebuild` memory file before running.

suppressPackageStartupMessages({
  library("MkPrime")
  library("ape")
  library("TreeTools")
})

# ---------------------------------------------------------------------------
# Knobs
# ---------------------------------------------------------------------------

GRID_NTIP  <- c(8L, 16L)
GRID_NCHAR <- c(16L, 32L)
N_REP      <- 4L          # replicate datasets per cell
N_ITER     <- 12000L
N_WARM     <- 4000L
KS_BAR     <- 0.01

OUT_DIR <- "dev/red-team/heavy-tests/marginal-k"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Forward simulation
# ---------------------------------------------------------------------------

simulate_dataset <- function(nTip, nChar, seed) {
  set.seed(seed)
  # Random unrooted binary tree with uniform branch lengths
  tree <- ape::rtree(nTip, br = function(n) runif(n, 0.02, 0.15))
  # All binary characters (kObs = 2) so the geometric arm is the
  # appropriate prior.
  mat <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                nrow = nTip, ncol = nChar,
                dimnames = list(tree$tip.label, NULL))
  list(tree = tree, mkd = MkPrimeData(MatrixToPhyDat(mat)))
}

# ---------------------------------------------------------------------------
# One-cell run: both modes on the same dataset, same seed.
# ---------------------------------------------------------------------------

run_cell <- function(nTip, nChar, rep) {
  seed <- as.integer(1e3 * rep + 10 * nTip + nChar)
  sim  <- simulate_dataset(nTip, nChar, seed)

  cell_name <- sprintf("n%02d_c%02d_r%02d", nTip, nChar, rep)
  message("[T-OVL] running ", cell_name, " ...")

  run_one <- function(mode) {
    model <- MkPrimeModel(kPrimePrior  = "geometric",
                          likelihoodMode = mode,
                          coding = "variable")
    RunMkPrime(
      sim$tree, sim$mkd, model,
      mcmc = list(nIter   = N_ITER,
                  warmup  = N_WARM,
                  thin    = "auto",
                  nChains = 1L),
      seed = seed
    )
  }

  fit_s <- run_one("sampled_k")
  fit_m <- run_one("marginal_k")

  list(cell    = cell_name,
       sampled = fit_s,
       marginal = fit_m,
       seed    = seed)
}

# ---------------------------------------------------------------------------
# KS test per parameter
# ---------------------------------------------------------------------------

ks_per_param <- function(fit_s, fit_m, params = c("tree_length",
                                                    "rate_log_sd",
                                                    "p")) {
  out <- list()
  for (pname in params) {
    xs <- fit_s$samples[, pname]
    xm <- fit_m$samples[, pname]
    ks <- suppressWarnings(stats::ks.test(xs, xm))
    out[[pname]] <- list(statistic = unname(ks$statistic),
                         p.value   = unname(ks$p.value))
  }
  out
}

# ---------------------------------------------------------------------------
# Grid sweep
# ---------------------------------------------------------------------------

if (sys.nframe() == 0L) {
  results <- list()
  for (nTip in GRID_NTIP)
    for (nChar in GRID_NCHAR)
      for (rep in seq_len(N_REP)) {
        cell <- run_cell(nTip, nChar, rep)
        cell$ks <- ks_per_param(cell$sampled, cell$marginal)
        results[[cell$cell]] <- cell
        saveRDS(results, file.path(OUT_DIR, "T-OVL-results.rds"))
      }

  # Verdict
  all_pass <- TRUE
  for (cn in names(results)) {
    ks <- results[[cn]]$ks
    for (pn in names(ks)) {
      pv <- ks[[pn]]$p.value
      pass <- pv > KS_BAR
      all_pass <- all_pass && pass
      cat(sprintf("[T-OVL] %s : %s : KS p = %.4f : %s\n",
                  cn, pn, pv, if (pass) "PASS" else "FAIL"))
    }
  }
  cat(if (all_pass) "[T-OVL] ALL PASS\n" else "[T-OVL] SOME FAIL\n")
  writeLines(if (all_pass) "PASS" else "FAIL",
             file.path(OUT_DIR, "T-OVL-verdict.txt"))
}
