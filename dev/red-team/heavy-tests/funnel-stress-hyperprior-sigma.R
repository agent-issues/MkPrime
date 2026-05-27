#!/usr/bin/env Rscript
#
# Funnel-stress comparison: pooled half-normal hyperprior vs legacy
# independent Gamma prior on per-class σ_c = `class_rate_log_sd`.
#
# Fixture: 5-class partition, sizes (5, 50, 50, 50, 50). Class 1 carries
# only 5 characters and is the small class that drives the σ_c funnel
# under an independent-prior baseline; the larger classes anchor the
# population scale τ under the pooled hyperprior.
#
# What this script does:
#   1. Build the fixture (deterministic seed).
#   2. Run MCMC twice — once with the new hyperprior (default), once with
#      the legacy gamma_independent prior.
#   3. Compute per-class σ_c ESS with `coda::effectiveSize`.
#   4. Print a side-by-side table.
#
# Runs in ~ minutes on a workstation. Standalone Rscript, NOT wired into
# testthat (per repo convention for dev/red-team/heavy-tests/).
#
# Usage (from the package root):
#   Rscript dev/red-team/heavy-tests/funnel-stress-hyperprior-sigma.R

suppressMessages({
  # devtools::load_all so the just-built feature branch sources are used,
  # not whatever MkPrime is installed in the user library.
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION")) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(MkPrime)
  }
  library(TreeTools)
  library(coda)
})


# ----------------------------------------------------------------------------
# Fixture
# ----------------------------------------------------------------------------

build_fixture <- function(seed = 20260527L,
                          class_sizes = c(5L, 50L, 50L, 50L, 50L),
                          n_tip = 12L) {
  set.seed(seed)
  n_char <- sum(class_sizes)
  mat <- matrix(sample(0:1, n_tip * n_char, replace = TRUE),
                nrow = n_tip, ncol = n_char,
                dimnames = list(paste0("t", seq_len(n_tip)), NULL))
  # Force every character to be variable so MkPrimeData accepts them.
  for (j in seq_len(n_char)) {
    if (length(unique(mat[, j])) < 2L) mat[1L, j] <- 1L - mat[1L, j]
  }
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrime::MkPrimeData(pd)

  partition <- rep.int(seq_along(class_sizes), class_sizes)

  tree <- TreeTools::Preorder(
    TreeTools::NJTree(pd, edgeLengths = TRUE) %||%
    TreeTools::RandomTree(pd, root = TRUE)
  )
  if (is.null(tree$edge.length) || any(tree$edge.length <= 0)) {
    tree$edge.length <- rep(0.1, nrow(tree$edge))
  }

  list(mkd = mkd, tree = tree, partition = partition)
}


# ----------------------------------------------------------------------------
# Run one MCMC under a given prior
# ----------------------------------------------------------------------------

run_one <- function(fixture, prior_kind,
                    n_gen = 20000L, n_warmup = 4000L,
                    seed = 1L) {
  set.seed(seed)
  model <- MkPrimeModel(kPrimePrior = "geometric",
                        priorOnClassRateLogSd = prior_kind)
  mcmc <- MkPrimeMCMC(
    nIter      = n_gen,
    minWarmup  = n_warmup,
    maxWarmup  = n_warmup,
    nChains    = 1L,
    nRuns      = 1L,
    thin       = 1L,
    autoTune   = FALSE
  )
  result <- RunMkPrime(
    data      = fixture$mkd,
    tree      = fixture$tree,
    model     = model,
    mcmc      = mcmc,
    partition = fixture$partition,
    unlink    = "shape"
  )
  result
}


ess_per_class <- function(samples, n_classes) {
  cols <- paste0("class", seq_len(n_classes), "_rate_log_sd")
  vapply(cols, function(col) {
    x <- samples[, col]
    x <- x[is.finite(x)]
    if (length(x) < 10L) return(NA_real_)
    as.numeric(coda::effectiveSize(coda::as.mcmc(x)))
  }, numeric(1L))
}


# ----------------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------------

cat("Building fixture: 5 classes, sizes (5, 50, 50, 50, 50)\n")
fixture <- build_fixture()
n_classes <- length(unique(fixture$partition))
n_char    <- length(fixture$partition)
cat(sprintf("  nChar = %d, nTip = %d, partition table:\n",
            n_char, length(fixture$tree$tip.label)))
print(table(fixture$partition))

cat("\n--- Running with NEW prior: hyperprior_pooled (default) ---\n")
t0 <- Sys.time()
res_hyper <- run_one(fixture, "hyperprior_pooled")
cat(sprintf("  elapsed: %.1f s, nSamples: %d\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            res_hyper$nSamples))
ess_hyper <- ess_per_class(res_hyper$samples, n_classes)

cat("\n--- Running with LEGACY prior: gamma_independent ---\n")
t0 <- Sys.time()
res_gamma <- run_one(fixture, "gamma_independent")
cat(sprintf("  elapsed: %.1f s, nSamples: %d\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            res_gamma$nSamples))
ess_gamma <- ess_per_class(res_gamma$samples, n_classes)


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------

tbl <- data.frame(
  class     = seq_len(n_classes),
  size      = as.integer(table(fixture$partition)),
  ess_gamma_old = round(ess_gamma, 1L),
  ess_hyper_new = round(ess_hyper, 1L),
  ratio_new_over_old = round(ess_hyper / ess_gamma, 2L)
)

cat("\n================== ESS per class (σ_c) ==================\n")
print(tbl, row.names = FALSE)
cat("=========================================================\n\n")

# Headline summary (for NEWS.md / cherry-pick summary)
small_class <- which.min(tbl$size)
cat(sprintf(
  "Small-class (size %d) ESS: %s (gamma_independent) -> %s (hyperprior_pooled)\n",
  tbl$size[small_class],
  format(tbl$ess_gamma_old[small_class]),
  format(tbl$ess_hyper_new[small_class])
))
cat(sprintf("Improvement factor: %.1fx\n",
            tbl$ess_hyper_new[small_class] / tbl$ess_gamma_old[small_class]))

invisible(list(
  hyper = res_hyper,
  gamma = res_gamma,
  table = tbl
))
