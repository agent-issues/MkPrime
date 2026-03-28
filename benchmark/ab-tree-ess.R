# A/B benchmark: tree ESS pipeline (Ref vs Dev)
#
# Compares the full .TreeESS() and its sub-components between the
# reference build (coda-based MedianPseudoESS) and the dev build
# (C++ Geyer-based MedianPseudoESS with row subsampling).
#
# Usage:
#   Rscript benchmark/ab-tree-ess.R

ref_lib <- normalizePath(".builds/MkPrimeRef")
dev_lib <- normalizePath(".builds/MkPrimeDev")
.libPaths(c(dev_lib, ref_lib, .libPaths()))

library(MkPrimeRef, lib.loc = ref_lib)
library(MkPrimeDev, lib.loc = dev_lib)
library(ape)
library(TreeDist)
library(bench)

set.seed(2781)

# --- Generate realistic NNI chain ---
make_nni_chain <- function(n_trees, n_tips, seed = 4517L) {
  set.seed(seed)
  tr <- rtree(n_tips, rooted = FALSE)
  trees <- vector("list", n_trees)
  trees[[1]] <- tr
  for (i in 2:n_trees) {
    cand <- tryCatch(rNNI(trees[[i - 1]], moves = 1),
                     error = function(e) trees[[i - 1]])
    trees[[i]] <- cand
  }
  do.call(c, trees)
}

# --- Sub-component: median_pseudo_ess only ---
cat("=== MedianPseudoESS sub-component benchmark ===\n\n")

ref_median <- MkPrimeRef:::.MedianPseudoESS
dev_median <- MkPrimeDev:::.MedianPseudoESS

sizes <- c(100L, 250L, 500L, 1000L)
for (sz in sizes) {
  trees <- make_nni_chain(sz, 50L, seed = 4517L + sz)
  dmat <- as.matrix(TreeDist::RobinsonFoulds(trees))

  # Warmup
  for (i in 1:2) { ref_median(dmat); dev_median(dmat) }

  b <- bench::mark(
    ref = ref_median(dmat),
    dev = dev_median(dmat),
    min_iterations = 3,
    check = FALSE
  )
  cat(sprintf("  n=%4d: ref=%10s  dev=%10s  speedup=%.0fx\n",
              sz, format(b$median[1]), format(b$median[2]),
              as.numeric(b$median[1]) / as.numeric(b$median[2])))
}

# --- Full .TreeESS() pipeline ---
cat("\n=== Full .TreeESS() pipeline benchmark ===\n\n")

ref_tree_ess <- MkPrimeRef:::.TreeESS
dev_tree_ess <- MkPrimeDev:::.TreeESS

for (sz in sizes) {
  trees <- make_nni_chain(sz, 50L, seed = 4517L + sz)

  # Warmup
  for (i in 1:2) { ref_tree_ess(trees); dev_tree_ess(trees) }

  b <- bench::mark(
    ref = ref_tree_ess(trees),
    dev = dev_tree_ess(trees),
    min_iterations = 3,
    check = FALSE
  )
  cat(sprintf("  n=%4d: ref=%10s  dev=%10s  speedup=%.1fx\n",
              sz, format(b$median[1]), format(b$median[2]),
              as.numeric(b$median[1]) / as.numeric(b$median[2])))
}

# --- Canary: Frechet ESS (should NOT differ) ---
cat("\n=== Canary: FrechetCorrelationESS (unchanged code) ===\n")
ref_frechet <- MkPrimeRef:::frechet_correlation_ess_cpp
dev_frechet <- MkPrimeDev:::frechet_correlation_ess_cpp

trees <- make_nni_chain(500L, 50L, seed = 9023L)
dmat <- as.matrix(TreeDist::RobinsonFoulds(trees))
dmat_sq <- dmat * dmat

for (i in 1:3) { ref_frechet(dmat_sq, 5L); dev_frechet(dmat_sq, 5L) }

bc <- bench::mark(
  ref = ref_frechet(dmat_sq, 5L),
  dev = dev_frechet(dmat_sq, 5L),
  min_iterations = 50,
  check = FALSE
)
cat(sprintf("  canary: ref=%8s  dev=%8s  ratio=%.2f\n",
            format(bc$median[1]), format(bc$median[2]),
            as.numeric(bc$median[1]) / as.numeric(bc$median[2])))
cat("  (Canary should be ~1.0; >1.03 or <0.97 suggests environmental noise)\n")
