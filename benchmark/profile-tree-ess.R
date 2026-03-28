# End-to-end profiling of the tree ESS pipeline
#
# Breaks .TreeESS() into its component stages and times each independently
# across a grid of (n_trees, n_tips) combinations.
#
# Stages:
#   RF        — TreeDist::RobinsonFoulds() pairwise distance matrix
#   asmat     — as.matrix() conversion from dist
#   sq        — element-wise squaring (dmat * dmat) for Fréchet ESS
#   frechet   — C++ Fréchet correlation ESS inner loop
#   median    — C++ Geyer initial-monotone-sequence median pseudo-ESS
#
# Usage (from mkp/ root, using a release-installed build):
#   Rscript benchmark/profile-tree-ess.R
#
# Or from the tree-ess worktree:
#   Rscript -e "pkgbuild::compile_dll(debug=FALSE); devtools::load_all(quiet=TRUE); source('benchmark/profile-tree-ess.R')"

suppressPackageStartupMessages({
  library(ape)
  library(TreeDist)
  library(bench)
})

# If mkp is installed, load it; otherwise assume load_all() was done
if (requireNamespace("mkp", quietly = TRUE)) {
  library(mkp)
} else if (!exists("frechet_correlation_ess_cpp", mode = "function")) {
  # Loaded via devtools::load_all() in the wrapper call
  if (!exists(".FrechetCorrelationESS", mode = "function")) {
    stop("mkp package not loaded. Run via devtools::load_all() or install first.")
  }
}

# --- Tree chain generator ---
# Simulate an autocorrelated NNI chain (like an MCMC tree sample)
make_nni_chain <- function(n_trees, n_tips, seed = 3847L) {
  set.seed(seed)
  tr <- rtree(n_tips, rooted = FALSE)
  trees <- vector("list", n_trees)
  trees[[1]] <- tr
  for (i in 2:n_trees) {
    # Apply a random NNI move; if it fails, keep the same tree
    cand <- tryCatch(
      rNNI(trees[[i - 1]], moves = 1),
      error = function(e) trees[[i - 1]]
    )
    trees[[i]] <- cand
  }
  do.call(c, trees)
}

# --- Profiling grid ---
grid <- expand.grid(
  n_trees = c(100L, 250L, 500L, 1000L),
  n_tips  = c(20L, 50L, 100L)
)

cat(sprintf("Profiling %d configurations...\n\n", nrow(grid)))

results <- vector("list", nrow(grid))

for (idx in seq_len(nrow(grid))) {
  n <- grid$n_trees[idx]
  t <- grid$n_tips[idx]
  cat(sprintf("[%d/%d] n_trees=%d, n_tips=%d ... ", idx, nrow(grid), n, t))

  trees <- make_nni_chain(n, t, seed = 3847L + idx)

  # Warmup (get JIT/cache effects out of the way)
  invisible(RobinsonFoulds(trees[1:min(10, n)]))

  # Stage 1: RF distance computation
  t_rf <- bench::mark(
    rf = RobinsonFoulds(trees),
    min_iterations = 3, max_iterations = 10, check = FALSE
  )

  # Compute once for use in subsequent stages
  d_obj <- RobinsonFoulds(trees)
  
  # Stage 1b: dist → matrix conversion
  t_asmat <- bench::mark(
    asmat = as.matrix(d_obj),
    min_iterations = 5, max_iterations = 20, check = FALSE
  )

  dmat <- as.matrix(d_obj)

  # Stage 1c: element-wise squaring
  t_sq <- bench::mark(
    sq = dmat * dmat,
    min_iterations = 10, max_iterations = 50, check = FALSE
  )

  dmat_sq <- dmat * dmat

  # Stage 2: Frechet correlation ESS (C++ inner loop)
  t_frechet <- bench::mark(
    frechet = frechet_correlation_ess_cpp(dmat_sq, 5L),
    min_iterations = 10, max_iterations = 100, check = FALSE
  )

  # Stage 3: Median pseudo-ESS (C++ Geyer with row subsampling)
  t_median <- bench::mark(
    median_ess = median_pseudo_ess_cpp(dmat, 5L, 200L),
    min_iterations = 10, max_iterations = 100, check = FALSE
  )

  # Collect results
  ms <- function(b) as.numeric(b$median) * 1000  # median in milliseconds
  row <- data.frame(
    n_trees    = n,
    n_tips     = t,
    rf_ms      = ms(t_rf),
    asmat_ms   = ms(t_asmat),
    sq_ms      = ms(t_sq),
    frechet_ms = ms(t_frechet),
    median_ms  = ms(t_median)
  )
  row$total_ms <- with(row, rf_ms + asmat_ms + sq_ms + frechet_ms +
                          ifelse(is.na(median_ms), 0, median_ms))
  results[[idx]] <- row

  cat(sprintf("total=%.0fms (RF=%.0f, asmat=%.0f, sq=%.1f, frechet=%.1f, median=%.0f)\n",
              row$total_ms, row$rf_ms, row$asmat_ms, row$sq_ms,
              row$frechet_ms, row$median_ms))
}

res <- do.call(rbind, results)

# --- Summary table ---
cat("\n=== Stage Breakdown (median ms) ===\n\n")
print(res[, c("n_trees", "n_tips", "rf_ms", "asmat_ms", "sq_ms",
              "frechet_ms", "median_ms", "total_ms")], digits = 3, row.names = FALSE)

# --- Percentage breakdown ---
cat("\n=== Percentage of Total ===\n\n")
pct <- res
for (col in c("rf_ms", "asmat_ms", "sq_ms", "frechet_ms", "median_ms")) {
  pct[[sub("_ms", "_pct", col)]] <- round(100 * res[[col]] / res$total_ms, 1)
}
print(pct[, c("n_trees", "n_tips", "rf_pct", "asmat_pct", "sq_pct",
              "frechet_pct", "median_pct")], row.names = FALSE)
