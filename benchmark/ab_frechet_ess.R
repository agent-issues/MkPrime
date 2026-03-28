# A/B benchmark: O(n^2+nL) (Ref) vs O(n^2+L) (Dev) Frechet ESS
#
# Runs both implementations in the same process using renamed packages.
# The "canary" is a function from a shared .cpp file that should NOT
# differ between Ref and Dev (pruning_jc_acrv).

setTimeLimit(elapsed = 55, transient = FALSE)

ref_lib <- normalizePath("../.builds/MkPrime-Ref")
dev_lib <- normalizePath("../.builds/MkPrime-Dev")
.libPaths(c(dev_lib, ref_lib, .libPaths()))

library(MkPrime.Ref, lib.loc = ref_lib)
library(MkPrime.Dev, lib.loc = dev_lib)
library(bench)

# Access the C++ functions directly
ref_fn <- MkPrime.Ref:::frechet_correlation_ess_cpp
dev_fn <- MkPrime.Dev:::frechet_correlation_ess_cpp

set.seed(7491)

# --- Generate test matrices of various sizes ---
# Simulate autocorrelated distance matrices (realistic MCMC scenario)
make_autocorrelated_dmat_sq <- function(n, rho = 0.95) {
  # AR(1) positions in 1D; pairwise Euclidean distances -> squared
  x <- numeric(n)
  x[1] <- rnorm(1)
  for (i in 2:n) x[i] <- rho * x[i - 1] + sqrt(1 - rho^2) * rnorm(1)
  d <- as.matrix(dist(x))
  d * d
}

make_random_dmat_sq <- function(n) {
  x <- rnorm(n)
  d <- as.matrix(dist(x))
  d * d
}

sizes <- c(50L, 100L, 200L, 500L, 1000L)
min_ns <- 5L

cat("\n=== Autocorrelated chains (rho=0.95) ===\n\n")
for (sz in sizes) {
  dmat_sq <- make_autocorrelated_dmat_sq(sz)
  
  # Verify agreement
  r <- ref_fn(dmat_sq, min_ns)
  v <- dev_fn(dmat_sq, min_ns)
  if (abs(r - v) / max(abs(r), 1) > 1e-6) {
    cat(sprintf("  n=%d: MISMATCH ref=%.6f dev=%.6f\n", sz, r, v))
  }
  
  # Warmup
  for (i in 1:3) { ref_fn(dmat_sq, min_ns); dev_fn(dmat_sq, min_ns) }
  
  b <- bench::mark(
    ref = ref_fn(dmat_sq, min_ns),
    dev = dev_fn(dmat_sq, min_ns),
    min_iterations = 20,
    check = FALSE
  )
  cat(sprintf("  n=%4d: ref=%8s  dev=%8s  speedup=%.2fx\n",
              sz,
              format(b$median[1]),
              format(b$median[2]),
              as.numeric(b$median[1]) / as.numeric(b$median[2])))
}

cat("\n=== Random (low autocorrelation) ===\n\n")
for (sz in sizes) {
  dmat_sq <- make_random_dmat_sq(sz)
  
  r <- ref_fn(dmat_sq, min_ns)
  v <- dev_fn(dmat_sq, min_ns)
  if (abs(r - v) / max(abs(r), 1) > 1e-6) {
    cat(sprintf("  n=%d: MISMATCH ref=%.6f dev=%.6f\n", sz, r, v))
  }
  
  for (i in 1:3) { ref_fn(dmat_sq, min_ns); dev_fn(dmat_sq, min_ns) }
  
  b <- bench::mark(
    ref = ref_fn(dmat_sq, min_ns),
    dev = dev_fn(dmat_sq, min_ns),
    min_iterations = 20,
    check = FALSE
  )
  cat(sprintf("  n=%4d: ref=%8s  dev=%8s  speedup=%.2fx\n",
              sz,
              format(b$median[1]),
              format(b$median[2]),
              as.numeric(b$median[1]) / as.numeric(b$median[2])))
}

# --- Canary: pruning_jc_acrv (shared code, no changes expected) ---
cat("\n=== Canary: pruning_jc_acrv ===\n")
parent <- c(5L, 5L, 6L, 6L)
child <- c(1L, 2L, 3L, 4L)
el <- c(0.1, 0.2, 0.15, 0.25)
tip <- matrix(c(1L, 2L, 1L, 2L), nrow = 4, ncol = 1)
rf <- c(0.5, 0.5)
rm_vec <- 1.0

ref_canary <- MkPrime.Ref:::pruning_jc_acrv
dev_canary <- MkPrime.Dev:::pruning_jc_acrv

for (i in 1:3) {
  ref_canary(parent, child, el, tip, 2L, rf, rm_vec)
  dev_canary(parent, child, el, tip, 2L, rf, rm_vec)
}

bc <- bench::mark(
  ref = ref_canary(parent, child, el, tip, 2L, rf, rm_vec),
  dev = dev_canary(parent, child, el, tip, 2L, rf, rm_vec),
  min_iterations = 100,
  check = FALSE
)
cat(sprintf("  canary: ref=%8s  dev=%8s  ratio=%.2f\n",
            format(bc$median[1]),
            format(bc$median[2]),
            as.numeric(bc$median[1]) / as.numeric(bc$median[2])))
cat("  (Canary should be ~1.0; >1.03 or <0.97 suggests environmental noise)\n")
