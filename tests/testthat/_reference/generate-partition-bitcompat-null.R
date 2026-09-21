#!/usr/bin/env Rscript
# Regenerate the §7a partition-bitcompat-null reference RDS.
#
# Run from the package root, against the SAME installed build the test suite
# runs against (build-agent.sh / test-agent.sh; see ../AGENTS.md):
#
#   MKP_REF_LIB=.builds/MkPrime-<id> MKP_REF_PKG=MkPrime.<id> \
#     Rscript tests/testthat/_reference/generate-partition-bitcompat-null.R
#
# The installed build is not a convenience: devtools::load_all() compiles
# without -O2, so a reference generated that way can disagree in the last
# bits with the -O2 build the test then compares it against. It also
# collides across concurrent sessions sharing the checkout. The load_all
# fallback below exists only for a single interactive session.
#
# Only re-run this script after an INTENTIONAL change to the legacy code
# path (or to the helper's MCMC config). Every value in result$samples must be
# reproducible bit-for-bit under a re-run with the same seed.

refLib <- Sys.getenv("MKP_REF_LIB", "")
refPkg <- Sys.getenv("MKP_REF_PKG", "")

if (nzchar(refLib) && nzchar(refPkg)) {
  .libPaths(c(normalizePath(refLib), .libPaths()))
  suppressPackageStartupMessages({
    library(refPkg, character.only = TRUE)
    library("TreeTools")
  })
  # Namespace shim: register the renamed build under the original name so
  # MkPrime::: calls resolve (same mechanism as test-agent.sh).
  if (is.null(.Internal(getRegisteredNamespace("MkPrime")))) {
    .Internal(registerNamespace("MkPrime", asNamespace(refPkg)))
  }
  env <- new.env(parent = asNamespace(refPkg))
} else {
  suppressPackageStartupMessages({
    devtools::load_all(".", quiet = TRUE)
    library("TreeTools")
  })
  env <- globalenv()
}

sys.source("tests/testthat/helper-partition-ref.R", envir = env)

r <- env$.RunPartitionBitcompatReference()

out <- list(
  samples       = r$samples,
  samples_cols  = colnames(r$samples),
  samples_dim   = dim(r$samples),
  nSamples      = r$nSamples,
  stop_reason   = r$stop_reason,
  actual_iter   = r$actual_iter,
  treeThin      = r$treeThin
)

ref_path <- env$.PartitionBitcompatReferencePath()
dir.create(dirname(ref_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(out, ref_path)

cat("Wrote:", ref_path, "\n")
cat("  built from:", if (nzchar(refPkg)) refPkg else "devtools::load_all", "\n")
cat("  samples dim:", paste(dim(r$samples), collapse = " x "), "\n")
cat("  final logP:", r$samples[nrow(r$samples), "log_posterior"], "\n")
cat("  size bytes:", file.size(ref_path), "\n")
