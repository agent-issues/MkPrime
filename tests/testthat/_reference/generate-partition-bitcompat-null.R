#!/usr/bin/env Rscript
# Regenerate the §7a partition-bitcompat-null reference RDS.
#
# Run from the package root:
#   Rscript tests/testthat/_reference/generate-partition-bitcompat-null.R
#
# Only re-run this script after an INTENTIONAL change to the legacy code
# path (or to the helper's MCMC config), with reviewer sign-off and a NEWS.md
# entry. Every value in result$samples must be reproducible bit-for-bit
# under a re-run with the same seed.

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(TreeTools)
})

source("tests/testthat/helper-partition-ref.R")

r <- .RunPartitionBitcompatReference()

out <- list(
  samples       = r$samples,
  samples_cols  = colnames(r$samples),
  samples_dim   = dim(r$samples),
  nSamples      = r$nSamples,
  stop_reason   = r$stop_reason,
  actual_iter   = r$actual_iter,
  treeThin      = r$treeThin
)

ref_path <- .PartitionBitcompatReferencePath()
dir.create(dirname(ref_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(out, ref_path)

cat("Wrote:", ref_path, "\n")
cat("  samples dim:", paste(dim(r$samples), collapse = " x "), "\n")
cat("  final logP:", r$samples[nrow(r$samples), "log_posterior"], "\n")
cat("  size bytes:", file.size(ref_path), "\n")
