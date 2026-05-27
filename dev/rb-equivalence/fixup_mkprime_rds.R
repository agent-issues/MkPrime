#!/usr/bin/env Rscript
# Repair an existing mkprime_<pid>_<model>.rds in-place by adding
# per_run_scalars (from the streamed *.log files) and per_run_trees (from the
# posterior$per_run trees field).
#
# Usage:
#   Rscript dev/rb-equivalence/fixup_mkprime_rds.R <pid> <model> \
#     [--out-dir=dev/rb-equivalence/out]

suppressPackageStartupMessages({
  library(MkPrime)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: fixup_mkprime_rds.R <pid> <model> [opts]")
pid <- args[[1]]; model <- args[[2]]

script_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- a[grep("^--file=", a)]
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})()

opt <- list(out_dir = file.path(script_dir, "out"))
for (a in args[-(1:2)]) {
  kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  if (length(kv) != 2L) stop("Bad option: ", a)
  opt[[gsub("-", "_", kv[[1]])]] <- kv[[2]]
}

rds_path <- file.path(opt$out_dir,
                     sprintf("mkprime_%s_%s.rds", pid, model))
stopifnot(file.exists(rds_path))
result <- readRDS(rds_path)

# Locate per-run log files (mkprime_<pid>_<model>_<r>.log)
base <- file.path(opt$out_dir, sprintf("mkprime_%s_%s", pid, model))
per_run_log_files <- sprintf("%s_%d.log", base, seq_len(2L))
per_run_scalars <- lapply(per_run_log_files, function(f) {
  if (!file.exists(f)) {
    cli::cli_warn("Missing log: {f}"); return(NULL)
  }
  as.data.frame(ReadMkLog(f), stringsAsFactors = FALSE)
})

per_run_trees <- lapply(result$posterior$per_run, function(pr) pr$trees)

result$per_run_scalars <- per_run_scalars
result$per_run_trees <- per_run_trees

saveRDS(result, rds_path)
cat(sprintf("[fixup] %s patched: nRows per scalar run = %s; nTrees per run = %s\n",
            rds_path,
            paste(vapply(per_run_scalars, nrow, integer(1)), collapse = ","),
            paste(vapply(per_run_trees, length, integer(1)), collapse = ",")))
