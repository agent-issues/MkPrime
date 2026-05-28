#!/usr/bin/env Rscript
#
# casali-ess-aggregate.R — read all per-(cell, prior) RDS outputs
# written by casali-ess-hyperprior-vs-gamma.R and emit a side-by-side
# comparison table.
#
# Usage:
#   Rscript dev/red-team/heavy-tests/casali-ess-aggregate.R \
#     --results-dir /nobackup/pjjg18/mkp-hyperprior-bench/results \
#     [--out-md     dev/red-team/heavy-tests/casali-ess-results.md] \
#     [--out-csv    dev/red-team/heavy-tests/casali-ess-results.csv]
#
# Output: a wide table per (matrix, treatment) with min, mean, and
# per-class σ_c ESS under each prior, and the ratio.

suppressPackageStartupMessages({
  library(optparse)
})

`%||%` <- function(a, b) if (is.null(a)) b else a


.parseArgs <- function(argv = commandArgs(trailingOnly = TRUE)) {
  optList <- list(
    make_option("--results-dir", type = "character", default = NULL),
    make_option("--out-md",      type = "character", default = NULL),
    make_option("--out-csv",     type = "character", default = NULL)
  )
  opt <- parse_args(OptionParser(option_list = optList), args = argv)
  if (is.null(opt[["results-dir"]])) {
    stop("--results-dir is required", call. = FALSE)
  }
  opt
}


.LoadAll <- function(resultsDir) {
  files <- list.files(resultsDir, pattern = "\\.rds$", full.names = TRUE)
  if (!length(files)) {
    stop("No .rds files in ", resultsDir, call. = FALSE)
  }
  lapply(files, function(f) {
    x <- readRDS(f)
    x$.file <- basename(f)
    x
  })
}


.Build <- function(records) {
  rows <- lapply(records, function(r) {
    sig <- r$ess_sigma
    data.frame(
      matrix     = r$cell$matrix,
      treatment  = r$cell$treatment,
      nClasses   = r$cell$nClasses,
      prior      = r$prior,
      ess_sig_min = if (length(sig)) min(sig, na.rm = TRUE) else NA_real_,
      ess_sig_mean = if (length(sig)) mean(sig, na.rm = TRUE) else NA_real_,
      ess_sig_max  = if (length(sig)) max(sig, na.rm = TRUE) else NA_real_,
      ess_tree_len = r$ess_tree_len %||% NA_real_,
      ess_log_post = r$ess_log_post %||% NA_real_,
      ess_hyper_tau = r$ess_hyper_tau %||% NA_real_,
      nSamples   = r$nSamples,
      wallSecs   = r$wallSecs,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


.Pivot <- function(long) {
  hyper <- long[long$prior == "hyperprior_pooled", ]
  gamma <- long[long$prior == "gamma_independent", ]
  key <- function(d) paste(d$matrix, d$treatment, sep = "__")
  merged <- merge(
    gamma[, c("matrix", "treatment", "nClasses",
              "ess_sig_min", "ess_sig_mean", "ess_tree_len", "wallSecs")],
    hyper[, c("matrix", "treatment",
              "ess_sig_min", "ess_sig_mean", "ess_tree_len",
              "ess_hyper_tau", "wallSecs")],
    by = c("matrix", "treatment"),
    suffixes = c("_gamma", "_hyper"),
    all = TRUE
  )
  merged$ratio_min  <- merged$ess_sig_min_hyper  / merged$ess_sig_min_gamma
  merged$ratio_mean <- merged$ess_sig_mean_hyper / merged$ess_sig_mean_gamma
  merged[order(merged$ess_sig_min_gamma), ]
}


.Format <- function(wide) {
  out <- data.frame(
    cell = sprintf("%s/%s", wide$matrix, wide$treatment),
    K    = wide$nClasses,
    sig_min_gamma  = round(wide$ess_sig_min_gamma,  0),
    sig_min_hyper  = round(wide$ess_sig_min_hyper,  0),
    ratio_min      = round(wide$ratio_min,          2),
    sig_mean_gamma = round(wide$ess_sig_mean_gamma, 0),
    sig_mean_hyper = round(wide$ess_sig_mean_hyper, 0),
    ratio_mean     = round(wide$ratio_mean,         2),
    tl_gamma       = round(wide$ess_tree_len_gamma, 0),
    tl_hyper       = round(wide$ess_tree_len_hyper, 0),
    tau_ess        = round(wide$ess_hyper_tau,      0),
    wall_gamma_h   = round(wide$wallSecs_gamma / 3600, 2),
    wall_hyper_h   = round(wide$wallSecs_hyper / 3600, 2),
    stringsAsFactors = FALSE
  )
  out
}


.WriteMd <- function(formatted, path) {
  hd <- c("cell", "K",
          "σ_min ESS γ", "σ_min ESS hyper", "×",
          "σ_mean ESS γ", "σ_mean ESS hyper", "×",
          "TL γ", "TL hyper", "τ ESS",
          "wall(γ) h", "wall(hyper) h")
  body <- apply(formatted, 1L, function(row) paste(row, collapse = " | "))
  txt <- c(
    paste0("| ", paste(hd, collapse = " | "), " |"),
    paste0("| ", paste(rep("---", length(hd)), collapse = " | "), " |"),
    paste0("| ", body, " |")
  )
  writeLines(txt, path)
  message("Wrote markdown table: ", path)
}


.Main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  opt <- .parseArgs(argv)
  records <- .LoadAll(opt[["results-dir"]])
  message("Loaded ", length(records), " result files from ",
          opt[["results-dir"]])
  long <- .Build(records)
  wide <- .Pivot(long)
  fmt  <- .Format(wide)

  cat("\n================ Casali ESS: hyperprior_pooled vs gamma_independent ================\n")
  print(fmt, row.names = FALSE)
  cat("====================================================================================\n\n")

  if (!is.null(opt[["out-csv"]])) {
    write.csv(fmt, opt[["out-csv"]], row.names = FALSE)
    message("Wrote CSV: ", opt[["out-csv"]])
  }
  if (!is.null(opt[["out-md"]])) {
    .WriteMd(fmt, opt[["out-md"]])
  }
  invisible(0L)
}

if (sys.nframe() == 0L) {
  status <- tryCatch(.Main(), error = function(e) {
    message("ERROR: ", conditionMessage(e))
    1L
  })
  quit(status = if (is.numeric(status)) as.integer(status) else 0L)
}
