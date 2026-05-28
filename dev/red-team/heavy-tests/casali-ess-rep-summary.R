#!/usr/bin/env Rscript
# Quick per-rep diagnostic on Burns2011/T4 and CarranoSampson2008/T2a
# across batch-1 (rep1, paired-seed default) and batch-2 (rep2, rep3).
# Prints σ per class, σ-min, σ-mean, τ ESS, TL ESS per (cell, prior, rep).

dir <- "dev/red-team/heavy-tests/casali-ess-results"

show_cell <- function(cell_tag) {
  cat("\n=========================================================\n")
  cat(cell_tag, "\n")
  cat("=========================================================\n")
  files <- list.files(dir, pattern = paste0("^", cell_tag, "__"), full.names = TRUE)
  for (f in sort(files)) {
    x <- readRDS(f)
    parts <- strsplit(sub("[.]rds$", "", basename(f)), "__")[[1]]
    prior <- parts[3]
    tag   <- if (length(parts) >= 4) parts[4] else "rep1"
    cat(sprintf("\n  %-20s %-6s seed=%d\n", prior, tag, x$config$seed))
    cat(sprintf("    sig per class: %s\n",
                paste(sprintf("%5.0f", x$ess_sigma), collapse = " ")))
    cat(sprintf("    sig-min: %5.0f  sig-mean: %5.0f  tau: %s  TL: %5.0f\n",
                min(x$ess_sigma, na.rm = TRUE),
                mean(x$ess_sigma, na.rm = TRUE),
                if (is.na(x$ess_hyper_tau)) "    -" else sprintf("%5.0f", x$ess_hyper_tau),
                x$ess_tree_len))
  }
}

show_cell("Burns2011__T4")
show_cell("CarranoSampson2008__T2a")
