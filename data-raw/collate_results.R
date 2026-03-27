#!/usr/bin/env Rscript
# collate_results.R
# Aggregates per-dataset RDS files from the Hamilton MkPrime validation array job
# into two analysis-ready CSVs.
#
# Outputs written to <out_dir>:
#   summary.csv      – one row per dataset (tree × rep):
#                        CID statistics, sample counts, stop reasons
#   u_estimates.csv  – one row per variable character per dataset:
#                        kObs and posterior mean u
#   missing.txt      – list of expected files not found (if any)
#
# Usage (run from mkp/ root):
#   Rscript data-raw/collate_results.R [results_dir] [out_dir] [data_root]
#
#   results_dir  directory containing result_t??_r??.rds files
#                (default: data-raw/hamilton-results)
#   out_dir      where to write CSVs
#                (default: same as results_dir)
#   data_root    [optional] local copy of tree-inference data root;
#                when supplied, n_tips and total_tree_length are read
#                from tree_NN/tree.nwk and added to summary.csv

suppressPackageStartupMessages(library(ape))

# ---- Arguments ---------------------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1L) args[[1L]] else "data-raw/hamilton-results"
out_dir     <- if (length(args) >= 2L) args[[2L]] else results_dir
data_root   <- if (length(args) >= 3L) args[[3L]] else NULL

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== MkPrime validation: collate results ===\n")
cat(sprintf("  results_dir : %s\n", results_dir))
cat(sprintf("  out_dir     : %s\n", out_dir))
if (!is.null(data_root))
  cat(sprintf("  data_root   : %s (tree properties enabled)\n", data_root))

# ---- Discover RDS files ------------------------------------------------------
rds_files <- sort(list.files(results_dir,
                              pattern    = "^result_t[0-9]+_r[0-9]+\\.rds$",
                              full.names = TRUE))
cat(sprintf("\nFound %d result file(s) (expected 260)\n", length(rds_files)))

# Identify any missing datasets from the full 26 × 10 grid
expected_tags <- sprintf("result_t%02d_r%02d.rds",
                         rep(seq_len(26L), each = 10L),
                         rep(seq_len(10L), times = 26L))
found_tags    <- basename(rds_files)
missing_tags  <- setdiff(expected_tags, found_tags)
if (length(missing_tags) > 0L) {
  cat(sprintf("  Missing: %d dataset(s) — written to missing.txt\n",
              length(missing_tags)))
  writeLines(missing_tags, file.path(out_dir, "missing.txt"))
} else {
  cat("  All 260 files present.\n")
}

# ---- Load with error handling ------------------------------------------------
load_result <- function(f) {
  tryCatch(
    readRDS(f),
    error = function(e) {
      message(sprintf("  FAILED %s: %s", basename(f), conditionMessage(e)))
      NULL
    }
  )
}

results <- lapply(rds_files, load_result)
ok      <- !vapply(results, is.null, logical(1L))
cat(sprintf("  Loaded: %d  |  Failed to parse: %d\n", sum(ok), sum(!ok)))
results <- results[ok]

if (length(results) == 0L)
  stop("No results loaded — nothing to collate.")

# ---- Optional: tree properties from true trees --------------------------------
# Returns a list(n_tips, total_tree_length) for a given tree index.
get_tree_props <- function(root, tree_idx) {
  f <- file.path(root, sprintf("tree_%02d/tree.nwk", tree_idx))
  if (!file.exists(f))
    return(list(n_tips = NA_integer_, total_tree_length = NA_real_))
  tr <- ape::read.tree(f)
  list(
    n_tips           = ape::Ntip(tr),
    total_tree_length = sum(tr$edge.length)
  )
}

# ---- Build per-dataset summary -----------------------------------------------
summary_rows <- lapply(results, function(r) {

  tp <- if (!is.null(data_root)) {
    get_tree_props(data_root, r$tree_idx)
  } else {
    list(n_tips = NA_integer_, total_tree_length = NA_real_)
  }

  # mean_delta_cid > 0: Mk' further from true tree than Mk (Mk' worse)
  # mean_delta_cid < 0: Mk' closer to true tree than Mk (Mk' better)
  data.frame(
    tree_idx          = r$tree_idx,
    rep_idx           = r$rep_idx,
    n_char            = r$n_char,
    n_tips            = as.integer(tp$n_tips),
    total_tree_length = tp$total_tree_length,
    n_samples_mk      = r$n_samples_mk,
    n_samples_mkp     = r$n_samples_mkp,
    stop_reason_mk    = r$stop_reason_mk,
    stop_reason_mkp   = r$stop_reason_mkp,
    mean_cid_mk       = mean(r$cid_mk,  na.rm = TRUE),
    mean_cid_mkp      = mean(r$cid_mkp, na.rm = TRUE),
    median_cid_mk     = median(r$cid_mk,  na.rm = TRUE),
    median_cid_mkp    = median(r$cid_mkp, na.rm = TRUE),
    mean_delta_cid    = mean(r$cid_mkp, na.rm = TRUE) -
                          mean(r$cid_mk, na.rm = TRUE),
    mean_kObs         = mean(r$kObs,         na.rm = TRUE),
    mean_u_post       = mean(r$u_post_means, na.rm = TRUE),
    stringsAsFactors  = FALSE
  )
})

summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[order(summary_df$tree_idx, summary_df$rep_idx), ]
rownames(summary_df) <- NULL

# ---- Build per-character u-estimates -----------------------------------------
u_rows <- lapply(results, function(r) {
  n <- length(r$kObs)
  data.frame(
    tree_idx    = r$tree_idx,
    rep_idx     = r$rep_idx,
    char_idx    = seq_len(n),
    kObs        = r$kObs,
    u_post_mean = unname(r$u_post_means),
    stringsAsFactors = FALSE
  )
})

u_df <- do.call(rbind, u_rows)
u_df <- u_df[order(u_df$tree_idx, u_df$rep_idx, u_df$char_idx), ]
rownames(u_df) <- NULL

# ---- Quality checks ----------------------------------------------------------
cat("\n--- Quality checks ---\n")

low_mk  <- sum(summary_df$n_samples_mk  < 500L, na.rm = TRUE)
low_mkp <- sum(summary_df$n_samples_mkp < 500L, na.rm = TRUE)
if (low_mk > 0L || low_mkp > 0L) {
  cat(sprintf("WARNING: datasets with < 500 posterior trees  Mk: %d  Mk': %d\n",
              low_mk, low_mkp))
  flag <- (summary_df$n_samples_mk  < 500L |
           summary_df$n_samples_mkp < 500L)
  flag[is.na(flag)] <- FALSE
  print(summary_df[flag, c("tree_idx", "rep_idx",
                            "n_samples_mk", "n_samples_mkp",
                            "stop_reason_mk", "stop_reason_mkp")])
} else {
  cat("  All datasets have >= 500 posterior trees for both arms.\n")
}

# Flag any datasets where NA crept into CID
na_cid <- sum(is.na(summary_df$mean_cid_mk) | is.na(summary_df$mean_cid_mkp))
if (na_cid > 0L)
  cat(sprintf("WARNING: %d dataset(s) have NA CID values.\n", na_cid))

# ---- Write outputs -----------------------------------------------------------
write.csv(summary_df, file.path(out_dir, "summary.csv"),     row.names = FALSE)
write.csv(u_df,       file.path(out_dir, "u_estimates.csv"), row.names = FALSE)
cat(sprintf("\nWrote summary.csv      (%d rows × %d cols)\n",
            nrow(summary_df), ncol(summary_df)))
cat(sprintf("Wrote u_estimates.csv  (%d rows × %d cols)\n",
            nrow(u_df),       ncol(u_df)))

# ---- Quick summary stats -----------------------------------------------------
cat("\n--- CID summary (across all datasets) ---\n")
cat(sprintf("  Mk   mean CID: %.4f  (median %.4f)\n",
            mean(summary_df$mean_cid_mk,   na.rm = TRUE),
            median(summary_df$mean_cid_mk, na.rm = TRUE)))
cat(sprintf("  Mk'  mean CID: %.4f  (median %.4f)\n",
            mean(summary_df$mean_cid_mkp,   na.rm = TRUE),
            median(summary_df$mean_cid_mkp, na.rm = TRUE)))
cat(sprintf("  Mean delta (Mk' - Mk): %.4f  (positive = Mk worse at recovering true tree)\n",
            mean(summary_df$mean_delta_cid, na.rm = TRUE)))
n_valid <- sum(!is.na(summary_df$mean_delta_cid))
n_mkp_better <- sum(summary_df$mean_delta_cid < 0, na.rm = TRUE)
cat(sprintf("  Mk' better in %d / %d datasets (%.1f%%)\n",
            n_mkp_better, n_valid, 100 * n_mkp_better / n_valid))

cat("\n--- Sample counts ---\n")
cat(sprintf("  Mk   samples: min=%d  median=%d  max=%d\n",
            min(summary_df$n_samples_mk,    na.rm = TRUE),
            median(summary_df$n_samples_mk, na.rm = TRUE),
            max(summary_df$n_samples_mk,    na.rm = TRUE)))
cat(sprintf("  Mk'  samples: min=%d  median=%d  max=%d\n",
            min(summary_df$n_samples_mkp,    na.rm = TRUE),
            median(summary_df$n_samples_mkp, na.rm = TRUE),
            max(summary_df$n_samples_mkp,    na.rm = TRUE)))

cat("\n--- Stop reasons (Mk) ---\n")
print(table(summary_df$stop_reason_mk))
cat("\n--- Stop reasons (Mk') ---\n")
print(table(summary_df$stop_reason_mkp))

cat("\n--- Posterior u (Mk' arm) ---\n")
cat(sprintf("  Mean of per-dataset mean u: %.3f\n",
            mean(summary_df$mean_u_post, na.rm = TRUE)))
cat(sprintf("  Range: [%.3f, %.3f]\n",
            min(summary_df$mean_u_post,  na.rm = TRUE),
            max(summary_df$mean_u_post,  na.rm = TRUE)))

cat("\nDone.\n")
