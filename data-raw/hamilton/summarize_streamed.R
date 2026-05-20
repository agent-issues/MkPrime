#!/usr/bin/env Rscript
# Post-hoc summarisation of streamed MCMC output from a TIMEOUT'd
# RunMkPrime job (the in-R `saveRDS(partial, ...)` never fired because
# HARNESS-001 graceful-exit didn't trigger before SIGKILL at 8h wall).
#
# Reads the streamed .log + _trees.nwk files for ONE (tree, rep, arm)
# task and writes a small RDS summary to
#   <out_dir>/<arm>_t{NN}_r{MM}.rds
#
# Usage:
#   Rscript summarize_streamed.R <tree_idx> <rep_idx> <arm> <results_root> \
#                                <data_root> <out_dir> [n_thin]
#
# Arguments:
#   tree_idx      integer 1..26
#   rep_idx       integer 1..10
#   arm           "mk" or "mkp_geo"
#   results_root  /nobackup/pjjg18/mkp-study/results
#   data_root     /nobackup/pjjg18/mkprime-files/tree-inference  (for true tree)
#   out_dir       /nobackup/pjjg18/mkp-study/summary
#   n_thin        optional integer, default 1000 (target # thinned trees)

.libPaths(c("/nobackup/pjjg18/mkp-study/lib", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(TreeDist)
  library(MkPrime)
  library(TreeTools)
})

args         <- commandArgs(trailingOnly = TRUE)
tree_idx     <- as.integer(args[1])
rep_idx      <- as.integer(args[2])
arm          <- match.arg(args[3], c("mk", "mkp_geo", "mkp_eg", "mkp", "mk_kp1", "mk_kp2", "mk_k9", "mk_k15", "mk_k24", "mk_k40", "mk_ktrue", "mkp_highk", "mkp_logs"))
results_root <- args[4]
data_root    <- args[5]
out_dir      <- args[6]
n_thin       <- if (length(args) >= 7L) as.integer(args[7]) else 1000L

tag       <- sprintf("t%02d_r%02d", tree_idx, rep_idx)
task_dir  <- file.path(results_root, tag)
log_file  <- file.path(task_dir,  sprintf("%s_run_1.log", arm))
tree_file <- file.path(task_dir,  sprintf("%s_trees.nwk", arm))
true_tree_file <- file.path(data_root,
                            sprintf("tree_%02d/tree.nwk", tree_idx))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("%s_%s.rds", arm, tag))

cat(sprintf("[%s] tag=%s arm=%s\n", format(Sys.time(), "%H:%M:%S"), tag, arm))
cat("  log : ", log_file,  "\n")
cat("  tree: ", tree_file, "\n")
cat("  true: ", true_tree_file, "\n")
cat("  out : ", out_file,  "\n")

stopifnot(file.exists(log_file), file.exists(tree_file),
          file.exists(true_tree_file))

# ---- Load character data in MCMC lex order to get kObs alignment -----------
# IMPORTANT: run_one.R uses plain sort() on chr*.nex filenames, which gives
# lexicographic order (chr1, chr10, chr11, ..., chr2, ...) for unpadded names.
# kPrime_<i> in the log refers to position i in this lex-sorted character list.
# Downstream code MUST use this kObs vector to pair kp_means[i] with kObs[i].
# Using natural sort (chr1, chr2, ...) misaligns kObs for tasks with >= 10
# characters — that was the source of the per-kObs mean flatness bug.
dataset_dir <- file.path(data_root,
                          sprintf("tree_%02d/rep_%02d", tree_idx, rep_idx))
stopifnot(dir.exists(dataset_dir))
# Lex sort: matches sort() in run_one.R exactly.
nex_files_lex <- sort(list.files(dataset_dir, pattern = "^chr[0-9]+\\.nex$",
                                  full.names = TRUE))
mat_list <- lapply(nex_files_lex, TreeTools::ReadCharacters)
combined_mat <- do.call(cbind, mat_list)
n_char_raw <- ncol(combined_mat)
tmp_nex <- tempfile(fileext = ".nex")
ape::write.nexus.data(
  setNames(lapply(seq_len(nrow(combined_mat)), function(i) combined_mat[i, ]),
           rownames(combined_mat)),
  file = tmp_nex, format = "standard"
)
pd <- TreeTools::ReadAsPhyDat(tmp_nex)
file.remove(tmp_nex)
mkd_local <- MkPrimeData(pd)
# kObs in lex-sorted character order (matches kPrime_i log column positions).
kObs_lex <- as.integer(mkd_local$kObs)
cat(sprintf("  char data: %d raw chars -> %d variable (kObs range %d-%d)\n",
            n_char_raw, mkd_local$nChar,
            min(kObs_lex), max(kObs_lex)))
rm(mat_list, combined_mat, pd, mkd_local)
invisible(gc(verbose = FALSE))

# ---- Read log via fread, dropping #-comment rows ----------------------------
# The MkPrime streamed logger emits periodic "# Topology:" / "# Branches:" /
# "# Characters:" / "# Rates:" acceptance lines interleaved with data rows;
# fread with sep="\t" + skip header + comment.char="#" handles both.
t0 <- Sys.time()
dt <- fread(log_file, sep = "\t", header = TRUE,
            data.table = TRUE, showProgress = FALSE,
            fill = Inf)
# Drop any leaked #-rows (Sample column will be NA for them).
dt <- dt[!is.na(Sample)]
cat(sprintf("  fread: %.1fs, n_samples=%d, ncols=%d\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            nrow(dt), ncol(dt)))

n_samples <- nrow(dt)
stopifnot(n_samples > 0L)

# ---- Posterior means --------------------------------------------------------
scalar_cols <- intersect(c("log_posterior", "log_likelihood", "tree_length",
                           "rate_log_sd", "p"),
                         names(dt))
scalar_means <- vapply(scalar_cols,
                       function(cn) mean(dt[[cn]], na.rm = TRUE),
                       numeric(1L))

br_cols <- grep("^br_", names(dt), value = TRUE)
br_means <- vapply(br_cols,
                   function(cn) mean(dt[[cn]], na.rm = TRUE),
                   numeric(1L))

if (arm %in% c("mkp_geo", "mkp_eg", "mkp")) {
  kp_cols <- grep("^kPrime_", names(dt), value = TRUE)
  kp_means <- vapply(kp_cols,
                     function(cn) mean(dt[[cn]], na.rm = TRUE),
                     numeric(1L))
  p_mean   <- if ("p" %in% names(dt)) mean(dt$p, na.rm = TRUE) else NA_real_
  # Sanity: kp_means length must equal kObs_lex length (both are over
  # variable characters in lex order). Warn but do not abort.
  if (length(kp_means) != length(kObs_lex)) {
    cat(sprintf("  WARN: kp_means length (%d) != kObs_lex length (%d)\n",
                length(kp_means), length(kObs_lex)))
  }
} else {
  kp_cols <- character(0L)
  kp_means <- numeric(0L)
  p_mean   <- NA_real_
}

# Last-row sanity-check vector.
last_row <- as.list(dt[n_samples,
                       intersect(c(scalar_cols, "Sample"), names(dt)),
                       with = FALSE])

rm(dt); invisible(gc(verbose = FALSE))

# ---- Stream trees & thin uniformly ------------------------------------------
# Count tree lines without loading them all.
n_trees <- length(count.fields(tree_file, sep = "\n", quote = "")) # robust
cat(sprintf("  n_trees=%d\n", n_trees))

if (n_trees == 0L) {
  stop("No trees in ", tree_file)
}

# Thin: take the last min(n_trees, 5*n_thin) and uniformly sample n_thin
# from them so we lean on post-burn-in samples.
tail_n   <- min(n_trees, max(n_thin, 5L * n_thin))
start_at <- n_trees - tail_n + 1L
keep_idx <- if (tail_n <= n_thin) {
              seq.int(start_at, n_trees)
            } else {
              round(seq.int(start_at, n_trees, length.out = n_thin))
            }
keep_set <- unique(keep_idx)
n_keep   <- length(keep_set)

t1 <- Sys.time()
# Single pass: read the file, retain only the rows we need.
all_lines <- readLines(tree_file, warn = FALSE)
sel_lines <- all_lines[keep_set]
rm(all_lines); invisible(gc(verbose = FALSE))

# Drop malformed Newicks (e.g. truncated last line on SIGKILL): require
# trailing ';' and matched parens.
ok <- grepl(";\\s*$", sel_lines) &
      (lengths(regmatches(sel_lines, gregexpr("\\(", sel_lines))) ==
       lengths(regmatches(sel_lines, gregexpr("\\)", sel_lines))))
n_dropped <- sum(!ok)
if (n_dropped > 0L) {
  cat(sprintf("  WARN: dropped %d malformed Newick line(s)\n", n_dropped))
  sel_lines <- sel_lines[ok]
}
stopifnot(length(sel_lines) > 0L)

trees <- read.tree(text = sel_lines)
if (inherits(trees, "phylo")) trees <- structure(list(trees),
                                                  class = "multiPhylo")
cat(sprintf("  trees read+parsed: %.1fs, kept=%d\n",
            as.numeric(Sys.time() - t1, units = "secs"),
            length(trees)))

# ---- True tree --------------------------------------------------------------
true_tree <- read.tree(true_tree_file)
true_tree <- ape::reorder.phylo(true_tree, "cladewise")

reorder_multiPhylo <- function(tt) {
  out <- lapply(tt, ape::reorder.phylo, order = "cladewise")
  class(out) <- "multiPhylo"
  out
}

# ---- CID --------------------------------------------------------------------
t2 <- Sys.time()
cid <- as.numeric(ClusteringInfoDistance(reorder_multiPhylo(trees),
                                          true_tree, normalize = TRUE))
cat(sprintf("  CID: %.1fs, n=%d, range=[%.3f, %.3f]\n",
            as.numeric(Sys.time() - t2, units = "secs"),
            length(cid), min(cid), max(cid)))

# ---- Assemble summary -------------------------------------------------------
summary_list <- list(
  tree_idx      = tree_idx,
  rep_idx       = rep_idx,
  arm           = arm,
  tag           = tag,
  n_samples     = n_samples,
  n_trees_total = n_trees,
  n_trees_kept  = length(trees),
  thinned_trees = trees,
  thinned_idx   = keep_set,
  cid           = cid,
  scalar_means  = scalar_means,
  br_means      = br_means,
  kp_means      = kp_means,
  # kObs paired to kp_means: position i of kObs_lex corresponds to kp_means[i].
  # Both are in lex-sorted chr*.nex order (matching run_one.R's plain sort()).
  # Do NOT use natural-sort kObs vectors from raw nexus files to align kp_means.
  kObs          = kObs_lex,
  p_mean        = p_mean,
  last_row      = last_row,
  log_file      = log_file,
  tree_file     = tree_file
)

saveRDS(summary_list, out_file, compress = "gzip")
size_mb <- file.info(out_file)$size / 1024^2
cat(sprintf("  Wrote %s  (%.2f MB)\n", out_file, size_mb))
cat("Done.\n")
