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
#   arm           one of the fourteen arms match.arg'd at `arm <-` below
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
arm          <- match.arg(args[3], c("mk", "mkp_geo", "mkp_eg", "mkp", "mk_kp1", "mk_kp2", "mk_k9", "mk_k15", "mk_k24", "mk_k40", "mk_ktrue", "mkp_highk", "mkp_logs", "mk_tlshrink"))
results_root <- args[4]
data_root    <- args[5]
out_dir      <- args[6]
n_thin       <- if (length(args) >= 7L) as.integer(args[7]) else 1000L

tag       <- sprintf("t%02d_r%02d", tree_idx, rep_idx)
task_dir  <- file.path(results_root, tag)

# Per-run files, ordered by run index.
#
# `list.files()` sorts lexically, which puts run 10 between 1 and 2; and the
# glob used to look for a `_run_` infix (`{arm}_trees_run_1.nwk`) that
# `.TreeFilePaths()` has never written -- it writes `{arm}_trees_1.nwk`. So the
# summariser aborted at the stopifnot below on every current-layout task, and
# no summary it produced can have come from one (#99, #104).
# `.log` files may have been gzipped by a previous pass (see the cleanup at the
# foot of this script); fread reads .gz directly.
.run_files <- function(suffix, ext) {
  found <- list.files(
    task_dir,
    pattern = sprintf("^%s_%s_([0-9]+)\\.%s(\\.gz)?$", arm, suffix, ext),
    full.names = TRUE
  )
  if (!length(found)) return(character(0L))
  idx <- as.integer(sub(sprintf("^.*_%s_([0-9]+)\\.%s(\\.gz)?$", suffix, ext),
                        "\\1", basename(found)))
  found[order(idx)]
}

# Pre-2026-05-20 runs passed one shared path to every run, so the legacy single
# file is `[run 1][run 2]` CONCATENATED, not interleaved -- runs execute
# sequentially. It is kept readable, but its run boundary is unrecoverable, so
# it is summarised as a single chain and flagged as such.
legacy_log   <- file.path(task_dir, sprintf("%s_run.log", arm))
legacy_tree  <- file.path(task_dir, sprintf("%s_trees.nwk", arm))
log_files    <- .run_files("run", "log")
if (!length(log_files) && file.exists(legacy_log)) log_files <- legacy_log
tree_files   <- .run_files("trees", "nwk")
if (!length(tree_files) && file.exists(legacy_tree)) tree_files <- legacy_tree
legacy_layout <- identical(tree_files, legacy_tree)
true_tree_file <- file.path(data_root,
                            sprintf("tree_%02d/tree.nwk", tree_idx))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("%s_%s.rds", arm, tag))

cat(sprintf("[%s] tag=%s arm=%s%s\n", format(Sys.time(), "%H:%M:%S"), tag, arm,
            if (legacy_layout) "  [legacy single-file layout]" else ""))
cat("  log  : ", paste(log_files,  collapse = ", "), "\n")
cat("  tree : ", paste(tree_files, collapse = ", "), "\n")
cat("  true : ", true_tree_file, "\n")
cat("  out  : ", out_file,  "\n")

stopifnot(length(log_files) > 0L, all(file.exists(log_files)),
          length(tree_files) > 0L,
          all(file.exists(tree_files)), file.exists(true_tree_file))

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
# cbind joins positionally and keeps the first matrix's rownames (#104).
stopifnot(all(vapply(mat_list,
                     function(m) identical(rownames(m), rownames(mat_list[[1L]])),
                     logical(1L))))
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

# ---- Read every run's log, drop #-comment rows, drop per-run burn-in --------
# The MkPrime streamed logger writes a "# Topology:" / "# Branches:" /
# "# Characters:" / "# Rates:" acceptance block once per run, at the
# Warmup -> Sample transition (R/RunMkPrime.R:1360, :1430) -- not periodically,
# as this comment used to claim. `fread` has no comment.char, so those rows
# arrived as data and the `!is.na(Sample)` guard did not remove them: Sample is
# coerced to character, so the comment text is a non-NA string. Every mean was
# unaffected (all pass na.rm = TRUE) but `n_samples` over-reported by exactly
# the comment-line count (#104).
#
# Reading `_run_1` alone and then deleting `_run_2` unread, as this script did,
# made the summary internally inconsistent -- scalars from one run, topology
# from both -- and destroyed the cross-run material that `nRuns = 2L,
# maxRhat = 1.1` is paid for (#102).
.BURNIN_FRAC <- 0.25

t0 <- Sys.time()
per_run <- lapply(log_files, function(f) {
  d <- fread(f, sep = "\t", header = TRUE,
             data.table = TRUE, showProgress = FALSE, fill = Inf)
  d <- d[!grepl("^#", as.character(d[[1L]]))]
  # The comment block has no tabs, so it lands entirely in column 1 and coerces
  # `Sample` to character while leaving every other column numeric. Put it back.
  d[, Sample := as.numeric(as.character(Sample))]
  d[!is.na(Sample)]
})
n_raw_per_run <- vapply(per_run, nrow, integer(1L))

# Burn-in is dropped per run and the remainder pooled. A tail taken after
# pooling can consist entirely of the last run (#99).
per_run <- lapply(per_run, function(d) {
  n <- nrow(d)
  if (n == 0L) return(d)
  d[seq.int(floor(n * .BURNIN_FRAC) + 1L, n)]
})
n_kept_per_run <- vapply(per_run, nrow, integer(1L))

dt <- data.table::rbindlist(per_run, use.names = TRUE, fill = TRUE)
cat(sprintf("  fread: %.1fs, runs=%d, rows/run=%s, kept/run=%s, ncols=%d\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            length(log_files),
            paste(n_raw_per_run, collapse = "+"),
            paste(n_kept_per_run, collapse = "+"),
            ncol(dt)))

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

# Driven off what the log actually contains, not off a hard-coded arm list that
# every new arm has to remember to join. `mkp_highk` -- whose entire purpose is
# a Beta(1, 20) hyperprior on `p` -- was absent from that list, so its twelve
# kPrime_* columns and its `p` column were discarded unread and then deleted
# (#99). `mkp_logs` legitimately has no `p`: .ParamNames omits it for the
# logseries prior (R/RunMkPrime.R:4251-4255), which the `%in% names(dt)` test
# below already handles.
kp_cols  <- grep("^kPrime_", names(dt), value = TRUE)
kp_means <- vapply(kp_cols,
                   function(cn) mean(dt[[cn]], na.rm = TRUE),
                   numeric(1L))
p_mean   <- if ("p" %in% names(dt)) mean(dt$p, na.rm = TRUE) else NA_real_
# Sanity: kp_means length must equal kObs_lex length (both are over
# variable characters in lex order). Warn but do not abort.
if (length(kp_means) && length(kp_means) != length(kObs_lex)) {
  cat(sprintf("  WARN: kp_means length (%d) != kObs_lex length (%d)\n",
              length(kp_means), length(kObs_lex)))
}

# Last-row sanity-check vector.
last_row <- as.list(dt[n_samples,
                       intersect(c(scalar_cols, "Sample"), names(dt)),
                       with = FALSE])

rm(dt); invisible(gc(verbose = FALSE))

# ---- Stream trees & thin uniformly ------------------------------------------
# Concatenate lines across all per-run tree files (or just the single
# legacy file). For CID, all trees go into one bag; for future R-hat on
# tree statistics, the per-run files remain separate on disk.
t1 <- Sys.time()
run_lines     <- lapply(tree_files, readLines, warn = FALSE)
n_trees_run   <- vapply(run_lines, length, integer(1L))
n_trees       <- sum(n_trees_run)
cat(sprintf("  n_trees=%d (%s across %d file(s))\n",
            n_trees, paste(n_trees_run, collapse = "+"), length(tree_files)))

if (n_trees == 0L) {
  stop("No trees in ", paste(tree_files, collapse = ", "))
}

# Thin PER RUN, then pool.
#
# The old code pooled first and kept the last 5 * n_thin of the concatenation.
# Because the runs are concatenated rather than interleaved, that silently
# discarded whole runs: with 10000 trees in run 1 and 5000 in run 2, the window
# [10001, 15000] took nothing at all from run 1 (#99). Each run now contributes
# its own tail, so the bag is balanced and `n_trees_per_run` records the split
# that the saved object previously could not recover.
.thin_one <- function(lines, n_target) {
  n <- length(lines)
  if (n == 0L) return(integer(0L))
  tail_n   <- min(n, max(n_target, 5L * n_target))
  start_at <- n - tail_n + 1L
  idx <- if (tail_n <= n_target) {
    seq.int(start_at, n)
  } else {
    round(seq.int(start_at, n, length.out = n_target))
  }
  unique(idx)
}

n_thin_run <- max(1L, n_thin %/% max(1L, length(run_lines)))
keep_per_run <- lapply(run_lines, .thin_one, n_target = n_thin_run)

# Offsets make the kept indices global, so `thinned_idx` still points into the
# concatenated stream the way a reader expects.
offsets  <- cumsum(c(0L, n_trees_run[-length(n_trees_run)]))
keep_set <- unlist(Map(function(idx, off) idx + off, keep_per_run, offsets),
                   use.names = FALSE)
sel_lines <- unlist(Map(function(lines, idx) lines[idx], run_lines, keep_per_run),
                    use.names = FALSE)
rm(run_lines); invisible(gc(verbose = FALSE))

# Drop malformed Newicks (e.g. truncated last line on SIGKILL): require
# trailing ';' and matched parens. `keep_set` is filtered alongside, because it
# used to be fixed BEFORE this step -- so any drop desynchronised thinned_idx
# from thinned_trees with no record of which entries went (#104).
ok <- grepl(";\\s*$", sel_lines) &
      (lengths(regmatches(sel_lines, gregexpr("\\(", sel_lines))) ==
       lengths(regmatches(sel_lines, gregexpr("\\)", sel_lines))))
n_dropped <- sum(!ok)
if (n_dropped > 0L) {
  cat(sprintf("  WARN: dropped %d malformed Newick line(s)\n", n_dropped))
  sel_lines <- sel_lines[ok]
  keep_set  <- keep_set[ok]
}
stopifnot(length(sel_lines) > 0L, length(keep_set) == length(sel_lines))

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
  # Per-run counts, so the split stays recoverable after the raws are gone.
  n_samples_raw_per_run  = n_raw_per_run,
  n_samples_kept_per_run = n_kept_per_run,
  n_trees_per_run        = n_trees_run,
  burnin_frac            = .BURNIN_FRAC,
  legacy_layout          = legacy_layout,
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
  # The numeric character index behind each lex position, so a consumer can
  # join on char_idx rather than assume an order. Pairing lex-ordered posteriors
  # against numeric-ordered ground_truth.csv is EG-003 (#54).
  char_idx      = as.integer(sub("^chr([0-9]+)\\.nex$", "\\1",
                                 basename(nex_files_lex))),
  p_mean        = p_mean,
  last_row      = last_row,
  log_files     = log_files,
  tree_files    = tree_files
)

saveRDS(summary_list, out_file, compress = "gzip")
size_mb <- file.info(out_file)$size / 1024^2
cat(sprintf("  Wrote %s  (%.2f MB)\n", out_file, size_mb))

# Self-cleaning. The quota pressure is real -- 12 arms x 260 tasks accumulates
# ~400 GB and exhausted /nobackup on 2026-05-16 and 2026-05-20 -- but this
# script used to answer it by deleting the logs in the same pass that
# summarised them, including `_run_2.log`, which it had never read. If the
# summary is wrong, the evidence for saying so is already gone; that is exactly
# how EG-003 (#54) would have become uncorrectable (#99, #102).
#
# The bulk is the tree streams, so those are removed; the logs, which carry
# every scalar and every kPrime_ column, are gzipped instead (~10x, and fread
# reads .gz directly). The checkpoint is kept either way, so a run can still be
# extended with a longer walltime (see feedback_resumable_runs.md).
if (file.exists(out_file) && file.info(out_file)$size > 1024L) {
  gone <- tree_files[file.exists(tree_files)]
  if (length(gone) > 0L) unlink(gone)

  kept <- log_files[file.exists(log_files) & !grepl("\\.gz$", log_files)]
  for (f in kept) {
    ok <- tryCatch({
      system2("gzip", c("-f", shQuote(f)), stdout = FALSE, stderr = FALSE) == 0L
    }, error = function(e) FALSE)
    if (!ok) cat(sprintf("  WARN: could not gzip %s; left uncompressed\n", f))
  }
  cat(sprintf("  Cleaned %d tree stream(s), gzipped %d log(s) for %s/%s%s\n",
              length(gone), length(kept), tag, arm, " (checkpoint kept)"))
}
cat("Done.\n")
