#!/usr/bin/env Rscript
# Cross-sampler equivalence test for one (pid, model) cell.
#
# Usage:
#   Rscript dev/rb-equivalence/compare.R <pid> <model> \
#     [--target-rhat=1.025] [--target-ess=256] [--tolerance=0.25] [--alpha=0.05]
#     [--max-trees=500] [--out-dir=dev/rb-equivalence/out]
#
# --max-trees caps the trees per source used to pick the pooled median tree
# and to report tree ESS; both are quadratic in the pool. CID to that median
# is still gated on every post-burn-in tree.
#
# Loads mkprime_<pid>_<model>.rds and rb_<pid>_<model>.rds. For each scalar
# parameter and for CID to the pooled median tree, GateParam() (R/gate.R)
# checks each sampler's own convergence, cross-sampler rank-normalised R-hat,
# and a standardised mean-difference test (family-wise `alpha`); a cell that
# could miss a `tolerance`-SD shift is UNDERPOWERED, not PASS. Per-source tree
# ESS is reported only.
# Writes one row per (param) to out/summary.csv and one row per (source)
# to out/summary_tree.csv. Exit status: 0 PASS, 1 FAIL, 2 UNDERPOWERED.

suppressPackageStartupMessages({
  library(ape)
  library(MkPrime)
  library(TreeDist)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: compare.R <pid> <model> [opts]")
pid <- args[[1]]
model <- args[[2]]

script_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- a[grep("^--file=", a)]
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})()
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "gate.R"))

opt <- list(target_rhat = 1.025, target_ess = 256, tolerance = 0.25,
            alpha = 0.05, max_trees = 500,
            out_dir = file.path(script_dir, "out"))
# target_rhat = 1.025 is the canonical equivalence threshold. Cross-sampler
# rhat above this currently flags two known source-level asymmetries (see
# dev/rb-equivalence/notes/cross-sampler-rhat-investigation.md) which are
# treated as bugs to fix, not noise to absorb.
for (a in args[-(1:2)]) {
  kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  if (length(kv) != 2L) stop("Bad option: ", a)
  opt[[gsub("-", "_", kv[[1]])]] <- kv[[2]]
}
opt$target_rhat <- as.numeric(opt$target_rhat)
opt$target_ess <- as.numeric(opt$target_ess)
opt$tolerance <- as.numeric(opt$tolerance)
opt$alpha <- as.numeric(opt$alpha)
opt$max_trees <- as.integer(opt$max_trees)

# --- Load both rds ---------------------------------------------------------

mk_path <- file.path(opt$out_dir, sprintf("mkprime_%s_%s.rds", pid, model))
rb_path <- file.path(opt$out_dir, sprintf("rb_%s_%s.rds", pid, model))
stopifnot(file.exists(mk_path), file.exists(rb_path))
mk <- readRDS(mk_path)
rb <- readRDS(rb_path)

# --- Sanity asserts --------------------------------------------------------

# model_spec is read from the MkPrimeModel object and from the Rev script RB
# ran, not hard-coded (RB-109).
AssertSameSetup(mk, rb)
AssertMkRunTargets(mk$posterior, opt$target_rhat, opt$target_ess)
if (!identical(mk$host, rb$host)) {
  warning(sprintf("HOST MISMATCH (wall-clock not comparable): mk=%s rb=%s",
                  mk$host, rb$host))
}

# agent-issues/MkPrime#214: refuse to compare cells that weren't provably fed
# the same k/nChar/taxa (both sides get this from the shared cellinfo_*.rds
# cache written by render_rev.R / run_mkprime.R, so a match here means the
# equality check already passed at run time -- this just refuses to compare
# if either rds predates that fix and has no cell_info at all).
if (is.null(mk$cell_info) || is.null(rb$cell_info)) {
  stop(
    "cell_info missing from mkprime and/or rb rds -- rerun run_mkprime.R / ",
    "render_rev.R + post_rb.R (agent-issues/MkPrime#214); refusing to ",
    "compare a cell with no data-provenance guarantee."
  )
}
if (!identical(mk$cell_info, rb$cell_info)) {
  stop("cell_info mismatch between mkprime and rb rds -- the two samplers ",
       "were not fed the same data/k/taxa for this cell (agent-issues/MkPrime#214).")
}

# agent-issues/MkPrime#215 (RB-107): a template edit must invalidate a
# previously-rendered rb_*.rds even though nothing here recomputes mtimes.
mkTH <- mk$provenance$templateHash
rbTH <- rb$provenance$templateHash
if (is.na(mkTH) || is.na(rbTH)) {
  warning("templateHash unavailable on one or both rds (older run, or no ",
          "git checkout at run time) -- cannot verify the Rev templates ",
          "haven't drifted since this cell was rendered.")
} else if (!identical(mkTH, rbTH)) {
  stop("templateHash mismatch between mkprime and rb rds -- the RB templates ",
       "changed between the two runs (agent-issues/MkPrime#215); rerun both.")
}

# agent-issues/MkPrime#215: a resumed MkPrime run's wall_ratio is not
# comparable (wall_total only times the resumed segment).
if (isTRUE(mk$resumed)) {
  warning("mk$resumed is TRUE: wall_ratio for this cell is not meaningful ",
          "(wall_total only covers the resumed segment).")
}

# --- Extract per-source scalar samples (MkP) -------------------------------

mk_scalar_per_run <- mk$per_run_scalars
if (is.null(mk_scalar_per_run) || !length(mk_scalar_per_run)) {
  stop("mk$per_run_scalars missing; rerun run_mkprime.R with the streamed-log fix")
}
AssertScalarsPresent(mk_scalar_per_run, "mkprime")
AssertScalarsPresent(rb$per_run_scalars, "rb")
scalar_cols <- ExpectedScalars

# --- Subsample trees to common N per source (bounds the median's O(N^2)) ---

mk_trees <- mk$per_run_trees
if (is.null(mk_trees) || !length(mk_trees)) {
  stop("mk$per_run_trees missing; rerun run_mkprime.R")
}
source_lengths <- lengths(c(mk_trees, rb$per_run_trees))
n_common <- min(source_lengths, opt$max_trees) # nolint: object_name_linter. # file-wide snake_case
cat(sprintf("[compare] %d tree sources; common N = %d\n",
            length(source_lengths), n_common))
source_labels <- c(paste0("mkp_", seq_along(mk_trees)),
                   paste0("rb_", seq_along(rb$per_run_trees)))

mk_trees_sub <- lapply(mk_trees, Subsample, n = n_common)
rb_trees_sub <- lapply(rb$per_run_trees, Subsample, n = n_common)
all_trees_list <- c(mk_trees_sub, rb_trees_sub)
names(all_trees_list) <- source_labels

# Flatten into one multiPhylo
pooled <- do.call(c, lapply(all_trees_list, function(x) {
  if (inherits(x, "multiPhylo")) x else structure(x, class = "multiPhylo")
}))

# Median = pool tree minimising sum CID to all others
dmat <- as.matrix(TreeDist::ClusteringInfoDistance(pooled, normalize = TRUE))
median_idx <- which.min(rowSums(dmat))
median_tree <- pooled[[median_idx]]
cat(sprintf("[compare] Pooled median tree: index %d / %d (sum-CID = %.3f)\n",
            median_idx, length(pooled), sum(dmat[median_idx, ])))

# Per-source CID-to-median, on every post-burn-in tree
CidToMedian <- function(trs) {
  as.numeric(TreeDist::ClusteringInfoDistance(trs, median_tree, normalize = TRUE))
}
mkCid <- lapply(mk_trees, CidToMedian)
rbCid <- lapply(rb$per_run_trees, CidToMedian)

# --- Gate: per-sampler convergence, cross R-hat, mean shift ---------------

gate <- do.call(rbind, c(
  lapply(scalar_cols, function(p) {
    GateParam(lapply(mk_scalar_per_run, `[[`, p),
              lapply(rb$per_run_scalars, `[[`, p),
              opt$target_rhat, opt$target_ess,
              logScale = p %in% LogScaleScalars, alpha = opt$alpha,
              nTests = length(scalar_cols) + 1L, tolerance = opt$tolerance)
  }),
  list(GateParam(mkCid, rbCid, opt$target_rhat, opt$target_ess,
                 alpha = opt$alpha, nTests = length(scalar_cols) + 1L,
                 tolerance = opt$tolerance))
))

# --- Tree ESS per source (reporting only) ---------------------------------

tree_ess_per_source <- vapply(all_trees_list, function(trs) {
  e <- try(
    MkPrime::TreeESS(trs, frechet = TRUE)[["frechetCorrelationESS"]],
    silent = TRUE
  )
  if (inherits(e, "try-error") || !is.finite(e)) NA_real_ else as.numeric(e)
}, numeric(1))

# --- Assemble long summary ------------------------------------------------

# agent-issues/MkPrime#215 (RB-107, last bullet): summary.csv previously
# de-duplicated purely on (pid, model, param), so a stale row could not be
# told apart from a current one. cellinfo_hash + mkp_git_sha let a reader
# (or a future automated check) see whether two rows for the same cell came
# from the same data/code, without changing the de-dup key itself -- the
# newest run for a cell is still what should win.
cellinfo_hash <- rlang::hash(mk$cell_info)

# RevBayes always runs to srMaxTime, so its time to target is pro-rated from
# its ESS (RB-112) to the ESS the MkPrime run targeted, so both walls are to
# the same target; MkPrime's is NA unless it stopped on convergence.
mkRunEss <- mk$posterior$mcmc$minEss %||% opt$target_ess
wallRbEst <- rb$wall_total * mkRunEss / rb$diag_scalar$minEss
rows <- data.frame(
  pid = pid, model = model,
  param = c(scalar_cols, "cid_to_median"),
  gate,
  wall_mkp = mk$wall_to_target %||% NA_real_,
  wall_mkp_tree_target_est = mk$wall_to_tree_target_estimated %||% NA_real_,
  wall_rb_est = wallRbEst,
  wall_ratio = (mk$wall_to_target %||% NA_real_) / wallRbEst,
  cellinfo_hash = cellinfo_hash,
  mkp_git_sha = mk$provenance$gitSha %||% NA_character_,
  mkp_resumed = isTRUE(mk$resumed),
  stringsAsFactors = FALSE
)

tree_rows <- data.frame(
  pid = pid, model = model,
  source = names(tree_ess_per_source),
  tree_ess = as.numeric(tree_ess_per_source),
  cellinfo_hash = cellinfo_hash,
  stringsAsFactors = FALSE
)

cat(sprintf(paste0(
  "\n=== summary: tolerance %.2f SD, family-wise alpha %.3f over %d params, ",
  "target R-hat < %.3f, ESS > %g ===\n"),
  opt$tolerance, opt$alpha, nrow(rows), opt$target_rhat, opt$target_ess))
print(rows, row.names = FALSE, digits = 4)
cat("\n=== tree ESS per source ===\n")
print(tree_rows, row.names = FALSE, digits = 4)

# --- Append to global summary CSVs ----------------------------------------

AppendCsv(rows, file.path(opt$out_dir, "summary.csv"))
AppendCsv(tree_rows, file.path(opt$out_dir, "summary_tree.csv"))

cat(sprintf("\n[compare] Updated %s and summary_tree.csv\n",
            file.path(opt$out_dir, "summary.csv")))

# --- Exit status: 0 PASS, 1 FAIL, 2 UNDERPOWERED; no exemptions -----------

failed <- rows[rows$verdict == "FAIL", ]
if (nrow(failed)) {
  cat(sprintf("\n[compare] FAIL: %s\n",
              paste0(failed$param, " (", failed$reasons, ")", collapse = ", ")))
  quit(status = 1L)
}
weak <- rows[rows$verdict == "UNDERPOWERED", ]
if (nrow(weak)) {
  cat(sprintf(paste0(
    "\n[compare] UNDERPOWERED: %s cannot detect a %.2f-SD shift ",
    "(smallest detectable: %s); run both samplers longer.\n"),
    paste(weak$param, collapse = ", "), opt$tolerance,
    paste(sprintf("%.2f", weak$mdd_sd), collapse = ", ")))
  quit(status = 2L)
}
cat(sprintf("\n[compare] PASS: no shift of %.2f SD or more on any parameter\n",
            opt$tolerance))
