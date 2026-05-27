#!/usr/bin/env Rscript
# Post-process a completed RB cell into the rds schema compatible with the
# MkPrime driver output.  Expects to be run in the RB working directory
# (where <model>_run_1.log, <model>_run_2.log, etc. live).
#
# Usage:
#   Rscript dev/rb-equivalence/post_rb.R <pid> <model> --wall=<sec> \
#     [--out-dir=dev/rb-equivalence/out] [--rb-dir=.]

suppressPackageStartupMessages({
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: post_rb.R <pid> <model> --wall=<sec> [opts]")
pid <- args[[1]]
model <- args[[2]]

script_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- a[grep("^--file=", a)]
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})()

opt <- list(wall = NA_real_,
            out_dir = file.path(script_dir, "out"),
            rb_dir = ".")
for (a in args[-(1:2)]) {
  kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  if (length(kv) != 2L) stop("Bad option: ", a)
  opt[[gsub("-", "_", kv[[1]])]] <- kv[[2]]
}
opt$wall <- as.numeric(opt$wall)
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load per-run scalar logs ----------------------------------------------

read_rb_log <- function(f) {
  if (!file.exists(f)) return(NULL)
  d <- read.delim(f, stringsAsFactors = FALSE, comment.char = "")
  # Drop Iteration column; keep all scalar params
  d
}

# RB writes <model>_run_<N>.log (Tracer format)
log_files <- list.files(opt$rb_dir, pattern = sprintf("^%s_run_\\d+\\.log$", model),
                        full.names = TRUE)
log_files <- sort(log_files)
if (!length(log_files)) {
  stop("No <model>_run_*.log files under ", opt$rb_dir)
}
per_run_scalars <- lapply(log_files, read_rb_log)
nRuns <- length(per_run_scalars)

# --- Load per-run trees ----------------------------------------------------

read_rb_trees <- function(f) {
  if (!file.exists(f)) return(NULL)
  tab <- read.delim(f, stringsAsFactors = FALSE)
  trs <- lapply(tab$phylogeny, function(s) ape::read.tree(text = s))
  structure(trs, class = "multiPhylo")
}

tree_files <- list.files(opt$rb_dir, pattern = sprintf("^%s_run_\\d+\\.trees$", model),
                         full.names = TRUE)
tree_files <- sort(tree_files)
per_run_trees <- lapply(tree_files, read_rb_trees)

# --- Burn-in: drop first 25% ----------------------------------------------

drop_burnin <- function(x, frac = 0.25) {
  n <- if (is.data.frame(x)) nrow(x) else length(x)
  if (n < 4) return(x)
  start <- max(2L, floor(n * frac))
  if (is.data.frame(x)) x[start:n, , drop = FALSE] else x[start:n]
}

per_run_scalars_pb <- lapply(per_run_scalars, drop_burnin)
per_run_trees_pb <- lapply(per_run_trees, drop_burnin)

# --- Diagnostics: scalar ESS + R-hat ---------------------------------------

scalar_cols <- intersect(
  c("tree_length", "rate_log_sd", "rate_loss", "rate_neo"),
  Reduce(intersect, lapply(per_run_scalars_pb, names))
)

ess_per_run <- lapply(per_run_scalars_pb, function(d) {
  vapply(scalar_cols, function(p) coda::effectiveSize(d[[p]]), numeric(1))
})
ess <- Reduce(`+`, ess_per_run)  # summed across runs
names(ess) <- scalar_cols

# R-hat via posterior::rhat_basic when available
rhat <- if (requireNamespace("posterior", quietly = TRUE) && nRuns >= 2L) {
  vapply(scalar_cols, function(p) {
    mat <- do.call(cbind, lapply(per_run_scalars_pb, function(d) d[[p]]))
    # Truncate to common length
    nMin <- min(vapply(per_run_scalars_pb, function(d) length(d[[p]]), integer(1)))
    mat <- do.call(cbind, lapply(per_run_scalars_pb,
                                 function(d) d[[p]][seq_len(nMin)]))
    posterior::rhat_basic(mat)
  }, numeric(1))
} else {
  setNames(rep(NA_real_, length(scalar_cols)), scalar_cols)
}

cat(sprintf("[post_rb] scalar params: %s\n", paste(scalar_cols, collapse = ", ")))
cat("[post_rb] ESS per param:\n"); print(round(ess, 1))
cat("[post_rb] R-hat per param:\n"); print(round(rhat, 4))

# --- Save -----------------------------------------------------------------

result <- list(
  pid = pid,
  model = model,
  wall_to_target = opt$wall,
  wall_total = opt$wall,
  wall_to_tree_target_estimated = NA_real_,  # RB doesn't expose tree ESS
  diag_scalar = list(ess = ess, rhat = rhat,
                     minEss = min(ess, na.rm = TRUE),
                     maxRhat = max(rhat, na.rm = TRUE)),
  diag_tree = NULL,  # filled by compare.R using MkPrime::FrechetESS for fairness
  per_run_scalars = per_run_scalars_pb,
  per_run_trees = per_run_trees_pb,
  host = Sys.info()[["nodename"]],
  sysinfo = list(
    sysname = Sys.info()[["sysname"]],
    release = Sys.info()[["release"]],
    machine = Sys.info()[["machine"]],
    rversion = R.version.string,
    started = NA_character_, finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ),
  coding = "variable",
  nCat = 6L,
  prior_spec = list(
    treeLengthShape = 2, treeLengthRate = 2,
    rateLogSdShape = 1, rateLogSdRate = 1,
    rateLossMeanlog = 0, rateLossSdlog = 2,
    rateNeoMeanlog = 0, rateNeoSdlog = 2
  )
)
out_path <- file.path(opt$out_dir, sprintf("rb_%s_%s.rds", pid, model))
saveRDS(result, out_path)
cat(sprintf("[post_rb] Saved: %s\n", out_path))
