#!/usr/bin/env Rscript
# Render RB cell scripts from templates for a given (pid, model).
#
# Usage:
#   Rscript dev/rb-equivalence/render_rev.R <pid> <model> \
#     [--ess=128] [--max-time=3600] \
#     [--neo-path=PATH] [--trans-path=PATH] \
#     [--out-dir=dev/rb-equivalence/rev/<pid>] [--cellinfo-dir=dev/rb-equivalence/out]
#
# --neo-path/--trans-path name the *source* nexus files (the neotrans
# project splits); this script reads them, restricts both to the taxa they
# have in common, and writes the result as <out-dir>/{trans,neo}_filtered.nex
# for RB to read (agent-issues/MkPrime#214, RB-119: RB's template previously
# read the unfiltered splits and clamped to `neo.names()` alone, which can
# differ from the harness's trans/neo taxon intersection). k (state count)
# and nChar are also derived here, from the same PrepareCellInfo() that
# run_mkprime.R calls, and rendered into the model template so both samplers
# see matching values; WriteOrCheckCellInfo() asserts the two scripts agree,
# whichever runs first.
#
# Writes <out-dir>/{<model>.Rev, long_<model>.Rev, trans_filtered.nex,
# neo_filtered.nex} ready for `rb long_<model>.Rev` run from <out-dir>.

suppressPackageStartupMessages({
  library(MkPrime)
  library(TreeTools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: render_rev.R <pid> <model> [opts]")
pid <- args[[1]]
model <- args[[2]]
stopifnot(model %in% c("by_nt_9v", "by_nt_kv"))

script_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- a[grep("^--file=", a)]
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})()
source(file.path(script_dir, "R", "utils.R"))

opt <- list(ess = 128, max_time = 3600,
            neo_path = sprintf("project%s.neo.nex", pid),
            trans_path = sprintf("project%s.trans.nex", pid),
            out_dir = file.path(script_dir, "rev", pid),
            cellinfo_dir = file.path(script_dir, "out"))
for (a in args[-(1:2)]) {
  kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  if (length(kv) != 2L) stop("Bad option: ", a)
  opt[[gsub("-", "_", kv[[1]])]] <- kv[[2]]
}
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

read_template <- function(name) {
  paste(readLines(file.path(script_dir, "templates", name)), collapse = "\n")
}

# --- Derive the canonical cell info (shared with run_mkprime.R) ------------

trans_mat <- ReadSplit(opt$trans_path)
neo_mat <- ReadSplit(opt$neo_path)
cellInfo <- PrepareCellInfo(trans_mat, neo_mat, model)

cellinfo_path <- file.path(opt$cellinfo_dir,
                           sprintf("cellinfo_%s_%s.rds", pid, model))
WriteOrCheckCellInfo(cellinfo_path, cellInfo)

# --- Write the taxon-matched, polymorphism-sanitised data RB will actually
# read -- the exact matrices MkPrimeData's own input was built from ---------

WriteFilteredNexus(cellInfo$transMatCanonical,
                   file.path(opt$out_dir, "trans_filtered.nex"))
WriteFilteredNexus(cellInfo$neoMatCanonical,
                   file.path(opt$out_dir, "neo_filtered.nex"))

# --- Model body: substitute the post-invariant-drop nChar -------------------

model_body <- read_template(sprintf("%s.template.Rev", model))
model_body <- gsub("\\$\\{NCHAR_NEO\\}", cellInfo$nNeoFinal, model_body, perl = TRUE)
model_body <- gsub("\\$\\{NCHAR_TRANS\\}", cellInfo$nTransFinal, model_body, perl = TRUE)
writeLines(model_body, file.path(opt$out_dir, sprintf("%s.Rev", model)))

# --- Substitute into the driver wrapper -------------------------------------

long <- read_template("long.template.Rev")
subs <- list(
  PID = pid, MODEL = model,
  NEO_PATH = "neo_filtered.nex", TRANS_PATH = "trans_filtered.nex",
  TARGET_ESS = format(as.integer(opt$ess)),
  MAX_TIME_SEC = format(as.integer(opt$max_time))
)
for (k in names(subs)) {
  long <- gsub(sprintf("\\$\\{%s\\}", k), subs[[k]], long, perl = TRUE)
}
writeLines(long, file.path(opt$out_dir, sprintf("long_%s.Rev", model)))

cat(sprintf(
  "[render_rev] Wrote %s/{%s.Rev, long_%s.Rev, trans_filtered.nex, neo_filtered.nex}\n",
  opt$out_dir, model, model
))
cat(sprintf("[render_rev] %d common taxa; nChar = (neo=%d, trans=%d)\n",
            length(cellInfo$commonTaxa), cellInfo$nNeoFinal, cellInfo$nTransFinal))
