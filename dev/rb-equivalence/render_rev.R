#!/usr/bin/env Rscript
# Render RB cell scripts from templates for a given (pid, model).
#
# Usage:
#   Rscript dev/rb-equivalence/render_rev.R <pid> <model> \
#     [--ess=128] [--max-time=3600] \
#     [--neo-path=PATH] [--trans-path=PATH] \
#     [--out-dir=dev/rb-equivalence/rev/<pid>]
#
# Writes <out-dir>/{<model>.Rev, long_<model>.Rev} ready for `rb long_<model>.Rev`.

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

opt <- list(ess = 128, max_time = 3600,
            neo_path = sprintf("project%s.neo.nex", pid),
            trans_path = sprintf("project%s.trans.nex", pid),
            out_dir = file.path(script_dir, "rev", pid))
for (a in args[-(1:2)]) {
  kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  if (length(kv) != 2L) stop("Bad option: ", a)
  opt[[gsub("-", "_", kv[[1]])]] <- kv[[2]]
}
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

read_template <- function(name) {
  paste(readLines(file.path(script_dir, "templates", name)), collapse = "\n")
}

# Copy the model body verbatim
model_body <- read_template(sprintf("%s.template.Rev", model))
writeLines(model_body, file.path(opt$out_dir, sprintf("%s.Rev", model)))

# Substitute into the driver wrapper
long <- read_template("long.template.Rev")
subs <- list(
  PID = pid, MODEL = model,
  NEO_PATH = opt$neo_path, TRANS_PATH = opt$trans_path,
  TARGET_ESS = format(as.integer(opt$ess)),
  MAX_TIME_SEC = format(as.integer(opt$max_time))
)
for (k in names(subs)) {
  long <- gsub(sprintf("\\$\\{%s\\}", k), subs[[k]], long, perl = TRUE)
}
writeLines(long, file.path(opt$out_dir, sprintf("long_%s.Rev", model)))

cat(sprintf("[render_rev] Wrote %s/{%s.Rev, long_%s.Rev}\n",
            opt$out_dir, model, model))
