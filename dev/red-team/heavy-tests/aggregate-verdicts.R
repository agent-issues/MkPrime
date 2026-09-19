#!/usr/bin/env Rscript
# Reduce the per-arm SBC verdicts into one run-level verdict.
#
# Each SLURM array task runs one arm and writes `verdict-<arm>.txt` into the
# shared results root (#16). This is the aggregation step that turns those into
# `verdict.txt` — separately, and after the array has finished, so that the
# aggregate is never a single task's output wearing the run-level name.
#
# Usage:
#   Rscript dev/red-team/heavy-tests/aggregate-verdicts.R <outRoot> [expected_arms...]
#
# With `expected_arms` given, a missing arm is an explicit INCOMPLETE rather
# than a silently shorter table — which is the whole point: a partial SBC
# verdict that looks complete is worse than no verdict.
#
# Exit status: 0 if every expected arm is present and PASS or EXEC_OK;
# 1 otherwise.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("usage: aggregate-verdicts.R <outRoot> [expected_arms...]")
}
outRoot  <- args[1L]
expected <- if (length(args) > 1L) args[-1L] else character(0L)

stopifnot(dir.exists(outRoot))

files <- list.files(outRoot, "^verdict-.*\\.txt$", full.names = TRUE)
if (!length(files)) {
  stop("No verdict-<arm>.txt files in ", outRoot)
}

# The per-arm file's arm lines look like "  <name>  <VERDICT> (good=N)".
ParseArmLines <- function(path) {
  txt <- readLines(path, warn = FALSE)
  hit <- grep("^\\s{2}\\S+\\s+(PASS|FAIL|ERROR|EXEC_OK|MARGINAL)\\b", txt,
              value = TRUE)
  if (!length(hit)) return(NULL)
  nm  <- sub("^\\s+(\\S+)\\s+.*$", "\\1", hit)
  vd  <- sub("^\\s+\\S+\\s+(\\S+).*$", "\\1", hit)
  data.frame(arm = nm, verdict = vd, source = basename(path),
             stringsAsFactors = FALSE)
}

per_arm <- do.call(rbind, lapply(files, ParseArmLines))
if (is.null(per_arm) || !nrow(per_arm)) {
  stop("No parseable arm lines in ", paste(basename(files), collapse = ", "))
}

dup <- per_arm$arm[duplicated(per_arm$arm)]
missing <- setdiff(expected, per_arm$arm)

lines <- c(
  "SBC run-level verdict (aggregated)",
  sprintf("aggregated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("sources:    %d file(s)", length(files)),
  "",
  "arms:"
)
for (i in seq_len(nrow(per_arm))) {
  lines <- c(lines, sprintf("  %-30s %-10s  [%s]",
                            per_arm$arm[i], per_arm$verdict[i],
                            per_arm$source[i]))
}
if (length(dup)) {
  lines <- c(lines, "", sprintf("DUPLICATE arm reported more than once: %s",
                                paste(unique(dup), collapse = ", ")))
}
if (length(missing)) {
  lines <- c(lines, "", sprintf("INCOMPLETE — expected arm(s) with no verdict: %s",
                                paste(missing, collapse = ", ")))
}

bad <- !per_arm$verdict %in% c("PASS", "EXEC_OK")
overall <- if (length(missing)) {
  "INCOMPLETE"
} else if (any(bad)) {
  "FAIL"
} else {
  "PASS"
}
lines <- c(lines, "", sprintf("overall: %s", overall))

writeLines(lines, file.path(outRoot, "verdict.txt"))
cat(paste0(lines, "\n"), sep = "")

quit(save = "no", status = if (identical(overall, "PASS")) 0L else 1L)
