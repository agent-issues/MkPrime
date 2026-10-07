#!/usr/bin/env Rscript
# Reduce the per-arm SBC verdicts into one run-level verdict.
#
# Each SLURM array task runs one arm and writes `verdict-<arm>.txt` into the
# shared results root (#16). This is the aggregation step that turns those into
# `verdict.txt` — separately, and after the array has finished, so that the
# aggregate is never a single task's output wearing the run-level name.
#
# Usage:
#   Rscript dev/red-team/heavy-tests/aggregate-verdicts.R <outRoot> --run-id <id> [expected_arms...]
#
# With `expected_arms` given, a missing arm is an explicit INCOMPLETE rather
# than a silently shorter table — which is the whole point: a partial SBC
# verdict that looks complete is worse than no verdict.
#
# `outRoot` is reused across submissions, so a task that dies leaves the
# previous run's `verdict-<arm>.txt` in place. Only files stamped with this
# run's `run_id:` count; any other is reported STALE and its arm falls through
# to INCOMPLETE.
#
# Exit status: 0 if every expected arm is present and PASS or EXEC_OK;
# 1 otherwise.

usage <- "usage: aggregate-verdicts.R <outRoot> --run-id <id> [expected_arms...]"
args <- commandArgs(trailingOnly = TRUE)
runIx <- which(args == "--run-id")
if (length(runIx) != 1L || runIx == length(args) || !nzchar(args[runIx + 1L])) {
  stop(usage)
}
runId <- args[runIx + 1L]
args <- args[-c(runIx, runIx + 1L)]
if (!length(args)) {
  stop(usage)
}
outRoot  <- args[1L]
expected <- if (length(args) > 1L) args[-1L] else character(0L)

stopifnot(dir.exists(outRoot))

files <- list.files(outRoot, "^verdict-.*\\.txt$", full.names = TRUE)
if (!length(files)) {
  stop("No verdict-<arm>.txt files in ", outRoot)
}

FileRunId <- function(path) {
  hit <- grep("^run_id:", readLines(path, warn = FALSE), value = TRUE)
  if (length(hit) != 1L) return(NA_character_)
  trimws(sub("^run_id:", "", hit))
}
stale <- files[!vapply(files, function(f) identical(FileRunId(f), runId),
                       logical(1L))]
files <- setdiff(files, stale)

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
if (is.null(per_arm)) {
  if (length(files)) {
    stop("No parseable arm lines in ", paste(basename(files), collapse = ", "))
  }
  per_arm <- data.frame(arm = character(0L), verdict = character(0L),
                        source = character(0L), stringsAsFactors = FALSE)
}

dup <- per_arm$arm[duplicated(per_arm$arm)]
missing <- setdiff(expected, per_arm$arm)

lines <- c(
  "SBC run-level verdict (aggregated)",
  sprintf("aggregated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("run_id:     %s", runId),
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
if (length(stale)) {
  lines <- c(lines, "", sprintf("STALE — ignored, not stamped run_id %s: %s", runId,
                                paste(basename(stale), collapse = ", ")))
}
if (length(missing)) {
  lines <- c(lines, "", sprintf("INCOMPLETE — expected arm(s) with no verdict: %s",
                                paste(missing, collapse = ", ")))
}

bad <- !per_arm$verdict %in% c("PASS", "EXEC_OK")
overall <- if (length(missing) || !nrow(per_arm)) {
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
