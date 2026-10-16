#!/usr/bin/env Rscript
# Checks for aggregate-verdicts.R, against per-arm files written in exactly the
# format sbc.R and sbc-mixed.R produce.
#
# Run from the repository root:
#   Rscript dev/red-team/heavy-tests/test-aggregate-verdicts.R

script <- normalizePath("dev/red-team/heavy-tests/aggregate-verdicts.R",
                        winslash = "/", mustWork = TRUE)

ok <- function(label, cond) {
  cat(sprintf("%-58s %s\n", label, if (isTRUE(cond)) "PASS" else "*** FAIL"))
  if (!isTRUE(cond)) stop(label)
}

# Exactly sbc.R's layout: header block, blank line, "arms:", then
# sprintf("  %-30s %s", name, paste0(verdict, " (good=N)")).
WriteArmVerdict <- function(dir, arm, verdict, good = 200L, runId = RUN_ID) {
  lines <- c(
    "SBC harness top-level summary",
    sprintf("run_id:   %s", runId),
    "mode:     full",
    "seedBase: 20260528",
    "wall:     123.4s",
    "N_sim:    200  L:1000  thin:10",
    "",
    "arms:",
    sprintf("  %-30s %s", arm, paste0(verdict, sprintf(" (good=%d)", good)))
  )
  writeLines(lines, file.path(dir, sprintf("verdict-%s.txt", arm)))
}

RUN_ID <- "17330001"

Run <- function(dir, expected = character(0L), runId = RUN_ID) {
  runArg <- if (is.null(runId)) character(0L) else c("--run-id", runId)
  out <- system2("Rscript", c(shQuote(script), shQuote(dir), runArg, expected),
                 stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  list(text = paste(out, collapse = "\n"),
       status = if (is.null(status)) 0L else status)
}

ARMS <- c("MkNT_geometric", "Mkp_geometric", "Mkp_empirical_geometric",
          "Mkp_logseries", "MkNT_logseries")

## ---- All arms present and passing -------------------------------------------
d1 <- tempfile("sbc1"); dir.create(d1)
for (a in ARMS) WriteArmVerdict(d1, a, "PASS")
r1 <- Run(d1, ARMS)
ok("complete passing run exits 0", r1$status == 0L)
ok("complete passing run reports PASS", grepl("overall: PASS", r1$text))
agg1 <- readLines(file.path(d1, "verdict.txt"))
ok("every arm appears in the aggregate",
   all(vapply(ARMS, function(a) any(grepl(a, agg1, fixed = TRUE)), logical(1L))))

## ---- One arm never finished -- the case that used to be invisible -----------
d2 <- tempfile("sbc2"); dir.create(d2)
for (a in ARMS[-3L]) WriteArmVerdict(d2, a, "PASS")
r2 <- Run(d2, ARMS)
ok("missing arm exits non-zero", r2$status != 0L)
ok("missing arm is named", grepl(ARMS[3L], r2$text, fixed = TRUE))
ok("missing arm reports INCOMPLETE", grepl("overall: INCOMPLETE", r2$text))

## ---- A failing arm -----------------------------------------------------------
d3 <- tempfile("sbc3"); dir.create(d3)
for (a in ARMS) WriteArmVerdict(d3, a, if (a == "Mkp_logseries") "FAIL" else "PASS")
r3 <- Run(d3, ARMS)
ok("failing arm exits non-zero", r3$status != 0L)
ok("failing arm reports FAIL", grepl("overall: FAIL", r3$text))

## ---- EXEC_OK counts as passing (quick mode) ---------------------------------
d4 <- tempfile("sbc4"); dir.create(d4)
for (a in ARMS) WriteArmVerdict(d4, a, "EXEC_OK")
ok("EXEC_OK is not a failure", Run(d4, ARMS)$status == 0L)

## ---- sbc-mixed.R's own file is picked up ------------------------------------
d5 <- tempfile("sbc5"); dir.create(d5)
for (a in ARMS) WriteArmVerdict(d5, a, "PASS")
writeLines(c(
  "SBC mixed-partition summary",
  sprintf("run_id:   %s", RUN_ID),
  "mode:     full",
  "seedBase: 20260528",
  "N_sim:    200  L:1000  thin:10",
  "",
  sprintf("  %-30s %s", "MkNT_mixed", "PASS (good=200)")
), file.path(d5, "verdict-MkNT_mixed.txt"))
r5 <- Run(d5, c(ARMS, "MkNT_mixed"))
ok("sbc-mixed's verdict is aggregated too", r5$status == 0L)
ok("MkNT_mixed appears in the aggregate", grepl("MkNT_mixed", r5$text))

## ---- Two tasks reporting the same arm ---------------------------------------
d6 <- tempfile("sbc6"); dir.create(d6)
for (a in ARMS) WriteArmVerdict(d6, a, "PASS")
file.copy(file.path(d6, "verdict-Mkp_geometric.txt"),
          file.path(d6, "verdict-Mkp_geometric-rerun.txt"))
r6 <- Run(d6, ARMS)
ok("a duplicated arm is called out", grepl("DUPLICATE", r6$text))

## ---- A stale file from an earlier run must not stand in for a dead task ----
d7 <- tempfile("sbc7"); dir.create(d7)
for (a in ARMS) WriteArmVerdict(d7, a, "PASS", runId = "17220000")
for (a in ARMS[-3L]) WriteArmVerdict(d7, a, "PASS")
r7 <- Run(d7, ARMS)
ok("a previous run's verdict is not counted", r7$status != 0L)
ok("the arm whose task died this run is INCOMPLETE",
   grepl("overall: INCOMPLETE", r7$text) &&
     grepl(sprintf("no verdict: %s", ARMS[3L]), r7$text, fixed = TRUE))
ok("the stale file is named", grepl("STALE", r7$text) &&
     grepl(sprintf("verdict-%s.txt", ARMS[3L]), r7$text, fixed = TRUE))

## ---- A verdict with no run id at all is stale too ----------------------------
d8 <- tempfile("sbc8"); dir.create(d8)
for (a in ARMS) WriteArmVerdict(d8, a, "PASS")
txt <- readLines(file.path(d8, sprintf("verdict-%s.txt", ARMS[1L])))
writeLines(txt[!grepl("^run_id:", txt)],
           file.path(d8, sprintf("verdict-%s.txt", ARMS[1L])))
ok("an unstamped verdict is not counted", Run(d8, ARMS)$status != 0L)

## ---- The run id is required --------------------------------------------------
ok("aggregating without --run-id is refused", Run(d1, ARMS, runId = NULL)$status != 0L)

## ---- The aggregate never silently replaces a single task's output -----------
ok("the aggregate is a separate file from any verdict-<arm>.txt",
   file.exists(file.path(d1, "verdict.txt")) &&
     length(list.files(d1, "^verdict-.*\\.txt$")) == length(ARMS))

cat("\nALL CHECKS PASSED\n")
