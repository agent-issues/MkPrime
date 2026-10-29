#!/usr/bin/env Rscript
## Checks eg001_upost_compare.R's ground-truth join against fabricated data.
##
## Run from the repository root:
##   Rscript dev/pilots/2026-05-12-prior-validation/analysis/test_eg001_gt_for.R
##
## The helpers are lifted out of the analysis script by name, so this cannot
## drift from what it runs.

src    <- parse("dev/pilots/2026-05-12-prior-validation/analysis/eg001_upost_compare.R")
wanted <- c("gt_for", "check_lengths")
env    <- new.env(parent = globalenv())
env$GT_ROOT <- tempfile("gt")
found  <- character(0L)
for (e in src) {
  if (is.call(e) && identical(as.character(e[[1L]]), "<-") &&
      is.name(e[[2L]]) && as.character(e[[2L]]) %in% wanted) {
    eval(e, envir = env)
    found <- c(found, as.character(e[[2L]]))
  }
}
stopifnot(setequal(found, wanted))

ok <- function(label, cond) {
  cat(sprintf("%-58s %s\n", label, if (isTRUE(cond)) "PASS" else "*** FAIL"))
  if (!isTRUE(cond)) stop(label)
}
Fails <- function(expr) inherits(tryCatch(expr, error = identity), "error")

repDir <- file.path(env$GT_ROOT, "tree_01", "rep_01")
dir.create(repDir, recursive = TRUE)
nChar <- 11L
write.csv(data.frame(char_idx = seq_len(nChar), kObs = seq_len(nChar) + 1L,
                     k_true = seq_len(nChar) + 3L, u_true = 2L),
          file.path(repDir, "ground_truth.csv"), row.names = FALSE)
invisible(file.create(file.path(repDir, sprintf("chr%d.nex", seq_len(nChar)))))

ok("no recorded char_idx is refused, not assumed lexical",
   Fails(env$gt_for("t01_r01", NULL)))

variable <- c(1:4, 6:11)
gt <- env$gt_for("t01_r01", variable)
ok("a variable-only char_idx joins", identical(gt$char_idx, variable))
ok("rows follow char_idx, not file order",
   identical(env$gt_for("t01_r01", c(10L, 2L))$char_idx, c(10L, 2L)))
ok("an index absent from ground_truth.csv is refused",
   Fails(env$gt_for("t01_r01", c(1L, 12L))))
ok("a duplicated index is refused", Fails(env$gt_for("t01_r01", c(1L, 1L))))

ok("matching lengths pass", !Fails(env$check_lengths("t", rep(3, 10L), gt)))
ok("a length mismatch stops instead of dropping the task",
   Fails(env$check_lengths("t", rep(3, 11L), gt)))

cat("\nALL CHECKS PASSED\n")
