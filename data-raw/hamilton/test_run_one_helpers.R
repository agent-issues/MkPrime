#!/usr/bin/env Rscript
## Unit checks for the four helpers run_one.R gained in the area-4 pass.
##
## Run from the repository root:
##   Rscript data-raw/hamilton/test_run_one_helpers.R
##
## The helpers are lifted out of run_one.R by name rather than re-stated here,
## so this cannot drift from what the harness actually runs.

suppressPackageStartupMessages(library(MkPrime))

src   <- parse("data-raw/hamilton/run_one.R")
wanted <- c(".kObsFromMatrix", ".UPostMeans", ".BURNIN_FRAC", ".arm_own_files")
env   <- new.env(parent = globalenv())
## .arm_own_files() closes over the script's ckp_dir.
env$ckp_dir <- tempfile("ckp"); dir.create(env$ckp_dir)
found <- character(0L)
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

## ---- .kObsFromMatrix ---------------------------------------------------------
mat <- cbind(
  plain      = c("0", "1", "0", "1"),
  with_gaps  = c("0", "1", "?", "-"),
  polymorph  = c("0", "1", "{01}", "?"),
  paren_poly = c("0", "1", "(01)", "-"),
  invariant  = c("1", "1", "1", "?")
)
k <- env$.kObsFromMatrix(mat)
ok("plain binary counts 2", k[["plain"]] == 2L)
ok("? and - are not states", k[["with_gaps"]] == 2L)
ok("{01} is not a third state", k[["polymorph"]] == 2L)
ok("(01) is not a third state", k[["paren_poly"]] == 2L)
ok("invariant character counts 1", k[["invariant"]] == 1L)

## The old expression, for contrast: it is the inflation #104 describes.
old <- apply(mat, 2L, function(col) length(unique(col[!col %in% c("?", "-")])))
ok("old expression did inflate the polymorphic columns",
   old[["polymorph"]] == 3L && old[["paren_poly"]] == 3L)

## ---- .arm_own_files ----------------------------------------------------------
## The shared per-(tree, rep) directory as every arm leaves it.
siblings <- c(
  "mk_checkpoint.rds", "mk_run_1.log", "mk_run_2.log",
  "mk_trees_1.nwk", "mk_trees_2.nwk", "mk_trees.nwk", "mk_run.log",
  "mk_k9_checkpoint.rds", "mk_k9_run_1.log", "mk_k9_trees_1.nwk",
  "mk_k15_run_1.log", "mk_k24_run_1.log", "mk_k40_trees_2.nwk",
  "mk_kp1_checkpoint.rds", "mk_kp2_run_1.log", "mk_ktrue_trees_1.nwk",
  "mk_tlshrink_run_2.log",
  "mkp_checkpoint.rds", "mkp_run_1.log",
  "mkp_eg_run_1.log", "mkp_geo_checkpoint.rds", "mkp_highk_trees_2.nwk",
  "mkp_logs_run_1.log",
  ".slurm_job_id"
)
invisible(file.create(file.path(env$ckp_dir, siblings)))

own_mk  <- basename(env$.arm_own_files("mk"))
own_mkp <- basename(env$.arm_own_files("mkp"))

ok("mk claims its own seven files", setequal(own_mk, c(
  "mk_checkpoint.rds", "mk_run_1.log", "mk_run_2.log",
  "mk_trees_1.nwk", "mk_trees_2.nwk", "mk_trees.nwk", "mk_run.log")))
ok("mk claims nothing of mk_k9 / mk_kp1 / mk_ktrue / mk_tlshrink",
   !any(grepl("^mk_(k[0-9]|kp[0-9]|ktrue|tlshrink)", own_mk)))
ok("mkp claims its own two files",
   setequal(own_mkp, c("mkp_checkpoint.rds", "mkp_run_1.log")))
ok("mkp claims nothing of mkp_eg / mkp_geo / mkp_highk / mkp_logs",
   !any(grepl("^mkp_(eg|geo|highk|logs)", own_mkp)))
ok("the job sentinel is never purged",
   !(".slurm_job_id" %in% c(own_mk, own_mkp)))

## The old prefix glob, for contrast: this is the blast radius in #98.
old_mk <- grep("^mk_", siblings, value = TRUE)
ok("old prefix glob did reach ten foreign files",
   sum(!old_mk %in% own_mk) == 10L)

## ---- .UPostMeans -------------------------------------------------------------
## Two runs whose k' means differ sharply, so reading one run or skipping the
## burn-in both give visibly different answers from reading both correctly.
nChar <- 3L
mkd   <- list(nChar = nChar, kObs = rep(2L, nChar))

write_log <- function(path, n, kp_value) {
  hdr <- paste(c("Sample", "tree_length", paste0("kPrime_", seq_len(nChar))),
               collapse = "\t")
  rows <- vapply(seq_len(n), function(i) {
    ## First half of each run sits at 99 so a missing burn-in cut is obvious.
    kp <- if (i <= n / 2) 99 else kp_value
    paste(c(i, 1.0, rep(kp, nChar)), collapse = "\t")
  }, character(1L))
  writeLines(c(hdr, rows), path)
}

dir <- tempfile("upost"); dir.create(dir)
f1 <- file.path(dir, "a_run_1.log")
f2 <- file.path(dir, "a_run_2.log")
write_log(f1, 100L, 10)
write_log(f2, 100L, 20)

res <- list(logFile = c(f1, f2), samples = NULL)
up  <- env$.UPostMeans(res, mkd, burninFrac = 0.5)

ok("both runs read", up$nRuns == 2L)
ok("post-burn-in row count is 50 + 50", up$n == 100L)
ok("burn-in excludes the 99 block", all(abs(up$means - (15 - 2)) < 1e-9))
ok("burninFrac is reported back", identical(up$burninFrac, 0.5))

## Run 1 alone, no burn-in -- what the old code computed -- differs.
old_means <- colMeans(ReadMkLog(f1)[, paste0("kPrime_", seq_len(nChar)),
                                    drop = FALSE]) - mkd$kObs
ok("old run-1-only, no-burn-in estimate differs",
   all(abs(old_means - up$means) > 1))

## Degrades to res$samples when no log file is on disk.
res_nolog <- list(logFile = character(0L),
                  samples = ReadMkLog(f2))
ok("falls back to res$samples", env$.UPostMeans(res_nolog, mkd)$nRuns == 1L)

## No k' columns at all (the fixed-k arms): NA of the right length.
res_nokp <- list(logFile = character(0L),
                 samples = matrix(1, nrow = 4L, ncol = 1L,
                                  dimnames = list(NULL, "tree_length")))
ok("no kPrime_ columns gives NA per character",
   identical(env$.UPostMeans(res_nokp, mkd)$means, rep(NA_real_, nChar)))

ok("default burn-in fraction is 25%", identical(env$.BURNIN_FRAC, 0.25))

cat("\nALL CHECKS PASSED\n")
