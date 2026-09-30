#!/usr/bin/env Rscript
## End-to-end exercise of summarize_streamed.R against fabricated task
## directories shaped exactly like run_one.R / RunMkPrime's streaming output.
##
## Run from the repository root, with MkPrime installed (the summariser calls
## MkPrimeData()) alongside data.table, TreeDist, ape and TreeTools:
##   Rscript data-raw/hamilton/test_summarize_streamed.R
##
## This exists because #99 was a mismatch between what the summariser globbed
## for and what `.TreeFilePaths()` writes, and nothing anywhere compared the
## two. No cluster is needed.
##
## Checks the claims the fixes rest on:
##   1. the per-run tree glob matches the layout the package actually writes
##   2. every run's log is read, not just run 1
##   3. the #-comment block is excluded from n_samples
##   4. thinned_idx indexes thinned_trees after malformed lines are dropped
##   5. kPrime_* are summarised for an arm outside the old hard-coded list
##   6. kObs and char_idx follow <arm>_char_order.csv, cover the variable
##      characters only, and differ from what lexical order would give (#293)
##   7. a task without an up-to-date final .rds is refused unless --incomplete
##   8. two files for one run index are refused, not pooled (#289)
##   9. nothing is cleaned while a checkpoint remains (#289)

suppressPackageStartupMessages({library("ape"); library("TreeTools")})
set.seed(1)

script <- normalizePath("data-raw/hamilton/summarize_streamed.R", winslash = "/")
arm    <- "mkp_highk"
taxa   <- paste0("t", 1:8)
nChar  <- 11L                # >= 10, so lexical and numeric orders differ
invar  <- 5L                 # one invariant character, dropped by MkPrimeData
KObsOf <- function(i) if (i == invar) 1L else 2L + (i %% 5L)
varIdx <- setdiff(seq_len(nChar), invar)
kObsNumeric <- vapply(varIdx, KObsOf, integer(1L))
lexVar      <- as.integer(sub("chr", "", sort(paste0("chr", varIdx))))
kObsLexical <- vapply(lexVar, KObsOf, integer(1L))
n1 <- 300L; n2 <- 120L       # run 1 longer: the old pooled tail took none of it

ok <- function(label, cond) {
  cat(sprintf("%-62s %s\n", label, if (isTRUE(cond)) "PASS" else "*** FAIL"))
  if (!isTRUE(cond)) stop(label)
}
ok("fixture discriminates lexical from numeric kObs",
   !identical(kObsNumeric, kObsLexical))

## ---- Fixture ----------------------------------------------------------------
Fixture <- function(final = TRUE, checkpoint = FALSE) {
  root <- normalizePath(tempfile("harness"), winslash = "/", mustWork = FALSE)
  fx <- list(root = root,
             results = file.path(root, "results"),
             data = file.path(root, "data"),
             out = file.path(root, "summary"),
             task = file.path(root, "results", "t01_r01"))
  repDir <- file.path(fx$data, "tree_01/rep_01")
  dir.create(fx$task, recursive = TRUE)
  dir.create(repDir, recursive = TRUE)
  ape::write.tree(ape::rtree(8, tip.label = taxa, br = NULL),
                  file.path(fx$data, "tree_01/tree.nwk"))

  for (i in seq_len(nChar)) {
    states <- as.character(rep_len(seq_len(KObsOf(i)) - 1L, length(taxa)))
    ape::write.nexus.data(setNames(as.list(states), taxa),
                          file = file.path(repDir, sprintf("chr%d.nex", i)),
                          format = "standard")
  }
  write.csv(data.frame(position = seq_len(nChar),
                       file = sprintf("chr%d.nex", seq_len(nChar)),
                       char_idx = seq_len(nChar)),
            file.path(fx$task, sprintf("%s_char_order.csv", arm)),
            row.names = FALSE)

  MkLog(file.path(fx$task, sprintf("%s_run_1.log", arm)), n1, 0L)
  MkLog(file.path(fx$task, sprintf("%s_run_2.log", arm)), n2, 5L)
  MkTrees(file.path(fx$task, sprintf("%s_trees_1.nwk", arm)), n1, FALSE)
  MkTrees(file.path(fx$task, sprintf("%s_trees_2.nwk", arm)), n2, TRUE)
  if (checkpoint) {
    saveRDS(list(runs = list(1), mcmc = list(), iter = 1L),
            file.path(fx$task, sprintf("%s_checkpoint.rds", arm)))
  }
  fx$final <- file.path(fx$results, sprintf("%s_t01_r01.rds", arm))
  if (final) WriteFinal(fx)
  fx
}

WriteFinal <- function(fx, charIdx = varIdx) {
  saveRDS(list(stop_reason = "max_time", char_idx = charIdx), fx$final)
}

hdr <- paste(c("Sample", "log_posterior", "log_likelihood", "tree_length",
               "rate_log_sd", "p", paste0("br_", 1:13),
               paste0("kPrime_", seq_along(varIdx))), collapse = "\t")
MkLog <- function(path, n, offset) {
  rows <- vapply(seq_len(n), function(i) {
    paste(c(i * 10L + offset,
            round(rnorm(1, -60), 4), round(rnorm(1, -46), 4),
            round(runif(1, 1, 5), 4), round(runif(1, 0.5, 2), 4),
            round(runif(1, 0.2, 0.4), 4),
            round(runif(13, 0.01, 0.5), 4),
            kObsNumeric + sample(0:3, length(varIdx), replace = TRUE)),
          collapse = "\t")
  }, character(1L))
  ## The acceptance block: written once, directly under the header.
  comment <- c("# Topology: nni:13.8% gibbs_spr:5.4% tbr:5.4%",
               "# Branches: branch_lengths:7.1% tree_length:3.6%",
               "# Characters: mh_logit_p:10.7% gibbs_kPrime:7.1%",
               "# Rates: slice_rate_log_sd:5.4%")
  writeLines(c(hdr, comment, rows), path)
}

MkTrees <- function(path, n, truncateLast) {
  txt <- vapply(seq_len(n), function(i) {
    ape::write.tree(ape::rtree(8, tip.label = taxa))
  }, character(1L))
  if (truncateLast) txt[n] <- substr(txt[n], 1L, nchar(txt[n]) - 12L)
  writeLines(txt, path)
}

Summarise <- function(fx, ...) {
  out <- suppressWarnings(system2(
    "Rscript",
    c(shQuote(script), "1", "1", arm, shQuote(fx$results), shQuote(fx$data),
      shQuote(fx$out), "60", ...),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(out, "status")
  list(ok = is.null(status) || status == 0L, out = out,
       summary = file.path(fx$out, sprintf("%s_t01_r01.rds", arm)))
}

TaskFiles <- function(fx) sort(list.files(fx$task))

## ---- A finished task with no checkpoint -------------------------------------
fx <- Fixture()
rc <- Summarise(fx)
cat(paste(rc$out, collapse = "\n"), "\n")
ok("finished task summarised", rc$ok)
s <- readRDS(rc$summary)

ok("both tree files found", length(s$tree_files) == 2L)
ok("both logs read", length(s$log_files) == 2L)
ok("per-run raw counts are n1, n2",
   identical(as.integer(s$n_samples_raw_per_run), c(n1, n2)))
ok("comment block excluded from n_samples",
   s$n_samples == sum(s$n_samples_kept_per_run))
ok("burn-in dropped 25% of each run",
   identical(as.integer(s$n_samples_kept_per_run),
             as.integer(c(n1 - floor(n1 * 0.25), n2 - floor(n2 * 0.25)))))
ok("both runs contribute trees", identical(as.integer(s$n_trees_per_run),
                                           c(n1, n2)))
ok("truncated final Newick dropped", s$n_trees_kept < length(s$thinned_idx) + 1L)
ok("thinned_idx aligns with thinned_trees",
   length(s$thinned_idx) == length(s$thinned_trees))
ok("kept trees drawn from BOTH runs",
   any(s$thinned_idx <= n1) && any(s$thinned_idx > n1))
ok("kPrime_ columns summarised for mkp_highk",
   length(s$kp_means) == length(varIdx))
ok("p summarised for mkp_highk", !is.na(s$p_mean))
ok("burnin_frac recorded", identical(s$burnin_frac, 0.25))
ok("CID computed for every kept tree", length(s$cid) == s$n_trees_kept)
ok("kObs follows the recorded numeric order", identical(s$kObs, kObsNumeric))
ok("char_idx is the variable characters, numerically ordered",
   identical(s$char_idx, varIdx))
ok("char_idx pairs one-to-one with kObs", length(s$char_idx) == length(s$kObs))
ok("every posterior k' is at least its kObs", all(s$kp_means >= s$kObs))
ok("not flagged incomplete", identical(s$incomplete, FALSE))
ok("tree streams removed",
   !any(file.exists(file.path(fx$task, sprintf("%s_trees_%d.nwk", arm, 1:2)))))
ok("logs gzipped",
   all(file.exists(file.path(fx$task, sprintf("%s_run_%d.log.gz", arm, 1:2)))))

## A fresh chain beside the gzipped one: two files for run 1.
MkLog(file.path(fx$task, sprintf("%s_run_1.log", arm)), 50L, 0L)
MkTrees(file.path(fx$task, sprintf("%s_trees_1.nwk", arm)), 50L, FALSE)
MkTrees(file.path(fx$task, sprintf("%s_trees_2.nwk", arm)), 50L, FALSE)
WriteFinal(fx)
rc <- Summarise(fx)
ok(".log beside .log.gz for one run is refused",
   !rc$ok && any(grepl("More than one file for run 1", rc$out)))

## ---- A finished task whose checkpoint is kept -------------------------------
fx <- Fixture(checkpoint = TRUE)
before <- TaskFiles(fx)
rc <- Summarise(fx)
ok("task with checkpoint summarised", rc$ok && file.exists(rc$summary))
ok("nothing cleaned while the checkpoint remains",
   identical(TaskFiles(fx), before))

## ---- Tasks that may still be running ----------------------------------------
fx <- Fixture(final = FALSE)
before <- TaskFiles(fx)
rc <- Summarise(fx)
ok("no final .rds: refused", !rc$ok && !file.exists(rc$summary))
ok("refusal points at --incomplete", any(grepl("--incomplete", rc$out)))
ok("refusal touches nothing", identical(TaskFiles(fx), before))

rc <- Summarise(fx, "--incomplete")
ok("--incomplete summarises it", rc$ok && file.exists(rc$summary))
ok("and flags the summary incomplete", isTRUE(readRDS(rc$summary)$incomplete))

## A resubmitted task extending an earlier run: its old .rds predates the logs.
fx <- Fixture()
Sys.setFileTime(fx$final, Sys.time() - 3600)
before <- TaskFiles(fx)
rc <- Summarise(fx)
ok("final .rds older than the logs: refused",
   !rc$ok && identical(TaskFiles(fx), before))

## ---- Output from before the order record: lexical ---------------------------
fx <- Fixture()
invisible(file.remove(file.path(fx$task, sprintf("%s_char_order.csv", arm))))
WriteFinal(fx, charIdx = lexVar)
rc <- Summarise(fx)
s  <- readRDS(rc$summary)
ok("no order record: kObs and char_idx follow lexical order",
   rc$ok && identical(s$kObs, kObsLexical) && identical(s$char_idx, lexVar))

## ---- A final .rds recorded under another character order --------------------
fx <- Fixture()
WriteFinal(fx, charIdx = lexVar)
rc <- Summarise(fx)
ok("char_idx disagreeing with run_one.R's record is refused",
   !rc$ok && any(grepl("char_idx", rc$out)))

cat("\nALL CHECKS PASSED\n")
