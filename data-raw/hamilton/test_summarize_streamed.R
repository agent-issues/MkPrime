#!/usr/bin/env Rscript
## End-to-end exercise of summarize_streamed.R against a fabricated task
## directory shaped exactly like RunMkPrime's streaming output.
##
## Run from the repository root:
##   Rscript data-raw/hamilton/test_summarize_streamed.R
##
## This exists because #99 was a mismatch between what the summariser globbed
## for and what `.TreeFilePaths()` writes, and nothing anywhere compared the
## two. It needs no cluster and no MkPrime build -- only ape and TreeTools.
##
## Checks the claims the fix rests on:
##   1. the per-run tree glob matches the layout the package actually writes
##   2. every run's log is read, not just run 1
##   3. the #-comment block is excluded from n_samples
##   4. thinned_idx indexes thinned_trees after malformed lines are dropped
##   5. kPrime_* are summarised for an arm outside the old hard-coded list

suppressPackageStartupMessages({library(ape); library(TreeTools)})
set.seed(1)

root <- normalizePath(tempfile("harness"), winslash = "/", mustWork = FALSE)
dir.create(root, recursive = TRUE)
results_root <- file.path(root, "results")
data_root    <- file.path(root, "data")
out_dir      <- file.path(root, "summary")
task_dir     <- file.path(results_root, "t01_r01")
rep_dir      <- file.path(data_root, "tree_01/rep_01")
dir.create(task_dir, recursive = TRUE)
dir.create(rep_dir,  recursive = TRUE)

taxa <- paste0("t", 1:8)
true_tree <- ape::rtree(8, tip.label = taxa, br = NULL)
ape::write.tree(true_tree, file.path(data_root, "tree_01/tree.nwk"))

nChar <- 6L
for (i in seq_len(nChar)) {
  states <- as.character(c(0, 1, 0, 1, 0, 1, 0, 1))
  ape::write.nexus.data(setNames(as.list(states), taxa),
                        file = file.path(rep_dir, sprintf("chr%d.nex", i)),
                        format = "standard")
}

## ---- Fabricate two runs' streamed output -----------------------------------
brs <- paste0("br_", 1:13)
kps <- paste0("kPrime_", seq_len(nChar))
hdr <- paste(c("Sample", "log_posterior", "log_likelihood", "tree_length",
               "rate_log_sd", "p", brs, kps), collapse = "\t")

## Run 1 deliberately much longer than run 2: under the OLD pooled tail window
## run 1 would contribute nothing at all.
n1 <- 300L; n2 <- 120L
mk_log <- function(path, n, offset) {
  rows <- vapply(seq_len(n), function(i) {
    paste(c(i * 10L + offset,
            round(rnorm(1, -60), 4), round(rnorm(1, -46), 4),
            round(runif(1, 1, 5), 4), round(runif(1, 0.5, 2), 4),
            round(runif(1, 0.2, 0.4), 4),
            round(runif(13, 0.01, 0.5), 4),
            sample(2:6, nChar, replace = TRUE)),
          collapse = "\t")
  }, character(1L))
  ## The acceptance block: written once, directly under the header.
  comment <- c("# Topology: nni:13.8% gibbs_spr:5.4% tbr:5.4%",
               "# Branches: branch_lengths:7.1% tree_length:3.6%",
               "# Characters: mh_logit_p:10.7% gibbs_kPrime:7.1%",
               "# Rates: slice_rate_log_sd:5.4%")
  writeLines(c(hdr, comment, rows), path)
}
mk_log(file.path(task_dir, "mkp_highk_run_1.log"), n1, 0L)
mk_log(file.path(task_dir, "mkp_highk_run_2.log"), n2, 5L)

mk_trees <- function(path, n, truncate_last) {
  tt <- lapply(seq_len(n), function(i) ape::rtree(8, tip.label = taxa))
  txt <- vapply(tt, function(t) ape::write.tree(t), character(1L))
  if (truncate_last) txt[n] <- substr(txt[n], 1L, nchar(txt[n]) - 12L)
  writeLines(txt, path)
}
mk_trees(file.path(task_dir, "mkp_highk_trees_1.nwk"), n1, FALSE)
mk_trees(file.path(task_dir, "mkp_highk_trees_2.nwk"), n2, TRUE)  # SIGKILL tail

## ---- Run the summariser ------------------------------------------------------
script <- normalizePath("data-raw/hamilton/summarize_streamed.R", winslash = "/")
rc <- system2("Rscript",
              c(shQuote(script), "1", "1", "mkp_highk",
                shQuote(results_root), shQuote(data_root), shQuote(out_dir),
                "60"),
              stdout = TRUE, stderr = TRUE)
cat(paste(rc, collapse = "\n"), "\n")
stopifnot(!inherits(attr(rc, "status"), "integer") || attr(rc, "status") == 0L)

s <- readRDS(file.path(out_dir, "mkp_highk_t01_r01.rds"))

ok <- function(label, cond) {
  cat(sprintf("%-58s %s\n", label, if (isTRUE(cond)) "PASS" else "*** FAIL"))
  if (!isTRUE(cond)) stop(label)
}

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
ok("kPrime_ columns summarised for mkp_highk", length(s$kp_means) == nChar)
ok("p summarised for mkp_highk", !is.na(s$p_mean))
ok("burnin_frac recorded", identical(s$burnin_frac, 0.25))
ok("CID computed for every kept tree", length(s$cid) == s$n_trees_kept)

## ---- Cleanup behaviour -------------------------------------------------------
ok("tree streams removed",
   !any(file.exists(file.path(task_dir,
                              c("mkp_highk_trees_1.nwk", "mkp_highk_trees_2.nwk")))))
ok("logs retained (gzipped or plain)",
   all(file.exists(file.path(task_dir, "mkp_highk_run_1.log")) |
       file.exists(file.path(task_dir, "mkp_highk_run_1.log.gz"))))

cat("\nALL CHECKS PASSED\n")
