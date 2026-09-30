#!/usr/bin/env Rscript
# Did any Mk' task resume across the dee5668 character-order change (#286)?
#
# Paired correctly, no logged k' can fall below its character's kObs: the prior
# is -Inf there. So a row with kPrime_j < kObs_j, where kObs is taken in the
# numeric order run_one.R now uses, was written under the other (lexical)
# order. One row is proof; none means the numeric pairing held throughout.
#
# Usage:
#   Rscript check_kprime_order.R <results_root> <data_root> [arm ...]
#
# Writes kprime_order_check.csv to the working directory: one line per log
# file, with the rows read and the rows that violate k' >= kObs.

.libPaths(c("/nobackup/pjjg18/mkp-study/lib", .libPaths()))
suppressPackageStartupMessages({
  library("data.table")
  library("MkPrime")
  library("TreeTools")
})

args         <- commandArgs(trailingOnly = TRUE)
results_root <- args[1]
data_root    <- args[2]
arms         <- if (length(args) > 2L) args[-(1:2)] else
  c("mkp", "mkp_eg", "mkp_geo", "mkp_highk", "mkp_logs")

.NumericKObs <- function(repDir) {
  files <- list.files(repDir, "^chr[0-9]+\\.nex$", full.names = TRUE)
  files <- files[order(as.integer(sub("^chr([0-9]+)\\.nex$", "\\1",
                                      basename(files))))]
  mat <- do.call(cbind, lapply(files, ReadCharacters))
  tmp <- tempfile(fileext = ".nex")
  on.exit(unlink(tmp))
  ape::write.nexus.data(
    setNames(lapply(seq_len(nrow(mat)), function(i) mat[i, ]), rownames(mat)),
    file = tmp, format = "standard"
  )
  as.integer(suppressWarnings(MkPrimeData(ReadAsPhyDat(tmp)))$kObs)
}

tags <- list.files(results_root, "^t[0-9]{2}_r[0-9]{2}$")
rows <- list()
for (tag in tags) {
  taskDir <- file.path(results_root, tag)
  repDir  <- file.path(data_root,
                       sub("^t([0-9]{2})_r([0-9]{2})$", "tree_\\1/rep_\\2", tag))
  kObs    <- NULL
  for (arm in arms) {
    logs <- list.files(taskDir,
                       sprintf("^%s_run(_[0-9]+)?\\.log(\\.gz)?$", arm),
                       full.names = TRUE)
    for (f in logs) {
      if (is.null(kObs)) kObs <- .NumericKObs(repDir)
      d  <- fread(f, sep = "\t", fill = Inf, showProgress = FALSE)
      d  <- d[!grepl("^#", as.character(d[[1L]]))]
      kp <- grep("^kPrime_", names(d), value = TRUE)
      bad <- if (length(kp) == length(kObs)) {
        m <- as.matrix(d[, lapply(.SD, as.numeric), .SDcols = kp])
        sum(rowSums(sweep(m, 2L, kObs) < 0, na.rm = TRUE) > 0)
      } else {
        NA_integer_
      }
      rows[[length(rows) + 1L]] <- data.frame(
        tag = tag, arm = arm, log = basename(f), n_rows = nrow(d),
        n_kprime = length(kp), n_kobs = length(kObs), n_bad_rows = bad
      )
    }
  }
}

out <- do.call(rbind, rows)
write.csv(out, "kprime_order_check.csv", row.names = FALSE)
cat(sprintf(
  "%d log(s) checked; %d with a row below kObs; %d unpaired (length mismatch)\n",
  nrow(out), sum(out$n_bad_rows > 0L, na.rm = TRUE), sum(is.na(out$n_bad_rows))
))
