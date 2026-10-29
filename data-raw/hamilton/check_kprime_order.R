#!/usr/bin/env Rscript
# Did any Mk' task resume across the dee5668 character-order change (#286)?
#
# Paired correctly, no logged k' can fall below its character's kObs: the prior
# is -Inf there. kObs is taken in the order the summariser will pair the log
# with -- `<arm>_char_order.csv`, else lexical -- so a row with
# kPrime_j < kObs_j was written under another order. One row is proof. None is
# not: a permutation among characters of equal kObs cannot show up this way.
#
# Usage:
#   Rscript check_kprime_order.R <results_root> <data_root> [arm ...]
#
# Writes kprime_order_check.csv to the working directory: one line per log
# file, with the order used, the rows read and the rows that violate
# k' >= kObs (NA with an `error` if the log could not be read).

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

.KObsInOrder <- function(repDir, orderFile) {
  files <- if (file.exists(orderFile)) {
    file.path(repDir, read.csv(orderFile)$file)
  } else {
    sort(list.files(repDir, "^chr[0-9]+\\.nex$", full.names = TRUE))
  }
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
  for (arm in arms) {
    orderFile <- file.path(taskDir, sprintf("%s_char_order.csv", arm))
    logs <- list.files(taskDir,
                       sprintf("^%s_run(_[0-9]+)?\\.log(\\.gz)?$", arm),
                       full.names = TRUE)
    if (length(logs)) kObs <- .KObsInOrder(repDir, orderFile)
    for (f in logs) {
      d <- tryCatch(fread(f, sep = "\t", fill = Inf, showProgress = FALSE),
                    error = function(e) conditionMessage(e))
      if (is.character(d)) {
        rows[[length(rows) + 1L]] <- data.frame(
          tag = tag, arm = arm, log = basename(f), order = NA, n_rows = NA,
          n_kprime = NA, n_kobs = length(kObs), n_bad_rows = NA, error = d
        )
        next
      }
      d  <- d[!grepl("^#", as.character(d[[1L]]))]
      kp <- grep("^kPrime_", names(d), value = TRUE)
      bad <- if (length(kp) == length(kObs)) {
        m <- as.matrix(d[, lapply(.SD, as.numeric), .SDcols = kp])
        sum(rowSums(sweep(m, 2L, kObs) < 0, na.rm = TRUE) > 0)
      } else {
        NA_integer_
      }
      rows[[length(rows) + 1L]] <- data.frame(
        tag = tag, arm = arm, log = basename(f),
        order = if (file.exists(orderFile)) "recorded" else "lexical",
        n_rows = nrow(d), n_kprime = length(kp), n_kobs = length(kObs),
        n_bad_rows = bad, error = NA
      )
    }
  }
}

out <- do.call(rbind, rows)
write.csv(out, "kprime_order_check.csv", row.names = FALSE)
cat(sprintf(
  "%d log(s) checked; %d with a row below kObs; %d unreadable or unpaired\n",
  nrow(out), sum(out$n_bad_rows > 0L, na.rm = TRUE), sum(is.na(out$n_bad_rows))
))
# Exit status is what run_one.R's MKP_ACCEPT_UNSTAMPED_ORDER=1 relies on: only
# a clean check licenses resuming a task whose order record has no birth stamp.
if (any(out$n_bad_rows > 0L, na.rm = TRUE) || anyNA(out$n_bad_rows)) {
  quit(status = 1L)
}
