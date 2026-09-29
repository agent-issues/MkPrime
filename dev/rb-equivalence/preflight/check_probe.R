#!/usr/bin/env Rscript
# Compare MkPrime's likelihood with RevBayes' at the probe's fixed states.
#
# Usage:
#   Rscript check_probe.R <pid> <model> [--probe-dir=probe]
#     [--local-neotrans=../../../neotrans] [--tol=1e-4]
#
# Columns, in nats (MkPrime minus RevBayes, per state):
#   asRun    coding "variable" vs RB as the harness runs it (rescaled neo Q,
#            polymorphisms as partial ambiguity);
#   depoly   the same with polymorphic cells as '?' on the RB side, which is
#            how MkPrime reads them;
#   noAsc    coding "none" vs RB coding "all", polymorphisms as '?'.
# Exits non-zero if max |depoly| or |noAsc| exceeds --tol.
#
# by_nt_kv uses RevBayes' k (max state index + 1, polymorphisms included),
# not the harness's kObs; the two differ where a matrix has polymorphisms.

suppressPackageStartupMessages({
  library("TreeTools")
  library("MkPrime")
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: check_probe.R <pid> <model> [opts]")
pid <- args[[1]]
model <- args[[2]]
opt <- list(probe_dir = "probe", local_neotrans = "../../../neotrans", tol = "1e-4")
for (a in args[-(1:2)]) {
  kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  if (length(kv) != 2L) stop("Bad option: ", a)
  opt[[gsub("-", "_", kv[[1]])]] <- kv[[2]]
}
tol <- as.numeric(opt$tol)
scriptDir <- (function() {
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1]]))) else getwd()
})()
source(file.path(scriptDir, "..", "R", "utils.R"))
proj <- file.path(opt$local_neotrans, "inst", "projects")
cell <- sprintf("%s_%s", pid, model)

transMat <- ReadSplit(file.path(proj, sprintf("project%s.trans.nex", pid)))
comb <- CombineSplits(transMat, ReadSplit(file.path(proj, sprintf("project%s.neo.nex", pid))))
rbK <- apply(transMat, 2, function(col) {
  digits <- unlist(regmatches(col, gregexpr("[0-9]", col)))
  max(as.integer(digits)) + 1L
})
info <- list(nTrans = comb$nTrans, nNeo = comb$nNeo, kObs = rbK)
knownStates <- BuildKnownStates(info, model)
mkd <- suppressWarnings(MkPrimeData(MatrixToCombinedPhyDat(comb$matrix),
                                    neomorphic = BuildNeoIdx(info, model),
                                    knownStates = knownStates))

smp <- readRDS(file.path(opt$probe_dir, sprintf("sample_%s.rds", cell)))
lines <- readLines(file.path(opt$probe_dir, sprintf("out_%s.txt", cell)))
header <- strsplit(sub(".*HEADER ", "", grep("HEADER", lines, value = TRUE)[[1]]), " ")[[1]]
rb <- read.table(text = grep("^ROW", lines, value = TRUE),
                 col.names = c("tag", sub('"$', "", header)))
if (nrow(rb) != nrow(smp$par)) stop("Probe output has ", nrow(rb), " rows, expected ", nrow(smp$par))

res <- do.call(rbind, lapply(seq_len(nrow(rb)), function(i) {
  p <- smp$par[i, ]
  stopifnot(p$Iteration == rb$iter[[i]])
  # MkPrime pairs tips with data rows by position (#224)
  tr <- RenumberTips(smp$trees[[i]], mkd$taxon_names)
  lnL <- function(coding) {
    MkpLogLikelihood(tr, mkd, rate_loss = p$rate_loss, rate_log_sd = p$rate_log_sd,
                     nCat = 6L, coding = coding, rate_neo = p$rate_neo)
  }
  mkVar <- lnL("variable")
  data.frame(iter = p$Iteration, rateLoss = p$rate_loss,
             rbLoggedCheck = p$Likelihood - (rb$nTV[[i]] + rb$tV[[i]]),
             asRun = mkVar - (rb$nTV[[i]] + rb$tV[[i]]),
             depoly = mkVar - (rb$nTV_D[[i]] + rb$tV_D[[i]]),
             noAsc = lnL("none") - (rb$nTA_D[[i]] + rb$tA_D[[i]]))
}))
print(res, digits = 4, row.names = FALSE)
worst <- max(abs(c(res$depoly, res$noAsc)))
cat(sprintf("\n[check_probe] %s: max |depoly, noAsc| = %.2e nats (tol %.0e); max |asRun| = %.2e; probe reproduces RB's logged lnL to %.1e\n",
            cell, worst, tol, max(abs(res$asRun)), max(abs(res$rbLoggedCheck))))
if (worst > tol) {
  cat("[check_probe] FAIL\n")
  quit(status = 1L)
}
cat("[check_probe] PASS\n")
