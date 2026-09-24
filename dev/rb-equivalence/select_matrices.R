#!/usr/bin/env Rscript
# Build the frozen matrix selection for the RB-oracle benchmark.
#
# Picks the five smallest matrices from the `fast` set in
# `../neotrans/inst/1_PrepareMatrix.R` by (nChar + nTaxa), with the
# constraint that one matrix has nChar >= 30 and nChar >= nTaxa to
# stress the harness on a denser likelihood surface.
#
# Output: dev/rb-equivalence/matrices.rds, a list keyed by pID with:
#   $pid     character  e.g. "950"
#   $nChar   integer    from metadata
#   $nTaxa   integer    from metadata
#   $neoIdx  integer    char indices in the *combined* file
#   $transIdx integer   "
#   $kObs    integer    observed-state count for each trans char (order matches transIdx)
#
# Selection is hard-coded after a one-off audit; the script verifies that the
# audit still holds against current metadata.csv + fast list.
#
# NOTE (agent-issues/MkPrime#214): $kObs here is the distinct-non-ambiguous-
# token count, NOT RevBayes' state-count rule (max observed state index + 1,
# polymorphic bits included). render_rev.R and run_mkprime.R no longer read
# this rds for k/nChar/taxa -- they derive those from R/utils.R's
# PrepareCellInfo() directly from the split files, which both samplers are
# guaranteed to agree on. This script's output remains useful only for the
# PIDS selection audit above.

suppressPackageStartupMessages({
  library(TreeTools)
})

NEOTRANS <- normalizePath("../neotrans", mustWork = TRUE)

PIDS <- c("950", "6072", "635", "3832", "3408")

# --- Verify selection is still optimal against current metadata --------------

metaPath <- file.path(NEOTRANS, "inst", "matrices", "metadata.csv")
meta <- read.csv(metaPath, check.names = FALSE, fileEncoding = "UTF-8-BOM")
names(meta)[1] <- "project"  # strip any BOM artefact

prepPath <- file.path(NEOTRANS, "inst", "1_PrepareMatrix.R")
prepSrc <- readLines(prepPath)
fastLine <- grep("^\\s*fast\\s*<-", prepSrc)
if (!length(fastLine)) stop("Couldn't find `fast` definition in ", prepPath)
# Read continuation lines until the closing paren
fastTxt <- paste(prepSrc[fastLine:min(length(prepSrc), fastLine + 5)],
                 collapse = " ")
fastIds <- as.character(eval(parse(text = sub(".*fast\\s*<-\\s*", "",
                                              sub("\\).*", ")", fastTxt)))))

# Sanity: all selected pIDs are in the fast list, splits exist, and at least
# one pID meets the dense-likelihood constraint (>=30 chars, chars >= taxa).
metaSel <- meta[as.character(meta$project) %in% PIDS, ]
if (nrow(metaSel) != length(PIDS)) {
  stop("Selection includes pID(s) not in metadata.csv: ",
       paste(setdiff(PIDS, as.character(meta$project)), collapse = ", "))
}
notFast <- setdiff(PIDS, fastIds)
if (length(notFast)) {
  stop("Selection includes pID(s) not in fast list: ",
       paste(notFast, collapse = ", "))
}
denseOk <- any(metaSel$nChar >= 30 & metaSel$nChar >= metaSel$nTaxa)
if (!denseOk) {
  stop("Selection lacks any matrix with nChar >= 30 and nChar >= nTaxa")
}

# --- Read each split + record neo/trans indices in the combined file ---------

readPhy <- function(path) {
  stopifnot(file.exists(path))
  TreeTools::ReadCharacters(path)
}

build <- function(pid) {
  trans <- readPhy(file.path(NEOTRANS, "inst", "projects",
                             sprintf("project%s.trans.nex", pid)))
  neo <- readPhy(file.path(NEOTRANS, "inst", "projects",
                           sprintf("project%s.neo.nex", pid)))

  nTrans <- ncol(trans)
  nNeo <- ncol(neo)

  # Convention used downstream: trans chars come FIRST in the combined matrix
  transIdx <- seq_len(nTrans)
  neoIdx <- nTrans + seq_len(nNeo)

  # kObs per trans char — count distinct non-ambiguous tokens per column
  kObs <- apply(trans, 2, function(col) {
    states <- unique(col[!is.na(col) & col != "?" & col != "-"])
    length(states)
  })

  mrow <- meta[as.character(meta$project) == pid, ]
  list(
    pid = pid,
    nChar = as.integer(mrow$nChar),
    nTaxa = as.integer(mrow$nTaxa),
    nTrans = nTrans,
    nNeo = nNeo,
    transIdx = transIdx,
    neoIdx = neoIdx,
    kObs = as.integer(kObs),
    transTipLabels = rownames(trans),
    neoTipLabels = rownames(neo)
  )
}

matrices <- lapply(PIDS, build)
names(matrices) <- PIDS

# --- Sanity: taxon sets should match between neo and trans ------------------

for (m in matrices) {
  if (!setequal(m$transTipLabels, m$neoTipLabels)) {
    warning(sprintf("pid %s: trans (%d tips) and neo (%d tips) taxon sets differ",
                    m$pid, length(m$transTipLabels), length(m$neoTipLabels)))
  }
}

# --- Summary print ----------------------------------------------------------

summary_df <- do.call(rbind, lapply(matrices, function(m) {
  data.frame(pid = m$pid, nTaxa = length(m$transTipLabels),
             nTrans = m$nTrans, nNeo = m$nNeo,
             maxK = max(m$kObs), meanK = round(mean(m$kObs), 2),
             stringsAsFactors = FALSE)
}))
cat("=== RB-oracle matrix selection ===\n")
print(summary_df, row.names = FALSE)

outPath <- file.path("dev", "rb-equivalence", "matrices.rds")
saveRDS(matrices, outPath)
cat(sprintf("\nSaved: %s\n", outPath))
