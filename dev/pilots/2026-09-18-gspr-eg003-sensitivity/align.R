# Character-order and grid helpers. No package dependency, so test-align.R can source this alone.

CharIdxFromFiles <- function(files) {
  as.integer(sub("^chr([0-9]+)\\.nex$", "\\1", basename(files)))
}

# Rows of ground_truth.csv in the sampler's character order, joined on char_idx.
AlignGroundTruth <- function(gt, charIdx) {
  stopifnot(!anyDuplicated(gt$char_idx), !anyDuplicated(charIdx),
            all(charIdx %in% gt$char_idx))
  gt[match(charIdx, gt$char_idx), , drop = FALSE]
}

# Seeds per (tree, arm) when the on/off grid is complete; stops, naming the
# missing cells, otherwise. A cell killed by its time limit writes no .rds.
CheckCompleteGrid <- function(d, treeIdx, arms = c("on", "off")) {
  counts <- table(factor(d$tree, levels = treeIdx), factor(d$arm, levels = arms))
  nSeeds <- max(counts)
  short <- which(counts != nSeeds, arr.ind = TRUE)
  if (nSeeds == 0L || nrow(short)) {
    stop("Incomplete grid (expected ", nSeeds, " seeds per tree and arm): ",
         paste(sprintf("tree %s %s has %d", rownames(counts)[short[, 1]],
                       colnames(counts)[short[, 2]], counts[short]),
               collapse = "; "),
         ". Rerun 02b to fill the gaps.", call. = FALSE)
  }
  nSeeds
}
