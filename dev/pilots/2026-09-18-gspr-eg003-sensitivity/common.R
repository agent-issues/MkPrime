# Shared setup for the gibbs_spr / EG-003 sensitivity check.
.libPaths(c("C:/Users/pjjg18/GitHub/.builds/MkPrime-G3", .libPaths()))
suppressPackageStartupMessages({
  library(MkPrime.G3)
  library(ape); library(TreeTools)
})
# No namespace shim needed: nothing here calls MkPrime::: directly.

GT_ROOT <- "C:/Users/pjjg18/GitHub/mkprime/tree-inference"

# Replicates data-raw/hamilton/run_one.R exactly: sort(list.files(...)),
# cbind, write.nexus.data, ReadAsPhyDat.
LoadRep <- function(tree_idx, rep_idx = 1L) {
  d <- file.path(GT_ROOT, sprintf("tree_%02d/rep_%02d", tree_idx, rep_idx))
  nex <- sort(list.files(d, pattern = "^chr[0-9]+\\.nex$", full.names = TRUE))
  mats <- lapply(nex, TreeTools::ReadCharacters)
  m <- do.call(cbind, mats)
  tmp <- tempfile(fileext = ".nex")
  dl <- setNames(lapply(seq_len(nrow(m)), function(i) m[i, ]), rownames(m))
  ape::write.nexus.data(dl, file = tmp, format = "standard")
  pd <- TreeTools::ReadAsPhyDat(tmp)
  file.remove(tmp)
  gt <- read.csv(file.path(d, "ground_truth.csv"))
  list(pd = pd, gt = gt, mat = m, nexOrder = basename(nex),
       nChar = ncol(m), nTax = nrow(m))
}
