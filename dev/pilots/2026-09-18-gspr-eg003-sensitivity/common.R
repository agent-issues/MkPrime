# Shared setup for the gibbs_spr / EG-003 sensitivity check.
.libPaths(c("C:/Users/pjjg18/GitHub/.builds/MkPrime-G3", .libPaths()))
suppressPackageStartupMessages({
  library(MkPrime.G3)
  library(ape); library(TreeTools)
})
# No namespace shim needed: nothing here calls MkPrime::: directly.
source(file.path(Sys.getenv("GSPR_DIR", "."), "align.R"))

# Which build produced a result: package name, commit (RemoteSha, stamped by
# data-raw/hamilton/install_mkp.sh; NA for an unstamped build) and library path.
BuildInfo <- function(pkg = "MkPrime.G3") {
  d <- utils::packageDescription(pkg)
  list(package = pkg,
       remoteSha = if (is.null(d$RemoteSha)) NA_character_ else d$RemoteSha,
       version = d$Version,
       libPath = dirname(find.package(pkg)))
}

GT_ROOT <- "C:/Users/pjjg18/GitHub/mkprime/tree-inference"

# Replicates data-raw/hamilton/run_one.R: numeric file order, cbind,
# write.nexus.data, ReadAsPhyDat.
LoadRep <- function(tree_idx, rep_idx = 1L) {
  d <- file.path(GT_ROOT, sprintf("tree_%02d/rep_%02d", tree_idx, rep_idx))
  nex <- list.files(d, pattern = "^chr[0-9]+\\.nex$", full.names = TRUE)
  nex <- nex[order(as.integer(sub("^chr([0-9]+)\\.nex$", "\\1", basename(nex))))]
  mats <- lapply(nex, TreeTools::ReadCharacters)
  m <- do.call(cbind, mats)
  tmp <- tempfile(fileext = ".nex")
  dl <- setNames(lapply(seq_len(nrow(m)), function(i) m[i, ]), rownames(m))
  ape::write.nexus.data(dl, file = tmp, format = "standard")
  pd <- TreeTools::ReadAsPhyDat(tmp)
  file.remove(tmp)
  gt <- read.csv(file.path(d, "ground_truth.csv"))
  charIdx <- CharIdxFromFiles(nex)
  list(pd = pd, gt = gt, gtAligned = AlignGroundTruth(gt, charIdx),
       charIdx = charIdx, mat = m, nexOrder = basename(nex),
       nChar = ncol(m), nTax = nrow(m))
}
