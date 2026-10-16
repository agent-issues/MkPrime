# Q1: do real datasets produce single-topology posteriors?
# Short guarded single runs; counts distinct (unrooted) topologies in the
# post-warmup tree sample. Run:
#   timeout 900 Rscript dev/red-team/explore-195/q1-topology-counts.R
setTimeLimit(elapsed = 880, transient = FALSE)
source("dev/red-team/explore-195/common.R")

Subset <- function(dat, n) {
  # First n taxa, drop characters left uninformative
  dat <- dat[seq_len(n)]
  mat <- PhyDatToMatrix(dat)
  MatrixToPhyDat(mat)
}

sets <- list(
  "6 tips, 160 congruent (sim)" = CongruentData()$data,
  "Vinther2008 (23 x 50)" = inapplicable.phyData[["Vinther2008"]],
  "Longrich2010 (20 x 80)" = inapplicable.phyData[["Longrich2010"]],
  "Sansom2010 (23 x 97)" = inapplicable.phyData[["Sansom2010"]],
  "Longrich2010, first 8 taxa" = Subset(inapplicable.phyData[["Longrich2010"]], 8),
  "Zhu2013, first 8 taxa" = Subset(inapplicable.phyData[["Zhu2013"]], 8)
)

# Clean simulated 8-tip data under Mk(3), moderate branch lengths
set.seed(195)
simTree <- RandomTree(8, root = TRUE)
simTree$edge.length <- rep(0.15, nrow(simTree$edge))
sim <- phangorn::simSeq(simTree, l = 200, type = "USER",
                        levels = c("0", "1", "2"), rate = 1)
sets[["8 tips, 200 chars Mk(3), br 0.15 (sim)"]] <- sim

out <- lapply(names(sets), function(nm) {
  dat <- sets[[nm]]
  res <- tryCatch({
    r <- ShortRun(dat, AdditionTree(dat), nIter = 6000, minWarmup = 1000,
                  maxTime = 90, seed = 1)
    s <- TopoSummary(r$post$trees)
    data.frame(dataset = nm, tips = length(dat), chars = sum(attr(dat, "weight")),
               trees = s$n, distinct = s$distinct, mapFreq = round(s$mapFreq, 3),
               changes = s$changes, wall = round(r$wall, 1))
  }, error = function(e) {
    data.frame(dataset = nm, tips = length(dat), chars = NA, trees = NA,
               distinct = NA, mapFreq = NA, changes = NA, wall = NA,
               error = conditionMessage(e))
  })
  print(res)
  res
})
out <- do.call(rbind, lapply(out, function(x) x[, 1:8]))
print(out)
write.csv(out, "dev/red-team/explore-195/q1-results.csv", row.names = FALSE)
