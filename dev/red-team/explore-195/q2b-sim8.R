# Q2b: the assessment phase on a less trivial concentrated posterior, the
# clean 8-tip Mk(3) simulation from Q1 (one topology in 5000 samples), at
# two per-chain budgets: does too short a budget produce false "stuck"?
# Run: timeout 900 Rscript dev/red-team/explore-195/q2b-sim8.R
setTimeLimit(elapsed = 880, transient = FALSE)
source("dev/red-team/explore-195/common.R")

set.seed(195)
simTree <- RandomTree(8, root = TRUE)
simTree$edge.length <- rep(0.15, nrow(simTree$edge))
sim <- phangorn::simSeq(simTree, l = 200, type = "USER",
                        levels = c("0", "1", "2"), rate = 1)

mainIter <- 6000
rows <- list()
for (s in 1:2) {
  main <- ShortRun(sim, AdditionTree(sim), nIter = mainIter, minWarmup = 1000,
                   maxTime = 120, seed = s)
  stuckTree <- main$post$trees[[length(main$post$trees)]]
  for (budget in c(600, 1500)) {
    a <- AssessPhase(sim, stuckTree, k = 3, m = 3, budget = budget,
                     minWarmup = 300, seed = s, maxTime = 60)
    shares <- a$perChain$share
    row <- data.frame(
      seed = s, budget = budget,
      mainDistinct = TopoSummary(main$post$trees)$distinct,
      lenient = a$lenient, strict = a$strict,
      nniShare = paste(round(shares[a$perChain$kind == "nni"], 2),
                       collapse = "/"),
      randShare = paste(round(shares[a$perChain$kind == "random"], 2),
                        collapse = "/"),
      costIter = a$iter / mainIter, costWall = round(a$wall / main$wall, 2)
    )
    print(row, row.names = FALSE)
    rows[[length(rows) + 1]] <- row
  }
}
out <- do.call(rbind, rows)
print(out, row.names = FALSE)
write.csv(out, "dev/red-team/explore-195/q2b-results.csv", row.names = FALSE)
