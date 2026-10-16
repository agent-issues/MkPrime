# Q2: prototype of Martin's assessment phase (#195).
# From a chain's final topology, start k chains at NNI neighbours and m at
# random trees; "favoured" if every chain settles back on that topology.
# Run: timeout 1500 Rscript dev/red-team/explore-195/q2-assessment-phase.R
setTimeLimit(elapsed = 1480, transient = FALSE)
source("dev/red-team/explore-195/common.R")

args <- commandArgs(trailingOnly = TRUE)
hardName <- if (length(args)) args[1] else "Longrich2010"
mainIter <- 6000
budget <- 1500
seeds <- 1:3

cong <- CongruentData()$data
hard <- inapplicable.phyData[[hardName]]

cases <- list(
  # truth = what a correct classifier should say
  list(name = "(a) 6-tip congruent, normal moves", data = cong,
       start = function(s) AdditionTree(cong), weights = NULL,
       truth = "favoured"),
  list(name = paste0("(b1) ", hardName, ", topology moves ~0, random start"),
       data = hard,
       start = function(s) { set.seed(s); RandomTree(hard, root = TRUE) },
       weights = WeakTopology, truth = "stuck"),
  list(name = paste0("(b2) ", hardName, ", topology moves ~0, parsimony start"),
       data = hard, start = function(s) AdditionTree(hard),
       weights = WeakTopology, truth = "stuck"),
  list(name = paste0("(c) ", hardName, ", normal moves"), data = hard,
       start = function(s) AdditionTree(hard), weights = NULL,
       truth = "mixing")
)

rows <- list()
for (case in cases) for (s in seeds) {
  main <- ShortRun(case$data, case$start(s), nIter = mainIter,
                   minWarmup = 1000, maxTime = 120,
                   moveWeights = case$weights, seed = s)
  ms <- TopoSummary(main$post$trees)
  stuckTree <- main$post$trees[[length(main$post$trees)]]
  a <- AssessPhase(case$data, stuckTree, k = 3, m = 3, budget = budget,
                   minWarmup = 300, seed = s, maxTime = 60)
  shares <- a$perChain$share
  row <- data.frame(
    case = case$name, seed = s, truth = case$truth,
    mainDistinct = ms$distinct, mainWall = round(main$wall, 1),
    lenient = a$lenient, strict = a$strict,
    nniShare = paste(round(shares[a$perChain$kind == "nni"], 2),
                     collapse = "/"),
    randShare = paste(round(shares[a$perChain$kind == "random"], 2),
                      collapse = "/"),
    assessWall = round(a$wall, 1),
    costIter = a$iter / mainIter,
    costWall = round(a$wall / main$wall, 2)
  )
  print(row, row.names = FALSE)
  rows[[length(rows) + 1]] <- row
  write.csv(do.call(rbind, rows),
            paste0("dev/red-team/explore-195/q2-results-", hardName, ".csv"),
            row.names = FALSE)
}
print(do.call(rbind, rows), row.names = FALSE)
