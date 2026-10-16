# Q3: today's behaviour of the RT9-D configurations at the current main.
# 6 tips, 160 congruent characters, minTreeEss = 100, maxTime = 60.
# Run: timeout 400 Rscript dev/red-team/explore-195/q3-reachability.R
setTimeLimit(elapsed = 380, transient = FALSE)
source("dev/red-team/explore-195/common.R")
cd <- CongruentData()
start <- AdditionTree(cd$data)
start$edge.length <- rep(0.1, nrow(start$edge))

Try <- function(label, ...) {
  set.seed(1)
  t0 <- proc.time()[["elapsed"]]
  warns <- character(0)
  post <- withCallingHandlers(
    RunMkPrime(cd$data, tree = start, minTreeEss = 100, maxTime = 60,
               minWarmup = 500, verbosity = 0L, ...),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  stopReason <- post$stopReason %||% post$stop_reason %||%
    attr(post, "stopReason") %||% NA
  data.frame(config = label, wall = round(proc.time()[["elapsed"]] - t0),
             stopReason = paste(stopReason, collapse = ","),
             nTrees = length(post$trees),
             distinct = TopoSummary(post$trees)$distinct,
             stuckWarnings = sum(grepl("single tree topology", warns)),
             budgetWarning = any(grepl("reached during run", warns)))
}

out <- rbind(
  Try("nRuns = 1", nRuns = 1L),
  Try("serial nRuns = 2", nRuns = 2L, nCore = 1L),
  Try("parallel nRuns = 2", nRuns = 2L, nCore = 2L)
)
print(out, row.names = FALSE)
write.csv(out, "dev/red-team/explore-195/q3-results.csv", row.names = FALSE)
