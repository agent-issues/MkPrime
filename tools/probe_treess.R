#!/usr/bin/env Rscript
# Isolate the treess::treess() crash.
cat("PROBE_TREESS: starting\n")
suppressMessages(devtools::load_all(quiet = TRUE))
library(ape)

tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
               dimnames = list(paste0("t", 1:4), NULL))
pd   <- TreeTools::MatrixToPhyDat(mat)

set.seed(2204)
result <- RunMkPrime(pd, tree,
  mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 1000L, thin = 5L, warmup = 500L))

cat("PROBE_TREESS: RunMkPrime returned\n")

pb <- .PostBurninData(result)
cat(sprintf("PROBE_TREESS: per_run has %d runs\n", length(pb$per_run)))

perRunTrees <- lapply(pb$per_run, `[[`, "trees")
cat(sprintf("PROBE_TREESS: run 1 has %d trees, run 2 has %d trees\n",
  length(perRunTrees[[1]]), length(perRunTrees[[2]])))
cat(sprintf("PROBE_TREESS: tree class: %s\n",
  class(perRunTrees[[1]][[1]])))
cat(sprintf("PROBE_TREESS: tree tip labels: %s\n",
  paste(perRunTrees[[1]][[1]]$tip.label, collapse=", ")))
cat(sprintf("PROBE_TREESS: nTip: %d\n",
  length(perRunTrees[[1]][[1]]$tip.label)))

# Get treess method names first
cat("PROBE_TREESS: getting ESS methods\n")
methods <- treess::getESSMethods(TRUE)
cat(sprintf("PROBE_TREESS: methods: %s\n", paste(names(methods), collapse=", ")))

# Now try the call that crashes
cat("PROBE_TREESS: calling treess::treess()\n")
result_ess <- treess::treess(perRunTrees, TreeDist::RobinsonFoulds,
                              methods = methods)
cat("PROBE_TREESS: treess returned — RESULT=OK\n")
