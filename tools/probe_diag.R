#!/usr/bin/env Rscript
# Probe: which step of ConvergenceDiagnostics crashes?
# Run from the mkp package root directory.
# Uses load_all() so all internal functions are available directly.

cat("PROBE_DIAG: starting\n")
suppressMessages(devtools::load_all(quiet = TRUE))
library(ape)

tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
               dimnames = list(paste0("t", 1:4), NULL))
pd   <- TreeTools::MatrixToPhyDat(mat)

set.seed(2204)
result <- RunMkPrime(pd, tree,
  mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 1000L, thin = 5L, warmup = 500L))

cat("PROBE_DIAG: RunMkPrime returned\n")

# Is treess/TreeDist installed?
cat(sprintf("PROBE_DIAG: treess available: %s\n",
  requireNamespace("treess", quietly = TRUE)))
cat(sprintf("PROBE_DIAG: TreeDist available: %s\n",
  requireNamespace("TreeDist", quietly = TRUE)))

# Step through ConvergenceDiagnostics manually (functions from load_all)
cat("PROBE_DIAG: calling .PostBurninData\n")
pb <- .PostBurninData(result)
cat(sprintf("PROBE_DIAG: .PostBurninData done; nrow(samples)=%d, length(per_run)=%d\n",
  nrow(pb$samples), length(pb$per_run)))

cat("PROBE_DIAG: calling .KeyParamCols\n")
keyCols <- .KeyParamCols(pb$samples)
cat(sprintf("PROBE_DIAG: .KeyParamCols done: %d cols\n", length(keyCols)))

cat("PROBE_DIAG: calling .ComputeEss\n")
ess <- .ComputeEss(pb$samples[, keyCols, drop = FALSE])
cat("PROBE_DIAG: .ComputeEss done\n")

cat("PROBE_DIAG: calling .ComputePsrf\n")
psrf <- .ComputePsrf(pb$per_run, keyCols)
cat("PROBE_DIAG: .ComputePsrf done\n")

cat("PROBE_DIAG: calling .ComputeTreeEss\n")
treeEss <- .ComputeTreeEss(pb, TRUE)
cat(sprintf("PROBE_DIAG: .ComputeTreeEss done: %s\n",
  if (is.null(treeEss)) "NULL" else paste(names(treeEss), collapse = ", ")))

cat("PROBE_DIAG: calling ConvergenceDiagnostics\n")
diag <- ConvergenceDiagnostics(result)
cat("PROBE_DIAG: ConvergenceDiagnostics done\n")

cat("PROBE_DIAG: calling print(diag)\n")
print(diag)
cat("PROBE_DIAG: print(diag) done — RESULT=OK\n")
