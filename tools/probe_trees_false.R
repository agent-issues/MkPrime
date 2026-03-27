#!/usr/bin/env Rscript
# Confirm trees=FALSE bypasses the crash.
cat("PROBE_TF: starting\n")
suppressMessages(devtools::load_all(quiet = TRUE))
library(ape)

tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
mat  <- matrix(c(0L, 1L, 0L, 1L, 0L, 0L, 1L, 1L), 4, 2,
               dimnames = list(paste0("t", 1:4), NULL))
pd   <- TreeTools::MatrixToPhyDat(mat)

set.seed(2204)
result <- RunMkPrime(pd, tree,
  mcmc = MkPrimeMCMC(nRuns = 2L, nIter = 1000L, thin = 5L, warmup = 500L))

cat("PROBE_TF: calling ConvergenceDiagnostics(trees=FALSE)\n")
diag <- ConvergenceDiagnostics(result, trees = FALSE)
cat("PROBE_TF: ConvergenceDiagnostics done\n")
print(diag)
cat("PROBE_TF: RESULT=OK\n")
