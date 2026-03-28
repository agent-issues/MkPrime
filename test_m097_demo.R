# M-097 progress display demo
# Run this interactively in the console to see:
#   1. "warmup N" vs "iter N" labelling
#   2. Whether consecutive progress tables overwrite or stack
#
# Runs a short MCMC (6 taxa, 3 chars) with checkEvery = 500 so
# you'll see several table prints in quick succession.

library(ape)
library(TreeTools)

tree <- read.tree(text = "((t1:0.1,(t2:0.1,t3:0.1):0.1):0.1,(t4:0.1,(t5:0.1,t6:0.1):0.1):0.1);")
mat <- matrix(
  c(0,1,0,1,0,1, 0,0,1,1,0,1, 0,1,1,0,1,0), 6, 3,
  dimnames = list(paste0("t", 1:6), NULL)
)
pd <- MatrixToPhyDat(mat)

mkd <- MkPrimeData(pd)
model <- MkPrimeModel(coding = "variable")
mcmc <- MkPrimeMCMC(
  nIter      = 5000L,
  warmup     = 1000L,
  thin       = 1L,
  nRuns      = 1L,
  checkEvery = 500L
)

cat("\n=== Starting MCMC — watch for warmup/iter label and table overwriting ===\n\n")
post <- RunMkPrime(mkd, model, mcmc, tree = tree, fixTopology = TRUE)
cat("\n=== Done ===\n")
