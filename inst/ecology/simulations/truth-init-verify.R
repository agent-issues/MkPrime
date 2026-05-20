# truth-init-verify.R — recompute log_post at chain sample 1 to test
# whether chain.log values are accurate or stale by ~120 nats.
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/ecology/simulations/sim3-helpers.R")
source("inst/ecology/simulations/sim3-simulate.R")

# Rebuild data (matching truth_init.R, seed 20260602)
nEco <- 120L; nBase <- 360L; phi <- 4
stemBr <- 0.30; rootBr <- 0.15; tipBr <- 0.5
tree <- .BuildConvergentTree(tipBranch = tipBr,
                              stemBranch = stemBr, rootBranch = rootBr)
eco  <- .ConvergentEcology(tree)
edgeEc <- .AssignEdgeEcology(tree, eco)
set.seed(20260602L)
zFull <- matrix(0L, nrow = nEco + nBase, ncol = 2)
zFull[seq_len(nEco), 2] <- 1L
cFull <- c(rep("neomorphic", nEco), rep("transformational", nBase))
datSim <- .SimulateMkPrimeEcology(tree, edgeEc, zFull, phi = phi,
                                  type = cFull, baseRate = 1.0,
                                  normalize = TRUE,
                                  pi0 = 0.75, theta = 1.0,
                                  refEcology = 0L, rateLoss = 1)
pdSim <- MatrixToPhyDat(datSim)
mkdBlind <- MkPrimeData(pdSim, neomorphic = seq_len(nEco))
ecoTipVec <- eco[mkdBlind$taxon_names]
mkdAware <- MkPrimeData(pdSim, neomorphic = seq_len(nEco), ecology = ecoTipVec)

modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)

# Read sample 1's tree
trees <- ape::read.tree("inst/ecology/simulations/truth-init-results/truth-init_trees.nwk")
cat("nTrees:", length(trees), "\n")
treeS1 <- if (inherits(trees, "multiPhylo")) trees[[1]] else trees
treeS1 <- Preorder(treeS1)
cat("Tree 1 nTip:", TreeTools::NTip(treeS1), "  TL:", sum(treeS1$edge.length), "\n")

# Pull sample-1 parameters from log
logCsv <- read.table("inst/ecology/simulations/truth-init-results/truth-init.log",
                     sep = "\t", header = TRUE, comment.char = "#",
                     check.names = FALSE)
s1 <- logCsv[1, ]
cat("\n=== Sample 1 logged values ===\n")
print(s1[, c("Sample", "log_posterior", "log_likelihood", "tree_length",
             "rate_loss", "rate_log_sd", "rate_neo", "p",
             "phi", "pi0", "theta_1")])

# The log doesn't store z. We don't know z exactly. Use a reasonable proxy:
# pretend z = all-zero (impossible to match) OR use the saved chain state's
# final z. For sample 1 specifically, we don't have z. So compare a range.

# Read final chain z from the result rds if available
res <- readRDS("inst/ecology/simulations/truth-init-results/truth-init-result.rds")
zFinal <- res$z_samples
cat("\n=== z_samples structure ===\n")
str(zFinal, max.level = 2)
if (!is.null(zFinal) && is.list(zFinal) && length(zFinal) > 0) {
  z1 <- zFinal[[1]]
  cat("z first sample dim:", paste(dim(z1), collapse = "x"), "\n")
  cat("z first sample table:\n"); print(table(z1))

  zLast <- zFinal[[length(zFinal)]]
  cat("z last sample table:\n"); print(table(zLast))
}

# Compute log_lik at sample 1 with its z
if (!is.null(zFinal) && length(zFinal) > 0) {
  z1 <- zFinal[[1]]
  scaleTL <- sum(treeS1$edge.length)
  # Pull sample 1's actual kPrime values from the log. The chain has
  # been running gibbs_kPrime since init, so kObs is stale and using it
  # would inflate the recomputed gap by 100+ nats. kPrime columns are
  # named kPrime_118, kPrime_119, ..., one per character (transformational
  # only; neomorphic = kObs).
  kpCols <- grep("^kPrime_\\d+$", colnames(logCsv), value = TRUE)
  if (length(kpCols) == 0) stop("no kPrime columns in log")
  # logCsv kpCols cover transformational chars only. Build the full kPrime
  # vector by starting from kObs and overwriting trans positions.
  kpFull <- as.integer(mkdAware$kObs)
  charIdx <- as.integer(sub("kPrime_", "", kpCols))
  kpFull[charIdx] <- as.integer(unlist(s1[, kpCols]))
  cat("kPrime cols in log:", length(kpCols),
      " range charIdx:", min(charIdx), "..", max(charIdx), "\n")
  cat("nChar mkd:", mkdAware$nChar, " nNeo:",
      sum(mkdAware$type == "neomorphic"), "\n")
  cat("kpFull table (head):\n"); print(table(kpFull))
  cat("kObs   table (head):\n"); print(table(as.integer(mkdAware$kObs)))
  cat("kpFull vs kObs differ at:", sum(kpFull != as.integer(mkdAware$kObs)),
      "positions\n")
  ll <- .MkpEcologyLogLikelihood(
    treeS1, mkdAware,
    kPrime = kpFull,
    rate_loss = s1$rate_loss,
    rate_log_sd = s1$rate_log_sd,
    nCat = modelAware$nCat,
    rate_neo = s1$rate_neo,
    relabel = modelAware$relabel,
    phi = s1$phi, zMat = z1,
    magnitudeMode = modelAware$magnitudeMode,
    coding = modelAware$coding,
    refEcology = mkdAware$refEcology,
    theta = s1$theta_1, pi0 = s1$pi0)
  cat("\nRECOMPUTED log_lik at sample-1 (with z1):", ll, "\n")

  st <- list(
    tree = treeS1,
    tree_length = sum(treeS1$edge.length),
    rel_br_lengths = treeS1$edge.length / sum(treeS1$edge.length),
    rate_loss = s1$rate_loss, rate_log_sd = s1$rate_log_sd,
    rate_neo = s1$rate_neo, p = s1$p,
    kPrime = kpFull,
    phi = s1$phi, pi0 = s1$pi0, theta = s1$theta_1, z = z1
  )
  lp <- LogPrior(st, modelAware, mkdAware)
  cat("RECOMPUTED log_prior at sample-1:", lp, "\n")
  cat("RECOMPUTED log_post  at sample-1:", ll + lp, "\n")
  cat("LOGGED     log_post  at sample-1:", s1$log_posterior, "\n")
  cat("Difference (recomputed - logged):", (ll + lp) - s1$log_posterior, "nats\n")
}
