# sim3-pilot.R ----------------------------------------------------------------
# Pilot grid for the Sim 3 convergent-clades validation: across
# (nEco, nBase, phi), simulate 8 replicates each and check whether the
# convergent-misled topology beats the true topology under parsimony.
# Lower-bound demonstration: if parsimony already fails on a config,
# Bayesian Mk inference is very likely to fail too — and the ecology
# layer is the proposed remedy.
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/simulations/ecology/sim3-helpers.R")
source("inst/simulations/ecology/sim3-simulate.R")

treeAB <- Preorder(ape::read.tree(text =
  paste0("(((A1,A2),(A3,A4)),((B1,B2),(B3,B4)),",
         "(((C1,C2),(C3,C4)),((D1,D2),(D3,D4))));")))

# A single trial draws (tree, ecology) once and a config supplies the
# character counts / phi / branch-length scaling.
configs <- list(
  list(nEco = 30,  nBase = 90,  phi = 4,   stem = 0.05, root = 0.10),
  list(nEco = 30,  nBase = 180, phi = 4,   stem = 0.05, root = 0.10),
  list(nEco = 30,  nBase = 90,  phi = 4,   stem = 0.10, root = 0.15),
  list(nEco = 30,  nBase = 180, phi = 4,   stem = 0.10, root = 0.15),
  list(nEco = 60,  nBase = 180, phi = 4,   stem = 0.05, root = 0.10),
  list(nEco = 60,  nBase = 180, phi = 4,   stem = 0.10, root = 0.15)
)

out <- data.frame()
for (cfg in configs) {
  nEco  <- cfg$nEco; nBase <- cfg$nBase; phi <- cfg$phi
  stem  <- cfg$stem; root  <- cfg$root
  tree   <- .BuildConvergentTree(stemBranch = stem, rootBranch = root)
  eco    <- .ConvergentEcology(tree)
  edgeEc <- .AssignEdgeEcology(tree, eco)
  psDiff <- function(d) {
    pd <- MatrixToPhyDat(d)
    phangorn::parsimony(tree, pd) - phangorn::parsimony(treeAB, pd)
  }
  zFull <- matrix(0L, nrow = nEco + nBase, ncol = 2)
  zFull[seq_len(nEco), 2] <- 1L
  cFull <- c(rep("neomorphic", nEco),
             rep("transformational", nBase))
  zBase <- matrix(0L, nrow = nBase, ncol = 2)
  cBase <- rep("transformational", nBase)
  fullDiffs <- numeric(0); baseDiffs <- numeric(0)
  for (seed in 1:8) {
    set.seed(seed)
    dFull <- .SimulateMkPrimeEcology(tree, edgeEc, zFull, phi = phi,
                                     type = cFull, baseRate = 0.5,
                                     rateLoss = 1)
    fullDiffs <- c(fullDiffs, psDiff(dFull))
    set.seed(seed + 1000L)
    dBase <- .SimulateMkPrimeEcology(tree, edgeEc, zBase, phi = phi,
                                     type = cBase, baseRate = 0.5,
                                     rateLoss = 1)
    baseDiffs <- c(baseDiffs, psDiff(dBase))
  }
  out <- rbind(out, data.frame(
    nEco = nEco, nBase = nBase, phi = phi,
    stem = stem, root = root,
    medBase = median(baseDiffs), minBase = min(baseDiffs),
    nTrueBase = sum(baseDiffs < 0),
    medFull = median(fullDiffs), maxFull = max(fullDiffs),
    nFooled = sum(fullDiffs > 0)
  ))
}
print(out, row.names = FALSE)
cat("\nReading the table:\n")
cat("  medBase  median(true - wrong) parsimony on BASELINE-only data.\n")
cat("            Strongly negative = true tree wins on signal alone.\n")
cat("  nTrueBase / 8 baseline-only replicates where truth wins.\n")
cat("  medFull  median(true - wrong) parsimony on FULL data.\n")
cat("            Positive = blind Mk fooled.\n")
cat("  nFooled  / 8 full-data replicates where wrong tree wins.\n")
cat("\nGoldilocks: nTrueBase = 8 AND nFooled high (>= 6) AND medFull modest (3-7).\n")
