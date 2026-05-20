# Sim 3 v3 redesign — parsimony grid scoping
#
# Goal: find a (nChar, stem, root, tipBranch, phi) combination where:
#   - Baseline data (no ecology effect): truth tree wins on parsimony (8/8)
#   - Full data (ecology effect on): blind parsimony is fooled (nFooled high,
#     medFull modestly positive)
#   - Signal stronger than the existing Goldilocks (nEco=60, nBase=180,
#     phi=4, stem=0.10, root=0.15) so the Bayesian chain can actually
#     recover the truth topology.
#
# Strategy: scale chars up (2x and 4x current Goldilocks), lengthen
# stems (0.10, 0.20, 0.30) for more discriminating-edge signal, and try
# stronger phi (4, 6).
#
# Output: a table of (config, medBase, nTrueBase, medFull, nFooled).
# Goldilocks v3 should have nTrueBase = 8, nFooled >= 6, medFull in 3-12.
#
# Run from worktree root:
#   Rscript -e "setwd('C:/Users/pjjg18/GitHub/mkp/.claude/worktrees/ecology-aware'); source('dev/pilots/2026-05-14-sim3-v3-redesign/parsimony-grid.R')"

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
  library("phangorn")
})
source("inst/ecology/simulations/sim3-helpers.R")
source("inst/ecology/simulations/sim3-simulate.R")

# Wrong topology: ((A,B),(C,D)) — the convergent-driven false grouping
treeAB <- Preorder(ape::read.tree(text =
  paste0("(((A1,A2),(A3,A4)),((B1,B2),(B3,B4)),",
         "(((C1,C2),(C3,C4)),((D1,D2),(D3,D4))));")))

# v3 grid: scale up chars + lengthen stems; sample phi
configs <- list(
  # Reference (existing Goldilocks)
  list(name = "v2-ref",     nEco = 60,  nBase = 180, phi = 4, stem = 0.10, root = 0.15, tip = 0.5),
  # 2x chars
  list(name = "2x-chars",   nEco = 120, nBase = 360, phi = 4, stem = 0.10, root = 0.15, tip = 0.5),
  # 4x chars
  list(name = "4x-chars",   nEco = 240, nBase = 720, phi = 4, stem = 0.10, root = 0.15, tip = 0.5),
  # 2x chars + 2x stems
  list(name = "2x-2x-stem", nEco = 120, nBase = 360, phi = 4, stem = 0.20, root = 0.15, tip = 0.5),
  # 2x chars + 3x stems
  list(name = "2x-3x-stem", nEco = 120, nBase = 360, phi = 4, stem = 0.30, root = 0.15, tip = 0.5),
  # 4x chars + 2x stems
  list(name = "4x-2x-stem", nEco = 240, nBase = 720, phi = 4, stem = 0.20, root = 0.15, tip = 0.5),
  # Stronger phi
  list(name = "2x-phi6",    nEco = 120, nBase = 360, phi = 6, stem = 0.10, root = 0.15, tip = 0.5),
  # 4x chars only (alternative top candidate)
  list(name = "4x-phi6",    nEco = 240, nBase = 720, phi = 6, stem = 0.10, root = 0.15, tip = 0.5)
)

cat("Sim 3 v3 parsimony grid\n")
cat(sprintf("Configs tested: %d\n", length(configs)))
cat("Truth topo: ((A,C),(B,D)).  Wrong (convergent) topo: ((A,B),(C,D))\n")
cat("Δ = parsimony(true) − parsimony(wrong); negative = truth wins\n\n")

out <- data.frame()
for (cfg in configs) {
  nEco  <- cfg$nEco; nBase <- cfg$nBase; phi <- cfg$phi
  stem  <- cfg$stem; root  <- cfg$root; tipBr <- cfg$tip
  tree   <- .BuildConvergentTree(tipBranch = tipBr, stemBranch = stem,
                                  rootBranch = root)
  eco    <- .ConvergentEcology(tree)
  edgeEc <- .AssignEdgeEcology(tree, eco)
  psDiff <- function(d) {
    pd <- MatrixToPhyDat(d)
    phangorn::parsimony(tree, pd) - phangorn::parsimony(treeAB, pd)
  }
  zFull <- matrix(0L, nrow = nEco + nBase, ncol = 2)
  zFull[seq_len(nEco), 2] <- 1L
  cFull <- c(rep("neomorphic", nEco), rep("transformational", nBase))
  zBase <- matrix(0L, nrow = nBase, ncol = 2)
  cBase <- rep("transformational", nBase)
  fullDiffs <- numeric(0); baseDiffs <- numeric(0)
  for (seed in 1:8) {
    set.seed(seed)
    dFull <- .SimulateMkPrimeEcology(tree, edgeEc, zFull, phi = phi,
                                     type = cFull, baseRate = 1.0,
                                     normalize = TRUE,
                                     pi0 = 0.75, theta = 1.0,
                                     refEcology = 0L,
                                     rateLoss = 1)
    fullDiffs <- c(fullDiffs, psDiff(dFull))
    set.seed(seed + 1000L)
    dBase <- .SimulateMkPrimeEcology(tree, edgeEc, zBase, phi = phi,
                                     type = cBase, baseRate = 1.0,
                                     normalize = TRUE,
                                     pi0 = 0.75, theta = 1.0,
                                     refEcology = 0L,
                                     rateLoss = 1)
    baseDiffs <- c(baseDiffs, psDiff(dBase))
  }
  out <- rbind(out, data.frame(
    name = cfg$name,
    nEco = nEco, nBase = nBase, phi = phi,
    stem = stem, root = root, tip = tipBr,
    treeTL = sum(tree$edge.length),
    medBase = median(baseDiffs), nTrueBase = sum(baseDiffs < 0),
    medFull = median(fullDiffs), nFooled = sum(fullDiffs > 0),
    rangeFull = paste0(min(fullDiffs), "..", max(fullDiffs))
  ))
}
print(out, row.names = FALSE)
cat("\nGoldilocks v3 criteria:\n")
cat("  - nTrueBase = 8  AND  nFooled >= 6  AND  medFull in 3..12\n")
cat("  - Higher medBase magnitude better (stronger phylogenetic signal\n")
cat("    when ecology turned off — gives blind a chance to recover)\n")

# Save table
saveRDS(out, "dev/pilots/2026-05-14-sim3-v3-redesign/parsimony-grid.rds")
write.csv(out, "dev/pilots/2026-05-14-sim3-v3-redesign/parsimony-grid.csv",
          row.names = FALSE)
cat("\nSaved to dev/pilots/2026-05-14-sim3-v3-redesign/parsimony-grid.{rds,csv}\n")
