# Shared helpers for the #195 exploration. Not package code.
#
# Load the agent build first:
#   library("MkPrime.a05d67", lib.loc = "C:/Users/pjjg18/GitHub/.builds/MkPrime-a05d67")
# so that `mkp` below resolves to its namespace.

suppressPackageStartupMessages({
  library("TreeTools")
  library("TreeSearch")
})

mkpLib <- Sys.getenv("MKP_LIB", "C:/Users/pjjg18/GitHub/.builds/MkPrime-a05d67")
mkpPkg <- Sys.getenv("MKP_PKG", "MkPrime.a05d67")
suppressPackageStartupMessages(
  library(mkpPkg, lib.loc = mkpLib, character.only = TRUE)
)
mkp <- asNamespace(mkpPkg)

# Rooting-insensitive topology fingerprint: sorted split strings
TopoKey <- function(tree) {
  tree <- RenumberTips(tree, sort(tree$tip.label))
  sp <- as.Splits(tree)
  paste(sort(as.character(PolarizeSplits(sp, 1))), collapse = "|")
}

TopoKeys <- function(trees) vapply(trees, TopoKey, character(1))

# Summary of a tree sample: distinct topologies and MAP frequency
TopoSummary <- function(trees) {
  keys <- TopoKeys(trees)
  tab <- sort(table(keys), decreasing = TRUE)
  list(n = length(keys), distinct = length(tab),
       mapFreq = as.numeric(tab[1]) / length(keys),
       changes = sum(keys[-1] != keys[-length(keys)]),
       mapKey = names(tab)[1], keys = keys)
}

# 6 tips, 160 characters each congruent with one tree (the RT9-D case)
CongruentData <- function(nTip = 6, nChar = 160) {
  tree <- BalancedTree(nTip)
  splits <- as.logical(as.Splits(tree))
  mat <- t(vapply(seq_len(nChar), function(i) {
    as.integer(splits[(i - 1) %% nrow(splits) + 1, ])
  }, integer(nTip)))
  colnames(mat) <- TipLabels(tree)
  list(data = MatrixToPhyDat(t(mat)), tree = tree)
}

# Short run with tuning off. Returns posterior + wall time.
ShortRun <- function(data, tree, nIter = 4000, minWarmup = 500,
                     maxWarmup = max(minWarmup, 1000L), maxTime = 60,
                     moveWeights = NULL, autoTune = FALSE,
                     seed = 1, ...) {
  set.seed(seed)
  tree$edge.length <- rep(0.1, nrow(tree$edge))
  t0 <- proc.time()[["elapsed"]]
  post <- RunMkPrime(data, tree = tree, nIter = nIter, minWarmup = minWarmup,
                     maxWarmup = maxWarmup,
                     autoTune = autoTune, nRuns = 1L, maxTime = maxTime,
                     thin = 1L, moveWeights = moveWeights, verbosity = 0L,
                     ...)
  list(post = post, wall = proc.time()[["elapsed"]] - t0)
}

# Topology move weights cut to near zero: manufactures a stuck sampler
WeakTopology <- c(nni = 1e-4, spr = 1e-4, tbr = 1e-4, pspr = 1e-4,
                  gibbs_spr = 1e-4, gibbs_subtree_swap = 1e-4)

# Martin's assessment phase (prototype). From `stuckTree`, start `k` chains
# at NNI neighbours and `m` at random trees; each runs `budget` iterations.
# "favoured" when every chain's modal topology over its second half is the
# stuck topology (lenient) / when every second-half tree is (strict).
AssessPhase <- function(data, stuckTree, k = 3, m = 3, budget = 2000,
                        minWarmup = 300, seed = 1, maxTime = 30) {
  set.seed(seed)
  target <- TopoKey(stuckTree)
  starts <- c(
    lapply(seq_len(k), function(i) {
      # A NNI neighbour that is a different topology
      repeat {
        nb <- phangorn::rNNI(Preorder(stuckTree), moves = 1)
        if (TopoKey(nb) != target) break
      }
      nb
    }),
    lapply(seq_len(m), function(i) {
      RandomTree(stuckTree$tip.label, root = TRUE)
    })
  )
  kind <- rep(c("nni", "random"), c(k, m))
  t0 <- proc.time()[["elapsed"]]
  res <- lapply(seq_along(starts), function(i) {
    run <- ShortRun(data, starts[[i]], nIter = budget, minWarmup = minWarmup,
                    maxWarmup = minWarmup, maxTime = maxTime, seed = seed * 100 + i)
    trees <- run$post$trees
    second <- trees[seq.int(ceiling(length(trees) / 2), length(trees))]
    keys <- TopoKeys(second)
    tab <- sort(table(keys), decreasing = TRUE)
    c(modal = names(tab)[1] == target, strict = all(keys == target),
      share = mean(keys == target))
  })
  res <- do.call(rbind, res)
  list(
    lenient = if (all(res[, "modal"] == 1)) "favoured" else "stuck",
    strict = if (all(res[, "strict"] == 1)) "favoured" else "stuck",
    perChain = data.frame(kind, res),
    wall = proc.time()[["elapsed"]] - t0,
    iter = length(starts) * budget
  )
}
