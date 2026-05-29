# Root-invariance probe: does the MCMC tree-move machinery preserve the root
# bipartition, or does it let the root float? (Model-agnostic — tree moves are
# the same regardless of likelihood model, so a blind MkNT chain answers it.)
#
#   root bipartition VARIES across chain  -> moves re-root freely => EBE-LIVE
#                                            (non-stationary likelihood ill-posed
#                                             unless we impose a root rule).
#   root bipartition CONSTANT (topology still moves) -> moves preserve root => moot.
#
# Doubles as the post-fix regression test: after Phase 2 anchors the root, the
# "distinct root bipartitions" count must collapse to 1.
suppressMessages({devtools::load_all(".", quiet = TRUE); library(ape); library(TreeTools)})

set.seed(20260529L)
nTip  <- 12L
tips  <- paste0("t", seq_len(nTip))
nChar <- 60L
m <- matrix(0L, nTip, nChar, dimnames = list(tips, NULL))
for (cc in seq_len(nChar)) {
  repeat { col <- sample(0:1, nTip, TRUE); if (length(unique(col)) > 1L) break }
  m[, cc] <- col
}
pd  <- TreeTools::MatrixToPhyDat(m)
mkd <- MkPrimeData(pd, neomorphic = seq_len(nChar))

# Start tree rooted so t1 is the lone basal outgroup.
tree0 <- TreeTools::Preorder(ape::root(ape::rtree(nTip, tip.label = tips),
                                       outgroup = "t1", resolve.root = TRUE))

model <- MkPrimeModel()             # blind MkNT; exercises the same move set
logFile <- "dev/ecology/rootprobe.log"
if (file.exists(logFile)) file.remove(logFile)
mcmc <- MkPrimeMCMC(nIter = 4000L, nChains = 1L, nRuns = 1L, nCore = 1L,
                    thin = 25L, treeThin = 25L, minWarmup = 300L, maxWarmup = 1000L,
                    logFile = logFile, checkpointFile = NULL)
res <- RunMkPrime(mkd, tree = tree0, model = model, mcmc = mcmc)

tf <- sub("\\.log$", "_trees.nwk", logFile)
if (!file.exists(tf)) tf <- paste0(logFile, "_trees.nwk")
tr <- ape::read.tree(tf)
if (inherits(tr, "phylo")) tr <- list(tr)
class(tr) <- "multiPhylo"

rootBipart <- function(t) {
  nT <- ape::Ntip(t); rt <- nT + 1L
  kids <- t$edge[t$edge[, 1] == rt, 2]
  s <- vapply(kids, function(k)
    if (k <= nT) t$tip.label[k]
    else paste(sort(ape::extract.clade(t, k)$tip.label), collapse = ","),
    character(1))
  paste(sort(s), collapse = "  |  ")
}
t1IsBasal <- function(t) {
  nT <- ape::Ntip(t); rt <- nT + 1L
  match("t1", t$tip.label) %in% t$edge[t$edge[, 1] == rt, 2]
}

rb  <- vapply(tr, rootBipart, character(1))
t1b <- vapply(tr, t1IsBasal, logical(1))
rf  <- vapply(tr, function(t) ape::dist.topo(ape::unroot(t), ape::unroot(tree0)), numeric(1))

cat(sprintf("\n=== ROOT PROBE (%d logged trees) ===\n", length(tr)))
cat(sprintf("distinct root bipartitions across chain : %d\n", length(unique(rb))))
cat(sprintf("t1 still lone basal outgroup            : %d / %d (%.0f%%)\n",
            sum(t1b), length(t1b), 100 * mean(t1b)))
cat(sprintf("unrooted RF from start  : min=%g max=%g  (moves accepted iff max>0)\n",
            min(rf), max(rf)))
cat("\nmost frequent root bipartitions:\n"); print(head(sort(table(rb), decreasing = TRUE), 6))
cat(if (length(unique(rb)) > 1) "\nVERDICT: ROOT FLOATS -> EBE-LIVE (Phase 2 must anchor root)\n"
    else "\nVERDICT: root stable across moves -> moot\n")
