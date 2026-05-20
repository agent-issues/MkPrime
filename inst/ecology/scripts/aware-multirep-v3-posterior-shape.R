# aware-multirep-v3-posterior-shape.R ------------------------------------------
#
# Characterise the multirep-v3 aware-chain posterior to see whether aware
# is doing the methodologically RIGHT thing (broad honest uncertainty
# including truth) or the wrong thing (failing to find truth at all).
#
# For each rep x chain, this script:
#   1. Loads the saved post-MCMC trees, discards 25 % burn-in.
#   2. Computes per-tree
#        - canonical (unrooted) topology hash via TreeTools::as.Splits
#        - presence of trueSister_AC and falseSister_AB (root-invariant)
#        - normalised Clustering Information Distance to the true tree.
#   3. Summarises CID distribution to truth (mean, median, 5/25/50/75/95).
#   4. Builds the 95 % credibility set (smallest set of unique canonical
#      topologies whose cumulative posterior >= 0.95).
#   5. Tests whether the canonical hash of the TRUTH tree is in
#         (a) the posterior sample at all, and
#         (b) the 95 % credibility set.
#   6. Reports top-5 canonical topologies and their frequencies, with
#      annotation of which contain AC or AB.
#   7. Writes a tidy CSV of per-(rep,chain) summary stats plus a CSV of
#      per-(rep,chain,topology) credible-set membership.
#   8. Writes a markdown report and a boxplot figure comparing aware vs
#      blind CID-to-truth across reps.
#
# Inputs:
#   inst/ecology/simulations/multirep-v3-results/rep%02d/{blind,aware}-result.rds
#
# Outputs:
#   inst/ecology/scripts/aware-multirep-v3-posterior-shape.rds       (R list)
#   inst/ecology/scripts/aware-multirep-v3-posterior-shape.csv       (per chain summary)
#   inst/ecology/scripts/aware-multirep-v3-posterior-shape-credset.csv (per topo, in CS or not)
#   inst/ecology/scripts/aware-multirep-v3-posterior-shape-report.md
#   inst/ecology/scripts/aware-multirep-v3-posterior-shape-plot.pdf
#
# Run from the MkPrime repo root with:
#   Rscript inst/ecology/scripts/aware-multirep-v3-posterior-shape.R

suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKP_LIB", "")
  if (nzchar(libPath)) .libPaths(c(libPath, .libPaths()))
  library("TreeTools")
  library("TreeDist")
})

mkpRoot <- Sys.getenv("MKP_REPO_ROOT", getwd())
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-scoring.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-helpers.R"))

# ----- multirep-v3 truth config (matches inst/ecology/hamilton/sim3-multirep-v3/run_rep.R)
TIPBR  <- 0.5
STEMBR <- 0.30
ROOTBR <- 0.15
truthTree <- .BuildConvergentTree(tipBranch  = TIPBR,
                                  stemBranch = STEMBR,
                                  rootBranch = ROOTBR)
truthTL <- sum(truthTree$edge.length)

# ----- bipartitions of interest
biparts <- list(
  trueSister_AC  = c(paste0("A", 1:4), paste0("C", 1:4)),
  falseSister_AB = c(paste0("A", 1:4), paste0("B", 1:4))
)

# ----- canonical-splits hash (root-invariant) for topology identification
canonicalSplitsHash <- function(tr) {
  spl <- TreeTools::as.Splits(tr, tipLabels = tr$tip.label)
  m <- as.logical(spl)
  if (is.null(dim(m))) m <- matrix(m, nrow = 1L)
  rows <- apply(m, 1L, function(r) {
    k1 <- paste(which(r),  collapse = ",")
    k2 <- paste(which(!r), collapse = ",")
    if (k1 < k2) k1 else k2
  })
  paste(sort(rows), collapse = "|")
}

truthHash <- canonicalSplitsHash(truthTree)

# ----- per-rep x chain processor returning a list of summary values plus
#       the per-tree vectors needed for plotting.
processChain <- function(repId, chain, repDir) {
  rdsPath <- file.path(repDir, sprintf("%s-result.rds", chain))
  if (!file.exists(rdsPath)) return(NULL)
  res <- readRDS(rdsPath)
  trees <- res$trees
  if (is.null(trees) || length(trees) == 0L) return(NULL)
  nTreesTotal <- length(trees)

  trBurned <- .DiscardBurnin(trees)
  class(trBurned) <- "multiPhylo"
  n <- length(trBurned)
  if (n == 0L) return(NULL)

  hashes <- vapply(trBurned, canonicalSplitsHash, character(1))
  inAC <- HasBipartSplits(trBurned, biparts$trueSister_AC)
  inAB <- HasBipartSplits(trBurned, biparts$falseSister_AB)
  cidVec <- as.numeric(TreeDist::ClusteringInfoDistance(trBurned, truthTree,
                                                       normalize = TRUE))

  freqTab <- sort(table(hashes), decreasing = TRUE)
  freq <- as.integer(freqTab)
  names(freq) <- names(freqTab)
  prob <- freq / n
  cumProb <- cumsum(prob)
  inCS <- cumProb <= 0.95
  # Always include the topology that crosses the 0.95 line
  firstOut <- which(!inCS)
  if (length(firstOut)) inCS[firstOut[1]] <- TRUE
  csTopo <- names(prob)[inCS]
  nCS <- length(csTopo)

  # AC / AB indicator per UNIQUE topology
  uniqHashes <- names(prob)
  acIdx <- vapply(uniqHashes, function(h) {
    j <- which(hashes == h)[1]
    inAC[j]
  }, logical(1))
  abIdx <- vapply(uniqHashes, function(h) {
    j <- which(hashes == h)[1]
    inAB[j]
  }, logical(1))
  cidPerTopo <- vapply(uniqHashes, function(h) {
    j <- which(hashes == h)[1]
    cidVec[j]
  }, numeric(1))

  truthInPosterior <- truthHash %in% uniqHashes
  truthInCS <- truthHash %in% csTopo

  list(
    rep = repId, chain = chain,
    n = n,
    pAC = mean(inAC),
    pAB = mean(inAB),
    cidMean   = mean(cidVec),
    cidMedian = median(cidVec),
    cidP05    = unname(quantile(cidVec, 0.05, names = FALSE)),
    cidP25    = unname(quantile(cidVec, 0.25, names = FALSE)),
    cidP75    = unname(quantile(cidVec, 0.75, names = FALSE)),
    cidP95    = unname(quantile(cidVec, 0.95, names = FALSE)),
    cidMin    = min(cidVec),
    cidMax    = max(cidVec),
    nUnique = length(uniqHashes),
    nCS = nCS,
    csMass = if (nCS) cumProb[nCS] else NA_real_,
    truthInPosterior = truthInPosterior,
    truthInCS = truthInCS,
    truthFreq = if (truthInPosterior) freq[truthHash] else 0L,
    truthProb = if (truthInPosterior) prob[truthHash] else 0,
    # Top-5 topologies (for the report)
    top5 = utils::head(data.frame(
      hash = names(prob), freq = freq, prob = prob,
      hasAC = acIdx, hasAB = abIdx, cid = cidPerTopo,
      stringsAsFactors = FALSE,
      row.names = NULL), 5L),
    # Vectors (for plotting and the credset CSV)
    cidVec = cidVec,
    hashes = hashes,
    freq = freq, prob = prob, inCS = inCS,
    hasAC = acIdx, hasAB = abIdx, cidPerTopo = cidPerTopo
  )
}

cat("=== aware-multirep-v3-posterior-shape.R ===\n")
cat("mkpRoot:", mkpRoot, "\n")
cat(sprintf("Truth: TL = %.3f, true sister = (A,C),(B,D); false eco sister = (A,B)\n",
            truthTL))
cat("Truth canonical hash (first 60 chars): ", substr(truthHash, 1, 60), "...\n", sep = "")

resultsDir <- file.path(mkpRoot, "inst/ecology/simulations/multirep-v3-results")
reps <- sort(list.files(resultsDir, pattern = "^rep[0-9]+$"))
cat(sprintf("Found %d rep dirs under %s\n\n", length(reps), resultsDir))

all <- list()
for (rep in reps) {
  repDir <- file.path(resultsDir, rep)
  repNum <- as.integer(sub("rep", "", rep))
  for (chain in c("blind", "aware")) {
    t0 <- Sys.time()
    r <- tryCatch(processChain(repNum, chain, repDir),
                  error = function(e) {
                    cat(sprintf("[%s %s] ERROR: %s\n", rep, chain, conditionMessage(e)))
                    NULL
                  })
    if (is.null(r)) {
      cat(sprintf("[%s %s] SKIP\n", rep, chain)); next
    }
    dt <- format(round(Sys.time() - t0, 1))
    cat(sprintf("[%s %s] n=%4d  P(AC)=%.3f P(AB)=%.3f  CID mean=%.3f [5%%=%.3f, 95%%=%.3f]  nUniq=%d  |CS|=%d  truth in CS? %s  %s\n",
                rep, chain, r$n, r$pAC, r$pAB,
                r$cidMean, r$cidP05, r$cidP95,
                r$nUnique, r$nCS,
                if (r$truthInCS) "YES" else "no", dt))
    all[[length(all) + 1L]] <- r
  }
}

if (!length(all)) stop("No data — aborting.")

# ----- assemble per-(rep,chain) summary data frame --------------------------
sumCols <- c("rep", "chain", "n", "pAC", "pAB",
             "cidMean", "cidMedian", "cidP05", "cidP25", "cidP75", "cidP95",
             "cidMin", "cidMax",
             "nUnique", "nCS", "csMass",
             "truthInPosterior", "truthInCS", "truthFreq", "truthProb")
sumDf <- do.call(rbind, lapply(all, function(r) {
  as.data.frame(r[sumCols], stringsAsFactors = FALSE)
}))
sumDf <- sumDf[order(sumDf$rep, sumDf$chain), , drop = FALSE]

outDir <- file.path(mkpRoot, "inst/ecology/scripts")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(summary = sumDf, perChain = all, truthHash = truthHash,
             truthTL = truthTL),
        file.path(outDir, "aware-multirep-v3-posterior-shape.rds"))
write.csv(sumDf, file.path(outDir, "aware-multirep-v3-posterior-shape.csv"),
          row.names = FALSE)

# ----- per-(rep, chain, topo) credset CSV -----------------------------------
csRows <- list()
for (r in all) {
  hashes <- names(r$freq)
  csRows[[length(csRows) + 1L]] <- data.frame(
    rep = r$rep, chain = r$chain,
    rank = seq_along(hashes),
    prob = unname(r$prob),
    hasAC = r$hasAC, hasAB = r$hasAB,
    cid = r$cidPerTopo,
    inCS95 = r$inCS,
    isTruth = hashes == truthHash,
    hash = substr(hashes, 1, 40),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
csDf <- do.call(rbind, csRows)
write.csv(csDf, file.path(outDir, "aware-multirep-v3-posterior-shape-credset.csv"),
          row.names = FALSE)

# ----- compute paired blind vs aware comparison ----------------------------
blindDf <- sumDf[sumDf$chain == "blind", , drop = FALSE]
awareDf <- sumDf[sumDf$chain == "aware", , drop = FALSE]
blindDf <- blindDf[order(blindDf$rep), , drop = FALSE]
awareDf <- awareDf[order(awareDf$rep), , drop = FALSE]
nReps <- nrow(blindDf)

deltaCID <- awareDf$cidMean - blindDf$cidMean
# Paired Wilcoxon: under H0 of no shift, where is aware vs blind?
wtest <- tryCatch(
  wilcox.test(awareDf$cidMean, blindDf$cidMean, paired = TRUE,
              alternative = "two.sided", exact = FALSE),
  error = function(e) NULL)

# Stochastic dominance check: compare CID distributions per rep
# Define: aware stoch dominated by blind iff at every quantile q,
#   awareCID(q) >= blindCID(q)? -- we'll loosely check 25/50/75/95.
percRows <- list()
for (i in seq_len(nReps)) {
  bi <- which(vapply(all, function(r) r$rep == blindDf$rep[i] && r$chain == "blind",
                     logical(1)))[1]
  ai <- which(vapply(all, function(r) r$rep == awareDf$rep[i] && r$chain == "aware",
                     logical(1)))[1]
  percRows[[i]] <- data.frame(
    rep = blindDf$rep[i],
    blindMean   = blindDf$cidMean[i],
    awareMean   = awareDf$cidMean[i],
    deltaMean   = deltaCID[i],
    blindMedian = blindDf$cidMedian[i],
    awareMedian = awareDf$cidMedian[i],
    blindP05    = blindDf$cidP05[i],
    awareP05    = awareDf$cidP05[i],
    blindP95    = blindDf$cidP95[i],
    awareP95    = awareDf$cidP95[i],
    stringsAsFactors = FALSE
  )
}
percDf <- do.call(rbind, percRows)

# ----- markdown report ------------------------------------------------------
mdLines <- c()
add <- function(...) mdLines <<- c(mdLines, sprintf(...))
add("# Aware multirep-v3 posterior shape characterisation")
add("")
add("Date: %s", format(Sys.Date()))
add("")
add("Goal: is the multirep-v3 aware chain doing the methodologically RIGHT thing ")
add("(broad honest uncertainty, including truth in the credibility set) or the WRONG ")
add("thing (failing to find truth at all)?")
add("")
add("Setup: 16 tips, full-clade eco-1 (clades A and B as the 8-tip false eco clade), ")
add("truth = ((A,C),(B,D)); tipBr=%.2f, stemBr=%.2f, rootBr=%.2f, truth TL=%.3f. ",
    TIPBR, STEMBR, ROOTBR, truthTL)
add("Scoring is root-invariant via `inst/ecology/simulations/sim3-scoring.R`.")
add("")
add("## Per-rep x chain posterior summary")
add("")
add("|  rep | chain |    n | P(AC) | P(AB) | CID mean | CID 05 | CID 95 | nUniq | |CS95| | truth in CS? |")
add("|-----:|:------|----:|----:|----:|------:|-----:|-----:|----:|----:|:----|")
for (i in seq_len(nrow(sumDf))) {
  add("|  %2d | %-5s | %4d | %0.3f | %0.3f | %0.3f | %0.3f | %0.3f | %4d | %4d | %s |",
      sumDf$rep[i], sumDf$chain[i], sumDf$n[i],
      sumDf$pAC[i], sumDf$pAB[i],
      sumDf$cidMean[i], sumDf$cidP05[i], sumDf$cidP95[i],
      sumDf$nUnique[i], sumDf$nCS[i],
      if (sumDf$truthInCS[i]) "**YES**" else "no")
}
add("")
add("Notes:")
add("- `P(AC)` = posterior probability that the TRUE sister-of-clades bipartition `((A,C) vs rest)` is present.")
add("- `P(AB)` = posterior probability of the FALSE (eco-driven) bipartition `((A,B) vs rest)`.")
add("- `CID` = normalised Clustering Information Distance to truth, lower is better; 0 = identical.")
add("- `|CS95|` = number of unique canonical topologies in the 95%% credibility set.")
add("")
add("## Aggregate aware-vs-blind comparison")
add("")
add("Mean across the %d reps (corrected, root-invariant):", nReps)
add("")
add("| metric                       |   blind |   aware |  delta (aware - blind) |")
add("|:-----------------------------|--------:|--------:|----------------------:|")
add("| mean P(AC) true              | %0.3f | %0.3f |  %+0.3f |",
    mean(blindDf$pAC), mean(awareDf$pAC), mean(awareDf$pAC) - mean(blindDf$pAC))
add("| mean P(AB) false             | %0.3f | %0.3f |  %+0.3f |",
    mean(blindDf$pAB), mean(awareDf$pAB), mean(awareDf$pAB) - mean(blindDf$pAB))
add("| mean CID-to-truth            | %0.3f | %0.3f |  %+0.3f |",
    mean(blindDf$cidMean), mean(awareDf$cidMean),
    mean(awareDf$cidMean) - mean(blindDf$cidMean))
add("| mean nUnique topologies      | %5.1f | %5.1f |  %+0.1f |",
    mean(blindDf$nUnique), mean(awareDf$nUnique),
    mean(awareDf$nUnique) - mean(blindDf$nUnique))
add("| mean |CS95|                  | %5.1f | %5.1f |  %+0.1f |",
    mean(blindDf$nCS), mean(awareDf$nCS),
    mean(awareDf$nCS) - mean(blindDf$nCS))
add("| reps with truth in CS95      | %d / %d | %d / %d | |",
    sum(blindDf$truthInCS), nReps,
    sum(awareDf$truthInCS), nReps)
add("| reps with truth in posterior | %d / %d | %d / %d | |",
    sum(blindDf$truthInPosterior), nReps,
    sum(awareDf$truthInPosterior), nReps)
add("")
add("Paired comparison of mean CID-to-truth across reps:")
add("- Reps where aware's mean CID is LOWER than blind's: **%d / %d**",
    sum(deltaCID < 0), nReps)
add("- Reps where aware's mean CID is HIGHER than blind's: %d / %d",
    sum(deltaCID > 0), nReps)
if (!is.null(wtest)) {
  add("- Paired Wilcoxon (aware - blind, two-sided): V = %.1f, p = %0.4f",
      unname(wtest$statistic), wtest$p.value)
}
add("")
add("Per-rep CID-to-truth comparison:")
add("")
add("|  rep | blind mean | aware mean | delta | blind median | aware median | blind 95%% | aware 95%% |")
add("|----:|----:|----:|----:|----:|----:|----:|----:|")
for (i in seq_len(nReps)) {
  add("|  %2d | %0.3f | %0.3f | %+0.3f | %0.3f | %0.3f | %0.3f | %0.3f |",
      percDf$rep[i], percDf$blindMean[i], percDf$awareMean[i],
      percDf$deltaMean[i], percDf$blindMedian[i], percDf$awareMedian[i],
      percDf$blindP95[i], percDf$awareP95[i])
}
add("")
add("## Top-5 canonical topologies per rep x chain")
add("")
add("`hasAC`, `hasAB` are unrooted-bipartition presence flags for the top topology of that mass.")
add("")
for (r in all) {
  add("### rep %d / %s — top 5 (of %d unique; |CS95|=%d, CS mass=%.3f)",
      r$rep, r$chain, r$nUnique, r$nCS, r$csMass)
  add("")
  add("| rank | prob  | hasAC | hasAB |   CID |")
  add("|----:|-----:|:---|:---|----:|")
  k <- min(5L, nrow(r$top5))
  for (j in seq_len(k)) {
    add("|  %2d | %0.3f | %s | %s | %0.3f |",
        j, r$top5$prob[j],
        if (r$top5$hasAC[j]) "AC" else ".",
        if (r$top5$hasAB[j]) "AB" else ".",
        r$top5$cid[j])
  }
  add("")
}

add("## Truth-in-posterior summary")
add("")
add("Each chain's posterior is a finite sample; a canonical topology hash either appears in the post-burnin draws or does not.")
add("")
add("| rep | blind: truth in post? | blind: truth prob | aware: truth in post? | aware: truth prob |")
add("|----:|:---|----:|:---|----:|")
for (i in seq_len(nReps)) {
  add("|  %2d | %s | %0.3f | %s | %0.3f |",
      blindDf$rep[i],
      if (blindDf$truthInPosterior[i]) "yes" else "no",
      blindDf$truthProb[i],
      if (awareDf$truthInPosterior[i]) "yes" else "no",
      awareDf$truthProb[i])
}
add("")
add("## Verdict criteria")
add("")
add("Interpretation grid:")
add("- **Broad honest uncertainty**: aware's CID-to-truth distribution is LOWER (closer to truth) than blind's on most reps, AND the credibility set is broader (larger |CS95|), AND truth is occasionally inside CS95.")
add("- **Failing to find truth**: aware's posterior never includes AC and CID-to-truth is no better than blind's.")
add("- **Mode-trapping into wrong neighbourhood**: aware concentrates probability on a NON-AB, NON-AC topology that is itself far from truth, so neither bipartition flag fires but CID is poor.")
add("")
add("See report Section 'Aggregate aware-vs-blind comparison' above for the numeric verdict inputs.")
add("")
add("Outputs:")
add("- `inst/ecology/scripts/aware-multirep-v3-posterior-shape.{rds,csv}`")
add("- `inst/ecology/scripts/aware-multirep-v3-posterior-shape-credset.csv`")
add("- `inst/ecology/scripts/aware-multirep-v3-posterior-shape-plot.pdf`")

writeLines(mdLines, file.path(outDir, "aware-multirep-v3-posterior-shape-report.md"))
cat("\nWrote markdown report.\n")

# ----- plots ----------------------------------------------------------------
pdf(file.path(outDir, "aware-multirep-v3-posterior-shape-plot.pdf"),
    width = 9, height = 6)
op <- par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))

# Boxplot of CID-to-truth per (rep, chain).
# Build a list ordered: rep1-blind, rep1-aware, rep2-blind, rep2-aware, ...
cidList <- list()
labels  <- character()
cols    <- character()
for (i in seq_len(nReps)) {
  r <- blindDf$rep[i]
  bi <- which(vapply(all, function(x) x$rep == r && x$chain == "blind",
                     logical(1)))[1]
  ai <- which(vapply(all, function(x) x$rep == r && x$chain == "aware",
                     logical(1)))[1]
  cidList[[length(cidList) + 1L]] <- all[[bi]]$cidVec
  labels <- c(labels, sprintf("r%02d.B", r))
  cols   <- c(cols, "steelblue")
  cidList[[length(cidList) + 1L]] <- all[[ai]]$cidVec
  labels <- c(labels, sprintf("r%02d.A", r))
  cols   <- c(cols, "firebrick")
}
boxplot(cidList, names = labels, las = 2, col = cols,
        ylab = "Normalised CID to truth (post-burnin)",
        main = "multirep-v3 posterior CID-to-truth, blind (B, blue) vs aware (A, red)")
abline(h = 0, lty = 2, col = "grey50")
legend("topright", legend = c("blind", "aware"), fill = c("steelblue", "firebrick"),
       bty = "n")

# Second page: mean P(AC) and P(AB) per rep x chain
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
mat1 <- rbind(blind = blindDf$pAC, aware = awareDf$pAC)
colnames(mat1) <- paste0("rep", blindDf$rep)
barplot(mat1, beside = TRUE, ylim = c(0, 1),
        col = c("steelblue", "firebrick"),
        main = "P(AC) — TRUE sister, post-burnin",
        ylab = "posterior support",
        legend.text = TRUE, args.legend = list(x = "topright"))
abline(h = 0.5, lty = 2, col = "grey50")

mat2 <- rbind(blind = blindDf$pAB, aware = awareDf$pAB)
colnames(mat2) <- paste0("rep", blindDf$rep)
barplot(mat2, beside = TRUE, ylim = c(0, 1),
        col = c("steelblue", "firebrick"),
        main = "P(AB) — FALSE eco-driven bipartition, post-burnin",
        ylab = "posterior support",
        legend.text = TRUE, args.legend = list(x = "topright"))
abline(h = 0.5, lty = 2, col = "grey50")

par(op)
dev.off()
cat("Wrote plot.pdf\n")

# ----- echo headline numbers -----------------------------------------------
cat("\n=== Headline ===\n")
cat(sprintf("Reps: %d\n", nReps))
cat(sprintf("Blind mean P(AC)=%.3f  P(AB)=%.3f  CID=%.3f  |CS95|=%.1f  truthInCS=%d/%d\n",
            mean(blindDf$pAC), mean(blindDf$pAB), mean(blindDf$cidMean),
            mean(blindDf$nCS), sum(blindDf$truthInCS), nReps))
cat(sprintf("Aware mean P(AC)=%.3f  P(AB)=%.3f  CID=%.3f  |CS95|=%.1f  truthInCS=%d/%d\n",
            mean(awareDf$pAC), mean(awareDf$pAB), mean(awareDf$cidMean),
            mean(awareDf$nCS), sum(awareDf$truthInCS), nReps))
cat(sprintf("Reps where aware has LOWER mean CID-to-truth than blind: %d/%d\n",
            sum(deltaCID < 0), nReps))
if (!is.null(wtest)) {
  cat(sprintf("Paired Wilcoxon p (mean CID aware vs blind): %.4f\n", wtest$p.value))
}
cat(sprintf("P(AC) min/mean/max across aware reps: %.3f / %.3f / %.3f\n",
            min(awareDf$pAC), mean(awareDf$pAC), max(awareDf$pAC)))
cat("\nDone.\n")
