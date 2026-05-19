# rescore-multirep-v3.R --------------------------------------------------------
#
# Re-score every saved multirep-v3 (Sim 3 multirep) blind + aware MCMC
# result under the root-invariant scorer (HasBipartSplits /
# ScoreTreesUnrooted) in inst/simulations/ecology/sim3-scoring.R.
#
# Motivation. The legacy prop.part-based scorer used by run_rep.R is
# rooted: when an MCMC sample's root sits inside a target tip-set, the
# bipartition appears in prop.part's output as its complement and is
# missed. Earlier audits showed this artefact hides what is in fact a
# strong false-eco-clade signal: rep 05 BLIND legacy P(AB) = 0.000 vs
# corrected P(AB) = 0.978. This driver rescores all 8 reps x 2 chains so
# we can decide whether multirep-v3 already constitutes the publication-
# quality "model dichotomy" demo.
#
# Inputs (one per rep x chain):
#   inst/simulations/ecology/multirep-v3-results/rep%02d/{blind,aware}-result.rds
#
# Outputs:
#   inst/scripts/rescore-multirep-v3-results.rds  (tidy tibble of all scores)
#   inst/scripts/rescore-multirep-v3-results.csv  (same, human-readable)
#   inst/scripts/rescore-multirep-v3-summary.csv  (one row per rep x chain,
#                                                  pivoted for readability)
#   inst/scripts/rescore-multirep-v3-report.md    (narrative summary)
#   inst/scripts/rescore-multirep-v3-plot.pdf     (per-rep blind vs aware
#                                                  P(AC) + P(AB) bar chart)
#
# Run locally from the MkPrime repo root with:
#   Rscript inst/scripts/rescore-multirep-v3.R

suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKP_LIB", "")
  if (nzchar(libPath)) .libPaths(c(libPath, .libPaths()))
  library("TreeTools")
  library("TreeDist")
})

mkpRoot <- Sys.getenv("MKP_REPO_ROOT", getwd())
source(file.path(mkpRoot, "inst/simulations/ecology/sim3-scoring.R"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3-helpers.R"))

# ----- multirep-v3 truth config (matches inst/hamilton/sim3-multirep-v3/run_rep.R)
TIPBR  <- 0.5
STEMBR <- 0.30
ROOTBR <- 0.15
truthTree <- .BuildConvergentTree(tipBranch  = TIPBR,
                                  stemBranch = STEMBR,
                                  rootBranch = ROOTBR)
truthTL <- sum(truthTree$edge.length)

# ----- bipartitions to score (root-invariant; all unrooted)
biparts <- list(
  trueSister_AC  = c(paste0("A", 1:4), paste0("C", 1:4)),
  trueSister_BD  = c(paste0("B", 1:4), paste0("D", 1:4)),
  falseSister_AB = c(paste0("A", 1:4), paste0("B", 1:4)),
  falseSister_CD = c(paste0("C", 1:4), paste0("D", 1:4)),
  trueClade_A    = paste0("A", 1:4),
  trueClade_B    = paste0("B", 1:4),
  trueClade_C    = paste0("C", 1:4),
  trueClade_D    = paste0("D", 1:4)
)

# ----- canonical-splits hash (root-invariant) for unique-topology counts
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

# ----- per-rep x chain processor
processChain <- function(repId, chain, repDir) {
  rdsPath <- file.path(repDir, sprintf("%s-result.rds", chain))
  if (!file.exists(rdsPath)) {
    return(NULL)
  }
  res <- readRDS(rdsPath)
  trees <- res$trees
  if (is.null(trees) || length(trees) == 0L) return(NULL)
  nTreesTotal <- length(trees)

  # 25% burn-in (same as run_rep.R::discard)
  trBurned <- .DiscardBurnin(trees)
  class(trBurned) <- "multiPhylo"
  nKept <- length(trBurned)

  scores <- ScoreTreesUnrooted(trees, biparts = biparts,
                               refTree = truthTree, discardBurnin = TRUE)

  # canonical-topology stats
  canHashes <- vapply(trBurned, canonicalSplitsHash, character(1))
  nUnique <- length(unique(canHashes))
  topTab <- sort(table(canHashes), decreasing = TRUE)
  topShare <- if (length(topTab)) unname(topTab[1] / length(canHashes)) else NA_real_

  # TL (root-invariant). Prefer the samples table if available; otherwise
  # compute from the saved trees themselves (which is the same quantity).
  sdf <- tryCatch(as.data.frame(res$samples), error = function(e) NULL)
  if (!is.null(sdf) && "tree_length" %in% colnames(sdf) && nrow(sdf) > 0L) {
    tl <- sdf$tree_length
    keepFrom <- ceiling(length(tl) / 4) + 1L
    tlKeep <- if (keepFrom > length(tl)) numeric(0) else tl[keepFrom:length(tl)]
  } else {
    tlKeep <- vapply(trBurned, function(tr) sum(tr$edge.length), numeric(1))
  }
  tlMean <- if (length(tlKeep)) mean(tlKeep) else NA_real_
  tlMed  <- if (length(tlKeep)) median(tlKeep) else NA_real_
  tlMin  <- if (length(tlKeep)) min(tlKeep) else NA_real_
  tlMax  <- if (length(tlKeep)) max(tlKeep) else NA_real_

  scores$rep      <- repId
  scores$chain    <- chain
  scores$nTreesTotal <- nTreesTotal
  scores$nUnique  <- nUnique
  scores$topShare <- topShare
  scores$tlMean   <- tlMean
  scores$tlMed    <- tlMed
  scores$tlMin    <- tlMin
  scores$tlMax    <- tlMax
  scores$truthTL  <- truthTL
  scores
}

cat("=== rescore-multirep-v3.R ===\n")
cat("mkpRoot:", mkpRoot, "\n")
cat(sprintf("Truth: 16 tips, TL = %.3f, true sister = (A,C),(B,D); false (eco) sister = (A,B)\n",
            truthTL))

resultsDir <- file.path(mkpRoot, "inst/simulations/ecology/multirep-v3-results")
reps <- sort(list.files(resultsDir, pattern = "^rep[0-9]+$"))
cat(sprintf("Found %d rep dirs under %s\n", length(reps), resultsDir))

allRows <- list()
for (rep in reps) {
  repDir <- file.path(resultsDir, rep)
  repNum <- as.integer(sub("rep", "", rep))
  for (chain in c("blind", "aware")) {
    t0 <- Sys.time()
    row <- tryCatch(processChain(repNum, chain, repDir),
                    error = function(e) {
                      cat(sprintf("[%s %s] ERROR: %s\n", rep, chain, conditionMessage(e)))
                      NULL
                    })
    if (is.null(row)) {
      cat(sprintf("[%s %s] SKIP (missing or empty)\n", rep, chain))
      next
    }
    dt <- format(round(Sys.time() - t0, 1))
    cat(sprintf("[%s %s] OK  nTrees=%d  nKept=%d  P(AC)=%.3f  P(AB)=%.3f  TLmean=%.2f  %s\n",
                rep, chain, row$nTreesTotal[1], row$nTrees[1],
                row$value_corrected[row$metric == "trueSister_AC"],
                row$value_corrected[row$metric == "falseSister_AB"],
                row$tlMean[1], dt))
    allRows[[length(allRows) + 1L]] <- row
  }
}

if (!length(allRows)) {
  stop("No rep results found — aborting.")
}
scoresDf <- do.call(rbind, allRows)

# Pivot to one row per rep x chain (wide)
metrics <- unique(scoresDf$metric)
metaCols <- c("rep", "chain", "nTreesTotal", "nTrees", "nUnique", "topShare",
              "tlMean", "tlMed", "tlMin", "tlMax", "truthTL")
wideRows <- list()
keys <- unique(scoresDf[, c("rep", "chain")])
keys <- keys[order(keys$rep, keys$chain), ]
for (i in seq_len(nrow(keys))) {
  r <- keys$rep[i]; ch <- keys$chain[i]
  sub <- scoresDf[scoresDf$rep == r & scoresDf$chain == ch, , drop = FALSE]
  rowOut <- as.list(sub[1, metaCols, drop = FALSE])
  for (m in metrics) {
    mrow <- sub[sub$metric == m, , drop = FALSE]
    rowOut[[paste0("corr_", m)]] <- if (nrow(mrow)) mrow$value_corrected[1] else NA_real_
    rowOut[[paste0("leg_",  m)]] <- if (nrow(mrow)) mrow$value_legacy[1]    else NA_real_
  }
  wideRows[[i]] <- as.data.frame(rowOut, stringsAsFactors = FALSE)
}
wideDf <- do.call(rbind, wideRows)

outDir <- file.path(mkpRoot, "inst/scripts")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
saveRDS(scoresDf, file.path(outDir, "rescore-multirep-v3-results.rds"))
write.csv(scoresDf, file.path(outDir, "rescore-multirep-v3-results.csv"),
          row.names = FALSE)
write.csv(wideDf,   file.path(outDir, "rescore-multirep-v3-summary.csv"),
          row.names = FALSE)
cat("\nWrote results, csv, summary.\n")

# ----- per-rep blind-vs-aware report ------------------------------------------
mkSubset <- function(df, ch) df[df$chain == ch, , drop = FALSE]
blindW <- mkSubset(wideDf, "blind")
awareW <- mkSubset(wideDf, "aware")
blindW <- blindW[order(blindW$rep), , drop = FALSE]
awareW <- awareW[order(awareW$rep), , drop = FALSE]

pAC_blind  <- blindW$corr_trueSister_AC
pAC_aware  <- awareW$corr_trueSister_AC
pAB_blind  <- blindW$corr_falseSister_AB
pAB_aware  <- awareW$corr_falseSister_AB
pBD_blind  <- blindW$corr_trueSister_BD
pBD_aware  <- awareW$corr_trueSister_BD
pCD_blind  <- blindW$corr_falseSister_CD
pCD_aware  <- awareW$corr_falseSister_CD

repsVec <- blindW$rep
nReps   <- length(repsVec)

# Headline counts (interpretive thresholds borrowed from the audit brief)
nBlindAB_gt5 <- sum(pAB_blind > 0.5, na.rm = TRUE)
nBlindAB_gt9 <- sum(pAB_blind > 0.9, na.rm = TRUE)
nAwareAC_gt5 <- sum(pAC_aware > 0.5, na.rm = TRUE)
nAwareAC_gt9 <- sum(pAC_aware > 0.9, na.rm = TRUE)
nAwareModeTrap <- sum(pAC_aware < 0.5, na.rm = TRUE)
nDichotomy <- sum(pAB_blind > 0.5 & pAC_aware > 0.5, na.rm = TRUE)
nDichotomyStrong <- sum(pAB_blind > 0.9 & pAC_aware > 0.9, na.rm = TRUE)
# reversed pattern: aware false-supports while blind correct
nReversed <- sum(pAB_aware > 0.5 & pAC_blind > 0.5, na.rm = TRUE)

# Build a markdown report
mdLines <- c()
add <- function(...) mdLines <<- c(mdLines, sprintf(...))
add("# Rescore: multirep-v3 (Sim 3 multirep) — 8 reps under root-invariant scoring")
add("")
add("Date: %s", format(Sys.Date()))
add("")
add("Setup: 16 tips arranged as `((A,C),(B,D))`, each clade balanced quartet, ")
add("eco-1 = ALL of clades A and B (full-clade convergent ecology), 480 chars ")
add("(120 neomorphic + 360 transformational), phi = 4, pi0 = 0.75, theta = 1, ")
add("tipBr = %.2f, stemBr = %.2f, rootBr = %.2f, truth TL = %.3f.", TIPBR, STEMBR, ROOTBR, truthTL)
add("Scoring uses `inst/simulations/ecology/sim3-scoring.R::ScoreTreesUnrooted` ")
add("(root-invariant) on the post-burnin 75%% of each chain's saved trees.")
add("")
add("## Per-rep table — corrected support")
add("")
add("|  rep | chain |  nKept | P(AC)  | P(BD)  | P(AB)  | P(CD)  | TLmean | nUniq | topShare |")
add("|-----:|:------|------:|-------:|-------:|-------:|-------:|-------:|------:|--------:|")
for (i in seq_len(nReps)) {
  for (ch in c("blind", "aware")) {
    sub <- wideDf[wideDf$rep == repsVec[i] & wideDf$chain == ch, , drop = FALSE]
    if (!nrow(sub)) {
      add("|  %2d | %-5s |     - |     - |     - |     - |     - |     - |     - |     - |",
          repsVec[i], ch)
      next
    }
    add("|  %2d | %-5s |  %4d | %0.3f | %0.3f | **%0.3f** | %0.3f | %5.2f |  %4d | %0.3f |",
        repsVec[i], ch, sub$nTrees[1],
        sub$corr_trueSister_AC[1], sub$corr_trueSister_BD[1],
        sub$corr_falseSister_AB[1], sub$corr_falseSister_CD[1],
        sub$tlMean[1], sub$nUnique[1], sub$topShare[1])
  }
}
add("")

add("## Mean across reps (corrected)")
add("")
add("| metric          |  blind |  aware |")
add("|:----------------|------:|------:|")
add("| P(AC) true      | %0.3f | %0.3f |",
    mean(pAC_blind, na.rm = TRUE), mean(pAC_aware, na.rm = TRUE))
add("| P(BD) true      | %0.3f | %0.3f |",
    mean(pBD_blind, na.rm = TRUE), mean(pBD_aware, na.rm = TRUE))
add("| **P(AB) false** | **%0.3f** | **%0.3f** |",
    mean(pAB_blind, na.rm = TRUE), mean(pAB_aware, na.rm = TRUE))
add("| P(CD) false     | %0.3f | %0.3f |",
    mean(pCD_blind, na.rm = TRUE), mean(pCD_aware, na.rm = TRUE))
add("| TL mean         | %0.2f | %0.2f | (truth %0.2f)",
    mean(blindW$tlMean, na.rm = TRUE),
    mean(awareW$tlMean, na.rm = TRUE), truthTL)
add("| nUnique topo    | %0.1f | %0.1f |",
    mean(blindW$nUnique, na.rm = TRUE),
    mean(awareW$nUnique, na.rm = TRUE))
add("")

add("## Headline counts (out of %d reps)", nReps)
add("")
add("- Blind P(AB) > 0.5: **%d / %d**", nBlindAB_gt5, nReps)
add("- Blind P(AB) > 0.9: **%d / %d**", nBlindAB_gt9, nReps)
add("- Aware P(AC) > 0.5: **%d / %d**", nAwareAC_gt5, nReps)
add("- Aware P(AC) > 0.9: %d / %d", nAwareAC_gt9, nReps)
add("- Aware mode-trapped (P(AC) < 0.5): %d / %d", nAwareModeTrap, nReps)
add("- Model dichotomy (blind P(AB)>0.5 AND aware P(AC)>0.5): **%d / %d**",
    nDichotomy, nReps)
add("- Strong dichotomy (both > 0.9): %d / %d", nDichotomyStrong, nReps)
add("- Reversed (aware false-supports AB > 0.5 AND blind correct AC > 0.5): %d / %d",
    nReversed, nReps)
add("")

add("## Per-rep blind-vs-aware delta on the headline metrics")
add("")
add("|  rep | P(AC) blind | P(AC) aware | P(AB) blind | P(AB) aware |  TL blind |  TL aware | dichotomy? |")
add("|-----:|----:|----:|----:|----:|----:|----:|:---|")
for (i in seq_len(nReps)) {
  dic <- (pAB_blind[i] > 0.5) && (pAC_aware[i] > 0.5)
  add("|  %2d | %0.3f | %0.3f | %0.3f | %0.3f | %5.2f | %5.2f | %s |",
      repsVec[i], pAC_blind[i], pAC_aware[i], pAB_blind[i], pAB_aware[i],
      blindW$tlMean[i], awareW$tlMean[i],
      if (isTRUE(dic)) "**YES**" else "no")
}
add("")
add("(`dichotomy?` = blind supports the false eco-clade AB AND aware supports the true AC.)")
add("")
add("Outputs:")
add("- `inst/scripts/rescore-multirep-v3-results.{rds,csv}` (long form)")
add("- `inst/scripts/rescore-multirep-v3-summary.csv` (wide form)")
add("- `inst/scripts/rescore-multirep-v3-plot.pdf` (bar chart)")

writeLines(mdLines, file.path(outDir, "rescore-multirep-v3-report.md"))
cat("Wrote report.md\n")

# ----- bar chart -------------------------------------------------------------
pdf(file.path(outDir, "rescore-multirep-v3-plot.pdf"), width = 8, height = 5)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
barCols <- c("steelblue", "firebrick")
mat1 <- rbind(blind = pAC_blind, aware = pAC_aware)
colnames(mat1) <- paste0("rep", repsVec)
barplot(mat1, beside = TRUE, ylim = c(0, 1), col = barCols,
        main = "Corrected P(AC) — TRUE sister",
        xlab = "rep", ylab = "posterior support",
        legend.text = TRUE, args.legend = list(x = "topright"))
abline(h = 0.5, lty = 2, col = "grey50")

mat2 <- rbind(blind = pAB_blind, aware = pAB_aware)
colnames(mat2) <- paste0("rep", repsVec)
barplot(mat2, beside = TRUE, ylim = c(0, 1), col = barCols,
        main = "Corrected P(AB) — FALSE eco sister",
        xlab = "rep", ylab = "posterior support",
        legend.text = TRUE, args.legend = list(x = "topright"))
abline(h = 0.5, lty = 2, col = "grey50")
par(op)
dev.off()
cat("Wrote plot.pdf\n")

cat("\n=== Headline ===\n")
cat(sprintf("Blind  mean P(AB false) = %.3f   (>0.5 in %d/%d reps; >0.9 in %d/%d)\n",
            mean(pAB_blind, na.rm = TRUE), nBlindAB_gt5, nReps, nBlindAB_gt9, nReps))
cat(sprintf("Aware  mean P(AC true)  = %.3f   (>0.5 in %d/%d reps)\n",
            mean(pAC_aware, na.rm = TRUE), nAwareAC_gt5, nReps))
cat(sprintf("Model dichotomy demonstrated (blind AB>0.5 AND aware AC>0.5): %d / %d reps\n",
            nDichotomy, nReps))
cat("\nDone.\n")
