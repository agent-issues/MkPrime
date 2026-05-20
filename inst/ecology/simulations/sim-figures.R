# sim-figures.R ---------------------------------------------------------------
# Headline figures for the ecology-aware simulation study.  Generates
# three PDF panels:
#
#   sim3-headline.pdf  -- Single replicate Sim 3: posterior bipartition
#                         support and CID for blind vs aware.
#   sim3b-negctrl.pdf  -- Sim 3b negative control: blind vs aware on
#                         no-effect data.
#   sim3-multirep.pdf  -- Across-replicates: paired blind/aware
#                         P(true), P(wrong) and CID gap.
#
# Uses base R graphics only (per house style).  Each script expects the
# corresponding *-result.rds / *-mcmc-result.rds files to exist in
# inst/ecology/simulations/.
suppressPackageStartupMessages({
  library("TreeTools")
  library("TreeDist")
})

trueTreeText <- paste0(
  "((((A1:0.5,A2:0.5):0.5,(A3:0.5,A4:0.5):0.5):0.1,",
  "(((C1:0.5,C2:0.5):0.5,(C3:0.5,C4:0.5):0.5):0.1)):0.15,",
  "(((B1:0.5,B2:0.5):0.5,(B3:0.5,B4:0.5):0.5):0.1,",
  "(((D1:0.5,D2:0.5):0.5,(D3:0.5,D4:0.5):0.5):0.1)):0.15);")
wrongTree <- Preorder(ape::read.tree(text = paste0(
  "(((A1,A2),(A3,A4)),((B1,B2),(B3,B4)),",
  "(((C1,C2),(C3,C4)),((D1,D2),(D3,D4))));")))

trueSplitTips  <- c(paste0("A", 1:4), paste0("C", 1:4))
wrongSplitTips <- c(paste0("A", 1:4), paste0("B", 1:4))

# Use root-invariant HasBipartSplits() — see sim3-scoring.R. Old
# inline hasBipart() was root-dependent and biased every figure that
# read it. Note: .PlotMultiRep reads `res$agg` which is pre-baked by
# sim3-multirep.R; regenerate the multirep RDS to refresh those values.
source("inst/ecology/simulations/sim3-scoring.R")
hasBipart <- HasBipartSplits
discard <- function(x) x[seq.int(ceiling(length(x) / 4) + 1L, length(x))]

scoreOne <- function(trees, trueTree) {
  if (!length(trees)) return(NULL)
  class(trees) <- "multiPhylo"
  list(
    pTrue  = mean(hasBipart(trees, trueSplitTips)),
    pWrong = mean(hasBipart(trees, wrongSplitTips)),
    cidTrue  = as.numeric(TreeDist::ClusteringInfoDist(
      trees, trueTree, normalize = TRUE)),
    cidWrong = as.numeric(TreeDist::ClusteringInfoDist(
      trees, wrongTree, normalize = TRUE))
  )
}

#' Side-by-side blind/aware figure for a single replicate.
.PlotSingleRep <- function(resFile, pdfOut, mainTitle) {
  if (!file.exists(resFile)) {
    message("Skipping (no file): ", resFile); return(invisible())
  }
  res <- readRDS(resFile)
  trB <- discard(res$resBlind$trees)
  trA <- discard(res$resAware$trees)
  sB  <- scoreOne(trB, res$tree)
  sA  <- scoreOne(trA, res$tree)
  pdf(pdfOut, width = 9, height = 4)
  on.exit(dev.off())
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  # Panel 1: P(true) vs P(wrong) bars.
  m <- rbind(c(sB$pTrue, sB$pWrong), c(sA$pTrue, sA$pWrong))
  rownames(m) <- c("blind", "aware"); colnames(m) <- c("P(true)", "P(wrong)")
  barplot(t(m), beside = TRUE, ylim = c(0, 1), legend = TRUE,
          col = c("seagreen", "tomato"), ylab = "Posterior probability",
          main = "Bipartition support")
  # Panel 2: CID distributions.
  boxplot(list(blind.true  = sB$cidTrue,  blind.wrong = sB$cidWrong,
               aware.true  = sA$cidTrue,  aware.wrong = sA$cidWrong),
          las = 2, ylab = "CID (lower = closer)",
          col = c("seagreen", "tomato", "seagreen", "tomato"),
          main = "Tree distance distributions")
  mtext(mainTitle, outer = TRUE, cex = 1.2, font = 2)
}

#' Multi-replicate paired comparison figure.
.PlotMultiRep <- function(resFile, pdfOut, mainTitle) {
  if (!file.exists(resFile)) {
    message("Skipping (no file): ", resFile); return(invisible())
  }
  res <- readRDS(resFile)
  agg <- res$agg
  pdf(pdfOut, width = 10, height = 4)
  on.exit(dev.off())
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  pair <- function(b, a, ylab, main, ylim = c(0, 1)) {
    plot(NA, xlim = c(0.8, 2.2), ylim = ylim, xaxt = "n",
         xlab = "", ylab = ylab, main = main)
    axis(1, at = 1:2, labels = c("blind", "aware"))
    for (i in seq_along(b)) {
      lines(c(1, 2), c(b[i], a[i]), col = adjustcolor("gray40", 0.6))
      points(c(1, 2), c(b[i], a[i]), pch = 16,
             col = c("seagreen", "tomato"))
    }
  }
  pair(agg$pTrue_blind, agg$pTrue_aware, "P(true)", "Truth support")
  pair(agg$pWrong_blind, agg$pWrong_aware, "P(wrong)", "Wrong support")
  pair(agg$cidTrue_blind - agg$cidWrong_blind,
       agg$cidTrue_aware - agg$cidWrong_aware,
       "CID(true) - CID(wrong)", "Bias toward wrong",
       ylim = range(c(agg$cidTrue_blind - agg$cidWrong_blind,
                      agg$cidTrue_aware - agg$cidWrong_aware, 0)))
  abline(h = 0, lty = 2, col = "gray60")
  mtext(mainTitle, outer = TRUE, cex = 1.2, font = 2)
}

#' Two-panel v1-vs-v2 head-to-head logL distribution.
#' v2 should close the rate-time ridge: aware logL median should be
#' >= blind logL median (aware nests blind).
.PlotV1V2LogL <- function(v1File, v2File, pdfOut, mainTitle) {
  if (!file.exists(v1File) || !file.exists(v2File)) {
    message("Skipping (no file): ", v1File, " or ", v2File)
    return(invisible())
  }
  v1 <- readRDS(v1File); v2 <- readRDS(v2File)
  discardSamp <- function(x) x[seq.int(ceiling(nrow(x) / 4) + 1L, nrow(x)), ]
  l1B <- v1$resBlind$samples[, "log_likelihood"]
  l1A <- v1$resAware$samples[, "log_likelihood"]
  l2B <- v2$samplesBlind[, "log_likelihood"]
  l2A <- v2$samplesAware[, "log_likelihood"]
  pdf(pdfOut, width = 7, height = 4)
  on.exit(dev.off())
  par(mar = c(4, 4, 3, 1))
  boxplot(list(`v1 blind` = l1B, `v1 aware` = l1A,
               `v2 blind` = l2B, `v2 aware` = l2A),
          las = 2, ylab = "log-likelihood",
          col = c("seagreen", "tomato", "seagreen", "tomato"),
          main = mainTitle)
  abline(v = 2.5, lty = 2, col = "gray60")
  legend("bottomright",
         legend = c(sprintf("v1 aware-blind gap: %.0f", median(l1A) - median(l1B)),
                    sprintf("v2 aware-blind gap: %.0f", median(l2A) - median(l2B))),
         bty = "n", cex = 0.85)
}

outDir <- "inst/ecology/simulations"
.PlotSingleRep(file.path(outDir, "sim3-mcmc-result.rds"),
               file.path(outDir, "sim3-headline.pdf"),
               "Sim 3: convergent ecology, single replicate")
.PlotSingleRep(file.path(outDir, "sim3b-mcmc-result.rds"),
               file.path(outDir, "sim3b-negctrl.pdf"),
               "Sim 3b: negative control (no ecology effect)")
.PlotV1V2LogL(file.path(outDir, "sim3-mcmc-result.rds"),
              file.path(outDir, "sim3-v2-result.rds"),
              file.path(outDir, "sim3-v1v2-logL.pdf"),
              "Sim 3 log-likelihood: v1 vs v2")
for (n in c(4, 20)) {
  f <- file.path(outDir, sprintf("sim3-multirep-n%d.rds", n))
  if (file.exists(f)) {
    .PlotMultiRep(f, file.path(outDir, sprintf("sim3-multirep-n%d.pdf", n)),
                  sprintf("Sim 3: %d replicates, blind vs aware", n))
  }
}
cat("Done.\n")
