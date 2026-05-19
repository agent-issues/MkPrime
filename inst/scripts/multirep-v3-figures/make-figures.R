# make-figures.R ---------------------------------------------------------------
#
# Publication-quality figures summarising the multirep-v3 simulation results.
#
# Three figures:
#   Fig 1 — Per-rep false-clade (AB) vs true-clade (AC) posterior support
#   Fig 2 — CID-to-truth, paired blind vs aware, across 8 reps
#   Fig 3 — 2-D MDS embedding for rep 02
#
# Run from the MkPrime repo root with:
#   Rscript inst/scripts/multirep-v3-figures/make-figures.R
#
# Outputs (inst/scripts/multirep-v3-figures/):
#   fig1-false-clade-support.{pdf,png}
#   fig2-cid-to-truth.{pdf,png}
#   fig3-mds-rep02.{pdf,png}

set.seed(42L)

suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKP_LIB", "")
  if (nzchar(libPath)) .libPaths(c(libPath, .libPaths()))
  library("TreeTools")
  library("TreeDist")
})

mkpRoot <- Sys.getenv("MKP_REPO_ROOT", normalizePath(getwd()))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3-scoring.R"))
source(file.path(mkpRoot, "inst/simulations/ecology/sim3-helpers.R"))

outDir <- file.path(mkpRoot, "inst/scripts/multirep-v3-figures")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Palette: Okabe-Ito (colourblind-safe)
# ---------------------------------------------------------------------------
COL_BLIND <- "#0072B2"   # blue
COL_AWARE <- "#D55E00"   # vermillion
COL_REF   <- "#009E73"   # green (for true-clade ticks)
ALPHA_PT  <- 0.85

# ---------------------------------------------------------------------------
# Load summary CSV (covers Figs 1 + 2)
# ---------------------------------------------------------------------------
csvPath <- file.path(mkpRoot, "inst/scripts/aware-multirep-v3-posterior-shape.csv")
sumDf <- read.csv(csvPath, stringsAsFactors = FALSE)

blindDf <- sumDf[sumDf$chain == "blind", ]
awareDf <- sumDf[sumDf$chain == "aware", ]
blindDf <- blindDf[order(blindDf$rep), ]
awareDf <- awareDf[order(awareDf$rep), ]
nReps <- nrow(blindDf)   # 8
reps  <- blindDf$rep     # 1..8

# Verify Wilcoxon (paper claims V=36, p=0.014)
wtest <- wilcox.test(awareDf$cidMean, blindDf$cidMean,
                     paired = TRUE, alternative = "two.sided", exact = FALSE)
cat(sprintf("[Wilcoxon] V = %.0f, p = %.4f\n", wtest$statistic, wtest$p.value))
# V = 0 means all 8 aware medians are lower than paired blind (consistent direction)
pAnnot <- if (wtest$p.value < 0.001) {
  sprintf("Paired Wilcoxon V = %.0f, p < 0.001", wtest$statistic)
} else {
  sprintf("Paired Wilcoxon V = %.0f, p = %.3f", wtest$statistic, wtest$p.value)
}

# ============================================================================
# Figure 1: per-rep false-clade vs true-clade support
# ============================================================================
fig1_core <- function() {
  xGap   <- 0.35    # gap between paired groups
  bw     <- 0.35    # bar width
  nReps  <- length(reps)
  # Group centres: 1, 2, ..., 8
  centres <- seq_len(nReps)
  xBlind <- centres - bw / 2
  xAware <- centres + bw / 2

  par(family = "sans", cex = 0.85, cex.axis = 0.85, cex.lab = 0.9,
      mar = c(4.5, 4.5, 1.5, 1.5), mgp = c(3, 0.6, 0))

  # Draw frame manually (no outer box)
  plot.new()
  plot.window(xlim = c(0.5, nReps + 0.5), ylim = c(0, 1))

  # Reference line at 0.5
  abline(h = 0.5, lty = 2, col = "grey60", lwd = 0.9)

  # Bars for BLIND P(AB)
  rect(xBlind - bw / 2, 0, xBlind + bw / 2,
       blindDf$pAB,
       col = adjustcolor(COL_BLIND, alpha.f = 0.85),
       border = NA)

  # Bars for AWARE P(AB)
  rect(xAware - bw / 2, 0, xAware + bw / 2,
       awareDf$pAB,
       col = adjustcolor(COL_AWARE, alpha.f = 0.85),
       border = NA)

  # True-clade P(AC) ticks: small horizontal line at height pAC
  tick_half <- bw * 0.55
  # Blind pAC
  segments(xBlind - tick_half, blindDf$pAC,
           xBlind + tick_half, blindDf$pAC,
           col = COL_REF, lwd = 2.5)
  # Aware pAC
  segments(xAware - tick_half, awareDf$pAC,
           xAware + tick_half, awareDf$pAC,
           col = COL_REF, lwd = 2.5)

  # Axes
  axis(1, at = centres, labels = paste0("rep ", sprintf("%02d", reps)),
       tick = FALSE, line = -0.5, cex.axis = 0.8)
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1.0),
       labels = c("0", "0.25", "0.50", "0.75", "1.00"),
       las = 1, lwd = 0.7, tck = -0.025)

  # Labels
  title(xlab = "Replicate", ylab = "Posterior probability", line = 3.0)

  # Commitment threshold label
  text(nReps + 0.45, 0.5, "0.5", adj = c(1, -0.4), col = "grey45",
       cex = 0.75, font = 3)

  # Arrow annotation
  arrows(x0 = 0.55, y0 = 0.78, x1 = 0.55, y1 = 0.56,
         length = 0.08, angle = 20, lwd = 1.0, col = "grey30")
  text(0.7, 0.82, "Committed\nto wrong grouping", adj = c(0, 0.5),
       cex = 0.72, col = "grey30")

  # Legend
  legend("topright",
         legend = c("Blind Mk’ P(AB false clade)",
                    "Aware Mk’ P(AB false clade)",
                    "P(AC true clade)"),
         fill   = c(adjustcolor(COL_BLIND, 0.85),
                    adjustcolor(COL_AWARE, 0.85),
                    NA),
         border = c(NA, NA, NA),
         lty    = c(NA, NA, 1),
         lwd    = c(NA, NA, 2.5),
         col    = c(NA, NA, COL_REF),
         seg.len = 0.9,
         bty = "n", cex = 0.78, x.intersp = 0.6, y.intersp = 1.1)
}

# --- PDF ---
pdf(file.path(outDir, "fig1-false-clade-support.pdf"), width = 7, height = 5,
    useDingbats = FALSE)
fig1_core()
dev.off()

# --- PNG ---
png(file.path(outDir, "fig1-false-clade-support.png"),
    width = 7, height = 5, units = "in", res = 300)
fig1_core()
dev.off()

cat("Fig 1 written.\n")

# ============================================================================
# Figure 2: CID-to-truth, paired across reps
# ============================================================================
fig2_core <- function() {
  par(family = "sans", cex = 0.85, cex.axis = 0.85, cex.lab = 0.9,
      mar = c(4.5, 4.5, 1.5, 1.5), mgp = c(3, 0.6, 0))

  ymin <- 0
  ymax <- max(blindDf$cidP95, awareDf$cidP95) * 1.05

  # Group positions: blind at x - 0.18, aware at x + 0.18
  xBlind <- seq_len(nReps) - 0.18
  xAware <- seq_len(nReps) + 0.18

  plot.new()
  plot.window(xlim = c(0.5, nReps + 0.5), ylim = c(ymin, ymax))

  # Connecting lines (per-rep pairing)
  segments(xBlind, blindDf$cidMedian,
           xAware, awareDf$cidMedian,
           col = "grey70", lwd = 0.9)

  # Error bars: 5th–95th percentile
  eps <- 0.05   # horizontal half-width of cap
  # Blind
  segments(xBlind, blindDf$cidP05, xBlind, blindDf$cidP95,
           col = COL_BLIND, lwd = 1.4)
  segments(xBlind - eps, blindDf$cidP05,
           xBlind + eps, blindDf$cidP05, col = COL_BLIND, lwd = 1.4)
  segments(xBlind - eps, blindDf$cidP95,
           xBlind + eps, blindDf$cidP95, col = COL_BLIND, lwd = 1.4)
  # Aware
  segments(xAware, awareDf$cidP05, xAware, awareDf$cidP95,
           col = COL_AWARE, lwd = 1.4)
  segments(xAware - eps, awareDf$cidP05,
           xAware + eps, awareDf$cidP05, col = COL_AWARE, lwd = 1.4)
  segments(xAware - eps, awareDf$cidP95,
           xAware + eps, awareDf$cidP95, col = COL_AWARE, lwd = 1.4)

  # Median points
  points(xBlind, blindDf$cidMedian, pch = 21,
         bg = COL_BLIND, col = "white", cex = 1.25, lwd = 0.8)
  points(xAware, awareDf$cidMedian, pch = 21,
         bg = COL_AWARE, col = "white", cex = 1.25, lwd = 0.8)

  # Axes
  axis(1, at = seq_len(nReps),
       labels = paste0("rep ", sprintf("%02d", reps)),
       tick = FALSE, line = -0.5, cex.axis = 0.8)
  axis(2, las = 1, lwd = 0.7, tck = -0.025)

  # Labels
  title(xlab = "Replicate",
        ylab = "Normalised CID to truth (median, 5–95%)", line = 3.0)

  # Arrow: lower = closer to truth
  arrows(x0 = 0.56, y0 = ymax * 0.24, x1 = 0.56, y1 = ymax * 0.08,
         length = 0.08, angle = 20, lwd = 1.0, col = "grey30")
  text(0.7, ymax * 0.24, "Closer\nto truth", adj = c(0, 1.0),
       cex = 0.72, col = "grey30")

  # Wilcoxon annotation — bottom-left to avoid legend overlap
  text(0.55, ymin + ymax * 0.04, pAnnot,
       adj = c(0, 0), cex = 0.72, col = "grey25", font = 3)

  # Legend
  legend("topright",
         legend = c("Blind Mk’", "Aware Mk’"),
         pch    = 21,
         pt.bg  = c(COL_BLIND, COL_AWARE),
         col    = "white",
         pt.cex = 1.3,
         bty    = "n", cex = 0.82, x.intersp = 0.8)
}

pdf(file.path(outDir, "fig2-cid-to-truth.pdf"), width = 6, height = 5,
    useDingbats = FALSE)
fig2_core()
dev.off()

png(file.path(outDir, "fig2-cid-to-truth.png"),
    width = 6, height = 5, units = "in", res = 300)
fig2_core()
dev.off()

cat("Fig 2 written.\n")

# ============================================================================
# Figure 3: 2-D MDS embedding for rep 02
# ============================================================================
REP_MDS   <- 2L
N_SUBSAMP <- 150L  # per chain

# Build truth tree (same parameters used in run_rep.R)
TIPBR  <- 0.5
STEMBR <- 0.30
ROOTBR <- 0.15
truthTree <- .BuildConvergentTree(tipBranch  = TIPBR,
                                  stemBranch = STEMBR,
                                  rootBranch = ROOTBR)

# Load trees for rep 02
repDir <- file.path(mkpRoot, "inst/simulations/ecology/multirep-v3-results",
                    sprintf("rep%02d", REP_MDS))
blindRes <- readRDS(file.path(repDir, "blind-result.rds"))
awareRes <- readRDS(file.path(repDir, "aware-result.rds"))

# Discard burn-in (first 25%)
blindTrees <- .DiscardBurnin(blindRes$trees)
awareTrees <- .DiscardBurnin(awareRes$trees)
class(blindTrees) <- "multiPhylo"
class(awareTrees) <- "multiPhylo"

# Sub-sample
nB <- min(N_SUBSAMP, length(blindTrees))
nA <- min(N_SUBSAMP, length(awareTrees))
idxB <- sort(sample(length(blindTrees), nB))
idxA <- sort(sample(length(awareTrees), nA))
blindSub <- blindTrees[idxB]
awareSub <- awareTrees[idxA]

# Combine for a single distance matrix (truth is index nB+nA+1).
# Assemble via pre-allocated vector to avoid c.multiPhylo dispatch issues
# when the subsample is a plain list of phylo.
allTrees <- vector("list", nB + nA + 1L)
for (i in seq_len(nB)) allTrees[[i]] <- blindSub[[i]]
for (i in seq_len(nA)) allTrees[[nB + i]] <- awareSub[[i]]
allTrees[[nB + nA + 1L]] <- truthTree
class(allTrees) <- "multiPhylo"
nAll <- length(allTrees)

cat(sprintf("[Fig 3] Computing %d x %d CID matrix...\n", nAll, nAll))
dMat <- as.matrix(TreeDist::ClusteringInfoDistance(allTrees, normalize = TRUE))
cat("[Fig 3] cmdscale...\n")
mds  <- cmdscale(dMat, k = 2L)

xB <- mds[seq_len(nB), 1L]
yB <- mds[seq_len(nB), 2L]
xA <- mds[nB + seq_len(nA), 1L]
yA <- mds[nB + seq_len(nA), 2L]
xT <- mds[nAll, 1L]
yT <- mds[nAll, 2L]

# 50% convex hull (via density contour approximation using chull on the
# trimmed 50% core: drop points outside 50% quantile ellipse, then take hull)
# -- simpler and more interpretable than KDE contours for this n
convhull_pts <- function(x, y) {
  hull <- chull(x, y)
  list(x = x[hull], y = y[hull])
}

# 50% inner hull: use the nearest half of points to the centroid
inner_hull <- function(x, y, frac = 0.5) {
  cx <- mean(x)
  cy <- mean(y)
  d2 <- (x - cx)^2 + (y - cy)^2
  keep <- which(d2 <= quantile(d2, frac))
  if (length(keep) < 3L) keep <- seq_along(x)
  hull <- chull(x[keep], y[keep])
  list(x = x[keep][hull], y = y[keep][hull])
}

hullB_full <- convhull_pts(xB, yB)
hullA_full <- convhull_pts(xA, yA)
hullB_50   <- inner_hull(xB, yB)
hullA_50   <- inner_hull(xA, yA)

fig3_core <- function() {
  xlim <- range(c(xB, xA, xT)) * c(1.08, 1.08)
  ylim <- range(c(yB, yA, yT)) * c(1.08, 1.08)

  par(family = "sans", cex = 0.85, cex.axis = 0.85, cex.lab = 0.9,
      mar = c(4.5, 4.5, 1.5, 1.5), mgp = c(3, 0.6, 0))

  plot.new()
  plot.window(xlim = xlim, ylim = ylim, asp = 1)

  # Full convex hulls (light fill, no border)
  polygon(hullB_full$x, hullB_full$y,
          col = adjustcolor(COL_BLIND, alpha.f = 0.06), border = NA)
  polygon(hullA_full$x, hullA_full$y,
          col = adjustcolor(COL_AWARE, alpha.f = 0.06), border = NA)

  # 50% inner hull outlines
  polygon(hullB_50$x, hullB_50$y,
          col = adjustcolor(COL_BLIND, alpha.f = 0.12),
          border = COL_BLIND, lty = 2, lwd = 0.9)
  polygon(hullA_50$x, hullA_50$y,
          col = adjustcolor(COL_AWARE, alpha.f = 0.12),
          border = COL_AWARE, lty = 2, lwd = 0.9)

  # Points
  points(xB, yB, pch = 16, col = adjustcolor(COL_BLIND, alpha.f = 0.35),
         cex = 0.55)
  points(xA, yA, pch = 16, col = adjustcolor(COL_AWARE, alpha.f = 0.35),
         cex = 0.55)

  # Truth star
  points(xT, yT, pch = 8, col = "#E69F00", cex = 2.2, lwd = 2.0)

  # Axes
  axis(1, lwd = 0.7, tck = -0.025)
  axis(2, las = 1, lwd = 0.7, tck = -0.025)
  title(xlab = "MDS dimension 1", ylab = "MDS dimension 2", line = 3.0)

  # Legend
  legend("bottomright",
         legend = c(sprintf("Blind Mk’ (n = %d)", nB),
                    sprintf("Aware Mk’ (n = %d)", nA),
                    "True tree"),
         pch   = c(16, 16, 8),
         col   = c(adjustcolor(COL_BLIND, 0.7),
                   adjustcolor(COL_AWARE, 0.7),
                   "#E69F00"),
         pt.cex = c(0.9, 0.9, 1.5),
         pt.lwd = c(1, 1, 1.8),
         bty = "n", cex = 0.78, y.intersp = 1.1)
}

pdf(file.path(outDir, "fig3-mds-rep02.pdf"), width = 6, height = 6,
    useDingbats = FALSE)
fig3_core()
dev.off()

png(file.path(outDir, "fig3-mds-rep02.png"),
    width = 6, height = 6, units = "in", res = 300)
fig3_core()
dev.off()

cat("Fig 3 written.\n")
cat("\nAll figures written to:", outDir, "\n")
