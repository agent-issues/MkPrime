# aware-mode-plots.R --------------------------------------------------------
# Two diagnostic plots for the v4 aware-vs-blind situation:
#   1. 4-panel PaintTree: reference + aware run1 cons + aware run2 cons + blind cons
#   2. CID-MDS with aware run1, aware run2, blind clouds + reference point
# Uses the same PaintTree / consensus / CID infrastructure already in this dir.

suppressPackageStartupMessages({
  library(ape)
  library(TreeTools)
  library(TreeDist)
})

outDir <- "inst/ecology/scripts/rodent-comparison"
mapFile <- file.path(outDir, "taxon-mapping.tsv")
refFile <- file.path(outDir, "reference-tree.nwk")

# --- Load ------------------------------------------------------------------
ref <- ape::read.tree(refFile)
map <- read.table(mapFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                  quote = "", comment.char = "", fill = TRUE)

# Drop ref tips with no matrix taxa
used_fams <- unique(map$family)
unused <- setdiff(ref$tip.label, used_fams)
if (length(unused)) ref <- ape::drop.tip(ref, unused)

# Expand reference families to genus polytomies (mirroring rodent-reference-plot.R)
ref_expanded <- ref
for (fam in ref$tip.label) {
  genera <- map$tip[map$family == fam]
  if (length(genera) == 1L) {
    ref_expanded$tip.label[ref_expanded$tip.label == fam] <- genera
  } else {
    star <- ape::read.tree(text = paste0("(", paste(genera, collapse = ","), ");"))
    tip_idx <- which(ref_expanded$tip.label == fam)
    ref_expanded <- ape::bind.tree(ref_expanded, star, where = tip_idx)
    ref_expanded <- ape::drop.tip(ref_expanded, fam)
  }
}

aware1 <- ape::read.tree(file.path(outDir, "rodent-MkNT-aware_trees_1.nwk"))
aware2 <- ape::read.tree(file.path(outDir, "rodent-MkNT-aware_trees_2.nwk"))
blind  <- do.call(c, lapply(Sys.glob(file.path(outDir, "rodent-MkNT-blind_trees_*.nwk")),
                            ape::read.tree))
class(aware1) <- "multiPhylo"
class(aware2) <- "multiPhylo"
class(blind)  <- "multiPhylo"

discard <- function(L) L[seq.int(ceiling(length(L)/4) + 1L, length(L))]
a1_post <- discard(aware1)
a2_post <- discard(aware2)
b_post  <- discard(blind)

cat("Sample counts: aware1=", length(a1_post),
    " aware2=", length(a2_post),
    " blind=", length(b_post), "\n", sep = "")

cons_a1 <- ape::consensus(a1_post, p = 0.5, rooted = FALSE)
cons_a2 <- ape::consensus(a2_post, p = 0.5, rooted = FALSE)
cons_b  <- ape::consensus(b_post,  p = 0.5, rooted = FALSE)

# --- PaintTree colours from reference --------------------------------------
painted <- TreeTools::PaintTree(ref_expanded)
tip_col_by_label <- setNames(painted$tipCol, ref_expanded$tip.label)
col_for <- function(tree) tip_col_by_label[tree$tip.label]

# Per-edge colour: if all descendants share a tip colour use it, else grey
edge_cols_for <- function(tree) {
  n_tip <- length(tree$tip.label)
  edges <- tree$edge
  children_of <- split(edges[, 2], edges[, 1])
  cache <- new.env(parent = emptyenv())
  descend <- function(node) {
    key <- as.character(node)
    if (!is.null(cache[[key]])) return(cache[[key]])
    out <- if (node <= n_tip) tree$tip.label[node]
           else unlist(lapply(children_of[[key]], descend), use.names = FALSE)
    cache[[key]] <- out
    out
  }
  vapply(edges[, 2], function(node) {
    cols <- unique(tip_col_by_label[descend(node)])
    if (length(cols) == 1L) cols else "grey70"
  }, character(1))
}

# --- Plot 1: 4-panel tree --------------------------------------------------
pdf4 <- file.path(outDir, "aware-modes-4panel.pdf")
png4 <- file.path(outDir, "aware-modes-4panel.png")

draw_panel <- function(tree, main, edge_col) {
  plot(tree, type = "phylogram", cex = 0.55,
       tip.color = col_for(tree),
       edge.color = edge_col, edge.width = 1.5,
       main = main, xpd = NA, no.margin = TRUE)
}

draw4 <- function() {
  op <- par(mfrow = c(1, 4), mar = c(0, 0, 1.5, 0), oma = c(0, 0, 1, 0))
  on.exit(par(op))
  draw_panel(ref_expanded, "Reference\n(Fabre 2012 + WCT)",  painted$edgeCol)
  draw_panel(cons_a1, sprintf("AWARE run 1 (N=%d)\nMkNT v4 PT", length(a1_post)),
             edge_cols_for(cons_a1))
  draw_panel(cons_a2, sprintf("AWARE run 2 (N=%d)\nMkNT v4 PT", length(a2_post)),
             edge_cols_for(cons_a2))
  draw_panel(cons_b,  sprintf("BLIND (N=%d)\nMkNT v1 4x100k", length(b_post)),
             edge_cols_for(cons_b))
}
pdf(pdf4, width = 18, height = 10); draw4(); dev.off()
png(png4, width = 18, height = 10, units = "in", res = 200); draw4(); dev.off()
cat("Saved:", pdf4, "\n")

# --- Plot 2: CID-MDS with 4 clouds + reference -----------------------------
# Subsample for tractable pairwise CID
set.seed(20260521)
target_n <- 100L
take <- function(trees, n) trees[sort(sample(length(trees), min(length(trees), n)))]
sub_a1 <- take(a1_post, target_n)
sub_a2 <- take(a2_post, target_n)
sub_b  <- take(b_post, target_n)

# Add the reference (expanded) as a single point — make a "trivial multiPhylo"
ref_for_dist <- list(ref_expanded)
class(ref_for_dist) <- "multiPhylo"

combined <- c(lapply(sub_a1, ape::unroot),
              lapply(sub_a2, ape::unroot),
              lapply(sub_b,  ape::unroot),
              lapply(ref_for_dist, ape::unroot))
class(combined) <- "multiPhylo"
N_total <- length(combined)
n_a1 <- length(sub_a1); n_a2 <- length(sub_a2); n_b <- length(sub_b)
cat("CID input size:", N_total, "\n")

dist_mat <- TreeDist::ClusteringInfoDistance(combined)
mds <- cmdscale(dist_mat, k = 2, eig = TRUE)
pts <- mds$points
var_pct <- 100 * abs(mds$eig[1:2]) / sum(abs(mds$eig))

idx_a1 <- seq_len(n_a1)
idx_a2 <- (n_a1 + 1):(n_a1 + n_a2)
idx_b  <- (n_a1 + n_a2 + 1):(n_a1 + n_a2 + n_b)
idx_ref <- N_total

col_a1 <- "#D55E00"  # vermillion
col_a2 <- "#E69F00"  # orange
col_b  <- "#56B4E9"  # sky blue
col_r  <- "#000000"

pdf_mds <- file.path(outDir, "aware-modes-cid-mds.pdf")
png_mds <- file.path(outDir, "aware-modes-cid-mds.png")

do_mds_plot <- function() {
  xlim <- range(pts[, 1]) * 1.1
  ylim <- range(pts[, 2]) * 1.1
  plot(pts[idx_a1, 1], pts[idx_a1, 2],
       col = scales::alpha(col_a1, 0.7), pch = 19, cex = 0.8,
       xlim = xlim, ylim = ylim,
       xlab = sprintf("MDS Axis 1 (%.1f%%)", var_pct[1]),
       ylab = sprintf("MDS Axis 2 (%.1f%%)", var_pct[2]),
       main = "Posterior tree space: AWARE r1 / r2 vs BLIND vs REFERENCE",
       las = 1, cex.main = 1.0)
  points(pts[idx_a2, 1], pts[idx_a2, 2],
         col = scales::alpha(col_a2, 0.7), pch = 19, cex = 0.8)
  points(pts[idx_b, 1], pts[idx_b, 2],
         col = scales::alpha(col_b, 0.7), pch = 19, cex = 0.8)

  # Hulls
  for (i_set in list(list(idx = idx_a1, col = col_a1),
                     list(idx = idx_a2, col = col_a2),
                     list(idx = idx_b,  col = col_b))) {
    h <- chull(pts[i_set$idx, ])
    polygon(pts[i_set$idx[h], 1], pts[i_set$idx[h], 2],
            border = i_set$col, lwd = 1.5, lty = 2, col = NA)
  }

  # Reference as a star
  points(pts[idx_ref, 1], pts[idx_ref, 2], col = col_r, pch = 8, cex = 2.5, lwd = 2)

  # Centroids
  for (i_set in list(list(idx = idx_a1, col = col_a1),
                     list(idx = idx_a2, col = col_a2),
                     list(idx = idx_b,  col = col_b))) {
    cent <- colMeans(pts[i_set$idx, ])
    points(cent[1], cent[2], col = i_set$col, pch = 8, cex = 1.5, lwd = 2)
  }

  legend("topright",
         legend = c(sprintf("Aware run 1 (N=%d)", n_a1),
                    sprintf("Aware run 2 (N=%d)", n_a2),
                    sprintf("Blind (N=%d)", n_b),
                    "Reference (Fabre)",
                    "Cloud centroid"),
         col    = c(col_a1, col_a2, col_b, col_r, "grey40"),
         pch    = c(19, 19, 19, 8, 8),
         pt.cex = c(0.9, 0.9, 0.9, 1.5, 1.0),
         pt.lwd = c(1, 1, 1, 2, 2),
         bty = "n", cex = 0.85)
}
pdf(pdf_mds, width = 7, height = 7); do_mds_plot(); dev.off()
png(png_mds, width = 7, height = 7, units = "in", res = 300); do_mds_plot(); dev.off()
cat("Saved:", pdf_mds, "\n")

# Summary of inter-cloud distances
cat("\n=== Inter-cloud mean CID ===\n")
ld <- as.matrix(dist_mat)
mean_d <- function(i, j) mean(ld[i, j])
cat(sprintf("aware1<->aware2 : %.3f\n", mean_d(idx_a1, idx_a2)))
cat(sprintf("aware1<->blind  : %.3f\n", mean_d(idx_a1, idx_b)))
cat(sprintf("aware2<->blind  : %.3f\n", mean_d(idx_a2, idx_b)))
cat(sprintf("aware1<->ref    : %.3f\n", mean(ld[idx_a1, idx_ref])))
cat(sprintf("aware2<->ref    : %.3f\n", mean(ld[idx_a2, idx_ref])))
cat(sprintf("blind <->ref    : %.3f\n", mean(ld[idx_b,  idx_ref])))
cat(sprintf("within aware1   : %.3f\n", mean_d(idx_a1, idx_a1)))
cat(sprintf("within aware2   : %.3f\n", mean_d(idx_a2, idx_a2)))
cat(sprintf("within blind    : %.3f\n", mean_d(idx_b,  idx_b)))
