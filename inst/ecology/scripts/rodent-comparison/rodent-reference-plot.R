# rodent-reference-plot.R ----------------------------------------------------
# Three-panel comparison of a reference rodent phylogeny against the AWARE
# and BLIND MkNT MCMC consensus trees.
#
# Data: MkNT v1 multirun chains (rodent-MkNT-v1, 2026-05-21):
#   - BLIND: 4 runs x 100k iter, treeThin=100, 4 chains PT, validated MkNT
#     likelihood (k = kObs per char; F81 for neomorphic). No Mk' (kPrime_*)
#     columns in the log. Per-run tree files rodent-MkNT-blind_trees_{1..4}.nwk
#     concatenated into a single multiPhylo.
#   - AWARE: same config + ecology-aware layer. Per-run tree files
#     rodent-MkNT-aware_trees_{1..4}.nwk.
#
# If MkNT aware tree files are not yet present (job still running) the script
# falls back to the older Mk'-based rodent-aware-v2_trees.nwk and stamps a
# prominent warning into the output panel title. Re-run after MkNT aware
# completes to refresh.
#
# Reference topology: Fabre et al. 2012 (rodent + lagomorph families,
# transcribed from the en.wikipedia Rodentia/Glires clade markup) with an
# outgroup backbone derived from `../neotrans/inst/wct/wellCorroboratedTrees.nwk`
# (asher + hallidayAll).
#
# Tip colours: TreeTools::PaintTree (default palette) applied to the painted
# reference; the resulting per-tip colour is reused as the tip colour on the
# AWARE and BLIND consensus panels so clade fidelity can be eyeballed.

suppressPackageStartupMessages({
  library(ape)
  library(TreeTools)
})

outDir   <- "inst/ecology/scripts/rodent-comparison"
mapFile  <- file.path(outDir, "taxon-mapping.tsv")
refFile  <- file.path(outDir, "reference-tree.nwk")

# Concatenate per-run NWK files into one multiPhylo. NWK files store one
# tree per line, so simple concat is equivalent.
read_multirun_trees <- function(pattern) {
  files <- Sys.glob(file.path(outDir, pattern))
  files <- files[file.info(files)$size > 0]
  if (length(files) == 0L) return(NULL)
  trees <- do.call(c, lapply(files, function(f) {
    tr <- ape::read.tree(f)
    if (inherits(tr, "phylo")) list(tr) else as.list(tr)
  }))
  class(trees) <- "multiPhylo"
  attr(trees, "source_files") <- files
  trees
}

# Blind: MkNT (current).
blind_all  <- read_multirun_trees("rodent-MkNT-blind_trees_*.nwk")
blind_tag  <- "MkNT v1 (4 x 100k iter, PT)"
if (is.null(blind_all)) {
  message("MkNT blind trees not found; falling back to Mk' rodent-blind-v2.")
  blind_all <- ape::read.tree(file.path(outDir, "rodent-blind-v2-full_trees.nwk"))
  blind_tag <- "Mk' v2 fallback"
}

# Aware: MkNT preferred; fall back to Mk' v2 if not yet streamed.
aware_all <- read_multirun_trees("rodent-MkNT-aware_trees_*.nwk")
aware_tag <- "MkNT v1 (4 x 100k iter, PT)"
if (is.null(aware_all)) {
  warning("MkNT aware trees not yet available; using Mk' v2 placeholder. ",
          "Re-run once /nobackup/pjjg18/mkp-rodent-MkNT-v1/aware/ ",
          "rodent-MkNT-aware_trees_*.nwk are populated.")
  aware_all <- ape::read.tree(file.path(outDir, "rodent-aware-v2_trees.nwk"))
  aware_tag <- "Mk' v2 PLACEHOLDER (MkNT aware run not yet complete)"
}

# --- 1. Mapping --------------------------------------------------------------
map <- read.table(mapFile, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "",
                  comment.char = "", fill = TRUE)
stopifnot(all(c("tip", "family") %in% names(map)))

ambig <- map$family == "??" | is.na(map$family) | map$family == ""
if (any(ambig)) {
  stop("Unresolved family assignments in ", mapFile, " for: ",
       paste(map$tip[ambig], collapse = ", "),
       ". Edit the TSV and re-run.")
}

# --- 2. Reference ------------------------------------------------------------
ref <- ape::read.tree(refFile)
ref_tips <- ref$tip.label

# Drop reference tips with no matrix taxa.
used_fams <- unique(map$family)
unused <- setdiff(ref_tips, used_fams)
if (length(unused)) {
  cat("Dropping reference tips with no matrix taxa:",
      paste(unused, collapse = ", "), "\n")
  ref <- ape::drop.tip(ref, unused)
  ref_tips <- ref$tip.label
}

missing_ref <- setdiff(used_fams, ref_tips)
if (length(missing_ref)) {
  stop("Mapping references families not present in reference tree: ",
       paste(missing_ref, collapse = ", "))
}

# --- 3. Expand family tips into polytomies of matrix genera ------------------
# For each reference tip, replace it with a soft polytomy of its assigned
# matrix taxa.  If a family has a single matrix taxon, simply rename the tip.
ref_expanded <- ref
for (fam in ref_tips) {
  genera <- map$tip[map$family == fam]
  if (length(genera) == 1L) {
    ref_expanded$tip.label[ref_expanded$tip.label == fam] <- genera
  } else {
    # Build a star polytomy of `genera` and graft at the focal tip.
    star <- ape::read.tree(text = paste0("(", paste(genera, collapse = ","), ");"))
    tip_idx <- which(ref_expanded$tip.label == fam)
    ref_expanded <- ape::bind.tree(ref_expanded, star, where = tip_idx)
    # bind.tree leaves the original tip in place; drop it.
    ref_expanded <- ape::drop.tip(ref_expanded, fam)
  }
}

# --- 4. Validate against matrix tip set --------------------------------------
if (inherits(aware_all, "phylo")) aware_all <- list(aware_all)
if (inherits(blind_all, "phylo")) blind_all <- list(blind_all)
class(aware_all) <- "multiPhylo"
class(blind_all) <- "multiPhylo"
cat(sprintf("Loaded: %d aware trees (%s), %d blind trees (%s)\n",
            length(aware_all), aware_tag, length(blind_all), blind_tag))

matrix_tips <- sort(aware_all[[1]]$tip.label)
ref_tips_x  <- sort(ref_expanded$tip.label)

if (!identical(matrix_tips, ref_tips_x)) {
  only_ref <- setdiff(ref_tips_x, matrix_tips)
  only_mat <- setdiff(matrix_tips, ref_tips_x)
  stop("Tip sets do not match.\n  Only in reference: ",
       paste(only_ref, collapse = ", "),
       "\n  Only in matrix: ",
       paste(only_mat, collapse = ", "))
}
cat(sprintf("Reference tip set matches matrix tip set (%d tips).\n",
            length(matrix_tips)))

# --- 5. Consensus trees (25% burnin, MR rule, unrooted) ----------------------
burnin_frac <- 0.25
post_burnin <- function(trees) {
  n <- length(trees)
  trees[seq.int(ceiling(n * burnin_frac) + 1L, n)]
}
aware_post <- post_burnin(aware_all)
blind_post <- post_burnin(blind_all)
cons_aware <- ape::consensus(aware_post, p = 0.5, rooted = FALSE)
cons_blind <- ape::consensus(blind_post, p = 0.5, rooted = FALSE)

# --- 6. Paint reference & build per-tip colour vector ------------------------
painted <- TreeTools::PaintTree(ref_expanded)
tip_col_by_label <- setNames(painted$tipCol, ref_expanded$tip.label)

col_for <- function(tree) tip_col_by_label[tree$tip.label]

# Edge colours for the reference (already from PaintTree).
ref_edge_col <- painted$edgeCol

# For aware/blind, derive a per-edge colour from descendant tip colours: if
# all descendants share a colour use it, otherwise grey (clade mismatch).
edge_cols_for <- function(tree) {
  n_tip <- length(tree$tip.label)
  edges <- tree$edge
  children_of <- split(edges[, 2], edges[, 1])
  desc_cache <- new.env(parent = emptyenv())
  descend <- function(node) {
    key <- as.character(node)
    if (!is.null(desc_cache[[key]])) return(desc_cache[[key]])
    out <- if (node <= n_tip) {
      tree$tip.label[node]
    } else {
      unlist(lapply(children_of[[key]], descend), use.names = FALSE)
    }
    desc_cache[[key]] <- out
    out
  }
  vapply(edges[, 2], function(node) {
    cols <- unique(tip_col_by_label[descend(node)])
    if (length(cols) == 1L) cols else "grey70"
  }, character(1))
}

# --- 7. Plot 3-panel ---------------------------------------------------------
pdfFile <- file.path(outDir, "rodent-3panel.pdf")
pngFile <- file.path(outDir, "rodent-3panel.png")

draw_panel <- function(tree, main, edge_col) {
  plot(tree, type = "phylogram", cex = 0.8,
       tip.color = col_for(tree),
       edge.color = edge_col,
       edge.width = 1.5,
       main = main,
       xpd = NA,
       no.margin = TRUE)
}

draw_all <- function() {
  op <- par(mfrow = c(1, 3), mar = c(0, 0, 0, 0), oma = c(0, 0, 1, 0))
  on.exit(par(op))
  draw_panel(ref_expanded, "Reference (Fabre 2012 + WCT outgroups)",
             ref_edge_col)
  draw_panel(cons_aware,
             sprintf("AWARE MR consensus (N=%d)\n%s", length(aware_post), aware_tag),
             edge_cols_for(cons_aware))
  draw_panel(cons_blind,
             sprintf("BLIND MR consensus (N=%d)\n%s", length(blind_post), blind_tag),
             edge_cols_for(cons_blind))
}

pdf(pdfFile, width = 8, height = 4)
draw_all()
dev.off()
cat("Saved:", pdfFile, "\n")

png(pngFile, width = 18, height = 9, units = "in", res = 200)
draw_all()
dev.off()
cat("Saved:", pngFile, "\n")
