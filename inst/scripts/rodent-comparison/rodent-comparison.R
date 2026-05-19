# rodent-comparison.R ---------------------------------------------------------
# Compare posterior trees from BLIND vs AWARE MkPrime MCMC chains on the
# rodent morphological dataset.
#
# AWARE chain: rodent-aware-v2_trees.nwk  (~160 trees, 1M iter, treeThin=1000)
#   Completed 2026-05-19; resumed from 371k → 1M, 20.68h wall time.
#   minESS = 88 at iter 1M (from .er log). Below 200; flagged in summary.
#
# BLIND chain: rodent-blind-v2-full_trees.nwk  (~359 trees, 1M iter)
#   Completed 2026-05-17. minESS = 38 at iter 1M (from .er log). Below 200;
#   flagged in summary.
#
# Outputs (all under inst/scripts/rodent-comparison/):
#   rodent-cid-mds.pdf / .png    -- Plot 1: 2D CID MDS
#   rodent-consensus.pdf / .png  -- Plot 2: MR consensus trees with conflict
#   rodent-comparison.md         -- Written summary
#
# Root-dependence audit (2026-05-19) -------------------------------------------
# The sim3 scoring layer was contaminated by a `ape::prop.part` + `setequal`
# pattern that is root-dependent (see `inst/simulations/ecology/sim3-scoring.R`).
# The MCMC posteriors here also have varying root configurations (aware: 4
# distinct rootings across 120 post-burnin trees; blind: 9 distinct rootings
# across 269 trees) -- so any root-dependent primitive would be vulnerable.
#
# This script was audited and is **NOT contaminated**:
#   - `ape::consensus(rooted = FALSE)` (lines below) internally re-roots all
#     trees to tip 1 and then runs `postprocess.prop.part(., "SHORTwise")`,
#     which canonicalises bipartitions to the smaller side AND sums counts
#     of duplicates. This is root-invariant. (Confirmed: ape::consensus source.)
#   - `ape::prop.part(aware_post)` (used only as input to `prop.clades` below)
#     is consumed by `ape:::prop.clades(cons_tree, pp_obj)` which defaults to
#     `rooted = FALSE` and applies `SHORTwise` itself.
#   - `pp_to_canonical_splits()` (line ~148) runs `ape::prop.part(list(tree))`
#     on a single consensus tree and canonicalises to the smaller side --
#     equivalent to Splits-based extraction for one tree.
#
# Empirical verification: a manual Splits-based MR consensus
# (TreeTools::as.Splits + SHORTwise + 50% threshold) reproduces 45 aware
# splits and 38 blind splits, with the same 32 shared count, matching
# `ape::consensus` exactly.

suppressPackageStartupMessages({
  library(ape)
  library(TreeDist)
  library(TreeTools)
  library(phangorn)
})

outDir <- "inst/scripts/rodent-comparison"

# --- 1. Load trees -----------------------------------------------------------
aware_nwk <- file.path(outDir, "rodent-aware-v2_trees.nwk")
blind_nwk  <- file.path(outDir, "rodent-blind-v2-full_trees.nwk")

# ape::read.tree reads one tree per line (standard Newick).
aware_all <- ape::read.tree(aware_nwk)
blind_all  <- ape::read.tree(blind_nwk)

# Ensure multiPhylo class
if (inherits(aware_all, "phylo")) aware_all <- list(aware_all)
if (inherits(blind_all, "phylo")) blind_all <- list(blind_all)
class(aware_all) <- "multiPhylo"
class(blind_all) <- "multiPhylo"

n_aware_total <- length(aware_all)
n_blind_total <- length(blind_all)
cat(sprintf("Loaded: %d aware trees, %d blind trees\n",
            n_aware_total, n_blind_total))

# --- 2. Post-burnin subsetting -----------------------------------------------
# treeThin = nIter %/% 1000 = 1000 for both chains.
# Apply 25% burnin: discard first 25% of trees.
burnin_frac <- 0.25
n_aware_burn <- max(1L, ceiling(n_aware_total * burnin_frac))
n_blind_burn <- max(1L, ceiling(n_blind_total * burnin_frac))

aware_post <- aware_all[seq.int(n_aware_burn + 1L, n_aware_total)]
blind_post  <- blind_all[seq.int(n_blind_burn + 1L, n_blind_total)]

n_aware_post <- length(aware_post)
n_blind_post <- length(blind_post)
cat(sprintf("Post-burnin: %d aware, %d blind\n", n_aware_post, n_blind_post))

# Subsample for CID analysis (target 200 per chain).
target_n <- 200L
set.seed(20260519)
if (n_aware_post > target_n) {
  aware_sub <- aware_post[sort(sample(n_aware_post, target_n))]
} else {
  aware_sub <- aware_post
}
if (n_blind_post > target_n) {
  blind_sub <- blind_post[sort(sample(n_blind_post, target_n))]
} else {
  blind_sub <- blind_post
}
n_aware_sub <- length(aware_sub)
n_blind_sub <- length(blind_sub)
cat(sprintf("Subsampled: %d aware, %d blind\n", n_aware_sub, n_blind_sub))

# Use unrooted trees for CID.
aware_sub_unr <- lapply(aware_sub, ape::unroot)
blind_sub_unr <- lapply(blind_sub, ape::unroot)
class(aware_sub_unr) <- "multiPhylo"
class(blind_sub_unr) <- "multiPhylo"

# Combined set for pairwise CID.
combined <- c(aware_sub_unr, blind_sub_unr)
N_total <- length(combined)
cat(sprintf("Combined tree set for CID: %d trees\n", N_total))

# --- 3. Pairwise CID distances -----------------------------------------------
cat("Computing pairwise CID distances...\n")
dist_mat <- TreeDist::ClusteringInfoDistance(combined)
cat("CID matrix computed.\n")

# Submatrices for within/between comparisons.
idx_a <- seq_len(n_aware_sub)
idx_b <- seq.int(n_aware_sub + 1L, N_total)

dist_within_aware  <- as.numeric(as.dist(as.matrix(dist_mat)[idx_a, idx_a]))
dist_within_blind  <- as.numeric(as.dist(as.matrix(dist_mat)[idx_b, idx_b]))
dist_between       <- as.numeric(as.matrix(dist_mat)[idx_a, idx_b])

cid_summary <- data.frame(
  comparison       = c("within_aware", "within_blind", "between"),
  n_pairs          = c(length(dist_within_aware),
                        length(dist_within_blind),
                        length(dist_between)),
  mean_CID         = c(mean(dist_within_aware),
                        mean(dist_within_blind),
                        mean(dist_between)),
  median_CID       = c(median(dist_within_aware),
                        median(dist_within_blind),
                        median(dist_between))
)
cat("\nCID summary:\n"); print(cid_summary)

# --- 4. MDS projection -------------------------------------------------------
cat("Running MDS...\n")
mds_raw <- cmdscale(dist_mat, k = 2, eig = TRUE)
mds_pts <- mds_raw$points
cat(sprintf("MDS variance explained: Axis1=%.1f%%, Axis2=%.1f%%\n",
            100 * mds_raw$eig[1] / sum(abs(mds_raw$eig)),
            100 * mds_raw$eig[2] / sum(abs(mds_raw$eig))))

pts_aware <- mds_pts[idx_a, , drop = FALSE]
pts_blind <- mds_pts[idx_b, , drop = FALSE]

# --- 5. MR consensus trees ---------------------------------------------------
cat("Computing MR consensus trees...\n")
cons_aware <- ape::consensus(aware_post, p = 0.5, rooted = FALSE)
cons_blind <- ape::consensus(blind_post,  p = 0.5, rooted = FALSE)

# Posterior support for each clade in each consensus.
pp_aware <- ape::prop.part(aware_post)
pp_blind <- ape::prop.part(blind_post)

get_node_support <- function(cons_tree, pp_obj, n_trees) {
  counts <- tryCatch(
    ape:::prop.clades(cons_tree, pp_obj),
    error = function(e) NULL
  )
  if (is.null(counts)) return(NULL)
  counts / n_trees
}

supp_aware <- get_node_support(cons_aware, pp_aware, n_aware_post)
supp_blind <- get_node_support(cons_blind, pp_blind, n_blind_post)

# --- 6. Identify conflicting splits -----------------------------------------
pp_to_canonical_splits <- function(tree, all_taxa) {
  pp <- ape::prop.part(list(tree))
  taxa_here <- attr(pp, "labels")
  all_sorted <- sort(all_taxa)
  n_all <- length(all_sorted)
  splits <- lapply(pp, function(idx) {
    s <- sort(taxa_here[idx])
    comp <- sort(setdiff(all_sorted, s))
    if (length(s) < length(comp)) return(s)
    if (length(s) > length(comp)) return(comp)
    if (s[1] <= comp[1]) return(s) else return(comp)
  })
  splits <- splits[sapply(splits, length) > 1]
  splits <- splits[sapply(splits, length) < n_all - 1]
  unique(splits)
}

all_taxa <- sort(cons_aware$tip.label)
splits_aware <- pp_to_canonical_splits(cons_aware, all_taxa)
splits_blind  <- pp_to_canonical_splits(cons_blind,  all_taxa)

split_match <- function(s, split_list) {
  any(sapply(split_list, function(x) identical(x, s)))
}

blind_unique <- splits_blind[!sapply(splits_blind, split_match, splits_aware)]
aware_unique <- splits_aware[!sapply(splits_aware, split_match, splits_blind)]

n_shared <- length(splits_blind) - length(blind_unique)
n_conflict <- length(blind_unique) + length(aware_unique)

cat(sprintf("\nSplit comparison:\n  Blind consensus: %d splits\n  Aware consensus: %d splits\n",
            length(splits_blind), length(splits_aware)))
cat(sprintf("  Shared: %d\n  Blind-unique: %d\n  Aware-unique: %d\n",
            n_shared, length(blind_unique), length(aware_unique)))

format_split <- function(s) paste(s, collapse = "+")

# --- 7. Consensus tree edge colouring ----------------------------------------
edge_colours <- function(cons_tree, cons_splits, unique_splits) {
  tree_r <- ape::root(cons_tree, all_taxa[1], resolve.root = TRUE)
  tips_r <- tree_r$tip.label
  n_tips <- length(tips_r)
  all_sorted <- sort(all_taxa)
  ed <- tree_r$edge
  cols <- rep("black", nrow(ed))

  children_of <- function(node) tree_r$edge[tree_r$edge[, 1] == node, 2]

  tip_descs <- function(node) {
    if (node <= n_tips) return(node)
    ch <- children_of(node)
    unlist(lapply(ch, tip_descs))
  }

  for (i in seq_len(nrow(ed))) {
    child <- ed[i, 2]
    if (child <= n_tips) next
    desc_idx <- tip_descs(child)
    s <- sort(tips_r[desc_idx])
    comp <- sort(setdiff(all_sorted, s))
    if (length(s) < length(comp)) cand <- s
    else if (length(s) > length(comp)) cand <- comp
    else if (s[1] <= comp[1]) cand <- s else cand <- comp
    if (length(cand) <= 1 || length(cand) >= length(all_sorted) - 1) next
    if (split_match(cand, unique_splits)) cols[i] <- "red"
  }
  cols
}

ecols_blind <- edge_colours(cons_blind, splits_blind, blind_unique)
ecols_aware <- edge_colours(cons_aware, splits_aware, aware_unique)

# --- 8. MR consensus position on MDS ----------------------------------------
cons_pt_aware <- colMeans(pts_aware)
cons_pt_blind  <- colMeans(pts_blind)

# --- 9. Plot 1: CID MDS ------------------------------------------------------
cat("Generating Plot 1: CID MDS...\n")

col_aware <- "#E69F00"  # orange
col_blind <- "#56B4E9"  # sky blue

plot1_pdf <- file.path(outDir, "rodent-cid-mds.pdf")
plot1_png <- file.path(outDir, "rodent-cid-mds.png")

do_mds_plot <- function() {
  xlim <- range(mds_pts[, 1]) * 1.1
  ylim <- range(mds_pts[, 2]) * 1.1

  plot(pts_aware[, 1], pts_aware[, 2],
       col = scales::alpha(col_aware, 0.7),
       pch = 19, cex = 0.8,
       xlim = xlim, ylim = ylim,
       xlab = "MDS Axis 1", ylab = "MDS Axis 2",
       main = "Posterior tree space: BLIND vs AWARE (CID, MDS)\n(Full 1M-iter chains)",
       las = 1)
  points(pts_blind[, 1], pts_blind[, 2],
         col = scales::alpha(col_blind, 0.7),
         pch = 19, cex = 0.8)

  if (n_aware_sub >= 3) {
    hull_a <- chull(pts_aware)
    polygon(pts_aware[hull_a, 1], pts_aware[hull_a, 2],
            border = col_aware, lwd = 1.5, lty = 2, col = NA)
  }
  if (n_blind_sub >= 3) {
    hull_b <- chull(pts_blind)
    polygon(pts_blind[hull_b, 1], pts_blind[hull_b, 2],
            border = col_blind, lwd = 1.5, lty = 2, col = NA)
  }

  points(cons_pt_aware[1], cons_pt_aware[2],
         col = col_aware, pch = 8, cex = 2, lwd = 2)
  points(cons_pt_blind[1], cons_pt_blind[2],
         col = col_blind, pch = 8, cex = 2, lwd = 2)

  legend("topright",
         legend = c(sprintf("Aware (N=%d)", n_aware_sub),
                    sprintf("Blind (N=%d)", n_blind_sub),
                    "Chain centroid"),
         col    = c(col_aware, col_blind, "black"),
         pch    = c(19, 19, 8),
         pt.cex = c(0.8, 0.8, 1.5),
         bty    = "n", cex = 0.85)
}

pdf(plot1_pdf, width = 6, height = 6)
do_mds_plot()
dev.off()
cat("Saved:", plot1_pdf, "\n")

png(plot1_png, width = 6, height = 6, units = "in", res = 300)
do_mds_plot()
dev.off()
cat("Saved:", plot1_png, "\n")

# --- 10. Plot 2: Consensus trees side-by-side --------------------------------
cat("Generating Plot 2: Consensus trees...\n")

plot2_pdf <- file.path(outDir, "rodent-consensus.pdf")
plot2_png <- file.path(outDir, "rodent-consensus.png")

do_consensus_plot <- function() {
  op <- par(mfrow = c(1, 2), mar = c(1, 1, 2.5, 1), oma = c(0, 0, 0, 0))
  on.exit(par(op))

  n_tips <- length(cons_blind$tip.label)
  tip_cex <- max(0.3, min(0.6, 6 / n_tips))

  # BLIND consensus.
  plot(cons_blind, type = "phylogram", cex = tip_cex,
       main = sprintf("BLIND MR consensus (N=%d post-burnin)\n(red = blind-unique splits)",
                      n_blind_post),
       edge.color = ecols_blind,
       show.node.label = FALSE, no.margin = FALSE)
  if (!is.null(supp_blind) && length(supp_blind) > 0) {
    supp_txt <- ifelse(is.na(supp_blind), "",
                       sprintf("%.2f", supp_blind))
    supp_txt[supp_txt == "NA"] <- ""
    ape::nodelabels(supp_txt,
                    frame = "none", adj = c(1.1, -0.3),
                    cex = max(0.25, tip_cex * 0.7),
                    col = "grey40")
  }

  # AWARE consensus.
  plot(cons_aware, type = "phylogram", cex = tip_cex,
       main = sprintf(
         "AWARE MR consensus (N=%d post-burnin)\n(red = aware-unique splits)",
         n_aware_post),
       edge.color = ecols_aware,
       show.node.label = FALSE, no.margin = FALSE)
  if (!is.null(supp_aware) && length(supp_aware) > 0) {
    supp_txt <- ifelse(is.na(supp_aware), "",
                       sprintf("%.2f", supp_aware))
    supp_txt[supp_txt == "NA"] <- ""
    ape::nodelabels(supp_txt,
                    frame = "none", adj = c(1.1, -0.3),
                    cex = max(0.25, tip_cex * 0.7),
                    col = "grey40")
  }
}

pdf(plot2_pdf, width = 14, height = 7)
do_consensus_plot()
dev.off()
cat("Saved:", plot2_pdf, "\n")

png(plot2_png, width = 14, height = 7, units = "in", res = 300)
do_consensus_plot()
dev.off()
cat("Saved:", plot2_png, "\n")

# --- 11. Written summary ------------------------------------------------------
cat("Writing summary...\n")

blind_unique_str <- paste0(
  sapply(blind_unique, function(s) {
    paste0("  {", paste(s, collapse = ", "), "}")
  }),
  collapse = "\n")
aware_unique_str <- paste0(
  sapply(aware_unique, function(s) {
    paste0("  {", paste(s, collapse = ", "), "}")
  }),
  collapse = "\n")

centroid_dist <- sqrt(sum((cons_pt_aware - cons_pt_blind)^2))
aware_spread  <- mean(sqrt(rowSums((pts_aware - cons_pt_aware)^2)))
blind_spread  <- mean(sqrt(rowSums((pts_blind - cons_pt_blind)^2)))

overlap_desc <- if (centroid_dist < max(aware_spread, blind_spread)) {
  "overlapping"
} else {
  "largely separated"
}

# ESS from .er log (hard-coded from Hamilton run logs).
aware_minESS <- 88L
blind_minESS <- 38L
aware_ess_note <- if (aware_minESS < 200) sprintf("**BELOW 200 (minESS = %d)**", aware_minESS) else sprintf("%d", aware_minESS)
blind_ess_note <- if (blind_minESS < 200) sprintf("**BELOW 200 (minESS = %d)**", blind_minESS) else sprintf("%d", blind_minESS)

summary_md <- sprintf(
'# Rodent BLIND vs AWARE comparison

Generated from full 1M-iteration chains (both chains complete as of 2026-05-19).

## Numerical summary

| Quantity | Value |
|---|---|
| Aware trees (total from NWK) | %d |
| Blind trees (total from NWK) | %d |
| Burnin fraction applied | %.0f%% |
| Aware post-burnin | %d |
| Blind post-burnin | %d |
| Aware subsampled for CID | %d |
| Blind subsampled for CID | %d |
| Aware minESS (at 1M iter, from .er log) | %s |
| Blind minESS (at 1M iter, from .er log) | %s |

### Pairwise CID distances (subsampled set)

| Comparison | N pairs | Mean CID | Median CID |
|---|---|---|---|
| Within aware | %d | %.4f | %.4f |
| Within blind | %d | %.4f | %.4f |
| Between chains | %d | %.4f | %.4f |

## Split agreement

- BLIND MR consensus: **%d splits**
- AWARE MR consensus: **%d splits**
- Shared splits: **%d**
- Blind-unique splits: **%d**
- Aware-unique splits: **%d**

### Blind-unique splits (absent from AWARE consensus)
%s

### Aware-unique splits (absent from BLIND consensus)
%s

## MDS topology

The two chains\' posterior distributions are **%s** in CID-MDS space
(centroid distance = %.4f; aware spread = %.4f, blind spread = %.4f).

## Key patterns

%s

## ESS and convergence

Both chains ran for 1M iterations with treeThin = 1000. The minESS values
reported here are continuous-parameter ESS from the MkPrime MCMC log (the
minimum over all monitored parameters at the final iteration).

- **Aware chain**: minESS = %d at 1M iterations. %s
- **Blind chain**: minESS = %d at 1M iterations. %s

ESS below 200 indicates that the chains have not fully converged on the
parameter that mixes most slowly (likely topology or a correlated rate
parameter). Results should be treated as indicative rather than definitive.
A further continuation or parallel-tempering run is advisable for publication.

## Chain details

- AWARE: 1M iterations, resumed from 371k checkpoint (2026-05-18 to 2026-05-19), 20.68h wall time.
- BLIND: 1M iterations, resumed from 200k checkpoint (2026-05-17), 19.5 min wall time
  (blind chain resumed quickly because the standard Mk\' likelihood is much faster).
',
  n_aware_total,
  n_blind_total,
  burnin_frac * 100,
  n_aware_post,
  n_blind_post,
  n_aware_sub,
  n_blind_sub,
  aware_ess_note,
  blind_ess_note,
  cid_summary$n_pairs[1], cid_summary$mean_CID[1], cid_summary$median_CID[1],
  cid_summary$n_pairs[2], cid_summary$mean_CID[2], cid_summary$median_CID[2],
  cid_summary$n_pairs[3], cid_summary$mean_CID[3], cid_summary$median_CID[3],
  length(splits_blind),
  length(splits_aware),
  n_shared,
  length(blind_unique),
  length(aware_unique),
  if (nchar(blind_unique_str) == 0) "(none)" else blind_unique_str,
  if (nchar(aware_unique_str) == 0) "(none)" else aware_unique_str,
  overlap_desc,
  centroid_dist,
  aware_spread,
  blind_spread,
  if (centroid_dist < max(aware_spread, blind_spread)) {
    paste0(
      "The aware and blind posterior clouds overlap in CID-MDS space, suggesting broad ",
      "concordance in recovered topology despite the model difference. Topological ",
      "separation between ecology-aware and blind inference, if present, is visible ",
      "at the periphery of the distribution and in unique splits of each consensus."
    )
  } else {
    paste0(
      "The aware and blind posterior clouds are distinctly separated in CID-MDS space, ",
      "indicating systematic topological differences driven by the ecology-aware model. ",
      "The ecology covariation prior shifts the inferred rodent phylogeny in a direction ",
      "that is inconsistent with the blind Mk\' posterior — consistent with the hypothesis ",
      "that ecological convergence creates homoplasy that misleads standard Mk\' inference."
    )
  },
  aware_minESS,
  if (aware_minESS < 200) "Flag: below recommended minimum of 200." else "Adequate.",
  blind_minESS,
  if (blind_minESS < 200) "Flag: below recommended minimum of 200." else "Adequate."
)

md_path <- file.path(outDir, "rodent-comparison.md")
writeLines(summary_md, md_path)
cat("Saved:", md_path, "\n")

cat("\n=== DONE ===\n")
cat("Outputs:\n")
cat("  ", plot1_pdf, "\n")
cat("  ", plot1_png, "\n")
cat("  ", plot2_pdf, "\n")
cat("  ", plot2_png, "\n")
cat("  ", md_path, "\n")
