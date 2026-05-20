# audit-multirep-pathologies-2026-05-19.R --------------------------------------
#
# Re-audit the two queued mixing-pathology claims against root-invariant
# scoring, after the discovery that the bipartition / topology-hash
# infrastructure used in the original diagnosis was root-dependent.
#
# Claims:
#   1. Rep 01 -- TL collapse. Chain TL collapses to ~3 vs truth 13.5.
#      TL is root-invariant; expected to survive.
#   2. Rep 05 -- topology mode-trap / valley. TL fine but truth tree
#      logPost +106 nats over chain MAP. Diagnosis based on (a) low
#      AC-bipartition support (legacy prop.part; root-dependent) and
#      (b) low unique-topology count via topo_hash (FNV-1a of parent
#      vector; root-dependent).
#
# For each rep we compute:
#   - TL trajectory: min/med/max/last-100-mean and full trace (saved)
#   - AC bipartition support: LEGACY (prop.part) vs CORRECTED (Splits)
#   - Unique topology count: LEGACY (topo_hash, root-dependent) vs
#     CORRECTED (canonical-splits hash, root-invariant)
#   - Root-configuration distribution (RootConfigCounts)
#
# Output:
#   - audit-multirep-pathologies-2026-05-19.txt (per-rep tables)
#   - audit-multirep-pathologies-2026-05-19.rds (full numeric results)

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/ecology/simulations/sim3-scoring.R")

OUT_DIR <- "dev/red-team"
TXT <- file.path(OUT_DIR, "audit-multirep-pathologies-2026-05-19.txt")
RDS <- file.path(OUT_DIR, "audit-multirep-pathologies-2026-05-19.rds")
sink(TXT, split = TRUE)

cat("=== Audit: multirep-v3 rep01 + rep05 against root-invariant scoring ===\n")
cat(sprintf("Date: %s\n", Sys.time()))

# Truth setup matches sim3-multirep.R
nEco <- 120L; nBase <- 360L; phi <- 4
stemBr <- 0.30; rootBr <- 0.15; tipBr <- 0.5
source("inst/ecology/simulations/sim3-helpers.R")
truthTree <- .BuildConvergentTree(tipBranch = tipBr,
                                   stemBranch = stemBr, rootBranch = rootBr)
truthTL <- sum(truthTree$edge.length)
trueSplit  <- c(paste0("A", 1:4), paste0("C", 1:4))
wrongSplit <- c(paste0("A", 1:4), paste0("B", 1:4))
cat(sprintf("Truth tree: %d tips, TL = %.3f, AC split = {%s}\n",
            length(truthTree$tip.label), truthTL,
            paste(trueSplit, collapse = ",")))

# Root-invariant canonical-splits hash for unrooted topology.
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

auditOne <- function(repId) {
  d <- sprintf("inst/ecology/simulations/multirep-v3-results/rep%s", repId)
  cat(sprintf("\n========== REP %s ==========\n", repId))
  resA <- readRDS(file.path(d, "aware-result.rds"))
  resA <- RelabelEcology(resA)
  resB <- readRDS(file.path(d, "blind-result.rds"))

  out <- list(rep = repId)

  for (chain in list(list(tag = "AWARE", res = resA),
                     list(tag = "BLIND", res = resB))) {
    tag <- chain$tag
    res <- chain$res
    sdf <- as.data.frame(res$samples)
    trees <- res$trees
    if (length(trees) == 0L || nrow(sdf) == 0L) {
      cat(sprintf("\n[%s] empty sample table -- skipping\n", tag))
      next
    }
    cat(sprintf("\n--- [%s] %d samples, %d saved trees ---\n",
                tag, nrow(sdf), length(trees)))

    # TL trajectory (root-invariant).
    tl <- sdf$tree_length
    cat(sprintf("  TL: min=%.3f q25=%.3f med=%.3f q75=%.3f max=%.3f last100mean=%.3f  truth=%.3f\n",
                min(tl), quantile(tl, 0.25), median(tl),
                quantile(tl, 0.75), max(tl), mean(tail(tl, 100)), truthTL))
    # Discard burn-in (first 25%) for support summaries
    keepFrom <- ceiling(length(trees) / 4) + 1L
    keepIdx <- if (keepFrom > length(trees)) integer(0)
               else seq.int(keepFrom, length(trees))
    trKeep <- trees[keepIdx]
    class(trKeep) <- "multiPhylo"
    cat(sprintf("  post-burnin trees retained: %d / %d\n",
                length(trKeep), length(trees)))

    # AC support: legacy vs corrected
    hasAC_leg  <- HasBipartLegacy(trKeep, trueSplit)
    hasAC_new  <- HasBipartSplits(trKeep, trueSplit)
    hasAB_leg  <- HasBipartLegacy(trKeep, wrongSplit)
    hasAB_new  <- HasBipartSplits(trKeep, wrongSplit)
    cat(sprintf("  P(AC) legacy=%.4f  corrected=%.4f  (delta = %+0.4f)\n",
                mean(hasAC_leg), mean(hasAC_new),
                mean(hasAC_new) - mean(hasAC_leg)))
    cat(sprintf("  P(AB) legacy=%.4f  corrected=%.4f  (delta = %+0.4f)\n",
                mean(hasAB_leg), mean(hasAB_new),
                mean(hasAB_new) - mean(hasAB_leg)))

    # Unique topology counts: legacy topo_hash vs canonical splits hash
    if ("topo_hash" %in% colnames(sdf)) {
      # The topo_hash column is per-sample. We have nrow(sdf) sampled rows
      # but only length(trees) saved trees -- the column samples are
      # at every thin step, trees only at treeThin step. Use the
      # subset of topo_hash that aligns with saved trees if possible.
      th_full <- sdf$topo_hash
      # Drop burn-in proportion on the sample-row axis
      th_keep <- th_full[seq.int(ceiling(length(th_full) / 4) + 1L,
                                  length(th_full))]
      nUniq_leg <- length(unique(th_keep))
    } else {
      nUniq_leg <- NA_integer_
    }
    canHashes <- vapply(trKeep, canonicalSplitsHash, character(1))
    nUniq_new <- length(unique(canHashes))
    cat(sprintf("  Unique topo (samples)   legacy topo_hash = %s  (out of %d)\n",
                if (is.na(nUniq_leg)) "NA" else format(nUniq_leg),
                if (is.na(nUniq_leg)) NA_integer_ else
                  length(seq.int(ceiling(nrow(sdf) / 4) + 1L, nrow(sdf)))))
    cat(sprintf("  Unique topo (saved trees) canonical splits hash = %d (out of %d)\n",
                nUniq_new, length(trKeep)))
    # Most-frequent canonical topology share
    mfTab <- sort(table(canHashes), decreasing = TRUE)
    cat(sprintf("  Most-frequent canonical topology share: %.4f\n",
                mfTab[1] / length(canHashes)))
    if (length(mfTab) >= 3) {
      cat(sprintf("  Top 3 canonical topology shares: %.4f  %.4f  %.4f\n",
                  mfTab[1] / length(canHashes),
                  mfTab[2] / length(canHashes),
                  mfTab[3] / length(canHashes)))
    }

    # Root-config distribution (how much root-shuffling does the chain do?)
    rc <- RootConfigCounts(trKeep)
    cat(sprintf("  Root configurations seen (distinct positions): %d\n",
                length(rc)))
    rcSorted <- sort(rc, decreasing = TRUE)
    cat(sprintf("  Top root-config share: %.4f\n",
                rcSorted[1] / sum(rcSorted)))
    if (length(rcSorted) >= 3) {
      cat(sprintf("  Top 3 root-config shares: %.4f  %.4f  %.4f\n",
                  rcSorted[1] / sum(rcSorted),
                  rcSorted[2] / sum(rcSorted),
                  rcSorted[3] / sum(rcSorted)))
    }
    # How often is the root inside the AC tip set?
    rootIn <- RootInsideSet(trKeep, trueSplit)
    cat(sprintf("  Root inside AC tipset: %d / %d (%.3f%%)  NA/crosses: %d\n",
                sum(rootIn, na.rm = TRUE), length(rootIn),
                100 * mean(rootIn, na.rm = TRUE), sum(is.na(rootIn))))

    out[[tag]] <- list(
      tl_summary = c(min = min(tl), median = median(tl), max = max(tl),
                     last100mean = mean(tail(tl, 100)), truth = truthTL),
      tl_trace = tl,
      pAC_legacy    = mean(hasAC_leg),
      pAC_corrected = mean(hasAC_new),
      pAB_legacy    = mean(hasAB_leg),
      pAB_corrected = mean(hasAB_new),
      nUniq_legacy_topohash    = nUniq_leg,
      nUniq_corrected_splits   = nUniq_new,
      nTreesKept = length(trKeep),
      nSamplesKept_topohash =
        if (is.na(nUniq_leg)) NA_integer_
        else length(seq.int(ceiling(nrow(sdf) / 4) + 1L, nrow(sdf))),
      top1_canonical_share = unname(mfTab[1] / length(canHashes)),
      rootConfig_distinct = length(rc),
      rootConfig_top1_share = unname(rcSorted[1] / sum(rcSorted)),
      rootInside_AC_frac = mean(rootIn, na.rm = TRUE),
      rootInside_AC_NA = sum(is.na(rootIn))
    )
  }
  out
}

results <- list()
for (repId in c("01", "05")) {
  results[[repId]] <- auditOne(repId)
}

saveRDS(results, RDS)
cat(sprintf("\nSaved RDS: %s\n", RDS))
cat("=== done ===\n")
sink()
cat(sprintf("Wrote %s\n", TXT))
