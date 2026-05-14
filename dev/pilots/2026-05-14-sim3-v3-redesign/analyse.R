# Analyse Sim 3 v3 Bayesian pilot
#
# Compares blind vs aware chains on the v3 redesigned dataset.
# Reports:
#   1. Topology recovery: P(true AC), P(wrong AB), CID-to-truth
#   2. Mixing: ESS, modal-topology fraction, chain logLik vs truth
#   3. Parameter recovery (aware only): phi, pi0, theta, tl, rate_neo
#
# Run after bayesian-pilot.R completes:
#   Rscript -e "setwd('C:/Users/pjjg18/GitHub/mkp/.claude/worktrees/ecology-aware'); source('dev/pilots/2026-05-14-sim3-v3-redesign/analyse.R')"

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(TreeTools)
  library(TreeDist)
})

PILOT_DIR <- "dev/pilots/2026-05-14-sim3-v3-redesign"
saved <- readRDS(file.path(PILOT_DIR, "result.rds"))

trueTree <- saved$tree
truthTL  <- saved$truthTL
phi_true <- saved$config$phi
nEco <- saved$config$nEco

cat("=== Sim 3 v3 pilot analysis ===\n")
cat(sprintf("Config: nEco=%d nBase=%d phi=%g stem=%.2f root=%.2f\n",
            nEco, saved$config$nBase, phi_true,
            saved$config$stemBr, saved$config$rootBr))
cat(sprintf("Truth tl = %.2f\n", truthTL))

# Helper: tip-set membership
hasSplit <- function(tr, tips) {
  cl <- ape::prop.part(tr)
  labs <- attr(cl, "labels")
  target <- which(labs %in% tips)
  any(sapply(cl, function(p) setequal(p, target)))
}
trueAC <- c(paste0("A", 1:4), paste0("C", 1:4))
wrongAB <- c(paste0("A", 1:4), paste0("B", 1:4))
trueAB_wrong <- wrongAB  # alias

# Process one chain
process_chain <- function(res, label, awareLog = FALSE) {
  cat(sprintf("\n--- %s chain ---\n", label))
  trees <- res$trees
  burnin <- floor(0.25 * length(trees))
  trees_pb <- trees[(burnin + 1):length(trees)]
  cat(sprintf("Trees: %d (post-burnin)\n", length(trees_pb)))

  # Topology
  cidn <- as.numeric(TreeDist::ClusteringInfoDistance(trees_pb, trueTree,
                                                      normalize = TRUE))
  cat(sprintf("Normalised CID to truth: min=%.3f  median=%.3f  max=%.3f\n",
              min(cidn), median(cidn), max(cidn)))
  cat(sprintf("Exact matches (CID=0): %d / %d (%.1f%%)\n",
              sum(cidn == 0), length(cidn),
              100 * sum(cidn == 0) / length(cidn)))

  pAC <- mean(sapply(trees_pb, hasSplit, tips = trueAC))
  pAB <- mean(sapply(trees_pb, hasSplit, tips = wrongAB))
  cat(sprintf("P(true bipartition AC together)  = %.3f\n", pAC))
  cat(sprintf("P(wrong bipartition AB together) = %.3f\n", pAB))

  # Modal topology
  rfHash <- sapply(trees_pb, function(tr) {
    paste(sort(unlist(lapply(ape::prop.part(tr), function(p) {
      paste(sort(p), collapse = "_")
    }))), collapse = "|")
  })
  nUnique <- length(unique(rfHash))
  topCount <- max(table(rfHash))
  cat(sprintf("Unique topologies: %d / %d, modal holds %d (%.1f%%)\n",
              nUnique, length(rfHash), topCount,
              100 * topCount / length(rfHash)))

  # ESS via batch means
  samp <- res$samples
  if (NROW(samp) == 0L && !is.null(res$logFile)) {
    samp <- ReadMkLog(res$logFile)
  }
  samp_pb <- samp[(burnin + 1):nrow(samp), , drop = FALSE]
  cat(sprintf("Log samples: %d post-burnin\n", nrow(samp_pb)))

  ess_bm <- function(x) {
    x <- x[!is.na(x)]
    b <- max(1L, floor(sqrt(length(x))))
    nb <- floor(length(x) / b)
    if (nb < 2) return(NA_real_)
    bm <- sapply(seq_len(nb), function(i) mean(x[((i-1)*b+1):(i*b)]))
    s2 <- var(bm) * b
    if (s2 <= 0 || !is.finite(s2)) return(length(x))
    min(length(x) * var(x) / (length(x) * s2), length(x))
  }
  cat(sprintf("ESS tl:     %.1f\n", ess_bm(samp_pb[, "tree_length"])))
  cat(sprintf("ESS rate_neo: %.1f\n", ess_bm(samp_pb[, "rate_neo"])))
  if (awareLog) {
    cat(sprintf("ESS phi:    %.1f\n", ess_bm(samp_pb[, "phi"])))
    cat(sprintf("ESS pi0:    %.1f\n", ess_bm(samp_pb[, "pi0"])))
    if ("theta_1" %in% colnames(samp_pb)) {
      cat(sprintf("ESS theta:  %.1f\n", ess_bm(samp_pb[, "theta_1"])))
    }
  }

  # Parameter medians
  cat(sprintf("Posterior tl: median=%.2f IQR [%.2f, %.2f]  (truth %.2f)\n",
              median(samp_pb[, "tree_length"]),
              quantile(samp_pb[, "tree_length"], 0.25),
              quantile(samp_pb[, "tree_length"], 0.75), truthTL))
  cat(sprintf("Posterior rate_neo: median=%.3f IQR [%.3f, %.3f] (truth 1.0)\n",
              median(samp_pb[, "rate_neo"]),
              quantile(samp_pb[, "rate_neo"], 0.25),
              quantile(samp_pb[, "rate_neo"], 0.75)))

  if (awareLog) {
    cat(sprintf("Posterior phi: median=%.3f IQR [%.3f, %.3f] (truth %g)\n",
                median(samp_pb[, "phi"]),
                quantile(samp_pb[, "phi"], 0.25),
                quantile(samp_pb[, "phi"], 0.75), phi_true))
    cat(sprintf("Posterior pi0: median=%.3f IQR [%.3f, %.3f] (truth 0.75)\n",
                median(samp_pb[, "pi0"]),
                quantile(samp_pb[, "pi0"], 0.25),
                quantile(samp_pb[, "pi0"], 0.75)))
    if ("theta_1" %in% colnames(samp_pb)) {
      cat(sprintf("Posterior theta: median=%.3f IQR [%.3f, %.3f] (truth 1.0)\n",
                  median(samp_pb[, "theta_1"]),
                  quantile(samp_pb[, "theta_1"], 0.25),
                  quantile(samp_pb[, "theta_1"], 0.75)))
    }
  }

  cat(sprintf("Chain logLik: median=%.2f  max=%.2f\n",
              median(samp_pb[, "log_likelihood"]),
              max(samp_pb[, "log_likelihood"])))

  list(cidn = cidn, pAC = pAC, pAB = pAB, modalFrac = topCount / length(rfHash),
       samp_pb = samp_pb)
}

resB <- process_chain(saved$resBlind, "BLIND", awareLog = FALSE)
resA <- process_chain(saved$resAware, "AWARE", awareLog = TRUE)

# ---- Side-by-side comparison ----
cat("\n=== AWARE vs BLIND comparison ===\n")
cat(sprintf("%-30s %-10s %-10s %-10s\n", "Metric", "Blind", "Aware", "Δ (A−B)"))
cat(sprintf("%-30s %-10.3f %-10.3f %+.3f\n", "P(true AC bipartition)",
            resB$pAC, resA$pAC, resA$pAC - resB$pAC))
cat(sprintf("%-30s %-10.3f %-10.3f %+.3f\n", "P(wrong AB bipartition)",
            resB$pAB, resA$pAB, resA$pAB - resB$pAB))
cat(sprintf("%-30s %-10.3f %-10.3f %+.3f\n", "min CID to truth",
            min(resB$cidn), min(resA$cidn), min(resA$cidn) - min(resB$cidn)))
cat(sprintf("%-30s %-10.3f %-10.3f %+.3f\n", "median CID to truth",
            median(resB$cidn), median(resA$cidn), median(resA$cidn) - median(resB$cidn)))
cat(sprintf("%-30s %-10.3f %-10.3f %+.3f\n", "modal topology fraction",
            resB$modalFrac, resA$modalFrac, resA$modalFrac - resB$modalFrac))

cat("\n=== Verdict ===\n")
delta_pAC <- resA$pAC - resB$pAC
delta_med_cid <- median(resA$cidn) - median(resB$cidn)
if (delta_pAC > 0.1 && delta_med_cid < -0.05) {
  cat("AWARE outperforms BLIND meaningfully — paper-worthy contrast.\n")
} else if (delta_pAC > 0.05) {
  cat("Aware improvement modest. Consider stronger v3 (4x chars or 4x-phi6).\n")
} else if (resA$modalFrac < 0.05) {
  cat("Aware still failing to mix (modal topology < 5%). Need mixing fixes.\n")
} else {
  cat("Aware ≈ blind. Either ecology layer not helping or v3 redesign too easy.\n")
}
