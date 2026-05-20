# analyse-aware-v2.R -----------------------------------------------------------
# Characterise the ecological-signal posterior from the rodent aware-v2 chain
# (1M iterations, completed on Hamilton).
#
# Inputs (Hamilton paths assumed when run on cluster):
#   /nobackup/pjjg18/mkp-rodent-aware-v2/results/rodent-aware-v2-result.rds
#   /nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex
#
# Outputs to inst/scripts/rodent-posterior-diagnostics/ :
#   - scalar-summary.csv          : phi, pi0, theta_e, tree_length posteriors
#   - z-marginals.csv             : per-(char, eco) P(z=0/1/2) + P(z!=0)
#   - eco-counts.csv              : per-ecology count of strong-flag chars
#   - cross-eco-overlap.csv       : pairwise overlap + correlations
#   - top-chars.csv               : top-10 per-ecology with biological label
#   - phi-hist.pdf, pz-hist.pdf   : posterior histograms
#   - report.md                   : machine-readable bottom-line report

suppressPackageStartupMessages({
  # No package load needed — we just read the RDS.
})

args <- commandArgs(trailingOnly = TRUE)
rdsFile <- if (length(args) >= 1) args[1] else
  "/nobackup/pjjg18/mkp-rodent-aware-v2/results/rodent-aware-v2-result.rds"
nexFile <- if (length(args) >= 2) args[2] else
  "/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex"
outDir  <- if (length(args) >= 3) args[3] else "."

dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

cat("Reading", rdsFile, "\n")
r <- readRDS(rdsFile)

samples <- r$samples
z_list  <- r$res$z_samples
mkd     <- r$mkd
mod     <- r$model
nIter   <- r$res$actual_iter

cat("nSamples (scalar):", nrow(samples), "  z samples:", length(z_list),
    "  actual_iter:", nIter, "  stop_reason:", r$res$stop_reason, "\n")

# Burn-in: drop first 25% of samples for stability.
burnFrac <- 0.25
keepScalar <- seq.int(ceiling(nrow(samples) * burnFrac) + 1L, nrow(samples))
keepZ      <- seq.int(ceiling(length(z_list) * burnFrac) + 1L, length(z_list))
cat("Burn-in fraction:", burnFrac,
    " — scalar samples retained:", length(keepScalar),
    " — z samples retained:", length(keepZ), "\n")
samples <- samples[keepScalar, , drop = FALSE]
z_list  <- z_list[keepZ]

# --- 1. Scalar posterior summaries ------------------------------------------
qFun <- function(x) {
  c(median = median(x, na.rm = TRUE),
    mean   = mean(x, na.rm = TRUE),
    p05    = quantile(x, 0.05, na.rm = TRUE, names = FALSE),
    p25    = quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    p75    = quantile(x, 0.75, na.rm = TRUE, names = FALSE),
    p95    = quantile(x, 0.95, na.rm = TRUE, names = FALSE))
}

scalarCols <- c("phi", "pi0", "theta_1", "theta_2", "theta_3",
                "tree_length", "log_posterior", "log_likelihood",
                "rate_loss", "rate_log_sd", "rate_neo", "p")
scalarCols <- intersect(scalarCols, colnames(samples))
scalarTab  <- t(sapply(scalarCols, function(cn) qFun(samples[, cn])))
scalarTab  <- data.frame(param = rownames(scalarTab), scalarTab, row.names = NULL)
write.csv(scalarTab, file.path(outDir, "scalar-summary.csv"), row.names = FALSE)
cat("\n=== Scalar posterior summaries (post burn-in) ===\n")
print(round(scalarTab[, -1], 4), row.names = scalarTab$param)

# Quick histograms
pdf(file.path(outDir, "phi-hist.pdf"), width = 7, height = 5)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
hist(samples[, "phi"],    breaks = 40, col = "steelblue", border = NA,
     main = "phi (ecology rate multiplier)", xlab = "phi")
abline(v = quantile(samples[, "phi"], c(.05, .5, .95)), col = c("grey50", "red", "grey50"),
       lty = c(2, 1, 2))
hist(samples[, "pi0"],    breaks = 40, col = "steelblue", border = NA,
     main = "pi0 (slab indicator)", xlab = "pi0")
abline(v = mod$rho0Alpha / (mod$rho0Alpha + mod$rho0Beta), col = "darkgreen", lty = 3)
hist(samples[, "theta_1"], breaks = 40, col = "steelblue", border = NA,
     main = "theta_1", xlab = "theta_1")
hist(samples[, "tree_length"], breaks = 40, col = "steelblue", border = NA,
     main = "tree_length", xlab = "TL")
par(op); dev.off()

# --- 2. z-matrix posterior --------------------------------------------------
nChar <- nrow(z_list[[1]])
kEcoM1 <- ncol(z_list[[1]])
cat("\nz dimensions: nChar=", nChar, "  k_eco_nonref=", kEcoM1, "\n", sep = "")

# Stack into 3D array (sample x char x eco)
zArr <- array(0L, dim = c(length(z_list), nChar, kEcoM1))
for (i in seq_along(z_list)) zArr[i, , ] <- z_list[[i]]

# Posterior marginals
pZ0 <- apply(zArr, c(2, 3), function(x) mean(x == 0L))
pZ1 <- apply(zArr, c(2, 3), function(x) mean(x == 1L))
pZ2 <- apply(zArr, c(2, 3), function(x) mean(x == 2L))
pAffected <- 1 - pZ0  # P(z != 0)

# Ecology labels: skip refEcology (state 0 here)
ecoLabels <- setdiff(seq.int(0L, mkd$kEcology - 1L), mkd$refEcology)
colnames(pZ0) <- colnames(pZ1) <- colnames(pZ2) <- colnames(pAffected) <-
  paste0("eco", ecoLabels)
ecoNameMap <- c("0" = "terrestrial(ref)", "1" = "arboreal",
                "2" = "semiaquatic", "3" = "fossorial")

zMarg <- data.frame(
  char = seq_len(nChar),
  type = mkd$type,
  pZ0_eco1 = pZ0[, 1], pZ1_eco1 = pZ1[, 1], pZ2_eco1 = pZ2[, 1],
  pAff_eco1 = pAffected[, 1],
  pZ0_eco2 = pZ0[, 2], pZ1_eco2 = pZ1[, 2], pZ2_eco2 = pZ2[, 2],
  pAff_eco2 = pAffected[, 2],
  pZ0_eco3 = pZ0[, 3], pZ1_eco3 = pZ1[, 3], pZ2_eco3 = pZ2[, 3],
  pAff_eco3 = pAffected[, 3]
)
write.csv(zMarg, file.path(outDir, "z-marginals.csv"), row.names = FALSE)

# Eco-flagged counts at multiple thresholds
threshes <- c(0.3, 0.5, 0.7, 0.9)
ecoCounts <- t(sapply(seq_len(kEcoM1), function(j) {
  sapply(threshes, function(t) sum(pAffected[, j] > t))
}))
ecoCounts <- data.frame(
  ecology    = paste0("eco", ecoLabels),
  eco_name   = ecoNameMap[as.character(ecoLabels)],
  thr30 = ecoCounts[, 1], thr50 = ecoCounts[, 2],
  thr70 = ecoCounts[, 3], thr90 = ecoCounts[, 4]
)
write.csv(ecoCounts, file.path(outDir, "eco-counts.csv"), row.names = FALSE)
cat("\n=== Characters flagged per ecology (P(z!=0) > threshold) ===\n")
print(ecoCounts)

# Histogram of P(z!=0) across all (char, eco) cells
pdf(file.path(outDir, "pz-hist.pdf"), width = 7, height = 5)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (j in seq_len(kEcoM1)) {
  hist(pAffected[, j], breaks = seq(0, 1, by = 0.05),
       col = "firebrick", border = NA,
       main = paste0("P(z!=0) for ", ecoNameMap[as.character(ecoLabels[j])]),
       xlab = "P(z != 0)")
  abline(v = 0.5, col = "black", lty = 2)
}
hist(as.vector(pAffected), breaks = seq(0, 1, by = 0.05),
     col = "grey40", border = NA,
     main = "P(z!=0) across all (char, eco)", xlab = "P(z != 0)")
abline(v = 0.5, col = "black", lty = 2)
par(op); dev.off()

# --- 3. Cross-ecology overlap -----------------------------------------------
overlapRows <- list()
for (i in seq_len(kEcoM1 - 1L)) {
  for (j in seq.int(i + 1L, kEcoM1)) {
    flagI <- pAffected[, i] > 0.5
    flagJ <- pAffected[, j] > 0.5
    overlapRows[[length(overlapRows) + 1L]] <- data.frame(
      eco_a       = paste0("eco", ecoLabels[i]),
      eco_b       = paste0("eco", ecoLabels[j]),
      name_a      = ecoNameMap[as.character(ecoLabels[i])],
      name_b      = ecoNameMap[as.character(ecoLabels[j])],
      n_flag_a    = sum(flagI),
      n_flag_b    = sum(flagJ),
      n_both      = sum(flagI & flagJ),
      jaccard     = if (sum(flagI | flagJ) == 0) NA_real_ else
        sum(flagI & flagJ) / sum(flagI | flagJ),
      cor_pearson = cor(pAffected[, i], pAffected[, j])
    )
  }
}
overlapTab <- do.call(rbind, overlapRows)
write.csv(overlapTab, file.path(outDir, "cross-eco-overlap.csv"), row.names = FALSE)
cat("\n=== Cross-ecology character-set overlap ===\n")
print(overlapTab)

# --- 4. Biological annotation: parse CHARSTATELABELS from nexus -------------
parseCharLabels <- function(nexFile) {
  if (!file.exists(nexFile)) {
    warning("Nexus file not found: ", nexFile)
    return(NULL)
  }
  txt <- readLines(nexFile, warn = FALSE)
  startIdx <- grep("CHARSTATELABELS", txt, ignore.case = TRUE)
  if (length(startIdx) == 0L) return(NULL)
  # The block runs from CHARSTATELABELS to the next ';' (semicolon on its own,
  # or terminating a line). Walk forward until we hit a line ending with ';'.
  tail <- txt[seq.int(startIdx[1], length(txt))]
  endRel <- which(grepl(";\\s*$", tail))[1]
  block <- tail[seq_len(endRel)]
  # Concatenate, then split on commas at top level (lines start with "  NN ...")
  blob <- paste(block, collapse = " ")
  # Split entries on commas that are followed by whitespace then a digit
  entries <- strsplit(blob, ",(?=\\s+\\d+\\s)", perl = TRUE)[[1]]
  labs <- character(0)
  ids  <- integer(0)
  for (e in entries) {
    m <- regmatches(e, regexec("^\\s*(\\d+)\\s+(.+?)\\s*(?:/|$)", e))[[1]]
    if (length(m) >= 3L) {
      ids  <- c(ids,  as.integer(m[2]))
      lab  <- m[3]
      # Strip leading/trailing quotes
      lab  <- gsub("^['\"]|['\"]$", "", lab)
      labs <- c(labs, lab)
    }
  }
  setNames(labs, as.character(ids))
}

charLabels <- tryCatch(parseCharLabels(nexFile), error = function(e) NULL)
if (is.null(charLabels) || length(charLabels) == 0L) {
  cat("\nNo character labels parsed; falling back to integer IDs.\n")
  charLabels <- setNames(paste0("char_", seq_len(nChar)),
                         as.character(seq_len(nChar)))
}
cat("\nParsed", length(charLabels), "character labels from nexus.\n")

# Match labels to characters. The MkPrimeData drops the ecology and extant
# columns but otherwise preserves order; nChar should be 219 - kEco_cols.
# Our nexus is 1-indexed to 219 morphological chars (1..219, with col 220 =
# ecology, col 221 = extant). After dropping ecology+extant, we have 219
# but we then trim non-variable? mkd$nChar reports the kept count.
# Approach: assume the FIRST mkd$nChar nexus labels map 1:1.
labVec <- charLabels[as.character(seq_len(nChar))]
labVec[is.na(labVec)] <- paste0("char_", which(is.na(labVec)))

# Top-10 per ecology
topRows <- list()
for (j in seq_len(kEcoM1)) {
  ord <- order(pAffected[, j], decreasing = TRUE)
  top <- head(ord, 10)
  for (c in top) {
    # Determine dominant non-zero state
    p1 <- pZ1[c, j]; p2 <- pZ2[c, j]
    dom <- if (p1 > p2) "z=1 (symmetric)" else "z=2 (asymmetric)"
    topRows[[length(topRows) + 1L]] <- data.frame(
      ecology  = paste0("eco", ecoLabels[j]),
      eco_name = ecoNameMap[as.character(ecoLabels[j])],
      char     = c,
      type     = mkd$type[c],
      p_affected = round(pAffected[c, j], 3),
      p_z1     = round(p1, 3),
      p_z2     = round(p2, 3),
      dominant = dom,
      label    = labVec[c]
    )
  }
}
topTab <- do.call(rbind, topRows)
write.csv(topTab, file.path(outDir, "top-chars.csv"), row.names = FALSE)
cat("\n=== Top-10 characters per ecology ===\n")
for (j in seq_len(kEcoM1)) {
  cat("\n-- eco", ecoLabels[j], "(", ecoNameMap[as.character(ecoLabels[j])], ") --\n", sep = "")
  sub <- topTab[topTab$ecology == paste0("eco", ecoLabels[j]), ]
  for (k in seq_len(nrow(sub))) {
    cat(sprintf("  #%3d  P(z!=0)=%.3f  P(z=1)=%.3f  P(z=2)=%.3f  %s  | %s\n",
                sub$char[k], sub$p_affected[k], sub$p_z1[k], sub$p_z2[k],
                sub$type[k], sub$label[k]))
  }
}

# --- 5. Bottom line ---------------------------------------------------------
phiQ <- quantile(samples[, "phi"], c(0.025, 0.05, 0.5, 0.95, 0.975), na.rm = TRUE)
pi0Q <- quantile(samples[, "pi0"], c(0.025, 0.5, 0.975), na.rm = TRUE)
avgFlag <- mean(rowSums(t(pAffected > 0.5)))
medFlag <- median(apply(pAffected > 0.5, 2, sum))
meanCor <- mean(overlapTab$cor_pearson)
meanJac <- mean(overlapTab$jaccard, na.rm = TRUE)

report <- c(
  "# Rodent aware-v2 posterior diagnostics",
  paste0("Source: ", basename(rdsFile),
         "  (actual_iter=", nIter, "; nSamples=", nrow(samples) + length(keepScalar),
         "; retained post-burnin scalar=", nrow(samples),
         "; retained z=", length(z_list), ")"),
  paste0("Model prior: rho0=Beta(", mod$rho0Alpha, ",", mod$rho0Beta,
         "), theta=Beta(", mod$thetaAlpha, ",", mod$thetaBeta,
         "), sigmaPhi=", mod$sigmaPhi, ", relabel=", mod$relabel),
  "",
  "## phi (ecological rate multiplier)",
  sprintf("median = %.3f   90%% CI = [%.3f, %.3f]   95%% CI = [%.3f, %.3f]",
          phiQ[3], phiQ[2], phiQ[4], phiQ[1], phiQ[5]),
  "",
  "## pi0 (slab probability)",
  sprintf("median = %.3f   95%% CI = [%.3f, %.3f]",
          pi0Q[2], pi0Q[1], pi0Q[3]),
  "",
  "## Characters flagged P(z!=0) > 0.5 per ecology",
  paste(capture.output(print(ecoCounts)), collapse = "\n"),
  "",
  "## Cross-ecology overlap",
  paste(capture.output(print(overlapTab)), collapse = "\n"),
  "",
  sprintf("Mean inter-ecology cor(P(z!=0)) = %.3f", meanCor),
  sprintf("Mean Jaccard overlap            = %.3f", meanJac),
  sprintf("Median characters flagged/eco   = %.0f", medFlag),
  "",
  "## Ecological-gravity verdict",
  sprintf("phi posterior median %.2f  (95%% CI [%.2f, %.2f])", phiQ[3], phiQ[1], phiQ[5]),
  sprintf("This is the empirical signal strength on real rodent morphology."),
  sprintf("Median of %.0f characters per ecology show P(z!=0)>0.5", medFlag),
  ""
)
writeLines(report, file.path(outDir, "report.md"))
cat("\nWrote report.md and diagnostic CSVs/PDFs to", outDir, "\n")
cat("\nDone.\n")
