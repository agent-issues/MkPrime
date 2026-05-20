# analyse.R -------------------------------------------------------------------
# Post-processing for Step 3b mode-persistence pilot.
# Reads chain.log, computes all diagnostics, writes REPORT.md and phi-trace.png.
#
# Run from the worktree root:
#   Rscript -e "setwd('C:/Users/pjjg18/GitHub/mkp/.claude/worktrees/ecology-aware'); source('dev/pilots/2026-05-13-step3b-mode-persistence-A/analyse.R')"

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("coda")
})
`%||%` <- function(a, b) if (!is.null(a)) a else b

PILOT_DIR <- "dev/pilots/2026-05-13-step3b-mode-persistence-A"
logFile   <- file.path(PILOT_DIR, "chain.log")

if (!file.exists(logFile)) stop("Log file not found: ", logFile)

samp_all <- ReadMkLog(logFile)
cat("Total rows in log:", nrow(samp_all), "\n")
cat("Columns:", paste(colnames(samp_all), collapse = ", "), "\n\n")

# Discard first 25% as burnin
n        <- nrow(samp_all)
discard_n <- ceiling(n / 4)
samp      <- samp_all[seq.int(discard_n + 1L, n), , drop = FALSE]
cat("Retained (post 25% discard):", nrow(samp), "rows\n\n")

qs <- function(x, p) stats::quantile(x, p, names = FALSE)

phi_col <- if ("phi"     %in% colnames(samp)) samp[, "phi"]     else NULL
pi0_col <- if ("pi0"     %in% colnames(samp)) samp[, "pi0"]     else NULL
th_col  <- if ("theta_1" %in% colnames(samp)) samp[, "theta_1"] else NULL
tl      <- samp[, "tree_length"]

# ---- 1. Mode-crossing check -------------------------------------------------
cat("=== 1. Mode-crossing check ===\n")
if (!is.null(phi_col)) {
  log_phi     <- log(phi_col)
  n_below1    <- sum(phi_col < 1)
  n_above1    <- sum(phi_col >= 1)
  # Zero-crossings: consecutive sign changes of log(phi)
  signs       <- sign(log_phi)
  crossings   <- sum(diff(signs) != 0)
  first_phi   <- phi_col[1]
  last_phi    <- phi_col[length(phi_col)]
  first_th    <- if (!is.null(th_col)) th_col[1] else NA
  last_th     <- if (!is.null(th_col)) th_col[length(th_col)] else NA
  min_phi     <- min(phi_col)
  max_phi     <- max(phi_col)

  cat(sprintf("  phi < 1  (reflected mode): %d samples\n", n_below1))
  cat(sprintf("  phi >= 1 (true mode)     : %d samples\n", n_above1))
  cat(sprintf("  Zero-crossings of log(phi): %d\n", crossings))
  cat(sprintf("  First sample: phi=%.4f, theta_1=%.4f\n", first_phi,
              ifelse(is.na(first_th), NaN, first_th)))
  cat(sprintf("  Last  sample: phi=%.4f, theta_1=%.4f\n", last_phi,
              ifelse(is.na(last_th), NaN, last_th)))
  cat(sprintf("  phi range: [%.4f, %.4f]\n\n", min_phi, max_phi))

  # Trace plot: save PNG
  pngFile <- file.path(PILOT_DIR, "phi-trace.png")
  png(pngFile, width = 900, height = 400)
  plot(log_phi, type = "l", col = "steelblue",
       xlab = "Retained sample index", ylab = "log(phi)",
       main = "Step 3b: log(phi) trace — 100k chain (post-burnin)")
  abline(h = 0, col = "red", lty = 2)
  legend("topright",
         legend = sprintf("Crossings: %d  |  phi<1: %d  |  phi>=1: %d",
                          crossings, n_below1, n_above1),
         bty = "n", cex = 0.85)
  dev.off()
  cat("  Trace plot saved to", pngFile, "\n\n")
} else {
  cat("  phi column missing — cannot assess mode-crossing.\n\n")
  n_below1 <- NA; n_above1 <- NA; crossings <- NA
  min_phi <- NA; max_phi <- NA; first_phi <- NA; last_phi <- NA
  first_th <- NA; last_th <- NA
}

# ---- 2. Settled mode --------------------------------------------------------
cat("=== 2. Settled mode ===\n")
if (!is.null(phi_col) && !is.null(th_col)) {
  med_phi <- median(phi_col)
  med_th  <- median(th_col)
  cat(sprintf("  Median phi    = %.4f  (truth: 4, reflected: 0.25)\n", med_phi))
  cat(sprintf("  Median theta_1= %.4f  (truth: ~0.97, reflected: ~0.03)\n", med_th))
  if (med_phi < 1) {
    cat("  --> Chain settled in the REFLECTED mode (phi < 1, theta_1 small)\n\n")
  } else {
    cat("  --> Chain settled in the TRUE mode (phi >= 1, theta_1 large)\n\n")
  }
} else {
  cat("  phi or theta_1 missing.\n\n")
}

# ---- 3. Mixing health: ESS --------------------------------------------------
cat("=== 3. Mixing health (ESS) ===\n")
ess_cols <- c("phi", "pi0", "theta_1", "tree_length")
ess_vals <- setNames(rep(NA_real_, length(ess_cols)), ess_cols)
for (nm in ess_cols) {
  if (nm %in% colnames(samp)) {
    ess_val      <- coda::effectiveSize(samp[, nm])
    ess_vals[nm] <- ess_val
    cat(sprintf("  %-14s ESS = %6.1f  (n=%d)  [20k pilot: phi=8.9, pi0=3.7, theta_1=10.2, tl=8.1 / 131]\n",
                nm, ess_val, nrow(samp)))
  } else {
    cat(sprintf("  %-14s MISSING\n", nm))
  }
}
cat("\n")

# ---- 4. Tail check: tree_length ---------------------------------------------
cat("=== 4. Tail check (tree_length) ===\n")
cat(sprintf("  median = %.3f  (20k pilot: ~29, truth: 12.7)\n",    median(tl)))
cat(sprintf("  IQR    = [%.3f, %.3f]  (20k pilot: [20.1, 42.3])\n",
            qs(tl, 0.25), qs(tl, 0.75)))
cat(sprintf("  max    = %.3f  (20k pilot: ~70)\n\n",               max(tl)))

# ---- 5. Comparison to 20k pilot ---------------------------------------------
cat("=== 5. Comparison to 20k pilot ===\n")
if (!is.na(crossings)) {
  if (crossings == 0) {
    cat("  100k chain: ZERO mode crossings — single-mode persistence confirmed.\n")
    cat("  Post-hoc relabelling (Option 1) is sufficient; prior symmetry-break not required by chain dynamics.\n\n")
  } else {
    cat(sprintf("  100k chain: %d mode crossings detected — chain is NOT mode-stable.\n", crossings))
    cat("  Option 2 (phi >= 1 constraint) or Option 3 (asymmetric Beta) should be considered.\n\n")
  }
}

# ---- Write REPORT.md --------------------------------------------------------
rpt <- file.path(PILOT_DIR, "REPORT.md")

# Collect pi0 summary if present
pi0_med <- if (!is.null(pi0_col)) median(pi0_col) else NA

report_text <- sprintf(
'# Step 3b Mode-Persistence Report — 100k Sim 3 chain

**Date:** %s
**Plan:** `dev/plans/2026-05-13-1100-v2-asymmetric-slab-prior-and-gamma-ridge.md`, Open question: label-switching

---

## Configuration

| Parameter | Value |
|---|---|
| Script | `dev/pilots/2026-05-13-step3b-mode-persistence-A/run.R` |
| Seed | 20260513 (one more than the 20k pilot — independent chain) |
| kEco | 2 (ecology 0 = C/D clades; ecology 1 = A/B clades) |
| nTheta | 1 (one non-reference ecology) |
| nEco | 60 neomorphic characters |
| nBase | 180 transformational characters |
| phi_truth | 4 |
| stem / root | 0.10 / 0.15 (Goldilocks) |
| nIter | 100,000 |
| thin | 40 |
| nChains / nRuns | 1 / 1 |
| minWarmup / maxWarmup | 2,000 / 5,000 |
| Initial state | default (.InitState) |
| model | ecologyAware=TRUE, magnitudeMode="global", kPrimePrior="geometric" |
| Log file | `dev/pilots/2026-05-13-step3b-mode-persistence-A/chain.log` |
| RDS | `dev/pilots/2026-05-13-step3b-mode-persistence-A/result.rds` |

**Post-burnin samples:** %d raw rows → %d retained after 25%% discard.

---

## 1. Mode-crossing check (critical)

| Metric | Value |
|---|---|
| phi < 1 (reflected mode) | %s samples |
| phi >= 1 (true mode) | %s samples |
| Zero-crossings of log(phi) | **%s** |
| First sample: phi / theta_1 | %.4f / %.4f |
| Last sample: phi / theta_1 | %.4f / %.4f |
| phi range: min / max | %.4f / %.4f |

%s

Trace plot: `phi-trace.png`

---

## 2. Settled mode

| Parameter | Median | Truth | Reflected |
|---|---|---|---|
| phi | %.4f | 4 | 0.25 |
| theta_1 | %.4f | ~0.97 | ~0.03 |

**%s**

---

## 3. Mixing health (ESS)

| Parameter | ESS | n retained | ESS/n | 20k pilot ESS/n |
|---|---|---|---|---|
| phi | %.1f | %d | %.3f | 8.9/131 = 0.068 |
| pi0 | %.1f | %d | %.3f | 3.7/131 = 0.028 |
| theta_1 | %.1f | %d | %.3f | 10.2/131 = 0.078 |
| tree_length | %.1f | %d | %.3f | 8.1/131 = 0.062 |

---

## 4. Tail check (tree_length)

| Metric | 100k chain | 20k pilot | Truth |
|---|---|---|---|
| median | %.3f | ~29 | 12.7 |
| IQR | [%.3f, %.3f] | [20.1, 42.3] | — |
| max | %.3f | ~70 | — |

---

## 5. Comparison to 20k pilot

%s

---

## Conclusion

%s
',
  format(Sys.Date()),
  n, nrow(samp),
  # mode crossing block
  format(n_below1, big.mark = ","), format(n_above1, big.mark = ","),
  ifelse(is.na(crossings), "NA", as.character(crossings)),
  ifelse(is.na(first_phi), NaN, first_phi),
  ifelse(is.na(first_th), NaN, first_th),
  ifelse(is.na(last_phi), NaN, last_phi),
  ifelse(is.na(last_th), NaN, last_th),
  ifelse(is.na(min_phi), NaN, min_phi),
  ifelse(is.na(max_phi), NaN, max_phi),
  if (!is.na(crossings) && crossings == 0)
    "**No mode crossings.** The chain never crossed log(phi) = 0."
  else if (!is.na(crossings))
    sprintf("**%d mode crossings detected.** The chain switched between phi modes.", crossings)
  else
    "Mode-crossing status unknown (phi column missing).",
  # settled mode
  ifelse(is.na(median(phi_col %||% NA_real_)), NaN, median(phi_col)),
  ifelse(is.na(median(th_col %||% NA_real_)), NaN, median(th_col)),
  if (!is.null(phi_col) && median(phi_col) < 1)
    "Chain settled in the **reflected mode** (phi < 1, theta_1 near zero)."
  else if (!is.null(phi_col))
    "Chain settled in the **true mode** (phi >= 1, theta_1 near 1)."
  else
    "Settled mode unknown.",
  # ESS block
  ess_vals["phi"], nrow(samp), ess_vals["phi"] / nrow(samp),
  ess_vals["pi0"], nrow(samp), ess_vals["pi0"] / nrow(samp),
  ess_vals["theta_1"], nrow(samp), ess_vals["theta_1"] / nrow(samp),
  ess_vals["tree_length"], nrow(samp), ess_vals["tree_length"] / nrow(samp),
  # tail check
  median(tl), qs(tl, 0.25), qs(tl, 0.75), max(tl),
  # comparison to 20k pilot
  if (!is.na(crossings) && crossings == 0)
    "100k chain: **zero mode crossings**. Single-chain mode persistence is confirmed. The 20k pilot result (phi inverted to reflected mode, theta_1 near zero) is replicated: the chain consistently finds one of the two symmetric modes and stays there. Post-hoc relabelling (Option 1) is sufficient."
  else if (!is.na(crossings))
    sprintf("100k chain: **%d mode crossings**. The chain does switch modes, meaning post-hoc relabelling is unreliable. Option 2 (phi >= 1 constraint) or Option 3 (asymmetric Beta) is needed.", crossings)
  else
    "Comparison unavailable (phi column missing).",
  # conclusion
  if (!is.na(crossings) && crossings == 0)
    "A single 100k chain shows zero log(phi) zero-crossings. Mode persistence is confirmed. Post-hoc relabelling (Option 1 from the plan) is sufficient to handle the phi <-> 1/phi symmetry. No prior symmetry-breaking is required for single-chain inference; multi-rep analyses should apply relabelling before cross-replicate summary."
  else if (!is.na(crossings))
    sprintf("%d zero-crossings detected. The chain does not stay in one mode. The (phi >= 1) prior constraint (Option 2) should be implemented before the HPC multi-rep run.", crossings)
  else
    "Conclusion unavailable — phi column missing from log."
)

writeLines(report_text, rpt)
cat("REPORT.md written to", rpt, "\n")
cat("Analysis complete.\n")
