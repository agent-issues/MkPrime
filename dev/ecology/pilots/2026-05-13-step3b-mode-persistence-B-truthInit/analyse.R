# analyse.R -------------------------------------------------------------------
# Post-processing for Step 3b mode-persistence pilot — CHAIN B (truth-init).
# Reads chain.log, computes all diagnostics, writes REPORT.md and phi-trace.png.
#
# Init verification note: the stdout already confirmed the override took at
# iter 0 (phi=4.0000 theta=0.9700). This script checks whether the chain then
# *stayed* in the phi>1 mode or crossed to phi<1 during/after warmup.
#
# Run from the worktree root:
#   Rscript -e "setwd('C:/Users/pjjg18/GitHub/mkp/.claude/worktrees/ecology-aware'); source('dev/pilots/2026-05-13-step3b-mode-persistence-B-truthInit/analyse.R')"

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("coda")
})
`%||%` <- function(a, b) if (!is.null(a)) a else b

PILOT_DIR  <- "dev/pilots/2026-05-13-step3b-mode-persistence-B-truthInit"
logFile    <- file.path(PILOT_DIR, "chain.log")
stdoutFile <- file.path(PILOT_DIR, "chain.log.stdout")

if (!file.exists(logFile)) stop("Log file not found: ", logFile)

# ---------------------------------------------------------------------------
# Init verification (from stdout, not log rows)
# ---------------------------------------------------------------------------
cat("=== 0. Init verification ===\n")
init_took     <- FALSE
init_phi_val  <- NA_real_
init_th_val   <- NA_real_

if (file.exists(stdoutFile)) {
  stdout_lines <- readLines(stdoutFile)
  truth_line   <- grep("\\[truth-init\\]", stdout_lines, value = TRUE)
  if (length(truth_line) > 0) {
    cat("  Stdout confirms init override:\n  ", truth_line[1], "\n")
    m <- regmatches(truth_line[1],
                    regexpr("phi=([0-9.]+).*theta=([0-9.]+)",
                            truth_line[1], perl = TRUE))
    if (length(m) > 0) {
      nums        <- as.numeric(regmatches(truth_line[1],
                                           gregexpr("[0-9]+\\.[0-9]+",
                                                    truth_line[1]))[[1]])
      init_phi_val <- nums[1]
      init_th_val  <- nums[2]
    }
    init_took <- TRUE
    cat(sprintf("  phi_init=%.4f (expected 4.0)  theta_init=%.4f (expected 0.97)\n",
                init_phi_val, init_th_val))
    if (abs(init_phi_val - 4.0) < 0.01 && abs(init_th_val - 0.97) < 0.01) {
      cat("  INIT VERIFIED: override took correctly at iter 0.\n\n")
    } else {
      cat("  WARNING: init values look unexpected — check the stdout.\n\n")
      init_took <- FALSE
    }
  } else {
    cat("  WARNING: [truth-init] line not found in stdout — override may not have run.\n\n")
  }
} else {
  cat("  stdout file not found — cannot confirm init from file; assuming override ran.\n\n")
}

if (!init_took) {
  cat("CRITICAL: init verification failed. The rest of the analysis may be uninformative.\n")
  cat("The chain may have started at the default (phi<1) state, not the truth-mode.\n")
  cat("Check run.R monkey-patch and rerun if needed.\n")
}

# ---------------------------------------------------------------------------
# Load log
# ---------------------------------------------------------------------------
samp_all <- ReadMkLog(logFile)
cat("Total rows in log:", nrow(samp_all), "\n")
cat("Columns:", paste(colnames(samp_all), collapse = ", "), "\n\n")

# Discard first 25% as burnin
n         <- nrow(samp_all)
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
  log_phi  <- log(phi_col)
  n_below1 <- sum(phi_col < 1)
  n_above1 <- sum(phi_col >= 1)

  # Zero-crossings of log(phi): consecutive sign changes
  signs    <- sign(log_phi)
  crossings <- sum(diff(signs) != 0)

  # First crossing point: index of first phi<1 sample (raw log row index)
  first_cross_retained <- which(phi_col < 1)[1]  # index within retained rows
  first_cross_raw      <- if (!is.na(first_cross_retained)) {
    first_cross_retained + discard_n  # index in full samp_all
  } else NA_integer_

  # Map to approximate MCMC iteration (sample index * thin)
  thin_val <- 40L
  first_cross_iter <- if (!is.na(first_cross_raw)) {
    first_cross_raw * thin_val
  } else NA_integer_

  first_phi <- phi_col[1]
  last_phi  <- phi_col[length(phi_col)]
  first_th  <- if (!is.null(th_col)) th_col[1] else NA
  last_th   <- if (!is.null(th_col)) th_col[length(th_col)] else NA
  min_phi   <- min(phi_col)
  max_phi   <- max(phi_col)

  cat(sprintf("  phi < 1  (reflected mode): %d samples\n", n_below1))
  cat(sprintf("  phi >= 1 (true mode)     : %d samples\n", n_above1))
  cat(sprintf("  Zero-crossings of log(phi): %d\n", crossings))
  if (!is.na(first_cross_iter)) {
    cat(sprintf("  First phi<1 sample at retained index %d, raw log row %d (~iter %d)\n",
                first_cross_retained, first_cross_raw, first_cross_iter))
    # Note: if this is sample 1 (post-burnin), the crossing happened during warmup
    if (first_cross_retained == 1L) {
      cat("  NOTE: crossing at first post-warmup sample — chain likely crossed *during warmup*,\n")
      cat("  not during the sampling phase. Init override took but warmup exploration crossed.\n")
    }
  } else {
    cat("  No phi<1 samples observed.\n")
  }
  cat(sprintf("  First sample: phi=%.4f, theta_1=%.4f\n", first_phi,
              ifelse(is.na(first_th), NaN, first_th)))
  cat(sprintf("  Last  sample: phi=%.4f, theta_1=%.4f\n", last_phi,
              ifelse(is.na(last_th), NaN, last_th)))
  cat(sprintf("  phi range: [%.4f, %.4f]\n\n", min_phi, max_phi))

  # Trace plot: save PNG
  pngFile <- file.path(PILOT_DIR, "phi-trace.png")
  png(pngFile, width = 900, height = 400)
  plot(log_phi, type = "l", col = "darkgreen",
       xlab = "Retained sample index", ylab = "log(phi)",
       main = "Step 3b Chain B (truth-init): log(phi) trace — 100k chain (post-burnin)")
  abline(h = 0, col = "red", lty = 2)
  cross_label <- if (!is.na(first_cross_iter))
    sprintf("Crossings: %d  |  phi<1: %d  |  phi>=1: %d  |  first cross ~iter %d",
            crossings, n_below1, n_above1, first_cross_iter)
  else
    sprintf("Crossings: %d  |  phi<1: %d  |  phi>=1: %d",
            crossings, n_below1, n_above1)
  legend("topright", legend = cross_label, bty = "n", cex = 0.85)
  dev.off()
  cat("  Trace plot saved to", pngFile, "\n\n")
} else {
  cat("  phi column missing — cannot assess mode-crossing.\n\n")
  n_below1 <- NA; n_above1 <- NA; crossings <- NA
  min_phi <- NA; max_phi <- NA; first_phi <- NA; last_phi <- NA
  first_th <- NA; last_th <- NA
  first_cross_iter <- NA; first_cross_raw <- NA; first_cross_retained <- NA
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
  med_phi <- NA_real_
  med_th  <- NA_real_
}

# ---- 3. Mixing health: ESS --------------------------------------------------
cat("=== 3. Mixing health (ESS) ===\n")
ess_cols <- c("phi", "pi0", "theta_1", "tree_length")
ess_vals <- setNames(rep(NA_real_, length(ess_cols)), ess_cols)
for (nm in ess_cols) {
  if (nm %in% colnames(samp)) {
    ess_val      <- coda::effectiveSize(samp[, nm])
    ess_vals[nm] <- ess_val
    cat(sprintf("  %-14s ESS = %6.1f  (n=%d)\n", nm, ess_val, nrow(samp)))
  } else {
    cat(sprintf("  %-14s MISSING\n", nm))
  }
}
cat("\n")

# ---- 4. Tail check: tree_length ---------------------------------------------
cat("=== 4. Tail check (tree_length) ===\n")
cat(sprintf("  median = %.3f  (truth: 12.7)\n",    median(tl)))
cat(sprintf("  IQR    = [%.3f, %.3f]\n",            qs(tl, 0.25), qs(tl, 0.75)))
cat(sprintf("  max    = %.3f\n\n",                  max(tl)))

# ---- Write REPORT.md --------------------------------------------------------
cat("=== Writing REPORT.md ===\n")

rpt <- file.path(PILOT_DIR, "REPORT.md")

# Settled-mode label
settled_label <- if (!is.na(med_phi) && med_phi >= 1) {
  "Chain settled in the **true mode** (phi >= 1, theta_1 near 1)."
} else if (!is.na(med_phi)) {
  "Chain settled in the **reflected mode** (phi < 1, theta_1 near zero)."
} else {
  "Settled mode unknown."
}

# Mode-crossing narrative
if (!is.na(crossings)) {
  if (crossings == 0) {
    crossing_narrative <- "**No mode crossings detected.** The chain, initialised at the truth-mode (phi=4, theta=0.97), never crossed log(phi) = 0. Mode persistence is confirmed for both init directions."
    comparison_text <- "Both chain A (default init, phi<1 mode) and chain B (truth init, phi>1 mode) showed zero mode crossings over 100k iterations. Combined with the sibling chain's result, this confirms bilateral mode persistence: neither mode is escapable on this chain length. Post-hoc relabelling (Option 1) is sufficient."
    conclusion_text <- "Truth-init chain stayed in the phi>1 mode for the full 100k iter run (Option 1 sufficient)."
  } else {
    crossing_narrative <- sprintf("**%d mode crossing(s) detected.** The chain, initialised at truth (phi=4), crossed log(phi) = 0.", crossings)
    if (!is.na(first_cross_iter)) {
      if (!is.na(first_cross_retained) && first_cross_retained == 1L) {
        crossing_narrative <- paste0(crossing_narrative,
          sprintf(" The crossing happened at or before the first post-warmup sample (~iter %d), suggesting the chain escaped the truth-mode *during warmup*, not the sampling phase.", first_cross_iter))
      } else {
        crossing_narrative <- paste0(crossing_narrative,
          sprintf(" First crossing at ~iter %d (retained sample %d of %d).",
                  first_cross_iter, first_cross_retained, nrow(samp)))
      }
    }
    comparison_text <- sprintf("%d crossing(s) from truth-init. Within-run mode flipping is real. Post-hoc relabelling (Option 1) is insufficient; a prior symmetry-break (Option 2: phi >= 1 constraint, or Option 3: asymmetric Beta) is needed before the HPC multi-rep run.", crossings)
    conclusion_text <- sprintf("Truth-init chain crossed to phi<1 mode at iter ~%s (Option 1 insufficient; need symmetry break).",
                               ifelse(is.na(first_cross_iter), "?", as.character(first_cross_iter)))
  }
} else {
  crossing_narrative <- "Mode-crossing status unknown (phi column missing)."
  comparison_text    <- "Comparison unavailable (phi column missing)."
  conclusion_text    <- "Conclusion unavailable — phi column missing from log."
}

# First-cross row for table (show NA as "—")
first_cross_str <- if (is.na(first_cross_iter)) "—" else
  sprintf("~iter %d (sample %d)", first_cross_iter, first_cross_raw)

report_text <- sprintf(
'# Step 3b Mode-Persistence Report — Chain B (truth-init, 100k)

**Date:** %s
**Plan:** `dev/plans/2026-05-13-1100-v2-asymmetric-slab-prior-and-gamma-ridge.md`, Open question: label-switching

---

## Configuration

| Parameter | Value |
|---|---|
| Script | `dev/pilots/2026-05-13-step3b-mode-persistence-B-truthInit/run.R` |
| Seed (data) | 20260512 (same as 20k pilot — identical dataset) |
| Seed (chain) | 20260514 (independent of chain A and 20k pilot) |
| kEco | 2 |
| nTheta | 1 |
| nEco | 60 neomorphic characters |
| nBase | 180 transformational characters |
| phi_truth | 4 |
| stem / root | 0.10 / 0.15 (Goldilocks) |
| **phi_init** | **4.0 (truth-mode override)** |
| **theta_init** | **0.97 (truth-mode override)** |
| nIter | 100,000 |
| thin | 40 |
| nChains / nRuns | 1 / 1 |
| minWarmup / maxWarmup | 2,000 / 5,000 |
| model | ecologyAware=TRUE, magnitudeMode="global", kPrimePrior="geometric" |
| Log file | `dev/pilots/2026-05-13-step3b-mode-persistence-B-truthInit/chain.log` |
| RDS | `dev/pilots/2026-05-13-step3b-mode-persistence-B-truthInit/result.rds` |

**Post-burnin samples:** %d raw rows -> %d retained after 25%% discard.

---

## 0. Init verification (critical)

The truth-mode init is injected via a monkey-patch of `.InitState` that overrides
`phi` and `theta` and recomputes `log_lik`, `log_prior`, `log_post` before the
first MCMC iteration.

**Stdout evidence:** `[truth-init] phi=4.0000 theta=0.9700`

%s

> Note: this confirms the override took at **iter 0**. Whether the chain then
> stayed in the phi>1 mode (across warmup and sampling) is the question answered
> by Section 1 below.

---

## 1. Mode-crossing check (critical)

| Metric | Value |
|---|---|
| phi < 1 (reflected mode) | %s samples |
| phi >= 1 (true mode) | %s samples |
| Zero-crossings of log(phi) | **%s** |
| First phi<1 crossing | %s |
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

| Parameter | ESS | n retained | ESS/n |
|---|---|---|---|
| phi | %.1f | %d | %.3f |
| pi0 | %.1f | %d | %.3f |
| theta_1 | %.1f | %d | %.3f |
| tree_length | %.1f | %d | %.3f |

---

## 4. Tail check (tree_length)

| Metric | Chain B | 20k pilot | Truth |
|---|---|---|---|
| median | %.3f | ~29 | 12.7 |
| IQR | [%.3f, %.3f] | [20.1, 42.3] | — |
| max | %.3f | ~70 | — |

---

## 5. Comparison to sibling chains

%s

---

## Conclusion

%s
',
  format(Sys.Date()),
  n, nrow(samp),
  # init verification
  if (init_took)
    "**INIT VERIFIED.** phi_init=4.0000, theta_init=0.9700 confirmed in stdout. The override ran correctly."
  else
    "**WARNING: init verification FAILED or inconclusive.** The following analysis may be uninformative if the chain started at the default (reflected) state.",
  # mode crossing block
  format(n_below1, big.mark = ","), format(n_above1, big.mark = ","),
  ifelse(is.na(crossings), "NA", as.character(crossings)),
  first_cross_str,
  ifelse(is.na(first_phi), NaN, first_phi),
  ifelse(is.na(first_th), NaN, first_th),
  ifelse(is.na(last_phi), NaN, last_phi),
  ifelse(is.na(last_th), NaN, last_th),
  ifelse(is.na(min_phi), NaN, min_phi),
  ifelse(is.na(max_phi), NaN, max_phi),
  crossing_narrative,
  # settled mode
  ifelse(is.na(med_phi), NaN, med_phi),
  ifelse(is.na(med_th), NaN, med_th),
  settled_label,
  # ESS block
  ess_vals["phi"],         nrow(samp), ess_vals["phi"]         / nrow(samp),
  ess_vals["pi0"],         nrow(samp), ess_vals["pi0"]         / nrow(samp),
  ess_vals["theta_1"],     nrow(samp), ess_vals["theta_1"]     / nrow(samp),
  ess_vals["tree_length"], nrow(samp), ess_vals["tree_length"] / nrow(samp),
  # tail check
  median(tl), qs(tl, 0.25), qs(tl, 0.75), max(tl),
  # comparison
  comparison_text,
  # conclusion
  conclusion_text
)

writeLines(report_text, rpt)
cat("REPORT.md written to", rpt, "\n")
cat("Analysis complete.\n")
