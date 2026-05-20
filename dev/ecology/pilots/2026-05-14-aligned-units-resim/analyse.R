# analyse.R --------------------------------------------------------------------
# Post-processing for the aligned-units re-sim.
# Reads result.rds + chain.log, applies RelabelEcology(), computes summary
# statistics, and writes REPORT.md.
#
# Run from the worktree root after run.R completes:
#   Rscript -e "setwd('C:/Users/pjjg18/GitHub/mkp/.claude/worktrees/ecology-aware'); source('dev/pilots/2026-05-14-aligned-units-resim/analyse.R')"

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})

PILOT_DIR <- "dev/pilots/2026-05-14-aligned-units-resim"
rdsFile   <- file.path(PILOT_DIR, "result.rds")
logFile   <- file.path(PILOT_DIR, "chain.log")
reportFile <- file.path(PILOT_DIR, "REPORT.md")

# ---- Load result --------------------------------------------------------------
saved <- readRDS(rdsFile)
if (isTRUE(saved$crashed) || is.null(saved$res)) stop("Chain crashed — no result to analyse.")
res <- saved$res

# ---- Truth values -------------------------------------------------------------
truth_tl    <- saved$truth_tl      # sum(tree$edge.length) in model units
truth_pi0   <- 0.75                 # 180 trans chars all z=0 -> 180/240 = 0.75
truth_phi   <- 4.0
truth_theta <- 1.0                  # 60 enc / (60+0) = 1

cat("Truth tree_length:", truth_tl, "\n")
cat("Truth pi0:", truth_pi0, "\n")
cat("Truth phi:", truth_phi, "\n")
cat("Truth theta_1:", truth_theta, "\n")

# ---- Apply RelabelEcology ----------------------------------------------------
cat("Applying RelabelEcology...\n")
res_rl <- RelabelEcology(res, magnitudeMode = "global")

# ---- Extract post-burnin samples from log file --------------------------------
# Read log and discard first 25% as burn-in
logData <- read.table(logFile, header = TRUE, sep = "\t", comment.char = "#")
nRaw    <- nrow(logData)
burnin  <- floor(0.25 * nRaw)
post    <- logData[(burnin + 1):nRaw, ]
n       <- nrow(post)
cat("Raw samples:", nRaw, "  Post-burnin:", n, "\n")

# Identify columns
hasTheta <- "theta_1" %in% names(post)
hasPhi   <- "phi"     %in% names(post)
hasPi0   <- "pi0"     %in% names(post)
hasTL    <- "tree_length" %in% names(post)

cat("Columns found — phi:", hasPhi, " pi0:", hasPi0,
    " theta_1:", hasTheta, " tree_length:", hasTL, "\n")

# ---- After relabelling, re-read log via RelabelEcology samples if available --
# RelabelEcology returns samples in $samples slot (streaming reads from logFile)
# Use whichever is available.
if (!is.null(res_rl$samples)) {
  samp_rl <- res_rl$samples
  # discard burnin
  n_rl <- nrow(samp_rl)
  bi_rl <- floor(0.25 * n_rl)
  samp  <- samp_rl[(bi_rl + 1):n_rl, ]
  cat("Using RelabelEcology samples (n=", nrow(samp), ")\n")
} else {
  # Fall back: apply relabelling logic directly to log columns
  samp <- post
  if (hasPhi && hasPi0 && hasTheta) {
    flip_idx <- samp$phi < 1
    if (any(flip_idx)) {
      samp$phi[flip_idx]     <- 1 / samp$phi[flip_idx]
      samp$theta_1[flip_idx] <- 1 - samp$theta_1[flip_idx]
      cat("Manual relabelling applied to", sum(flip_idx), "samples with phi<1\n")
    }
  }
}

# ---- Posterior summaries ------------------------------------------------------
iqr_str <- function(x) {
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.4f  IQR [%.4f, %.4f]  max %.4f", q[2], q[1], q[3], max(x, na.rm=TRUE))
}

phi_col    <- if ("phi"        %in% names(samp)) samp$phi        else NULL
pi0_col    <- if ("pi0"        %in% names(samp)) samp$pi0        else NULL
theta_col  <- if ("theta_1"    %in% names(samp)) samp$theta_1    else NULL
tl_col     <- if ("tree_length"%in% names(samp)) samp$tree_length else NULL

# gamma_e per sample: pi0 + (1-pi0)*(theta*phi + (1-theta)/phi)
if (!is.null(phi_col) && !is.null(pi0_col) && !is.null(theta_col)) {
  gE_col <- pi0_col + (1 - pi0_col) * (theta_col * phi_col + (1 - theta_col) / phi_col)
} else {
  gE_col <- NULL
}

# tl / gamma_e
if (!is.null(tl_col) && !is.null(gE_col)) {
  tl_gE_col <- tl_col / gE_col
} else {
  tl_gE_col <- NULL
}

cat("\n--- Posterior summaries (post-burnin, relabelled) ---\n")
if (!is.null(phi_col))   cat("phi:        ", iqr_str(phi_col),   "\n")
if (!is.null(pi0_col))   cat("pi0:        ", iqr_str(pi0_col),   "\n")
if (!is.null(theta_col)) cat("theta_1:    ", iqr_str(theta_col), "\n")
if (!is.null(tl_col))    cat("tree_length:", iqr_str(tl_col),    "\n")
if (!is.null(gE_col))    cat("gamma_e:    ", iqr_str(gE_col),    "\n")
if (!is.null(tl_gE_col)) cat("tl/gamma_e: ", iqr_str(tl_gE_col),"\n")

# Mode flag
if (!is.null(phi_col)) {
  n_phi_lt1 <- sum(phi_col < 1, na.rm = TRUE)
  mode_flag <- if (n_phi_lt1 < 0.05 * length(phi_col)) {
    sprintf("phi >= 1 (truth mode) — %d/%d samples < 1", n_phi_lt1, length(phi_col))
  } else if (n_phi_lt1 > 0.95 * length(phi_col)) {
    sprintf("phi < 1 (reflected mode) — %d/%d samples < 1", n_phi_lt1, length(phi_col))
  } else {
    sprintf("MIXED modes — %d/%d samples with phi < 1 (%.0f%%)",
            n_phi_lt1, length(phi_col), 100 * n_phi_lt1 / length(phi_col))
  }
  cat("Mode:", mode_flag, "\n")
}

# ESS
ess_scalar <- function(x, nm) {
  x <- x[!is.na(x)]
  if (length(x) < 5) return(invisible(NULL))
  # Simple batch-means ESS
  b  <- max(1L, floor(sqrt(length(x))))
  nb <- floor(length(x) / b)
  bm <- sapply(seq_len(nb), function(i) mean(x[((i-1)*b+1):(i*b)]))
  s2 <- var(bm) * b
  ess <- if (s2 > 0) length(x) * var(x) / (length(x) * s2) else Inf
  cat(sprintf("  ESS %-14s %.1f\n", nm, min(ess, length(x))))
  ess
}
cat("--- ESS ---\n")
ess_phi   <- if (!is.null(phi_col))   ess_scalar(phi_col,   "phi")   else NA
ess_pi0   <- if (!is.null(pi0_col))   ess_scalar(pi0_col,   "pi0")   else NA
ess_theta <- if (!is.null(theta_col)) ess_scalar(theta_col, "theta_1") else NA
ess_tl    <- if (!is.null(tl_col))    ess_scalar(tl_col,    "tree_length") else NA

# ---- Verdict ------------------------------------------------------------------
tl_med  <- if (!is.null(tl_col))  median(tl_col,  na.rm=TRUE) else NA
pi0_med <- if (!is.null(pi0_col)) median(pi0_col, na.rm=TRUE) else NA
phi_med <- if (!is.null(phi_col)) median(phi_col, na.rm=TRUE) else NA

tl_ok  <- !is.na(tl_med)  && abs(tl_med  - truth_tl)   / truth_tl   <= 0.20
pi0_ok <- !is.na(pi0_med) && abs(pi0_med - truth_pi0)  / truth_pi0  <= 0.20

verdict_tl  <- if (tl_ok)  "UNIT MISMATCH CONFIRMED"  else "UNIT MISMATCH NOT CONFIRMED"
verdict_pi0 <- if (pi0_ok) "pi0 bias RESOLVED"        else "pi0 bias PERSISTS"
headline    <- paste(verdict_tl, "+", verdict_pi0)
cat("\n=== VERDICT ===\n", headline, "\n")

# ---- Trace plots --------------------------------------------------------------
make_trace <- function(x, nm, truth, outfile) {
  tryCatch({
    png(outfile, width=900, height=350)
    plot(x, type="l", col="steelblue", xlab="Iteration (thinned)", ylab=nm,
         main=paste0(nm, " trace  (post-burnin, relabelled)"))
    abline(h=truth, col="red", lty=2, lwd=2)
    legend("topright", legend=paste0("truth=", truth), col="red", lty=2)
    dev.off()
    cat("Saved", outfile, "\n")
  }, error=function(e) cat("Trace plot failed:", conditionMessage(e), "\n"))
}

if (!is.null(phi_col))
  make_trace(phi_col,   "phi",         truth_phi,   file.path(PILOT_DIR, "phi-trace.png"))
if (!is.null(pi0_col))
  make_trace(pi0_col,   "pi0",         truth_pi0,   file.path(PILOT_DIR, "pi0-trace.png"))
if (!is.null(tl_col))
  make_trace(tl_col,    "tree_length", truth_tl,    file.path(PILOT_DIR, "tl-trace.png"))

# ---- Write REPORT.md ----------------------------------------------------------
fmt_iqr <- function(x, digits=4) {
  if (is.null(x)) return("N/A")
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm=TRUE)
  sprintf("%.*f  IQR [%.*f, %.*f]  max %.*f",
          digits, q[2], digits, q[1], digits, q[3], digits, max(x, na.rm=TRUE))
}
fmt_ess <- function(e) if (is.na(e)) "N/A" else sprintf("%.1f", e)
pct_dev <- function(est, truth) sprintf("%.1f%%", 100*(est-truth)/truth)

report <- c(
  "# REPORT: Aligned-units re-sim",
  "",
  paste0("**Date:** 2026-05-14"),
  paste0("**Run:** `dev/pilots/2026-05-14-aligned-units-resim/`"),
  paste0("**Seed (data):** 20260514 (generates new data with aligned sim params)"),
  paste0("**Seed (MCMC):** default init"),
  paste0("**nIter:** 100 000  thin=40  post-burnin n=", n),
  paste0("**Elapsed:** ", format(saved$elapsed)),
  "",
  "## Simulator alignment",
  "",
  "Previous chain A used `baseRate=0.5, normalize=FALSE` (defaults). The model's",
  "JC kernel has no separate `baseRate` — rate is fully absorbed into edge length.",
  "This makes simulator time units ≠ model substitution units, inflating posterior",
  "tree_length by ~2x.",
  "",
  "This run uses:",
  "- `baseRate = 1.0` — matches model parameterisation",
  "- `normalize = TRUE` — divides per-cell rate factor by gammaE so expected",
  "  substitutions equal the model's expectation at truth",
  "- `pi0 = 0.75, theta = 1.0` **explicit** — the simulator's default `pi0`",
  "  computes `mean(z==0)` over ALL columns including the reference column,",
  "  giving 420/480 = 0.875 rather than the model's 180/240 = 0.75. Explicit",
  "  values ensure simulator gammaE matches model gammaE at truth.",
  "",
  "## Truth values (aligned sim)",
  "",
  paste0("| Parameter | Truth | Notes |"),
  paste0("|-----------|-------|-------|"),
  paste0("| tree_length | ", round(truth_tl, 4),
         " | sum(tree$edge.length) in model substitution units |"),
  paste0("| pi0 | ", truth_pi0,
         " | 180 trans chars all z=0; 180/240 = 0.75 |"),
  paste0("| phi | ", truth_phi, " | simulation parameter |"),
  paste0("| theta_1 | ", truth_theta,
         " | 60 enc / (60 enc + 0 disc) = 1.0 |"),
  "",
  "## Posterior summaries (post-25%-burnin, RelabelEcology applied)",
  "",
  paste0("n = ", n, " samples"),
  "",
  paste0("| Parameter | Posterior | Truth | % dev |"),
  paste0("|-----------|-----------|-------|-------|"),
  paste0("| phi | ", fmt_iqr(phi_col, 3), " | ", truth_phi,
         " | ", if (!is.null(phi_col)) pct_dev(median(phi_col,na.rm=TRUE), truth_phi) else "N/A", " |"),
  paste0("| pi0 | ", fmt_iqr(pi0_col, 3), " | ", truth_pi0,
         " | ", if (!is.null(pi0_col)) pct_dev(pi0_med, truth_pi0) else "N/A", " |"),
  paste0("| theta_1 | ", fmt_iqr(theta_col, 3), " | ", truth_theta,
         " | ", if (!is.null(theta_col)) pct_dev(median(theta_col,na.rm=TRUE), truth_theta) else "N/A", " |"),
  paste0("| tree_length | ", fmt_iqr(tl_col, 3), " | ", round(truth_tl, 3),
         " | ", if (!is.null(tl_col)) pct_dev(tl_med, truth_tl) else "N/A", " |"),
  paste0("| gamma_e | ", fmt_iqr(gE_col, 3), " | ~1.75 (truth-pi0-aware) | — |"),
  paste0("| tl / gamma_e | ", fmt_iqr(tl_gE_col, 3), " | ~", round(truth_tl / 1.75, 2), " | — |"),
  "",
  paste0("**Mode:** ", if (!is.null(phi_col)) mode_flag else "N/A"),
  "",
  "### ESS",
  "",
  paste0("| Parameter | ESS |"),
  paste0("|-----------|-----|"),
  paste0("| phi | ", fmt_ess(ess_phi), " |"),
  paste0("| pi0 | ", fmt_ess(ess_pi0), " |"),
  paste0("| theta_1 | ", fmt_ess(ess_theta), " |"),
  paste0("| tree_length | ", fmt_ess(ess_tl), " |"),
  "",
  "## Verdict",
  "",
  paste0("**Headline: ", headline, "**"),
  "",
  paste0("- tree_length: posterior median = ", round(tl_med, 3),
         ", truth = ", round(truth_tl, 3),
         " (", pct_dev(tl_med, truth_tl), " deviation).  ",
         "Threshold ±20%.  → **", verdict_tl, "**"),
  paste0("- pi0: posterior median = ", round(pi0_med, 3),
         ", truth = ", truth_pi0,
         " (", pct_dev(pi0_med, truth_pi0), " deviation).  ",
         "Threshold ±20%.  → **", verdict_pi0, "**"),
  "",
  "## Implications for paper",
  ""
)

# Conditional implications
if (tl_ok && pi0_ok) {
  report <- c(report,
    "Unit mismatch confirmed AND pi0 resolved. Both biases were artefacts of the",
    "misaligned simulator parameters.",
    "",
    "**Action:** Dispatch Step 5 HPC multi-rep Sim 3. Update the multi-rep sim",
    "driver (`inst/ecology/simulations/sim3-simulate.R` call site) to use",
    "`baseRate=1.0, normalize=TRUE, pi0=<empirical>, theta=<empirical>` so all",
    "20 replicates are generated in aligned units.",
    "",
    "Document in vignette: `tree_length` is in model substitution-per-site units",
    "(not simulator time units). The aligned-sim convention is used throughout."
  )
} else if (tl_ok && !pi0_ok) {
  report <- c(report,
    "Unit mismatch confirmed (tree_length bias was a simulator/model unit issue).",
    "However, pi0 bias persists even with aligned simulator parameters.",
    "This is a real inference problem — the chain cannot recover the correct",
    "proportion of z=none characters even when the data are simulated in model units.",
    "",
    "**Action:** Do NOT dispatch HPC until pi0 bias is investigated. Possible",
    "causes: prior on pi0 is too weak, likelihood is weakly informative about",
    "pi0 when most chars are trans (z=0 looks like z=none under the symmetric",
    "slab), or gamma_e normalisation creates a pi0↔tl trade-off.",
    "",
    "Recommended next step: vary the rho0Alpha/rho0Beta prior hyperparameters",
    "and recheck. Alternatively, run a single chain on a dataset with larger",
    "nEco fraction to see if the pi0 bias shrinks."
  )
} else {
  report <- c(report,
    "Unit mismatch NOT confirmed. Even with aligned simulator parameters, posterior",
    "tree_length does not recover truth. This indicates a deeper problem:",
    "possibly gamma_e normalisation in the simulator does not match the model's",
    "gamma_e computation, or the model has an additional scaling we missed.",
    "",
    "**Action:** Do NOT dispatch HPC. Return to model diagnostics. Re-read",
    "`src/mcmc_ecology.cpp` gamma_e_compute and compare line-by-line with the",
    "simulator's normalize=TRUE path in `inst/ecology/simulations/sim3-simulate.R`."
  )
}

report <- c(report,
  "",
  "## Key files",
  "",
  "- `dev/pilots/2026-05-14-aligned-units-resim/run.R` — this run",
  "- `dev/pilots/2026-05-13-step3b-mode-persistence-A/run.R` — previous chain A (baseline)",
  "- `dev/pilots/2026-05-13-tl-bias-investigation/REPORT.md` — unit-mismatch diagnosis",
  "- `inst/ecology/simulations/sim3-simulate.R` lines 150–185, 210–249 — simulator",
  "- `src/mcmc_ecology.cpp` line 249 — model gamma_e_compute"
)

writeLines(report, reportFile)
cat("\nWrote", reportFile, "\n")
cat("=== DONE ===\n")
