# T-OVL: Posterior overlap, sampled-k vs marginal-k geometric arm
# ===============================================================
#
# Plan reference: dev/notes/2026-05-28-marginal-k-plan.md §7.2 (Rung 2);
# redesigned 2026-06-02 (FREEZE-003 follow-up) per advisor guidance.
#
# WHAT THIS CONFIRMS, AND WHY IT MATTERS
# --------------------------------------
# The deterministic RB identity (test-marginal-k-rb-free-topology.R) already
# proves marginal_k and sampled_k share the SAME target posterior to 1e-7 on
# arbitrary trees. That is a *likelihood-target* check. It does NOT exercise the
# move machinery: a topology/branch move can commit the correct logLik of
# wherever it landed (coherence-sweep gap = 0) yet still PROPOSE from a biased
# distribution. The committed-logLik sweep cannot see proposal bias. THIS harness
# can: if any enabled move samples the wrong target, the marginal_k posterior on
# (tree_length, rate_log_sd, p) drifts away from the sampled_k posterior beyond
# Monte-Carlo error. So this is the *acceptance gate* for re-enabling any move
# under marginal_k (FREEZE-003 "efficiency tomorrow"), not redundant confirmation.
#
# DESIGN (three fixes over the legacy version that false-FAILed)
# --------------------------------------------------------------
#  1. SIGNAL-BEARING data with GUARANTEED k' spread. The legacy harness used
#     sample(0:1) — signal-free, all kObs=2, so the latent-state marginalisation
#     was near-degenerate (k' pinned at kObs) and p was prior-dominated: the test
#     was blind to the very machinery it claimed to confirm. Here we forward-
#     simulate under JC with kTrue ~ 2 + Geo(p_true) (the inference rate
#     convention, SBC-HARNESS-003) so a real fraction of characters have
#     kObs < kTrue: the k' posterior has genuine spread and p is identifiable.
#  2. EQUIVALENCE-within-Monte-Carlo-error, not a difference test. KS assumes iid;
#     MCMC draws are autocorrelated, so KS treats effective-n as inflated and
#     over-rejects (the legacy false-FAIL). With enough draws KS rejects ANY two
#     finite chains. We instead test |mean_S - mean_M| (and |sd_S - sd_M|) against
#     an autocorrelation-corrected combined standard error, computed by BATCH
#     MEANS (self-contained; no coda dependency).
#  3. ESS FLOOR, pre-registered. An equivalence test with tiny ESS has huge bands
#     => vacuous pass. Cells below the floor report INCONCLUSIVE, never PASS.
#
# Matched seeds (common random numbers) start both modes from an identical RNG
# stream so the comparison isolates systematic move bias from seed noise.
#
# PARAMETERISED ON MOVE-CONFIG (Phase 2 reuse): set env MARGINAL_K_OVL_EXTRA to a
# comma list of MkPrimeMCMC flags to force ON under BOTH modes, e.g.
#   MARGINAL_K_OVL_EXTRA=weightedSpr Rscript T-OVL-sampled-vs-marginal.R
# With the var unset, both modes run the PRODUCTION gated schedule (what ships).
#
# DO NOT RUN THE FULL GRID INLINE — heavy. Submit to Hamilton via the companion
# submit-marginal-k-ovl.sh. For a local end-to-end check use MARGINAL_K_OVL_QUICK=1.

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library("ape")
  library("TreeTools")
})

# ---------------------------------------------------------------------------
# Knobs
# ---------------------------------------------------------------------------

# Knobs are env-overridable (so a runner/Hamilton task can resize without editing):
#   MARGINAL_K_OVL_NTIP, _NCHAR  : comma lists of ints   (e.g. "8,16")
#   MARGINAL_K_OVL_REP           : reps per cell
#   MARGINAL_K_OVL_ITER, _WARM   : chain length / warmup
#   MARGINAL_K_OVL_CELL          : 1-based index; run ONLY that cell (Hamilton array)
.env_ints <- function(key, default) {
  v <- Sys.getenv(key, unset = "")
  if (!nzchar(v)) return(default)
  as.integer(trimws(strsplit(v, ",", fixed = TRUE)[[1]]))
}
.env_int <- function(key, default) {
  v <- Sys.getenv(key, unset = "")
  if (!nzchar(v)) return(default)
  as.integer(v)
}
GRID_NTIP  <- .env_ints("MARGINAL_K_OVL_NTIP",  c(8L, 16L))
GRID_NCHAR <- .env_ints("MARGINAL_K_OVL_NCHAR", c(24L, 48L))
N_REP      <- .env_int("MARGINAL_K_OVL_REP",  4L)   # replicate datasets per cell
N_ITER     <- .env_int("MARGINAL_K_OVL_ITER", 60000L) # post-warmup; long => ESS over floor
N_WARM     <- .env_int("MARGINAL_K_OVL_WARM", 8000L)
P_TRUE     <- 0.25        # k' ~ 2 + Geo(0.25): mean ~5 states. Low p keeps a real
                          # fraction of characters with kObs<kTrue even at 16 tips
                          # (more tips observe more states), so k' spread + p
                          # identifiability survive => the marginalisation is exercised.
Z_BAR      <- 3.0         # equivalence band: |Delta| < Z_BAR * combined SE
ESS_FLOOR  <- 200L        # per param per mode; below => INCONCLUSIVE
PARAMS     <- c("tree_length", "rate_log_sd", "p")

# Extra MkPrimeMCMC flags to force ON in BOTH modes (Phase 2 move re-enable gate).
.parse_extra <- function() {
  raw <- Sys.getenv("MARGINAL_K_OVL_EXTRA", unset = "")
  if (!nzchar(raw)) return(list())
  flags <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  flags <- flags[nzchar(flags)]
  stats::setNames(as.list(rep(TRUE, length(flags))), flags)
}
EXTRA_MCMC <- .parse_extra()

# Quick smoke (env MARGINAL_K_OVL_QUICK=1): 1 small cell, short chains — verifies
# the driver runs end-to-end and the diagnostics are sane on the current binary.
# NOT a calibration check (ESS too small => INCONCLUSIVE by design).
QUICK <- identical(Sys.getenv("MARGINAL_K_OVL_QUICK", unset = ""), "1")
if (QUICK) {
  GRID_NTIP <- 8L; GRID_NCHAR <- 24L; N_REP <- 1L
  N_ITER <- 1500L; N_WARM <- 500L
  message("[T-OVL] QUICK smoke: 1 cell (8x24), N_ITER=1500 — end-to-end + diagnostics only")
}

OUT_DIR <- "dev/red-team/heavy-tests/marginal-k"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Forward simulation (signal-bearing; matches inference rate convention)
# Reuses the SBC-validated simulator (SBC-TL-MIX-001 / SBC-HARNESS-003).
# ---------------------------------------------------------------------------

EXPSTEPS_FIXED <- 50         # treeLengthRate = TREE_SHAPE / EXPSTEPS_FIXED
TREE_SHAPE     <- 2          # matches MkPrimeModel() default treeLengthShape

.simTree <- function(nTip) {
  tr <- ape::rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  nEdge <- nrow(tr$edge)
  tl <- stats::rgamma(1, shape = TREE_SHAPE, rate = TREE_SHAPE / EXPSTEPS_FIXED)
  w  <- stats::rexp(nEdge, rate = 1)               # Dirichlet(1,..,1) edge fractions
  tr$edge.length <- tl * w / sum(w)
  TreeTools::Preorder(tr)
}

# One character under JC(kTrue) by root->tip draws on a Preorder edge list.
# Rate arg matches inference (src/likelihood.cpp:85): -kStates*t/(kStates-1).
.simJCchar <- function(tree, kTrue) {
  nTip <- length(tree$tip.label)
  states <- integer(2L * nTip - 1L)
  rootIdx <- nTip + 1L
  states[rootIdx] <- sample.int(kTrue, 1L) - 1L
  edges <- tree$edge; el <- tree$edge.length
  for (e in seq_len(nrow(edges))) {
    pa <- edges[e, 1L]; ch <- edges[e, 2L]; t <- el[e]
    pSame <- 1 / kTrue + (1 - 1 / kTrue) * exp(-kTrue * t / (kTrue - 1))
    if (runif(1L) < pSame) {
      states[ch] <- states[pa]
    } else {
      states[ch] <- sample(setdiff(seq.int(0L, kTrue - 1L), states[pa]), 1L)
    }
  }
  states[seq_len(nTip)]
}

.canonicalise <- function(vec) {
  uvals <- sort(unique(vec))
  match(vec, uvals) - 1L
}

# Build a signal-bearing variable-character dataset with genuine latent spread.
# Returns the dataset, the true tree, and a k'-spread diagnostic.
simulate_dataset <- function(nTip, nChar, seed) {
  set.seed(seed)
  tree <- .simTree(nTip)
  cols <- vector("list", nChar)
  kTrueVec <- integer(nChar)
  kObsVec  <- integer(nChar)
  i <- 0L
  guard <- 0L
  while (i < nChar) {
    guard <- guard + 1L
    if (guard > 50L * nChar) stop("simulate_dataset: cannot fill variable chars")
    kTrue <- 2L + stats::rgeom(1L, P_TRUE)            # k' ~ 2 + Geo(p_true)
    raw   <- .simJCchar(tree, kTrue)
    can   <- .canonicalise(raw)
    kObs  <- length(unique(can))
    if (kObs < 2L) next                                # reject invariant chars
    i <- i + 1L
    cols[[i]]   <- can
    kTrueVec[i] <- kTrue
    kObsVec[i]  <- kObs
  }
  mat <- matrix(unlist(cols), nrow = nTip, ncol = nChar,
                dimnames = list(tree$tip.label, NULL))
  list(tree = tree,
       mkd  = MkPrimeData(MatrixToPhyDat(mat)),
       kTrue = kTrueVec, kObs = kObsVec,
       frac_latent = mean(kObsVec < kTrueVec))         # fraction with unobserved states
}

# ---------------------------------------------------------------------------
# Batch-means MCSE / ESS (autocorrelation-corrected, self-contained)
# MCSE(mean) = sd(batch means)/sqrt(a); ESS = var(x) / MCSE^2.
# ---------------------------------------------------------------------------

bm_stats <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 16L) return(list(mean = mean(x), sd = stats::sd(x),
                           mcse = NA_real_, ess = NA_real_, n = n))
  a <- floor(sqrt(n))                                  # number of batches
  b <- floor(n / a)                                    # batch size
  x <- x[(n - a * b + 1L):n]                           # keep last a*b (drop warmup remainder)
  bm <- colMeans(matrix(x, nrow = b, ncol = a))        # column k = batch k
  mcse <- stats::sd(bm) / sqrt(a)
  v    <- stats::var(x)
  ess  <- if (mcse > 0) v / (mcse^2) else NA_real_
  list(mean = mean(x), sd = sqrt(v), mcse = mcse, ess = ess, n = length(x))
}

# Equivalence verdict for one parameter, comparing two chains.
compare_param <- function(xs, xm) {
  s <- bm_stats(xs); m <- bm_stats(xm)
  ess_ok <- is.finite(s$ess) && is.finite(m$ess) &&
            s$ess >= ESS_FLOOR && m$ess >= ESS_FLOOR
  # standardised mean difference
  se_mean <- sqrt(s$mcse^2 + m$mcse^2)
  d_mean  <- if (is.finite(se_mean) && se_mean > 0)
               abs(s$mean - m$mean) / se_mean else NA_real_
  # standardised sd difference: MCSE(sd) ~ sd / sqrt(2*ESS)
  se_sd <- sqrt(s$sd^2 / (2 * s$ess) + m$sd^2 / (2 * m$ess))
  d_sd  <- if (is.finite(se_sd) && se_sd > 0)
             abs(s$sd - m$sd) / se_sd else NA_real_
  status <- if (!ess_ok) "INCONCLUSIVE"
            else if (is.finite(d_mean) && is.finite(d_sd) &&
                     d_mean < Z_BAR && d_sd < Z_BAR) "PASS"
            else "FAIL"
  list(status = status,
       mean_s = s$mean, mean_m = m$mean, d_mean = d_mean,
       sd_s = s$sd, sd_m = m$sd, d_sd = d_sd,
       ess_s = s$ess, ess_m = m$ess)
}

# ---------------------------------------------------------------------------
# One-cell run: both modes on the same dataset, matched seed, same schedule.
# ---------------------------------------------------------------------------

run_cell <- function(nTip, nChar, rep) {
  seed <- as.integer(1e3 * rep + 10 * nTip + nChar)
  sim  <- simulate_dataset(nTip, nChar, seed)
  cell_name <- sprintf("n%02d_c%02d_r%02d", nTip, nChar, rep)
  message(sprintf("[T-OVL] %s : frac_latent=%.2f (kObs<kTrue) ...",
                  cell_name, sim$frac_latent))

  run_one <- function(mode) {
    model <- suppressMessages(MkPrimeModel(kPrimePrior    = "geometric",
                                           likelihoodMode = mode,
                                           coding         = "variable"))
    # Pin thinning to store ~2000 draws so the batch-means ESS estimate is not
    # artificially capped below ESS_FLOOR by too few stored samples.
    thin_n <- max(1L, round(N_ITER / 2000L))
    mcmc_args <- c(list(nIter = N_ITER, thin = thin_n,
                        minWarmup = N_WARM, maxWarmup = N_WARM,
                        autoTune = FALSE, nRuns = 1L, nChains = 1L),
                   EXTRA_MCMC)
    mcmc <- do.call(MkPrimeMCMC, mcmc_args)
    set.seed(seed)   # common random numbers across modes
    suppressMessages(suppressWarnings(
      RunMkPrime(sim$mkd, sim$tree, model = model, mcmc = mcmc, overwrite = TRUE)
    ))
  }

  fit_s <- run_one("sampled_k")
  fit_m <- run_one("marginal_k")

  cmp <- list()
  for (pn in PARAMS) {
    if (pn %in% colnames(fit_s$samples) && pn %in% colnames(fit_m$samples)) {
      cmp[[pn]] <- compare_param(fit_s$samples[, pn], fit_m$samples[, pn])
    }
  }
  list(cell = cell_name, seed = seed, frac_latent = sim$frac_latent,
       extra = names(EXTRA_MCMC), cmp = cmp)
}

# ---------------------------------------------------------------------------
# Grid sweep
# ---------------------------------------------------------------------------

if (sys.nframe() == 0L) {
  # Enumerate the cell grid; MARGINAL_K_OVL_CELL selects one (Hamilton array).
  grid <- expand.grid(rep = seq_len(N_REP), nChar = GRID_NCHAR, nTip = GRID_NTIP,
                       KEEP.OUT.ATTRS = FALSE)
  grid <- grid[order(grid$nTip, grid$nChar, grid$rep), , drop = FALSE]
  cell_sel <- .env_int("MARGINAL_K_OVL_CELL", NA_integer_)
  tag <- if (length(EXTRA_MCMC)) paste(names(EXTRA_MCMC), collapse = "+") else "gated"
  if (!is.na(cell_sel)) {
    rows <- cell_sel
    out_rds <- file.path(OUT_DIR, sprintf("T-OVL-%s-cell%02d.rds", tag, cell_sel))
  } else {
    rows <- seq_len(nrow(grid))
    out_rds <- file.path(OUT_DIR, sprintf("T-OVL-%s-results.rds", tag))
  }
  results <- list()
  for (gi in rows) {
    cell <- run_cell(grid$nTip[gi], grid$nChar[gi], grid$rep[gi])
    results[[cell$cell]] <- cell
    saveRDS(results, out_rds)
  }

  # Verdict: PASS iff every conclusive cell/param PASSes and at least one is
  # conclusive. INCONCLUSIVE cells do not count as PASS or FAIL.
  any_fail <- FALSE; any_pass <- FALSE; n_incon <- 0L
  cat(sprintf("\n[T-OVL] schedule: %s | Z_BAR=%.1f ESS_FLOOR=%d\n",
              if (length(EXTRA_MCMC)) paste(names(EXTRA_MCMC), collapse = "+")
              else "PRODUCTION (gated)", Z_BAR, ESS_FLOOR))
  for (cn in names(results)) {
    fl <- results[[cn]]$frac_latent
    for (pn in names(results[[cn]]$cmp)) {
      r <- results[[cn]]$cmp[[pn]]
      if (identical(r$status, "FAIL")) any_fail <- TRUE
      if (identical(r$status, "PASS")) any_pass <- TRUE
      if (identical(r$status, "INCONCLUSIVE")) n_incon <- n_incon + 1L
      cat(sprintf("[T-OVL] %s : %-12s : %-12s d_mean=%5.2f d_sd=%5.2f ESS=(%.0f,%.0f)\n",
                  cn, pn, r$status,
                  ifelse(is.finite(r$d_mean), r$d_mean, NA),
                  ifelse(is.finite(r$d_sd), r$d_sd, NA),
                  ifelse(is.finite(r$ess_s), r$ess_s, NA),
                  ifelse(is.finite(r$ess_m), r$ess_m, NA)))
    }
  }
  verdict <- if (any_fail) "FAIL" else if (any_pass) "PASS" else "INCONCLUSIVE"
  cat(sprintf("[T-OVL] VERDICT: %s  (inconclusive params: %d)\n", verdict, n_incon))
  # Per-cell array tasks write a per-cell verdict; the gather step aggregates.
  vfile <- if (!is.na(cell_sel))
             sprintf("T-OVL-%s-cell%02d-verdict.txt", tag, cell_sel)
           else sprintf("T-OVL-%s-verdict.txt", tag)
  writeLines(verdict, file.path(OUT_DIR, vfile))
}
