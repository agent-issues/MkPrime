# dev/red-team/heavy-tests/subtree-swap-db.R
#
# Sub-harness D2(b): β=0 detailed-balance test for the two subtree-swap
#                    moves flagged OPEN by L3 (math-prover, Wave 1).
#
# CLAIMS CHALLENGED
#   SWAP-001 (HIGH, OPEN): `GibbsSubtreeSwap` (src/mcmc.cpp:1654-1878)
#     commits the chosen partner without an MH accept step. Detailed
#     balance requires Z(A;T) = Z(B;T'), which depends on partner-set
#     counts |partners_T(A)| = |partners_T'(B)|. These are NOT generally
#     equal across a swap.
#
#   SWAP-002 (HIGH, OPEN): `WeightedSubtreeSwap` (src/mcmc.cpp:2825-3009)
#     has an MH step, but its logHR (lines 2980-2983) omits the
#     log(|partners_T'(B)| / |partners_T(A)|) correction.
#
# THEORETICAL BASIS
#   With flat prior over labelled binary unrooted topologies (15 for n=5)
#   and a flat target (β=0 ⇒ all candidate logLs equal ⇒ weights uniform),
#   a correct MH chain MUST sample uniformly over the 15 topologies.
#   The 15-bucket chi-squared goodness-of-fit is the same test used by
#   `tests/testthat/test-tbr-detailed-balance.R`, which PASSES for NNI,
#   SPR, TBR.
#
# REDUCTION AT β=0
#   For Gibbs: exp(β·logLik) = 1 for every candidate, so the chosen
#     partner is drawn UNIFORMLY from {self, p_1, ..., p_K} where K =
#     |partners(A)|. We reproduce this exact draw in R.
#   For Weighted: at β=0 every per-bin weight is also 1; mCand[pi] = nBins
#     for every pi, mOrig = nBins. Topology selection is uniform.
#     The default-bin/new-bin Beta densities cancel except via the
#     Beta-density evaluations at f_default vs f_new — which are also
#     symmetric in expectation. The logHR submitted by the code is
#     therefore essentially zero (or close to it) and does NOT include
#     the partner-set size correction. We reproduce that exact logHR.
#
# IMPLEMENTATION
#   Pure R driver using two exported C++ helpers:
#     - get_valid_swap_partners_cpp(edge, nTip, pruneNode)
#     - swap_subtrees_cpp(edge, nTip, treeLength, relBrLengths, nodeA, nodeB)
#   These are the SAME functions called by gibbs_subtree_swap_impl and
#   weighted_subtree_swap_impl. Because the driver is pure R, the
#   "only this move is applied" check is trivial: no other proposals
#   are constructed.
#
# CONTROLS
#   Positive control: pure NNI on the same setup — must pass, replicates
#                     what test-tbr-detailed-balance.R does for SPR/TBR.
#
# PASS CRITERION
#   chi-squared p > 0.005 (matches the existing test threshold).
#   For 4 buckets (one positive control + two swap arms + one extra),
#   Bonferroni-adjusted alpha = 0.005 / 4 ≈ 0.00125 if we want strict
#   FWER. We use 0.005 here because the existing TBR test does so;
#   any p < 0.005 from the FULL run is reported as a likely real bug.
#
# Usage:
#   Rscript dev/red-team/heavy-tests/subtree-swap-db.R --quick
#   Rscript dev/red-team/heavy-tests/subtree-swap-db.R
#
# Output:
#   dev/red-team/heavy-tests/subtree-swap-db-results/summary.rds
#   dev/red-team/heavy-tests/subtree-swap-db-results/verdict.txt

suppressPackageStartupMessages({
  library(ape)
  library(TreeTools)
})

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
.nFromArgs <- {
  ix <- which(args == "--n")
  if (length(ix) && length(args) >= ix + 1L) as.integer(args[ix + 1L]) else 5L
}
N_TIP_CFG <- .nFromArgs
# Unrooted binary topologies on n tips = (2n - 5)!! = double-factorial of (2n-5)
.dfact_odd <- function(m) if (m <= 1L) 1L else prod(seq.int(m, 1L, by = -2L))
N_TOPOS_EXPECTED <- .dfact_odd(2L * N_TIP_CFG - 5L)
cat(sprintf("[D2b] n_tip=%d  expected_buckets=%d  quick=%s\n",
            N_TIP_CFG, N_TOPOS_EXPECTED, quick))

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(file.path(getwd()), quiet = TRUE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(file.path(getwd()), quiet = TRUE)
} else {
  stop("Need {pkgload} or {devtools}")
}

out_dir <- "dev/red-team/heavy-tests/subtree-swap-db-results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Topology key (bipartition-based, root-invariant; matches the helper used
# by tests/testthat/test-tbr-detailed-balance.R::.topo_key).
# ---------------------------------------------------------------------------

.topo_key <- function(tr) {
  splits <- TreeTools::as.Splits(tr)
  mat <- as.logical(splits)
  if (is.null(dim(mat))) {
    # Single split → mat may be a vector
    mat <- matrix(mat, nrow = 1L)
    tip_names <- tr$tip.label
  } else {
    tip_names <- colnames(mat)
  }
  if (is.null(tip_names)) tip_names <- tr$tip.label
  row_keys <- apply(mat, 1, function(r) {
    side_a <- sort(tip_names[r])
    side_b <- sort(tip_names[!r])
    if (side_a[1] < side_b[1]) paste(side_a, collapse = ",")
    else paste(side_b, collapse = ",")
  })
  paste(sort(row_keys), collapse = "|")
}

# ---------------------------------------------------------------------------
# Move drivers
# ---------------------------------------------------------------------------

# Pure-R re-implementation of gibbs_subtree_swap_impl at β=0.
# Inputs/outputs match the .flat_posterior_chain protocol so the test
# harness reuses .chi2_uniform below.
#
# At β=0 the Gibbs weights are: wOrig = 1, ws[pi] = 1 for every partner.
# So the sampler draws UNIFORMLY from {self, p_1, ..., p_K}.
# This is exactly what the C++ implementation does at β=0; no MH step.
.propose_gibbs_swap_beta0 <- function(tree, tree_length, rel_br_lengths) {
  edge   <- tree$edge
  n_tip  <- length(tree$tip.label)
  n_edge <- nrow(edge)
  child  <- edge[, 2L]

  # Pick a random edge child (matches gibbs_subtree_swap_impl line 1665).
  pick_idx <- sample.int(n_edge, 1L)
  nodeA <- child[pick_idx]

  partners <- MkPrime:::get_valid_swap_partners_cpp(edge, n_tip, nodeA)
  if (length(partners) == 0L) {
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                logHastings = -Inf))
  }

  # At β=0: uniform draw over {self} ∪ partners.
  K <- length(partners)
  draw <- sample.int(K + 1L, 1L)
  if (draw == 1L) {
    # self-draw → no-op (matches mcmc.cpp:1839)
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                logHastings = -Inf))
  }
  nodeB <- partners[draw - 1L]

  res <- MkPrime:::swap_subtrees_cpp(edge, n_tip, tree_length,
                                      rel_br_lengths, nodeA, nodeB)
  if (!is.finite(res$logHastings)) {
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                logHastings = -Inf))
  }

  newTree <- tree
  newTree$edge <- res$edge
  newTree$edge.length <- tree_length * res$rel_br_lengths

  # CRUCIAL: emulate "Gibbs commits without MH". We return logHastings = Inf
  # so that the .flat_posterior_chain driver ALWAYS accepts. This faithfully
  # mirrors gibbs_subtree_swap_impl's behaviour: it never rejects after the
  # uniform draw lands on a non-self partner.
  list(tree = newTree, rel_br_lengths = res$rel_br_lengths,
       logHastings = Inf)
}

# Pure-R re-implementation of weighted_subtree_swap_impl at β=0, with the
# code's actual (partner-count-omitting) logHR.
#
# At β=0:
#   wOrig          = 1
#   candW[pi][b]   = 1 for all (pi, b)
#   mCand[pi]      = nBins for all pi
#   sumM           = nBins * (nPartners + 1)
#
# Topology: probability of selecting candidate pi = nBins / sumM = 1/(K+1).
#   ⇒ uniform over {self, p_1, ..., p_K}.
# Bin: uniform over 1..nBins within the chosen topology.
# Beta(αNew, βNew) draw at chosenBin midpoint.
# Default fraction f_default = absLen[rowB] / (absLen[rowA] + absLen[rowB])
# decides defaultBin in 1..nBins.
#
# logHR = log(wOrig) + dbeta(f_default; α_old, β_old, log=TRUE)
#       - log(candW[chosen][chosenBin]) - dbeta(f_new; α_new, β_new, log=TRUE)
#       = dbeta(f_default; α_old, β_old, log=TRUE)
#       - dbeta(f_new;     α_new, β_new, log=TRUE)
#
# THIS IS THE FORMULA AS WRITTEN IN mcmc.cpp:2980-2983. By construction it
# omits the partner-count log-ratio. We submit the proposal to an MH step
# with exactly this logHR so the test reproduces the buggy code path
# faithfully.
.propose_weighted_swap_beta0 <- function(tree, tree_length, rel_br_lengths,
                                          n_bins = 5L, concentration = 100) {
  edge   <- tree$edge
  n_tip  <- length(tree$tip.label)
  n_edge <- nrow(edge)
  child  <- edge[, 2L]

  pick_idx <- sample.int(n_edge, 1L)
  nodeA <- child[pick_idx]

  partners <- MkPrime:::get_valid_swap_partners_cpp(edge, n_tip, nodeA)
  if (length(partners) == 0L) {
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                logHastings = -Inf))
  }

  rowA <- which(child == nodeA)[1L]

  # Uniform topology draw over {self, partners}.
  K <- length(partners)
  draw <- sample.int(K + 1L, 1L)
  if (draw == 1L) {
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                logHastings = -Inf))
  }
  nodeB <- partners[draw - 1L]
  rowB  <- which(child == nodeB)[1L]

  # Bin grid (uniform partition of (0,1), midpoints at (b - 0.5)/nBins).
  breaks <- seq(0, 1, length.out = n_bins + 1L)
  mids   <- (breaks[-1] + breaks[-length(breaks)]) / 2

  # Sample bin uniformly (at β=0 all bin weights are equal).
  chosen_bin <- sample.int(n_bins, 1L)
  chosen_mid <- mids[chosen_bin]
  alpha_new <- chosen_mid * concentration + 1
  beta_new  <- (1 - chosen_mid) * concentration + 1
  f_new <- rbeta(1, alpha_new, beta_new)
  f_new <- max(min(f_new, 1 - 1e-8), 1e-8)

  total <- rel_br_lengths[rowA] + rel_br_lengths[rowB]
  abs_total <- tree_length * total

  # f_default = current absLen[rowB] / (absLen[rowA] + absLen[rowB])
  f_default <- rel_br_lengths[rowB] / total
  default_bin <- max(1L, min(n_bins, findInterval(f_default, breaks,
                                                   rightmost.closed = TRUE,
                                                   all.inside = TRUE)))
  default_mid <- mids[default_bin]
  alpha_old <- default_mid * concentration + 1
  beta_old  <- (1 - default_mid) * concentration + 1

  # Build the swapped tree with new branch fractions
  new_rel <- rel_br_lengths
  new_rel[rowA] <- f_new * total
  new_rel[rowB] <- (1 - f_new) * total

  res <- MkPrime:::swap_subtrees_cpp(edge, n_tip, tree_length, new_rel,
                                      nodeA, nodeB)
  if (!is.finite(res$logHastings)) {
    return(list(tree = tree, rel_br_lengths = rel_br_lengths,
                logHastings = -Inf))
  }

  newTree <- tree
  newTree$edge <- res$edge
  newTree$edge.length <- tree_length * res$rel_br_lengths

  # logHR as written at mcmc.cpp:2980-2983 (partner-count log-ratio absent).
  log_hr <- dbeta(f_default, alpha_old, beta_old, log = TRUE) -
            dbeta(f_new,     alpha_new, beta_new, log = TRUE)

  list(tree = newTree, rel_br_lengths = res$rel_br_lengths,
       logHastings = log_hr)
}

# ---------------------------------------------------------------------------
# Driver: flat target, MH-style chain over topology keys
# ---------------------------------------------------------------------------

.flat_posterior_chain <- function(propose_fn, n_iter, thin, seed,
                                   n_tip = 5L) {
  set.seed(seed)
  tree <- ape::rtree(n_tip)
  tree <- ape::unroot(tree)
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)

  current <- tree
  current_rel <- tree$edge.length / tl
  n_samples <- n_iter %/% thin
  keys <- character(n_samples)

  for (i in seq_len(n_iter)) {
    prop <- propose_fn(current, tl, current_rel)
    if (is.finite(prop$logHastings) &&
        log(runif(1)) < prop$logHastings) {
      current <- prop$tree
      current_rel <- prop$rel_br_lengths
    } else if (identical(prop$logHastings, Inf)) {
      # Gibbs-commit path (no MH): always accept.
      current <- prop$tree
      current_rel <- prop$rel_br_lengths
    }
    if (i %% thin == 0L) {
      keys[i %/% thin] <- .topo_key(current)
    }
  }
  table(keys)
}

# Diagnostic: did this proposal type ever change the topology?
# (Sanity that we ARE applying the move under test, not silently rejecting.)
.acceptance_diag <- function(propose_fn, n_iter, seed, n_tip = 5L) {
  set.seed(seed)
  tree <- ape::rtree(n_tip)
  tree <- ape::unroot(tree)
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  current <- tree
  current_rel <- tree$edge.length / tl
  topology_changes <- 0L
  cur_key <- .topo_key(current)
  for (i in seq_len(n_iter)) {
    prop <- propose_fn(current, tl, current_rel)
    accept <- (is.finite(prop$logHastings) && log(runif(1)) < prop$logHastings) ||
              identical(prop$logHastings, Inf)
    if (accept) {
      current <- prop$tree
      current_rel <- prop$rel_br_lengths
      new_key <- .topo_key(current)
      if (new_key != cur_key) topology_changes <- topology_changes + 1L
      cur_key <- new_key
    }
  }
  topology_changes / n_iter
}

# ---------------------------------------------------------------------------
# Pass criterion
# ---------------------------------------------------------------------------

.expected_buckets <- N_TOPOS_EXPECTED
.alpha            <- 0.005

.assess_uniform <- function(freq, label) {
  n_buckets <- length(freq)
  chi2 <- if (n_buckets >= 2L) {
    suppressWarnings(chisq.test(as.numeric(freq)))
  } else {
    list(p.value = 0, statistic = NA_real_)
  }
  ratio <- if (n_buckets >= 2L && min(freq) > 0) {
    max(freq) / min(freq)
  } else {
    Inf
  }
  # SBC-HARNESS-001 fix (2026-05-28): max/min < 2 is too strict above ~50
  # buckets — at n=6 the expected count is ~48/bucket and Poisson noise
  # under H0 gives max/min ≈ 2.5 (job 17294906: NNI control chi² p=0.4514,
  # max/min=2.57 → spurious FAIL). chi² already encodes Poisson variance;
  # defer to it once binning is fine-grained.
  pass <- (n_buckets == .expected_buckets) &&
          (chi2$p.value > .alpha) &&
          (n_buckets > 50L || ratio < 2)
  list(label = label,
       freq = freq,
       n_buckets = n_buckets,
       p_value = chi2$p.value,
       chi2_stat = unname(chi2$statistic),
       max_over_min = ratio,
       pass = pass)
}

# ---------------------------------------------------------------------------
# Configurations
# ---------------------------------------------------------------------------

cfg <- if (quick) {
  list(n_iter = 50000L, thin = 25L, label = "quick")
} else {
  list(n_iter = 1000000L, thin = 200L, label = "full")
}

cat(sprintf("[D2b] Running detailed-balance chains: %s (n_iter=%d, thin=%d)\n",
            cfg$label, cfg$n_iter, cfg$thin))

# Each chain uses a distinct seed (independent test).
arms <- list(
  list(name = "nni_control",
       seed = 2025L,
       propose = function(tr, tl, rel) MkPrime:::ProposeNni(tr, tl, rel)),
  list(name = "gibbs_subtree_swap",
       seed = 3071L,
       propose = .propose_gibbs_swap_beta0),
  list(name = "weighted_subtree_swap",
       seed = 4099L,
       propose = .propose_weighted_swap_beta0)
)

results <- list()
for (arm in arms) {
  cat(sprintf("[D2b]   %-22s ", arm$name))
  t0 <- Sys.time()
  freq <- .flat_posterior_chain(arm$propose,
                                 n_iter = cfg$n_iter,
                                 thin   = cfg$thin,
                                 seed   = arm$seed,
                                 n_tip  = N_TIP_CFG)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  assess <- .assess_uniform(freq, arm$name)
  results[[arm$name]] <- list(arm = arm$name,
                              elapsed_sec = dt,
                              freq = freq,
                              assessment = assess)
  cat(sprintf("%.1fs  p=%.4g  n_topos=%d  max/min=%.2f  %s\n",
              dt, assess$p_value, assess$n_buckets,
              assess$max_over_min,
              if (assess$pass) "PASS" else "FAIL"))
}

# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------

control_pass <- results[["nni_control"]]$assessment$pass
gibbs_pass   <- results[["gibbs_subtree_swap"]]$assessment$pass
weighted_pass <- results[["weighted_subtree_swap"]]$assessment$pass

# Headline interpretation rules:
# - If control passes and a swap arm fails → CONFIRMS the corresponding
#   OPEN finding (SWAP-001 / SWAP-002).
# - If control fails → harness itself is suspect; report WARN.
# - If control passes and both swap arms pass → SWAP findings refuted
#   (at this sample size).

overall <- if (!control_pass) {
  "WARN"
} else if (gibbs_pass && weighted_pass) {
  "PASS"
} else {
  "FAIL"
}

summary <- list(
  cfg = cfg,
  control_arm = results[["nni_control"]],
  gibbs_arm   = results[["gibbs_subtree_swap"]],
  weighted_arm = results[["weighted_subtree_swap"]],
  control_pass = control_pass,
  gibbs_pass   = gibbs_pass,
  weighted_pass = weighted_pass,
  verdict = overall
)
saveRDS(summary, file.path(out_dir, "summary.rds"))

fmt_arm <- function(r) {
  a <- r$assessment
  paste(
    sprintf("  %s:", r$arm),
    sprintf("    elapsed         = %.1f s", r$elapsed_sec),
    sprintf("    distinct topos  = %d / %d", a$n_buckets, N_TOPOS_EXPECTED),
    sprintf("    chi2 statistic  = %.3f", a$chi2_stat),
    sprintf("    p-value         = %.4g", a$p_value),
    sprintf("    max/min         = %.2f", a$max_over_min),
    sprintf("    PASS?           = %s", a$pass),
    sep = "\n"
  )
}

verdict_lines <- c(
  sprintf("VERDICT: %s", overall),
  "",
  sprintf("Sub-harness:    D2(b) — subtree-swap β=0 detailed-balance"),
  sprintf("Mode:           %s", cfg$label),
  sprintf("n_iter:         %d", cfg$n_iter),
  sprintf("thin:           %d", cfg$thin),
  sprintf("alpha (chi2):   %.4f", .alpha),
  "",
  fmt_arm(results[["nni_control"]]),
  "",
  fmt_arm(results[["gibbs_subtree_swap"]]),
  "",
  fmt_arm(results[["weighted_subtree_swap"]]),
  "",
  "Interpretation:",
  if (overall == "PASS") {
    "  Control passes; both swap arms pass. SWAP-001/SWAP-002 not detected."
  } else if (overall == "FAIL") {
    paste(c(
      "  Control passes; one or more swap arms FAIL the chi-squared.",
      if (!gibbs_pass)   "  → CONFIRMS SWAP-001 (GibbsSubtreeSwap empirically biased).",
      if (!weighted_pass) "  → CONFIRMS SWAP-002 (WeightedSubtreeSwap empirically biased)."
    ), collapse = "\n")
  } else {
    "  WARN: NNI positive control failed — harness needs inspection."
  }
)
writeLines(verdict_lines, file.path(out_dir, "verdict.txt"))

cat("\n", paste(verdict_lines, collapse = "\n"), "\n", sep = "")
cat("\nWrote:", file.path(out_dir, "verdict.txt"), "\n")
cat("Wrote:", file.path(out_dir, "summary.rds"), "\n")

if (overall == "FAIL") {
  quit(status = 1L)
}
