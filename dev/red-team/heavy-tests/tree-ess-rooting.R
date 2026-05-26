# dev/red-team/heavy-tests/tree-ess-rooting.R
#
# Sub-harness D2(a): TreeESS root-dependence audit.
#
# CLAIM CHALLENGED
#   src/tree_ess.cpp (via R/treeESS.R::TreeESS) should produce identical
#   ESS estimates for a posterior of trees and for any per-tree re-rooted
#   version of that posterior, because Robinson-Foulds distance is an
#   unrooted-tree invariant (bipartitions are root-independent).
#
# PASS CRITERION
#   (i)  The full RF distance matrix is BIT-IDENTICAL between the original
#        and randomly re-rooted posterior. We assert identical(dmat_orig,
#        dmat_root) == TRUE.
#   (ii) medianPseudoESS and frechetCorrelationESS agree to within
#        floating-point tolerance (|a - b| < 1e-8).
#
#   If (i) holds but (ii) fails, the ESS routine itself has a non-deterministic
#   path (RNG, subsample seed). If (i) fails, RF computation is silently
#   root-dependent for these trees — a real bug.
#
# Project memory:
#   project_scoring_bug — ape::prop.part is root-dependent. TreeDist::RobinsonFoulds
#   computes splits via TreeTools::as.Splits which IS root-invariant, so the
#   expectation here is exact equality, not just Monte-Carlo equivalence.
#
# Usage:
#   Rscript dev/red-team/heavy-tests/tree-ess-rooting.R --quick
#   Rscript dev/red-team/heavy-tests/tree-ess-rooting.R
#
# Output:
#   dev/red-team/heavy-tests/tree-ess-rooting-results/summary.rds
#   dev/red-team/heavy-tests/tree-ess-rooting-results/verdict.txt

suppressPackageStartupMessages({
  library(ape)
  library(TreeTools)
  library(TreeDist)
})

quick <- "--quick" %in% commandArgs(trailingOnly = TRUE)

# Load package via devtools / pkgload — no install required.
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(file.path(getwd()), quiet = TRUE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(file.path(getwd()), quiet = TRUE)
} else {
  stop("Need {pkgload} or {devtools}")
}

out_dir <- "dev/red-team/heavy-tests/tree-ess-rooting-results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Generator
# ---------------------------------------------------------------------------

# Build a "posterior-like" sequence via an NNI random walk seeded at a random
# starting tree. This gives a correlated chain of trees on n tips, matching
# the kind of input TreeESS() is designed for.
nni_random_walk <- function(n_tip, n_trees, seed) {
  set.seed(seed)
  start <- ape::rtree(n_tip, br = NULL)
  start <- ape::unroot(start)
  start <- TreeTools::Preorder(start)
  tl <- if (is.null(start$edge.length)) {
    start$edge.length <- rep(1.0, nrow(start$edge))
    nrow(start$edge)
  } else {
    sum(start$edge.length)
  }
  rel <- start$edge.length / tl
  trees <- vector("list", n_trees)
  cur <- start
  cur_rel <- rel
  for (i in seq_len(n_trees)) {
    prop <- MkPrime:::ProposeNni(cur, tl, cur_rel)
    if (is.finite(prop$logHastings) && log(runif(1)) < prop$logHastings) {
      cur <- prop$tree
      cur_rel <- prop$rel_br_lengths
    }
    trees[[i]] <- cur
  }
  class(trees) <- "multiPhylo"
  trees
}

random_reroot <- function(trees, seed) {
  set.seed(seed)
  lapply(trees, function(tr) {
    tips <- tr$tip.label
    outgroup <- sample(tips, 1L)
    suppressWarnings(ape::root(tr, outgroup = outgroup, resolve.root = TRUE))
  }) -> rerooted
  class(rerooted) <- "multiPhylo"
  rerooted
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

cfg <- if (quick) {
  list(n_tip = 6L, n_trees = 200L, seed_walk = 4711L, seed_reroot = 9931L)
} else {
  list(n_tip = 8L, n_trees = 10000L, seed_walk = 4711L, seed_reroot = 9931L)
}

cat(sprintf("[D2a] Generating posterior: n_tip=%d, n_trees=%d %s\n",
            cfg$n_tip, cfg$n_trees, if (quick) "(QUICK)" else "(FULL)"))

t0 <- Sys.time()
posterior <- nni_random_walk(cfg$n_tip, cfg$n_trees, cfg$seed_walk)
rerooted  <- random_reroot(posterior, cfg$seed_reroot)
cat(sprintf("[D2a] Generated in %.2f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Sanity: trees aren't degenerate to a single topology.
n_unique <- length(unique(vapply(posterior, function(tr) {
  paste(sort(apply(as.logical(TreeTools::as.Splits(tr)), 1, function(r) {
    a <- sort(tr$tip.label[r])
    b <- sort(tr$tip.label[!r])
    if (a[1] < b[1]) paste(a, collapse = ",") else paste(b, collapse = ",")
  })), collapse = "|")
}, character(1))))
cat(sprintf("[D2a] Distinct topologies in posterior: %d\n", n_unique))

# (i) Distance-matrix invariance
cat("[D2a] Computing RF distance matrices...\n")
dmat_orig <- as.matrix(TreeDist::RobinsonFoulds(posterior))
dmat_root <- as.matrix(TreeDist::RobinsonFoulds(rerooted))

# Make order-independent: compare entries directly
dmat_identical  <- identical(dmat_orig, dmat_root)
dmat_max_absdiff <- max(abs(dmat_orig - dmat_root))

# (ii) ESS invariance
cat("[D2a] Computing TreeESS...\n")
ess_orig <- MkPrime::TreeESS(posterior, frechet = TRUE)
ess_root <- MkPrime::TreeESS(rerooted,  frechet = TRUE)

ess_diff_med  <- abs(ess_orig["medianPseudoESS"]    - ess_root["medianPseudoESS"])
ess_diff_frech <- abs(ess_orig["frechetCorrelationESS"] - ess_root["frechetCorrelationESS"])

# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------

tol_ess <- 1e-8
dmat_pass <- dmat_identical || dmat_max_absdiff < 1e-12
ess_pass  <- (ess_diff_med < tol_ess) && (ess_diff_frech < tol_ess)
overall   <- if (dmat_pass && ess_pass) "PASS" else "FAIL"

summary <- list(
  cfg               = cfg,
  n_unique_topos    = n_unique,
  dmat_identical    = dmat_identical,
  dmat_max_absdiff  = dmat_max_absdiff,
  ess_orig          = ess_orig,
  ess_root          = ess_root,
  ess_diff_med      = unname(ess_diff_med),
  ess_diff_frech    = unname(ess_diff_frech),
  tol_ess           = tol_ess,
  dmat_pass         = dmat_pass,
  ess_pass          = ess_pass,
  verdict           = overall
)
saveRDS(summary, file.path(out_dir, "summary.rds"))

verdict_lines <- c(
  sprintf("VERDICT: %s", overall),
  "",
  sprintf("Sub-harness:    D2(a) — TreeESS root-dependence audit"),
  sprintf("Mode:           %s", if (quick) "quick" else "full"),
  sprintf("n_tip:          %d", cfg$n_tip),
  sprintf("n_trees:        %d", cfg$n_trees),
  sprintf("distinct topos: %d", n_unique),
  "",
  "Distance matrix invariance:",
  sprintf("  identical(dmat_orig, dmat_root)     = %s", dmat_identical),
  sprintf("  max |dmat_orig - dmat_root|         = %.3e", dmat_max_absdiff),
  sprintf("  PASS (dmat invariant)?              = %s", dmat_pass),
  "",
  "ESS invariance:",
  sprintf("  ess_orig  medianPseudoESS           = %.6f",
          ess_orig["medianPseudoESS"]),
  sprintf("  ess_root  medianPseudoESS           = %.6f",
          ess_root["medianPseudoESS"]),
  sprintf("  |delta|   medianPseudoESS           = %.3e",
          unname(ess_diff_med)),
  sprintf("  ess_orig  frechetCorrelationESS     = %.6f",
          ess_orig["frechetCorrelationESS"]),
  sprintf("  ess_root  frechetCorrelationESS     = %.6f",
          ess_root["frechetCorrelationESS"]),
  sprintf("  |delta|   frechetCorrelationESS     = %.3e",
          unname(ess_diff_frech)),
  sprintf("  PASS (ESS within %.0e)?             = %s", tol_ess, ess_pass)
)
writeLines(verdict_lines, file.path(out_dir, "verdict.txt"))

cat("\n")
cat(paste(verdict_lines, collapse = "\n"), "\n")
cat("\nWrote:", file.path(out_dir, "verdict.txt"), "\n")
cat("Wrote:", file.path(out_dir, "summary.rds"), "\n")

if (overall == "FAIL") {
  quit(status = 1L)
}
