# Local smoke test for the empirical_geometric prior, mirroring the
# Hamilton pipeline at a small scale.
#
# Simulates characters under JC(k_true) on a small tree, runs the
# empirical_geometric and geometric Mk' arms, and reports CID + posterior
# u summaries.

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ape)
})

set.seed(2026)

# ---- Simulate a small dataset -------------------------------------------------
nTip <- 8L
nCharTotal <- 80L
true_tree <- rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
true_tree$edge.length <- runif(nrow(true_tree$edge), 0.05, 0.25)
true_tree <- TreeTools::Preorder(true_tree)
cat("True tree length:", sum(true_tree$edge.length), "\n")

# Mix of k' = 2, 3, 4, 5 — half the characters get k' > 2 so empirical
# prior can demonstrate its effect.
kTrueAll <- sample(c(2L, 3L, 4L, 5L), nCharTotal, replace = TRUE,
                    prob = c(0.40, 0.30, 0.20, 0.10))
sim_mat <- matrix(NA_integer_, nTip, nCharTotal,
                   dimnames = list(true_tree$tip.label, NULL))
for (ch in seq_len(nCharTotal)) {
  kTrue <- kTrueAll[ch]
  node_states <- integer(2 * nTip - 1)
  rootIdx <- nTip + 1L
  node_states[rootIdx] <- sample.int(kTrue, 1L) - 1L
  edges <- true_tree$edge
  el    <- true_tree$edge.length
  for (e in rev(seq_len(nrow(edges)))) {
    pa <- edges[e, 1]; ch2 <- edges[e, 2]; t <- el[e]
    pSame <- 1 / kTrue + (1 - 1 / kTrue) * exp(-kTrue * t)
    if (runif(1) < pSame) {
      node_states[ch2] <- node_states[pa]
    } else {
      node_states[ch2] <- sample(setdiff(seq.int(0, kTrue - 1L),
                                           node_states[pa]), 1L)
    }
  }
  sim_mat[, ch] <- node_states[seq_len(nTip)]
}
variable <- apply(sim_mat, 2, function(x) length(unique(x)) > 1)
sim_mat <- sim_mat[, variable, drop = FALSE]
kTrueAll <- kTrueAll[variable]

pd <- TreeTools::MatrixToPhyDat(sim_mat)
mkd <- MkPrimeData(pd)
cat(sprintf("Variable characters: %d\n", ncol(sim_mat)))
cat("kObs distribution: "); print(table(mkd$kObs))
cat("kTrue distribution: "); print(table(kTrueAll))
unseen <- sum(kTrueAll - mkd$kObs)
cat(sprintf("Total unseen states across all chars: %d\n", unseen))

# ---- NJ starting tree --------------------------------------------------------
start_tree <- TreeTools::NJTree(pd, edgeLengths = TRUE)

# ---- Run each arm ------------------------------------------------------------
runArm <- function(label, prior, seed) {
  cat(sprintf("\n=== Running arm: %s ===\n", label))
  set.seed(seed)
  t0 <- Sys.time()
  res <- RunMkPrime(
    mkd, start_tree,
    model = MkPrimeModel(coding = "variable", kPrimePrior = prior,
                          expSteps = sum(true_tree$edge.length)),
    mcmc = MkPrimeMCMC(nIter = 4000L, thin = 20L,
                       maxWarmup = 2000L, minWarmup = 2000L,
                       autoTune = FALSE,
                       nRuns = 2L, nChains = 2L)
  )
  cat(sprintf("  Wall time: %s\n",
              format(Sys.time() - t0, digits = 3)))

  # Per-move acceptance — catches the failure mode where a sampler is
  # crushed to its weight floor because its acceptance rate is too low to
  # contribute to ESS.  Under the empirical_geometric prior we expect
  # `mh_logit_p` to accept ~20–40%; anything below 5% indicates the move
  # is broken (e.g. multiplicative MH overshooting p = 1).
  if (!is.null(res$acceptance)) {
    cat("  Per-move acceptance:\n")
    print(round(res$acceptance, 3))
    if ("mh_logit_p" %in% names(res$acceptance)) {
      mhlpRate <- res$acceptance[["mh_logit_p"]]
      cat(sprintf("  mh_logit_p acceptance: %.3f\n", mhlpRate))
      if (is.finite(mhlpRate) && mhlpRate < 0.05) {
        warning("mh_logit_p acceptance below 5% — sampler likely broken.")
      }
    } else if ("mh_p" %in% names(res$acceptance)) {
      mhpRate <- res$acceptance[["mh_p"]]
      cat(sprintf("  mh_p acceptance: %.3f\n", mhpRate))
    }
  }

  # Posterior u per character
  kpCols <- grep("^kPrime_", colnames(res$samples), value = TRUE)
  if (length(kpCols) == 0L) {
    stop("No kPrime columns in samples")
  }
  kPostMean <- colMeans(res$samples[, kpCols, drop = FALSE])
  uPostMean <- kPostMean - mkd$kObs
  uPostMed <- apply(res$samples[, kpCols, drop = FALSE], 2,
                     function(x) median(x - mkd$kObs[match(NA, NA)]))
  # Reorder per char_idx (kpCols are kPrime_1..kPrime_n)
  kPrimeIdx <- as.integer(sub("kPrime_", "", kpCols))
  ord <- order(kPrimeIdx)
  uPostMean <- uPostMean[ord]

  cat(sprintf("  Posterior mean u (total): %.2f  (truth: %d)\n",
              sum(uPostMean), unseen))
  cat(sprintf("  Mean of per-char posterior u: %.3f\n", mean(uPostMean)))
  cat(sprintf("  Posterior u > 0 in %d / %d chars\n",
              sum(uPostMean > 0.1), length(uPostMean)))

  # CID against true tree
  if (length(res$trees) > 5L) {
    suppressPackageStartupMessages(library(TreeDist))
    trees <- lapply(res$trees, ape::reorder.phylo, order = "cladewise")
    class(trees) <- "multiPhylo"
    cid <- as.numeric(TreeDist::ClusteringInfoDistance(
      trees, ape::reorder.phylo(true_tree, "cladewise"), normalize = TRUE))
    cat(sprintf("  CID to true tree: mean=%.3f  median=%.3f  (n=%d trees)\n",
                mean(cid), median(cid), length(cid)))
  }

  list(uPostMean = uPostMean, samples = res$samples, trees = res$trees)
}

res_eg  <- runArm("empirical_geometric", "empirical_geometric", 101)
res_geo <- runArm("geometric",           "geometric",           101)

# ---- Comparison --------------------------------------------------------------
cat("\n=== Side-by-side ===\n")
cat(sprintf("Truth total unseen states:                  %4d\n", unseen))
cat(sprintf("Sum posterior u  (empirical_geometric):     %.2f\n",
            sum(res_eg$uPostMean)))
cat(sprintf("Sum posterior u  (geometric):               %.2f\n",
            sum(res_geo$uPostMean)))
cat(sprintf("Empirical recovers %.1f%% of unseen states\n",
            100 * sum(res_eg$uPostMean) / max(1, unseen)))
cat(sprintf("Geometric  recovers %.1f%% of unseen states\n",
            100 * sum(res_geo$uPostMean) / max(1, unseen)))

# Sanity checks
stopifnot(sum(res_eg$uPostMean) >= 0)
stopifnot(sum(res_geo$uPostMean) >= 0)

cat("\n--- Smoke test complete: no errors. ---\n")
