# Step 6b: replicate the Step 1 + acceptance findings on additional reps
# and compare eg vs non-eg arm on the SAME dataset.

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(TreeTools)
})

DATA_ROOT <- "C:/Users/pjjg18/GitHub/mkprime/tree-inference"

.LoadRep <- function(tree_idx, rep_idx) {
  rep_dir <- file.path(DATA_ROOT, sprintf("tree_%02d", tree_idx),
                       sprintf("rep_%02d", rep_idx))
  chr_files <- list.files(rep_dir, pattern = "\\.nex$", full.names = TRUE)
  mats <- lapply(chr_files, function(f) {
    ch <- TreeTools::ReadCharacters(f)
    if (is.matrix(ch)) ch else matrix(ch, ncol = 1L,
                                       dimnames = list(names(ch), NULL))
  })
  tips <- rownames(mats[[1]])
  mat <- do.call(cbind, lapply(mats, function(m) m[tips, , drop = FALSE]))
  mode(mat) <- "character"
  pd <- TreeTools::MatrixToPhyDat(mat)
  list(mkd = MkPrimeData(pd), tree = TreeTools::NJTree(pd, edgeLengths = TRUE))
}

.RunOne <- function(mkd, tree, prior, seed, n_iter = 5000L) {
  set.seed(seed)
  res <- RunMkPrime(
    mkd, tree,
    model = MkPrimeModel(coding = "variable", kPrimePrior = prior),
    mcmc = MkPrimeMCMC(nIter = n_iter, thin = 50L,
                       maxWarmup = 1000L, minWarmup = 1000L,
                       autoTune = FALSE,
                       nRuns = 1L, nChains = 2L,
                       progressFn = function(...) invisible())
  )
  lp <- res$samples[, "log_posterior"]
  ll <- res$samples[, "log_likelihood"]
  lpr <- lp - ll
  list(
    res = res,
    n = length(lp),
    logL_sd = sd(ll), logPr_sd = sd(lpr), logP_sd = sd(lp),
    cor_LP = cor(ll, lpr),
    accept = res$chain_acceptance[[1]],
    swap = res$swap_rates
  )
}

cat("=== Replicate check across 3 reps under empirical_geometric ===\n\n")
summary_rows <- list()
for (rep in 1:3) {
  cat(sprintf("--- tree_01/rep_%02d (eg) ---\n", rep))
  d <- .LoadRep(1L, rep)
  r <- .RunOne(d$mkd, d$tree, "empirical_geometric", seed = 100 + rep)
  cat(sprintf("  logL sd=%.2f  logPr sd=%.2f  logP sd=%.2f  cor(L, Pr)=%+0.3f  swap=%.4f\n",
              r$logL_sd, r$logPr_sd, r$logP_sd, r$cor_LP, r$swap))
  cat(sprintf("  spr=%.3f  tbr=%.3f  pspr=%.3f  nni=%.3f  mh_logit_p=%.3f  block_kP=%.3f\n",
              r$accept["spr"], r$accept["tbr"], r$accept["pspr"],
              r$accept["nni"], r$accept["mh_logit_p"], r$accept["block_kPrime"]))
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    rep = rep, arm = "eg", n = r$n,
    logL_sd = r$logL_sd, logPr_sd = r$logPr_sd, logP_sd = r$logP_sd,
    cor_LP = r$cor_LP, swap = r$swap,
    spr = r$accept["spr"], tbr = r$accept["tbr"], pspr = r$accept["pspr"],
    nni = r$accept["nni"], block_kP = r$accept["block_kPrime"]
  )
}

cat("\n=== Same datasets, non-eg arm (kPrimePrior=\"beta_geometric\") ===\n")
for (rep in 1:3) {
  cat(sprintf("--- tree_01/rep_%02d (bg) ---\n", rep))
  d <- .LoadRep(1L, rep)
  r <- .RunOne(d$mkd, d$tree, "beta_geometric", seed = 200 + rep)
  cat(sprintf("  logL sd=%.2f  logPr sd=%.2f  logP sd=%.2f  cor(L, Pr)=%+0.3f  swap=%.4f\n",
              r$logL_sd, r$logPr_sd, r$logP_sd, r$cor_LP, r$swap))
  cat(sprintf("  spr=%.3f  tbr=%.3f  pspr=%.3f  nni=%.3f\n",
              r$accept["spr"], r$accept["tbr"], r$accept["pspr"], r$accept["nni"]))
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    rep = rep, arm = "bg", n = r$n,
    logL_sd = r$logL_sd, logPr_sd = r$logPr_sd, logP_sd = r$logP_sd,
    cor_LP = r$cor_LP, swap = r$swap,
    spr = r$accept["spr"], tbr = r$accept["tbr"], pspr = r$accept["pspr"],
    nni = r$accept["nni"], block_kP = NA_real_
  )
}

cat("\n=== Summary table ===\n")
all_sum <- do.call(rbind, summary_rows)
print(all_sum, row.names = FALSE)
saveRDS(all_sum, "data-raw/step6b-replicate-summary.rds")
cat("\nSaved data-raw/step6b-replicate-summary.rds\n")
