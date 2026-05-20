# Real-data smoke test of ecology-aware MkPrime on the rodent matrix.
# Lessons embedded from the first pass (4.7 min, 5k iter, mbank_X24848):
#   * Extinct taxa (col 221 == "1") have no ecology coding (all 42 NA).
#     Dropping them gives 60 fully ecology-coded extant taxa, faster
#     mixing, and a cleaner posterior to interpret.
#   * The default 500-iter warmup cap was insufficient — logP was still
#     climbing.  Raise both warmup caps.
#   * gibbs_z weight 0.1% was a side-effect of the kPrime moves' raw
#     weights (nChar each) dominating the proposal budget.  Bumped in
#     .BuildMoves; expect ~1% now.
#   * Use ReadMkLog (streaming-mode log format), not read.table.
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})

nex <- "C:/Users/pjjg18/Downloads/mbank_X24848_2026-5-9-1135.nex"
stopifnot(file.exists(nex))

cat("Reading nexus...\n")
pd <- ReadAsPhyDat(nex)
mat <- as.matrix(pd)
cat("Full dataset: ", nrow(mat), " tips x ", ncol(mat), " chars\n", sep = "")

# Identify extant tips (col 221 == "0"); drop extinct.
extantTips <- rownames(mat)[mat[, 221] == "0"]
cat("Extant tips:", length(extantTips), " (dropping ",
    nrow(mat) - length(extantTips), " extinct)\n", sep = "")
pdExtant <- pd[extantTips]

# Build ecology-aware MkPrimeData.  Col 220 = ecology; the function
# drops it from the character matrix and exposes it as $ecology.
mkd <- MkPrimeData(pdExtant, ecology = 220L)
cat("kEcology: ", mkd$kEcology,
    "; nChar: ", mkd$nChar,
    "; nTip: ", mkd$nTip, "\n", sep = "")
cat("ecology table (post-extinct-drop):\n")
print(table(mkd$ecology, useNA = "ifany"))

# Starting tree: random; the MCMC will rearrange.
set.seed(1)
startTree <- Preorder(ape::rtree(mkd$nTip, tip.label = mkd$taxon_names))

model <- MkPrimeModel(
  ecologyAware  = TRUE,
  magnitudeMode = "global",
  kPrimePrior   = "geometric",
  coding        = "variable",
  expSteps      = 10
)

mcmc <- MkPrimeMCMC(
  nIter          = 10000L,
  nChains        = 1L,
  nRuns          = 1L,
  thin           = 100L,
  treeThin       = 500L,
  minWarmup      = 500L,
  maxWarmup      = 3000L,
  logFile        = "rodent-ecology-smoke.log",
  checkpointFile = NULL
)

for (f in c("rodent-ecology-smoke.ckp", "rodent-ecology-smoke.log")) {
  if (file.exists(f)) file.remove(f)
}
cat("Starting RunMkPrime...\n")
t0 <- Sys.time()
res <- RunMkPrime(mkd, tree = startTree, model = model, mcmc = mcmc)
cat("Elapsed:", format(Sys.time() - t0), "\n\n")

# Parse log via ReadMkLog (streaming-mode tabular format).
cat("== Posterior summary ==\n")
cat("logFile:", res$logFile, "\n")
samples <- ReadMkLog(res$logFile)
cat("Log rows: ", nrow(samples),
    "; columns: ", ncol(samples), "\n", sep = "")
cat("First few column names:",
    paste(head(colnames(samples), 8), collapse = ", "), "...\n")

summariseCol <- function(nm) {
  if (!nm %in% colnames(samples)) return(invisible())
  x <- samples[, nm]
  cat(sprintf("%-16s: mean=%.4g  sd=%.4g  q025=%.4g  q500=%.4g  q975=%.4g\n",
              nm, mean(x), stats::sd(x),
              stats::quantile(x, 0.025, names = FALSE),
              stats::quantile(x, 0.500, names = FALSE),
              stats::quantile(x, 0.975, names = FALSE)))
}
for (nm in c("log_likelihood", "log_posterior", "tree_length", "phi", "pi0",
             "rate_log_sd")) summariseCol(nm)

# z_samples summary: per-cell fraction with non-zero state.
if (!is.null(res$z_samples) && length(res$z_samples) > 0) {
  cat("\nz_samples: ", length(res$z_samples), " snapshots\n", sep = "")
  zArr <- simplify2array(res$z_samples)
  if (length(dim(zArr)) == 3) {
    pNonZero <- apply(zArr != 0, c(1, 2), mean)
    cat(sprintf("Per-cell P(z != 0): min=%.3f mean=%.3f max=%.3f\n",
                min(pNonZero), mean(pNonZero), max(pNonZero)))
    ord <- order(-pNonZero)[seq_len(min(20, length(pNonZero)))]
    idx <- arrayInd(ord, dim(pNonZero))
    df <- data.frame(char = idx[, 1], ecology = idx[, 2] - 1L,
                     pNonZero = round(pNonZero[ord], 3))
    cat("Top 20 (char, ecology) cells by P(z != 0):\n")
    print(df, row.names = FALSE)
  }
}
cat("\nDone.\n")
