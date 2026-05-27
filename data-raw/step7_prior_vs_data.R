# Step 7: prior–data mismatch sanity check.
#
# Compute kObs across the 260 datasets (26 trees x 10 reps) used in
# Hamilton array job 17140607 and compare to empiricalNObs body.
#
# If the empirical prior places appreciable mass at higher k' while the
# datasets' kObs are concentrated at low k', the prior fundamentally
# fights the data and the eg arm cannot mix well even with a perfect
# sampler — no MCMC fix will help.
#
# Output: data-raw/step7-prior-vs-data.rds + console summary.

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(TreeTools)
})

DATA_ROOT <- "C:/Users/pjjg18/GitHub/mkprime/tree-inference"
stopifnot(dir.exists(DATA_ROOT))

trees <- sprintf("tree_%02d", 1:26)
reps  <- sprintf("rep_%02d", 1:10)

# kObs of a single .nex file (one character).
.KObsFile <- function(path) {
  chars <- tryCatch(TreeTools::ReadCharacters(path), error = function(e) NULL)
  if (is.null(chars)) return(NA_integer_)
  # ReadCharacters returns a tip x char matrix; under one-char files it
  # may collapse — handle both shapes.
  if (is.matrix(chars)) {
    apply(chars, 2, function(col) {
      toks <- col[!col %in% c("?", "-", "")]
      if (!length(toks)) return(0L)
      isPoly <- grepl("[(){},]", toks)
      states <- character(0)
      if (any(!isPoly)) states <- c(states, toks[!isPoly])
      if (any(isPoly)) {
        states <- c(states, unlist(strsplit(gsub("[(){},]", " ",
                                                  toks[isPoly]), "\\s+")))
      }
      length(unique(states[nzchar(states)]))
    })
  } else {
    toks <- chars[!chars %in% c("?", "-", "")]
    if (!length(toks)) return(0L)
    length(unique(toks))
  }
}

cat("Computing kObs for 26 trees × 10 reps...\n")
results <- list()
for (tr in trees) {
  for (rp in reps) {
    rep_dir <- file.path(DATA_ROOT, tr, rp)
    if (!dir.exists(rep_dir)) {
      cat("  missing:", rep_dir, "\n"); next
    }
    chr_files <- list.files(rep_dir, pattern = "\\.nex$", full.names = TRUE)
    kObs <- unlist(lapply(chr_files, .KObsFile))
    kObs <- kObs[!is.na(kObs)]
    # Per pipeline: drop kObs <= 1 (invariant) and any kObs == 0.
    kObs <- kObs[kObs >= 2L]
    results[[paste(tr, rp, sep = "/")]] <- kObs
  }
}

allK <- unlist(results)
cat(sprintf("Total variable chars across 260 datasets: %d\n", length(allK)))
cat("Pooled kObs distribution:\n")
print(table(allK))

# Per-dataset summary.
medians <- vapply(results, median, numeric(1))
maxes   <- vapply(results, max,    numeric(1))
cat(sprintf("\nPer-dataset median kObs: median=%d, range=[%d, %d]\n",
            as.integer(median(medians)), as.integer(min(medians)),
            as.integer(max(medians))))
cat(sprintf("Per-dataset max kObs: median=%d, range=[%d, %d]\n",
            as.integer(median(maxes)), as.integer(min(maxes)),
            as.integer(max(maxes))))

# Compare to empiricalNObs.
data("empiricalNObs", package = "MkPrime")
cat("\nempiricalNObs structure:\n")
str(empiricalNObs)

body <- empiricalNObs$body
ks_body <- seq_along(body) + 1L  # body[i] is mass at k=i+1
cat("empiricalNObs body (P(k')):\n")
print(setNames(round(body, 3), paste0("k'=", ks_body)))

# Data pmf on same support.
dataPmf <- table(factor(allK, levels = ks_body)) / length(allK)
cat("\nPooled data pmf at same support:\n")
print(round(dataPmf, 3))

# Score: where is mass overlap good vs bad?
cat("\nPer-k comparison (prior body P vs data pmf):\n")
df <- data.frame(
  k = ks_body,
  prior = round(as.numeric(body), 3),
  data  = round(as.numeric(dataPmf), 3)
)
df$prior_minus_data <- round(df$prior - df$data, 3)
print(df)

# Per-dataset KL(data || prior) on the body support (additive eps).
.KL <- function(p, q, eps = 1e-6) {
  p <- p + eps; q <- q + eps
  p <- p / sum(p); q <- q / sum(q)
  sum(p * log(p / q))
}
kls <- vapply(results, function(k) {
  if (!length(k)) return(NA_real_)
  pmf <- table(factor(k, levels = ks_body)) / length(k)
  .KL(pmf, body)
}, numeric(1))
cat(sprintf("\nKL(data || empiricalNObs) per dataset: median=%.3f, max=%.3f\n",
            median(kls, na.rm = TRUE), max(kls, na.rm = TRUE)))

saveRDS(list(per_dataset = results, pooled_kObs = allK, kl = kls,
             empirical_body = body, comparison = df),
        "data-raw/step7-prior-vs-data.rds")
cat("\nSaved data-raw/step7-prior-vs-data.rds\n")
