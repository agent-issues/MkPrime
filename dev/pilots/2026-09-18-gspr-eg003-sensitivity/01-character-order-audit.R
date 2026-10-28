# EG-003 character-order audit
# =============================
#
# Before #54, `data-raw/hamilton/run_one.R` built the character matrix from
#   sort(list.files(dataset_dir, "^chr[0-9]+\\.nex$"))
# which is LEXICAL order (chr1, chr10, chr11, ..., chr19, chr2, ...), while
# `ground_truth.csv` is in NUMERIC order. `analysis/eg001_upost_compare.R`
# paired the sampler's kPrime_1..kPrime_50 against gt by position, so every
# per-character statistic in EG-003 from those summaries was a permuted pairing.
# Summaries written after #54 record `char_idx` (the sampler's order);
# older ones and the EG `post_means.rds` do not, and are assumed lexical.
#
# "published" is the positional pairing; "corrected" joins ground truth on
# the sampler's char_idx. For summaries that carry char_idx in numeric order the
# two coincide.
#
# Falsifier that needs no modelling assumption: k' >= kObs is enforced by
# construction, so k_post - kObs < 0 is impossible for a correctly aligned
# pair. This script counts those violations under each pairing and recomputes
# the EG-003 statistics both ways.
#
# Run:  Rscript dev/pilots/2026-09-18-gspr-eg003-sensitivity/01-character-order-audit.R
# Needs only the committed summary/ .rds files plus the ground-truth tree
# directory; no package build, no MCMC.

source("dev/pilots/2026-09-18-gspr-eg003-sensitivity/align.R")

SUMMARY_DIR <- Sys.getenv("AUDIT_SUMMARY_DIR", "dev/pilots/2026-05-12-prior-validation/summary")
EG_POST     <- Sys.getenv("AUDIT_EG_POST", "C:/Users/pjjg18/GitHub/mkprime/report-data/mkp-eg-260/post_means.rds")
GT_ROOT     <- Sys.getenv("AUDIT_GT_ROOT", "C:/Users/pjjg18/GitHub/mkprime/tree-inference")

# char_idx of the sampler's characters for a pre-#54 summary: lexical file order.
LexCharIdx <- function(tree, rep) {
  f <- list.files(file.path(GT_ROOT, sprintf("tree_%02d/rep_%02d", tree, rep)),
                  "^chr[0-9]+\\.nex$")
  CharIdxFromFiles(sort(f))
}

GroundTruth <- function(task, charIdx) {
  tree <- as.integer(sub("^t([0-9]+).*", "\\1", task))
  rep  <- as.integer(sub(".*_r([0-9]+)$", "\\1", task))
  gt <- read.csv(file.path(GT_ROOT,
                           sprintf("tree_%02d/rep_%02d/ground_truth.csv", tree, rep)))
  if (is.null(charIdx)) charIdx <- LexCharIdx(tree, rep)
  list(published = gt,
       corrected = AlignGroundTruth(gt, charIdx))
}

Stats <- function(kPost, gt) {
  u <- kPost - gt$kObs
  c(mean = mean(u), median = median(u), pLt05 = mean(u < 0.5),
    rho = suppressWarnings(cor(u, gt$u_true, method = "spearman")),
    violations = sum(u < -1e-9))
}

Collect <- function(entries) {
  do.call(rbind, lapply(entries, function(e) {
    g <- GroundTruth(e$task, e$charIdx)
    if (length(e$kPost) != nrow(g$corrected)) {
      warning("Skipping ", e$task, ": ", length(e$kPost), " k' values for ",
              nrow(g$corrected), " characters", call. = FALSE)
      return(NULL)
    }
    # Invariant characters are absent from the sampler, so the positional
    # pairing is defined only when lengths agree.
    pub <- if (length(e$kPost) == nrow(g$published)) {
      Stats(e$kPost, g$published)
    } else {
      setNames(rep(NA_real_, 5L), c("mean", "median", "pLt05", "rho", "violations"))
    }
    cor_ <- Stats(e$kPost, g$corrected)
    data.frame(task = e$task, prior = e$prior, nChar = length(e$kPost),
               t(setNames(pub,  paste0("pub_", names(pub)))),
               t(setNames(cor_, paste0("cor_", names(cor_)))))
  }))
}

geo <- lapply(list.files(SUMMARY_DIR, "^mkp_geo_.*\\.rds$", full.names = TRUE),
              function(f) {
                x <- readRDS(f)
                list(task = x$tag, prior = "geo", kPost = unname(x$kp_means),
                     charIdx = x$char_idx)
              })
egRaw <- readRDS(EG_POST)
eg <- Filter(function(e) grepl("_r01$", e$task), lapply(seq_along(egRaw), function(i) {
  r <- egRaw[[i]]
  list(task = if (!is.null(r$task)) r$task else names(egRaw)[i],
       prior = "EG", kPost = unname(r$k_post))
}))

res <- rbind(Collect(geo), Collect(eg))

for (p in unique(res$prior)) {
  d <- res[res$prior == p, ]
  n <- nrow(d)
  nPub <- sum(!is.na(d$pub_mean))
  cat("=== prior:", p, "(", n, "tasks ) ===\n")
  for (nm in c("mean", "median", "pLt05", "rho")) {
    a <- d[[paste0("pub_", nm)]]
    b <- d[[paste0("cor_", nm)]]
    cat(sprintf(
      "  %-7s published %+.4f +/- %.4f | corrected %+.4f +/- %.4f | paired diff %+.4f +/- %.4f\n",
      nm, mean(a, na.rm = TRUE), sd(a, na.rm = TRUE) / sqrt(nPub),
      mean(b), sd(b) / sqrt(n),
      mean(b - a, na.rm = TRUE), sd(b - a, na.rm = TRUE) / sqrt(nPub)))
  }
  cat(sprintf("  k' < kObs violations: published %d | corrected %d (of %d char-tasks; published pairing undefined for %d of %d tasks)\n\n",
              sum(d$pub_violations, na.rm = TRUE), sum(d$cor_violations), sum(d$nChar),
              n - nPub, n))
}

uTrue <- unlist(lapply(1:26, function(t) {
  read.csv(file.path(GT_ROOT, sprintf("tree_%02d/rep_01/ground_truth.csv", t)))$u_true
}))
cat(sprintf("reference scale -- u_true: mean %.3f, median %.0f, P(u < 0.5) %.3f (n = %d)\n",
            mean(uTrue), median(uTrue), mean(uTrue < 0.5), length(uTrue)))
