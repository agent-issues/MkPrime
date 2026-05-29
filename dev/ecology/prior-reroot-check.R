# Advisor's "T3 for the prior": does LogPrior depend on root PLACEMENT?
# Reading LogPrior (MkPrimeModel.R:579-708): the branch-shape term is a pure
# constant lfactorial(nEdge-1) (Dirichlet(1,...,1)); tree-length is Gamma on the
# preserved TOTAL; all else is parameters. So LogPrior should depend on branch
# lengths ONLY through (total, nEdge) -> invariant to distribution -> invariant
# to root placement (which preserves total and nEdge). Confirm empirically.
suppressMessages({devtools::load_all(".", quiet = TRUE); library(ape); library(TreeTools)})

set.seed(7L)
nTip <- 9L
tips <- paste0("t", seq_len(nTip))
base <- TreeTools::Preorder(ape::rtree(nTip, tip.label = tips))   # positive edges

m <- matrix(0L, nTip, 30L, dimnames = list(tips, NULL))
for (cc in seq_len(30L)) { repeat { col <- sample(0:1, nTip, TRUE)
  if (length(unique(col)) > 1L) break }; m[, cc] <- col }
mkd   <- MkPrimeData(TreeTools::MatrixToPhyDat(m), neomorphic = seq_len(30L))
model <- MkPrimeModel()
# MkPrimeModel() leaves data-dependent rates NULL (auto-set in RunMkPrime); fill
# them so LogPrior is finite. These values are irrelevant to root-invariance.
if (is.null(model$treeLengthRate))  model$treeLengthRate  <- 0.5
if (is.null(model$rateLogSdRate))   model$rateLogSdRate   <- 1
if (is.null(model$rateLossMeanlog)) model$rateLossMeanlog <- log(2)
if (is.null(model$rateLossSdlog))   model$rateLossSdlog   <- 0.73
if (is.null(model$rateNeoMeanlog))  model$rateNeoMeanlog  <- 0
if (is.null(model$rateNeoSdlog))    model$rateNeoSdlog    <- 1

mkState <- function(tl, rel) list(
  tree_length = tl, rel_br_lengths = rel,
  rate_loss = 1.3, rate_neo = 1.0, rate_log_sd = 0.5,
  kPrime = mkd$kObs, p = 0.5)

TL  <- sum(base$edge.length)
rel0 <- base$edge.length / TL
lp0  <- LogPrior(mkState(TL, rel0), model, mkd)
cat("=== PRIOR ROOT-INVARIANCE CHECK ===\n")
cat(sprintf("base: logPrior = %.12f  (nEdge=%d, TL=%.6f)\n", lp0, length(rel0), TL))

# (a) Same total, DIFFERENT positive branch distribution (what re-rooting does
#     to the prior inputs: changes the proportion vector, preserves total+count).
diffsA <- vapply(1:5, function(s) {
  set.seed(100 + s); p <- runif(length(rel0)); rel <- p / sum(p)
  abs(LogPrior(mkState(TL, rel), model, mkd) - lp0)
}, numeric(1))

# (b) Actual re-rooting, zero basal edge repaired to positive, renormalized to
#     the SAME total (a valid positive rooted representation of same unrooted top).
diffsB <- vapply(tips, function(og) {
  rr <- ape::root(base, outgroup = og, resolve.root = TRUE)
  el <- rr$edge.length; el[el <= 0] <- 1e-6
  abs(LogPrior(mkState(TL, el / sum(el)), model, mkd) - lp0)
}, numeric(1))

# control: a DIFFERENT total must change logPrior (test is sensitive, not const)
diffC <- abs(LogPrior(mkState(TL * 1.7, rel0), model, mkd) - lp0)

stopifnot(length(diffsA) == 5L, length(diffsB) == nTip)  # guard vacuous pass
cat("\n(a) redistribute diffs (n=5):\n"); print(signif(diffsA, 3))
cat("\n(b) reroot diffs (n=", nTip, "):\n", sep = ""); print(signif(diffsB, 3))
cat(sprintf("\ncontrol (total x1.7) diff = %.4f  (MUST be > 0)\n", diffC))
cat(sprintf("MAX |diff| invariance tests = %.3e ;  ran %d comparisons\n",
            max(diffsA, diffsB), length(diffsA) + length(diffsB)))
cat(if (max(diffsA, diffsB) < 1e-9 && diffC > 1e-6)
      "VERDICT: prior ROOT-INVARIANT (and test is sensitive) -> infer-root target-correct on prior\n"
    else "VERDICT: inconclusive/root-informative -> investigate\n")
