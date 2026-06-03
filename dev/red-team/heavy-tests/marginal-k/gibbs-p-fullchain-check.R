# DECISIVE end-to-end check: does gibbs_p_marginal (case 35), run through the full
# RunMkPrime pipeline on an n16_c48-like dataset (narrow near-boundary p, the
# production regime), lift marginal_k p-ESS toward sampled_k's ~2000 -- vs the
# mh_logit_p-only default's ~60? Also the advisor's end-to-end smoke for the
# opt-in flag (exercises moveTypeCodes / tuning / scalar-floor / pinned-weight
# plumbing that the unit tests bypass). Run in background; local free-topology
# marginal_k is slow (memory: long runs -> Hamilton), so N_ITER is modest.
suppressMessages(pkgload::load_all("C:/Users/pjjg18/GitHub/worktrees/mkp/marginal-k", quiet = TRUE))
library("ape"); library("TreeTools")

P_TRUE <- 0.25
.simJCchar <- function(tree, kTrue) {
  nTip <- length(tree$tip.label); states <- integer(2L*nTip-1L)
  rootIdx <- nTip + 1L; states[rootIdx] <- sample.int(kTrue, 1L) - 1L
  edges <- tree$edge; el <- tree$edge.length
  for (e in seq_len(nrow(edges))) {
    pa <- edges[e,1L]; ch <- edges[e,2L]; t <- el[e]
    pSame <- 1/kTrue + (1 - 1/kTrue) * exp(-kTrue * t / (kTrue - 1))
    states[ch] <- if (runif(1L) < pSame) states[pa]
                  else sample(setdiff(seq.int(0L, kTrue-1L), states[pa]), 1L)
  }
  states[seq_len(nTip)]
}
.canon <- function(v) { u <- sort(unique(v)); match(v, u) - 1L }
simulate_dataset <- function(nTip, nChar, seed) {
  set.seed(seed)
  tree <- rtree(nTip); tree$edge.length <- pmax(tree$edge.length, 0.02)
  tree$tip.label <- paste0("t", seq_len(nTip))
  cols <- vector("list", nChar); i <- 0L; g <- 0L
  while (i < nChar) {
    g <- g + 1L; if (g > 80L*nChar) stop("fill")
    kTrue <- 2L + rgeom(1L, P_TRUE); raw <- .simJCchar(tree, kTrue); can <- .canon(raw)
    if (length(unique(can)) < 2L) next
    i <- i + 1L; cols[[i]] <- can
  }
  mat <- matrix(unlist(cols), nrow = nTip, ncol = nChar, dimnames = list(tree$tip.label, NULL))
  list(tree = Preorder(tree), mkd = MkPrimeData(MatrixToPhyDat(mat)))
}
bm_ess <- function(x) {                       # non-overlapping batch-means ESS
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 20 || var(x) < 1e-14) return(0)
  a <- max(2L, floor(sqrt(n))); b <- floor(n / a)
  m <- colMeans(matrix(x[seq_len(a*b)], nrow = b))
  s2 <- b * var(m); if (s2 <= 0) return(NA_real_)
  n * var(x) / s2
}

N_ITER <- 12000L; thin_n <- max(1L, round(N_ITER / 3000L))
sim <- simulate_dataset(16L, 48L, 1648L)
# Four configs to locate the p-ESS ceiling and the role of move WEIGHT:
#  A mh-only        : marginal_k, gibbs off            (current default)
#  B gibbs-low      : marginal_k, gibbs on at weight 3 (current opt-in scheduling)
#  C gibbs-primary  : marginal_k, gibbs PINNED high + mh_logit_p low (proper use)
#  D sampled_k      : reference CEILING (conjugate Gibbs is the primary p-move)
run_cfg <- function(mode, gpm, mw = NULL) {
  model <- suppressMessages(MkPrimeModel(kPrimePrior = "geometric",
                                         likelihoodMode = mode, coding = "variable"))
  args <- list(nIter = N_ITER, thin = thin_n, minWarmup = 2000L, maxWarmup = 2000L,
               autoTune = FALSE, nRuns = 1L, nChains = 1L, gibbsPMarginal = gpm)
  if (!is.null(mw)) args$moveWeights <- mw
  mcmc <- suppressWarnings(do.call(MkPrimeMCMC, args))
  set.seed(99L)
  fit <- suppressMessages(suppressWarnings(
    RunMkPrime(sim$mkd, sim$tree, model = model, mcmc = mcmc, overwrite = TRUE)))
  s <- fit$samples
  pcol <- s[, "p"]; tcol <- if ("tree_length" %in% colnames(s)) s[, "tree_length"] else rep(NA, nrow(s))
  c(nsamp = nrow(s), ess_p = bm_ess(pcol), ess_tl = bm_ess(tcol), mean_p = mean(pcol), sd_p = sd(pcol))
}
A <- run_cfg("marginal_k", FALSE)              # mh_logit_p only (current default)
B <- run_cfg("marginal_k", TRUE)               # flag ON => gibbs_p_marginal PRIMARY
D <- run_cfg("sampled_k",  FALSE)              # conjugate-Gibbs ceiling reference
fmt <- function(tag, r) sprintf("%-34s nsamp=%.0f  p-ESS=%.0f  tl-ESS=%.0f  E[p]=%.3f sd=%.3f",
                                tag, r["nsamp"], r["ess_p"], r["ess_tl"], r["mean_p"], r["sd_p"])
out <- c(sprintf("n16_c48 free-topology, N_ITER=%d thin=%d  (p-posterior ~ narrow near-boundary)", N_ITER, thin_n),
         fmt("A marginal mh_logit_p only", A),
         fmt("B marginal gibbsPMarginal=TRUE", B),
         fmt("D sampled_k (ceiling ref)", D),
         sprintf("p-ESS: B/A = %.1fx (gibbs-primary vs RW) ; sampled ceiling D/A = %.1fx",
                 B["ess_p"]/max(A["ess_p"],1e-9), D["ess_p"]/max(A["ess_p"],1e-9)))
writeLines(out, "C:/Users/pjjg18/GitHub/worktrees/mkp/marginal-k/dev/red-team/heavy-tests/marginal-k/gibbs-p-fullchain-out.txt")
cat(out, sep = "\n")
