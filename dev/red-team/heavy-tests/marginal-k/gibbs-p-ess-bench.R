# Efficiency of gibbs_p_marginal (case 35) vs mh_logit_p (case 30), in BOTH
# p-regimes. The overlap p-posterior is NARROW + near the p->1 boundary
# (p~0.976, sd~0.02) because n16_c48 has ~48 chars mostly at u=0 => Beta(~49,1).
# A logit random-walk struggles near a boundary; a conjugate-style Gibbs draw
# does not. So this bench compares p-ESS at (i) an EASY broad mid-range regime
# (few chars) and (ii) a HARD narrow near-boundary regime (many chars), to see
# whether the Gibbs-p wins where it actually matters. Run after the regression
# (no concurrent load_all).
suppressMessages(pkgload::load_all("C:/Users/pjjg18/GitHub/worktrees/mkp/marginal-k", quiet = TRUE))

ess <- function(x) {                       # AR-spectral ESS (coda-free)
  x <- x[is.finite(x)]; n <- length(x)
  if (stats::var(x) < 1e-14) return(0)
  s <- spec.pgram(x, plot = FALSE, taper = 0, fast = FALSE, detrend = TRUE)
  n * stats::var(x) / s$spec[1]
}
make_tree <- function() ape::read.tree(text = paste0(
  "(((t1:0.05,t2:0.07):0.04,(t3:0.06,t4:0.05):0.03):0.05,",
  "((t5:0.04,t6:0.06):0.05,(t7:0.05,t8:0.04):0.06):0.04);"))
# nChar variable binary characters on 8 tips (deterministic given seed). Most
# random binary patterns favour u=0 (no hidden states), so more chars => higher,
# narrower p (Beta(a+nChar, b+Sum u) with Sum u small).
make_mkd <- function(nChar, seed) {
  set.seed(seed); tips <- paste0("t", 1:8)
  reps <- 0L
  cols <- vector("list", nChar)
  i <- 0L
  while (i < nChar) {
    v <- sample(0:1, 8, replace = TRUE)
    if (length(unique(v)) >= 2) { i <- i + 1L; cols[[i]] <- v }  # variable coding
  }
  m <- matrix(unlist(cols), nrow = 8, dimnames = list(tips, NULL))
  TreeTools::MatrixToPhyDat(m)
}
build <- function(nChar, seed) {
  tree  <- TreeTools::Preorder(make_tree()); mkd <- MkPrimeData(make_mkd(nChar, seed))
  model <- MkPrime:::.FinalizeModel(
    MkPrimeModel(kPrimePrior = "geometric", likelihoodMode = "marginal_k",
                 priorVariant = "conditional", coding = "variable"), tree, mkd)
  s0 <- MkPrime:::.InitState(tree, mkd, model)
  s0$p <- 0.7; s0$tree_length <- sum(tree$edge.length)
  s0$rel_br_lengths <- tree$edge.length / s0$tree_length
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model); statePtr <- MkPrime:::.InitMcmcChain(s0)
  fill_partition_cache(dataPtr, statePtr); invisible(eval_full_loglik_cpp(dataPtr, statePtr))
  list(dataPtr = dataPtr, statePtr = statePtr)
}
run <- function(fx, mt, nIter, st = 1.0) {
  ps <- numeric(nIter)
  for (i in seq_len(nIter)) {
    do_move_cpp(fx$dataPtr, fx$statePtr, moveType = mt, charIdx = 0L,
                scaleTuning = st, betaSimplexTuning = 0.5, intWalkWindow = 1L, beta = 1.0)
    ps[i] <- get_mcmc_state(fx$statePtr)$p
  }
  ps[-seq_len(nIter %/% 5L)]
}

N <- 40000L
out <- c("gibbs_p_marginal (case 35) vs mh_logit_p (case 30) -- p-ESS over 40000 iters",
         "regime varied by nChar (more chars => narrower, higher p, like the n16_c48 overlap)", "")
for (nc in c(4L, 24L, 48L)) {
  set.seed(7L)
  p35 <- run(build(nc, 101L), 35L, N)
  p30 <- run(build(nc, 101L), 30L, N, st = 1.0)
  out <- c(out, sprintf(
    "[nChar=%2d]  E[p]=%.3f sd=%.3f | ESS35=%.0f ESS30=%.0f ratio=%.2fx",
    nc, mean(p35), sd(p35), ess(p35), ess(p30), ess(p35) / max(ess(p30), 1e-9)))
}
writeLines(out, "C:/Users/pjjg18/GitHub/worktrees/mkp/marginal-k/dev/red-team/heavy-tests/marginal-k/gibbs-p-ess-bench-out.txt")
cat(out, sep = "\n")
