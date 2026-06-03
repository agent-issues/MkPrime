# Efficiency payoff for gibbs_p_marginal (case 35) vs mh_logit_p (case 30).
# Fires ONLY the p-move on a fixed-tree marginal_k fixture so both chains explore
# the same pi(p | theta, tree); reports p-ESS (per the same iteration budget).
# The point of the move is MIXING (independent draws vs random walk), so the
# ESS ratio is the headline. Run AFTER the regression (no concurrent load_all).
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
make_mkd <- function() {
  set.seed(42L); tips <- paste0("t", 1:8)
  m <- matrix(c(0,0,1,1,0,1,0,1, 0,1,1,0,1,0,0,1, 0,1,0,0,1,1,1,0, 1,1,0,1,0,0,1,0),
              nrow = 8, ncol = 4, dimnames = list(tips, NULL))
  TreeTools::MatrixToPhyDat(m)
}
build <- function(p, variant) {
  tree  <- TreeTools::Preorder(make_tree()); mkd <- MkPrimeData(make_mkd())
  model <- MkPrime:::.FinalizeModel(
    MkPrimeModel(kPrimePrior = "geometric", likelihoodMode = "marginal_k",
                 priorVariant = variant, coding = "variable"), tree, mkd)
  s0 <- MkPrime:::.InitState(tree, mkd, model)
  s0$p <- p; s0$tree_length <- sum(tree$edge.length)
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
res <- list()
for (variant in c("conditional", "unconditional")) {
  set.seed(7L)
  p35 <- run(build(0.5, variant), 35L, N)
  p30 <- run(build(0.5, variant), 30L, N, st = 1.0)
  res[[variant]] <- c(
    ess35 = ess(p35), ess30 = ess(p30), ratio = ess(p35) / max(ess(p30), 1e-9),
    mean35 = mean(p35), mean30 = mean(p30), sd35 = sd(p35), sd30 = sd(p30))
}
out <- capture.output({
  cat("gibbs_p_marginal (case 35) vs mh_logit_p (case 30) -- p-ESS over", N, "iters\n\n")
  for (v in names(res)) {
    r <- res[[v]]
    cat(sprintf("[%s]  ESS35=%.0f  ESS30=%.0f  ratio=%.1fx | E[p] 35=%.4f 30=%.4f | sd 35=%.4f 30=%.4f\n",
                v, r["ess35"], r["ess30"], r["ratio"], r["mean35"], r["mean30"], r["sd35"], r["sd30"]))
  }
})
writeLines(out, "C:/Users/pjjg18/GitHub/worktrees/mkp/marginal-k/dev/red-team/heavy-tests/marginal-k/gibbs-p-ess-bench-out.txt")
cat(out, sep = "\n")
