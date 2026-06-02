# Free-topology RB-equivalence (MARGINAL-K-FREEZE-003 validation).
#
# The Stage-2 RB identity (test-marginal-k-truncation.R) proves, on ONE fixed
# tree, that summing the sampled_k JOINT over k' reproduces the marginal_k joint:
#
#   lse_{k' in [kObs, K]} [ logPrior_S(k') + logLik_S(k') ] == logPrior_M + logLik_M
#
# The free-topology freeze fix changed WHICH topology the marginal evaluator
# reads (compute_per_kprime_log_lik / cpp_log_likelihood_marginal now honour the
# PASSED parent/child rather than state->parent/child). This test GENERALIZES the
# identity across many RANDOM topologies + varied p / branch lengths / rate
# heterogeneity, confirming the marginal evaluator is RB-correct on ARBITRARY
# trees -- the deterministic core of free-topology RB-equivalence, with no MCMC
# sampling noise. (The complementary "a topology MH proposal commits the correct
# marginal LL" facet is guarded behaviourally by test-marginal-k-free-topology.R:
# committed state->logLik == cold recompute after an accepted nni/spr/tbr/pspr.)
#
# One transformational character per tree keeps the k'-sum O(K); the
# multi-character closure is proven analytically (dev/red-team/proofs/
# marginal-k-sampled-rb-consistency.md) and tested empirically at n=2 in
# test-marginal-k-truncation.R. nCat = 4 with rls in {0, 0.6} exercises BOTH the
# homogeneous and the ACRV-heterogeneous pruning paths (the fix rerouted the
# parent/child argument in both).

library("ape")
library("TreeTools")

.rb_lse <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(-Inf)
  m <- max(x); m + log(sum(exp(x - m)))
}

# One evaluator closure per (mode, tree, data): build model + data once, then
# vary only the state (p / rate_log_sd / k') across calls.
.rb_eval_closure <- function(mode, trp, mkd, K) {
  mod  <- suppressMessages(MkPrimeModel(
    coding = "variable", nCat = 4L, kPrimePrior = "geometric",
    likelihoodMode = mode, priorVariant = "unconditional",
    kprimeTruncK = K, kprimeHyperA = 1, kprimeHyperB = 1, expSteps = 1.4))
  modf <- MkPrime:::.FinalizeModel(mod, trp, mkd)
  dp   <- MkPrime:::.InitMcmcData(mkd, modf)
  TL   <- sum(trp$edge.length); RBL <- trp$edge.length / TL
  function(mutate) {
    st <- MkPrime:::.InitState(trp, mkd, modf)
    st$tree_length <- TL; st$rel_br_lengths <- RBL
    st <- mutate(st)
    sp <- MkPrime:::.InitMcmcChain(st)
    fill_partition_cache(dp, sp)
    eval_log_prior_cpp(dp, sp) + eval_full_loglik_cpp(dp, sp)
  }
}

test_that("marginal-k RB identity holds across random topologies (FREEZE-003 validation)", {
  set.seed(909L)
  K <- 30L
  for (rep in seq_len(8L)) {
    ntip <- sample(6:11, 1L)
    tr <- ape::rtree(ntip, tip.label = paste0("t", seq_len(ntip)))
    tr$edge.length <- runif(nrow(tr$edge), 0.03, 0.5)   # varied (non-degenerate) branches
    trp <- TreeTools::Preorder(tr)
    repeat { v <- sample(0:1, ntip, replace = TRUE); if (length(unique(v)) == 2L) break }
    mkd  <- MkPrimeData(TreeTools::MatrixToPhyDat(
      matrix(v, ncol = 1, dimnames = list(tr$tip.label, "c1"))))
    kobs <- mkd$kObs[1]; ks <- kobs:K

    evalM <- .rb_eval_closure("marginal_k", trp, mkd, K)
    evalS <- .rb_eval_closure("sampled_k",  trp, mkd, K)

    for (p in c(0.1, 0.4)) for (rls in c(0, 0.6)) {
      M <- evalM(function(s) { s$p <- p; s$rate_log_sd <- rls; s })
      J <- vapply(ks, function(kp) evalS(function(s) {
        s$p <- p; s$rate_log_sd <- rls
        s$kPrime <- rep_len(as.integer(kp), length(s$kPrime)); s
      }), numeric(1))
      tag <- sprintf("rep%d ntip=%d p=%.1f rls=%.1f", rep, ntip, p, rls)
      # Non-vacuity: the joint genuinely varies with k' (k'-pinning is live).
      expect_gt(stats::sd(J), 1e-6, label = paste(tag, ": joint varies over k'"))
      # MAIN: summing the sampled_k joint over k' reproduces the marginal joint,
      # on THIS random topology.
      expect_equal(.rb_lse(J), M, tolerance = 1e-7,
                   label = paste(tag, ": lse(sampled joint over k') == marginal joint"))
    }
  }
})
