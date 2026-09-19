# Mechanistic bound: how sensitive is the per-character u marginal to the
# branch-length / topology configuration that gibbs_spr distorted?
#
# The defect pulled adjacent edge fractions toward equality (tau = 1/2 splits)
# and biased which topologies were accepted. Here we compute the EXACT
# conditional posterior mean of u_i given (T, b, p, rate_log_sd) -- no MCMC --
# and compare it across:
#   A. posterior-sampled (T, b)                       [what a correct sampler sees]
#   B. the same topology with ALL edge fractions equalised  [the defect's limit]
#   C. random NNI-perturbed topologies                 [topology sensitivity]
#
# Per-character log-likelihoods are recovered by additivity: the total logL is a
# sum over characters, so bumping k'_i alone isolates character i.
setTimeLimit(elapsed = 3000, transient = FALSE)
source(file.path(Sys.getenv("GSPR_DIR", "."), "common.R"))
SP <- Sys.getenv("GSPR_DIR", ".")

UMAX <- 6L   # truncation for the conditional sum; geometric tail beyond is tiny

CondU <- function(tree, mkd, kObs, p, rate_log_sd) {
  n <- length(kObs)
  base <- MkpLogLikelihood(tree, mkd, kPrime = kObs, rate_log_sd = rate_log_sd,
                           coding = "variable")
  dll <- matrix(0, n, UMAX + 1L)   # dll[i, u+1] = logL_i(kObs_i + u) - logL_i(kObs_i)
  for (u in seq_len(UMAX)) {
    for (i in seq_len(n)) {
      kp <- kObs; kp[i] <- kObs[i] + u
      dll[i, u + 1L] <- MkpLogLikelihood(tree, mkd, kPrime = kp,
                                         rate_log_sd = rate_log_sd,
                                         coding = "variable") - base
    }
  }
  lpri <- log(p) + (0:UMAX) * log1p(-p)          # Geometric(p) on u
  lw   <- sweep(dll, 2L, lpri, "+")
  w    <- exp(lw - apply(lw, 1L, max))
  w    <- w / rowSums(w)
  as.numeric(w %*% (0:UMAX))
}

Equalise <- function(tree) {                      # the tau = 1/2 limit
  tree$edge.length <- rep(sum(tree$edge.length) / nrow(tree$edge), nrow(tree$edge))
  tree
}

out <- list()
for (ti in 1:4) {
  r   <- LoadRep(ti, 1L)
  mkd <- MkPrimeData(r$pd)
  kObs <- as.integer(mkd$kObs)
  cell <- file.path(SP, "grid", sprintf("t%02d_on_s001.rds", ti))
  if (!file.exists(cell)) { cat("no cell for tree", ti, "- skipping\n"); next }
  x <- readRDS(cell)
  p <- x$p_mean; rls <- 0.5                      # p from the chain; rls at a typical value
  # posterior point tree: use the MCMC tree file is not stored, so rebuild a
  # plausible posterior tree by a short chain -- instead use the start tree with
  # the chain's posterior mean tree length and Dirichlet-drawn fractions.
  st <- TreeSearch::AdditionTree(r$pd)
  nE <- nrow(st$edge)
  set.seed(1000 + ti)
  fr <- rgamma(nE, 1); fr <- fr / sum(fr)
  tA <- st; tA$edge.length <- fr * x$tree_length              # A: uneven fractions
  tB <- Equalise(tA)                                          # B: fully equalised
  tC <- TreeTools::Preorder(phangorn::rNNI(st, moves = 5L))
  tC$edge.length <- fr * x$tree_length                        # C: perturbed topology

  uA <- CondU(tA, mkd, kObs, p, rls)
  uB <- CondU(tB, mkd, kObs, p, rls)
  uC <- CondU(tC, mkd, kObs, p, rls)
  ut <- r$gt[as.integer(sub("^chr([0-9]+)\\.nex$", "\\1",
        sort(list.files(file.path(GT_ROOT, sprintf("tree_%02d/rep_01", ti)),
                        "^chr[0-9]+\\.nex$")))), ]$u_true

  cat(sprintf("\ntree %02d  (p = %.3f, TL = %.3f)\n", ti, p, x$tree_length))
  cat(sprintf("  A uneven fractions : mean u = %+.4f  rho = %+.3f\n", mean(uA), cor(uA, ut, method="spearman")))
  cat(sprintf("  B equalised (tau=1/2 limit): mean u = %+.4f  rho = %+.3f   | B-A: mean %+.4f, max|du_i| %.4f\n",
              mean(uB), cor(uB, ut, method="spearman"), mean(uB - uA), max(abs(uB - uA))))
  cat(sprintf("  C  5 random NNIs   : mean u = %+.4f  rho = %+.3f   | C-A: mean %+.4f, max|du_i| %.4f\n",
              mean(uC), cor(uC, ut, method="spearman"), mean(uC - uA), max(abs(uC - uA))))
  out[[as.character(ti)]] <- list(uA = uA, uB = uB, uC = uC, u_true = ut, p = p)
}
saveRDS(out, file.path(SP, "conditional-u.rds"))
