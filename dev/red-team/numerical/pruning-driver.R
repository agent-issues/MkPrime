## pruning-driver.R — Lane N2 stress driver
##
## Audits numerical robustness of the JC Felsenstein pruning kernel and the
## Lewis ascertainment subtraction `log(1 - P_const)` (and the singleton
## extension `log(1 - P_const - P_singleton)`).
##
## Reference: an independent pure-R Felsenstein pruning that uses the JC
## analytic transition formula on the same tree.  Both implementations live
## in double precision, so the reference does *not* improve precision; it
## only catches algorithmic / index bugs.  The numerical questions
## (underflow, log-cancellation) are addressed by:
##
##   (a) direct construction of stress patterns (caterpillar + conflicting
##       tip states) that drive internal partial CL toward 0;
##   (b) extreme-p experiments for `log(1 - p)` using long-string
##       arithmetic to compute the exact `log1p(-p)` and comparing it to
##       `log(1 - p)`.
##
## The pruning kernels in src/likelihood.cpp, src/mcmc_likelihood.cpp, and
## src/ascertainment.cpp do **not** rescale partial likelihoods.  This
## driver quantifies the regime where that becomes silently wrong.
##
## Usage:
##   Rscript dev/red-team/numerical/pruning-driver.R           # full
##   Rscript dev/red-team/numerical/pruning-driver.R --quick   # <=2 min

quick <- any(commandArgs(trailingOnly = TRUE) == "--quick")

suppressPackageStartupMessages({
  library(MkPrime)   # installed version; provides pruning_jc, *_site_prob_*
  library(ape)
})

## Pruning kernel functions are non-exported; pull from the namespace.
mk_pruning_jc <- get("pruning_jc", envir = asNamespace("MkPrime"))
mk_const_jc   <- get("constant_site_prob_jc", envir = asNamespace("MkPrime"))
mk_sing_jc    <- get("singleton_site_prob_jc", envir = asNamespace("MkPrime"))

out_dir <- file.path("dev", "red-team", "numerical", "pruning-results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- helpers ---------------------------------------------------------------

## Build a caterpillar (ladder) tree on n tips.  Returns ape::phylo, with
## all edges set to `bl`.  Caterpillar => max depth = n - 1.
make_caterpillar <- function(n, bl = 0.1) {
  if (n < 2) stop("n >= 2")
  tip <- paste0("t", seq_len(n))
  ## Newick: (t1,(t2,(t3,(...(t_{n-1},t_n)...))))
  build <- function(i) {
    if (i == n - 1) sprintf("(%s:%g,%s:%g)", tip[i], bl, tip[n], bl)
    else sprintf("(%s:%g,%s:%g)", tip[i], bl, build(i + 1), bl)
  }
  newick <- paste0(build(1), ";")
  tr <- read.tree(text = newick)
  ## Set ALL edges (including internals) to bl
  tr$edge.length <- rep(bl, nrow(tr$edge))
  tr
}

## Build a balanced binary tree on n = 2^depth tips with all edges bl.
make_balanced <- function(depth, bl = 0.1) {
  n <- 2L^depth
  tr <- ape::stree(n, type = "balanced")
  tr$edge.length <- rep(bl, nrow(tr$edge))
  tr
}

## Pure-R reference Felsenstein on JC(k).  Returns site log-likelihoods.
## tip_states: nTip x nChar integer matrix (0-indexed; -1 missing).
ref_pruning_jc <- function(tree, tip_states, k) {
  ## Canonicalise edges in preorder by depth-first traversal.
  ## ape's edge ordering is already root-to-tip if we postorder by edge index
  ## from the end; we mimic the C++ kernel exactly.
  parent <- tree$edge[, 1]
  child  <- tree$edge[, 2]
  el     <- tree$edge.length
  nEdge  <- length(parent)
  nTip   <- length(tree$tip.label)
  nChar  <- ncol(tip_states)
  maxNode <- max(parent, child)

  ## CL[node][char, state]
  CL <- vector("list", maxNode + 1L)
  init <- logical(maxNode + 1L)
  for (n in 1:(maxNode + 1L)) CL[[n]] <- matrix(0, nChar, k)

  for (tip in seq_len(nTip)) {
    for (c in seq_len(nChar)) {
      s <- tip_states[tip, c]
      if (s < 0) CL[[tip]][c, ] <- 1.0
      else       CL[[tip]][c, s + 1L] <- 1.0
    }
    init[tip] <- TRUE
  }

  inv_k <- 1.0 / k
  for (e in nEdge:1) {
    p <- parent[e]; ch <- child[e]; t <- el[e]
    et <- exp(-k * t / (k - 1))
    p_same <- inv_k + (1 - inv_k) * et
    p_diff <- inv_k - inv_k * et
    ## per-char update
    for (c in seq_len(nChar)) {
      clCh <- CL[[ch]][c, ]
      s <- sum(clCh)
      newCL <- p_diff * s + (p_same - p_diff) * clCh
      if (!init[p]) CL[[p]][c, ] <- newCL
      else          CL[[p]][c, ] <- CL[[p]][c, ] * newCL
    }
    init[p] <- TRUE
  }

  root <- nTip + 1L
  rf <- rep(1/k, k)
  ll <- 0
  for (c in seq_len(nChar)) {
    sl <- sum(rf * CL[[root]][c, ])
    if (sl <= 0) return(-Inf)
    ll <- ll + log(sl)
  }
  ll
}

## Pure-R reference for constant-site probability (sum over k constant
## patterns) under JC(k) with optional rate categories.
ref_const_prob_jc <- function(tree, k, rates = 1) {
  ## P(all tips state s) for one s, then × k.
  nTip <- length(tree$tip.label)
  tip_states <- matrix(0L, nrow = nTip, ncol = 1)  ## all tips = state 0
  cp <- 0
  for (r in rates) {
    tt <- tree
    tt$edge.length <- tree$edge.length * r
    cp <- cp + exp(ref_pruning_jc(tt, tip_states, k))
  }
  cp / length(rates) * k
}

## Compute exact log1p(-p) and naive log(1-p) and return both.
## Uses base double precision; the *measured* error of log(1 - p) is
## (log1p(-p) - log(1 - p)) / log1p(-p).
log1m_compare <- function(p) {
  naive <- suppressWarnings(log(1 - p))
  safe  <- suppressWarnings(log1p(-p))
  rel <- if (!is.finite(safe) || safe == 0) NA_real_
         else abs(naive - safe) / abs(safe)
  list(naive = naive, safe = safe,
       abs_err = abs(naive - safe), rel_err = rel)
}

write_csv <- function(df, fname) {
  path <- file.path(out_dir, fname)
  write.csv(df, path, row.names = FALSE)
  cat("  wrote", path, "\n")
}

## ---- experiment 1: deep tree underflow -------------------------------------

cat("=== Experiment 1: caterpillar underflow on conflicting characters ===\n")

set.seed(42)

cater_depths <- if (quick) c(8L, 16L) else c(8L, 16L, 32L, 64L, 128L, 256L, 512L)
bl_grid <- c(1e-12, 1e-6, 1e-3, 0.1, 1.0, 10.0, 100.0)

cater_res <- data.frame()

for (n in cater_depths) {
  tr <- make_caterpillar(n, bl = 0.1)
  ## Conflicting pattern: alternating states along the ladder
  k <- 4L
  ts <- matrix((seq_len(n) %% k), nrow = n, ncol = 1)  ## 0,1,2,3,0,1,...
  ## To probe more, also include an "all distinct" pattern when n <= k
  for (bl in bl_grid) {
    tr$edge.length <- rep(bl, nrow(tr$edge))
    ll_kernel <- mk_pruning_jc(
      as.integer(tr$edge[, 1]),
      as.integer(tr$edge[, 2]),
      as.numeric(tr$edge.length),
      ts,
      kStates = k,
      root_freqs = rep(1/k, k)
    )
    ll_ref <- ref_pruning_jc(tr, ts, k)
    rel <- if (is.finite(ll_ref) && ll_ref != 0)
      abs(ll_kernel - ll_ref) / abs(ll_ref) else NA_real_
    cater_res <- rbind(cater_res, data.frame(
      n_tip = n, branch_len = bl,
      ll_kernel = ll_kernel, ll_ref = ll_ref,
      abs_err = abs(ll_kernel - ll_ref),
      rel_err = rel,
      kernel_underflow = !is.finite(ll_kernel) || ll_kernel == -Inf
    ))
  }
}

print(cater_res)
write_csv(cater_res, "01-caterpillar-underflow.csv")

## ---- experiment 2: deep balanced tree, equilibrium branches ----------------

cat("\n=== Experiment 2: balanced tree, branches near equilibrium ===\n")

bal_depths <- if (quick) c(4L, 6L) else c(4L, 6L, 8L, 10L)
## "near equilibrium" = bl ~ (k-1)/k * a few units, so exp(-k bl /(k-1)) is
## small.  Probe bl = 1, 10, 100.

bal_res <- data.frame()
for (d in bal_depths) {
  tr <- make_balanced(d, bl = 0.1)
  n <- length(tr$tip.label)
  k <- 4L
  ## "All tips in distinct states" pattern (rotates through k)
  ts <- matrix((seq_len(n) - 1L) %% k, nrow = n, ncol = 1)
  for (bl in c(0.5, 1, 5, 20)) {
    tr$edge.length <- rep(bl, nrow(tr$edge))
    ll_kernel <- mk_pruning_jc(
      as.integer(tr$edge[, 1]),
      as.integer(tr$edge[, 2]),
      as.numeric(tr$edge.length),
      ts,
      kStates = k,
      root_freqs = rep(1/k, k)
    )
    ll_ref <- ref_pruning_jc(tr, ts, k)
    rel <- if (is.finite(ll_ref) && ll_ref != 0)
      abs(ll_kernel - ll_ref) / abs(ll_ref) else NA_real_
    bal_res <- rbind(bal_res, data.frame(
      depth = d, n_tip = n, branch_len = bl,
      ll_kernel = ll_kernel, ll_ref = ll_ref,
      abs_err = abs(ll_kernel - ll_ref), rel_err = rel
    ))
  }
}

print(bal_res)
write_csv(bal_res, "02-balanced-equilibrium.csv")

## ---- experiment 3: tiny tree, near-zero branches => p_const -> 1 -----------

cat("\n=== Experiment 3: tiny tree, near-zero branches (p_const -> 1) ===\n")

## 3-tip tree with very short branches; check P(constant) and log(1 - p).

asc_res <- data.frame()
n_tips <- 3L
tr3 <- read.tree(text = "((t1:BL,t2:BL):BL,t3:BL);")
for (k in c(2L, 4L, 8L, 12L)) {
  for (bl in c(1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 0.1, 1.0, 10.0)) {
    tr3$edge.length <- rep(bl, nrow(tr3$edge))
    p <- mk_const_jc(
      as.integer(tr3$edge[, 1]),
      as.integer(tr3$edge[, 2]),
      as.numeric(tr3$edge.length),
      nTip = n_tips, kStates = k,
      root_freqs = rep(1/k, k),
      rate_multipliers = 1
    )
    p_ref <- ref_const_prob_jc(tr3, k)
    cmp <- log1m_compare(p)
    asc_res <- rbind(asc_res, data.frame(
      n_tip = n_tips, k_states = k, branch_len = bl,
      p_const = p, p_const_ref = p_ref,
      one_minus_p = 1 - p,
      log_naive = cmp$naive,
      log_safe  = cmp$safe,
      abs_err   = cmp$abs_err,
      rel_err   = cmp$rel_err
    ))
  }
}

print(asc_res)
write_csv(asc_res, "03-ascertainment-cancellation.csv")

## ---- experiment 4: tiny tree, very long branches => p_const -> 1/k ---------

cat("\n=== Experiment 4: tiny tree, long branches (p_const -> 1/k) ===\n")

long_res <- data.frame()
for (k in c(2L, 4L, 8L)) {
  for (bl in c(1, 10, 100, 1000)) {
    tr3$edge.length <- rep(bl, nrow(tr3$edge))
    p <- mk_const_jc(
      as.integer(tr3$edge[, 1]),
      as.integer(tr3$edge[, 2]),
      as.numeric(tr3$edge.length),
      nTip = n_tips, kStates = k,
      root_freqs = rep(1/k, k),
      rate_multipliers = 1
    )
    cmp <- log1m_compare(p)
    long_res <- rbind(long_res, data.frame(
      n_tip = n_tips, k_states = k, branch_len = bl,
      p_const = p, one_minus_p = 1 - p,
      log_naive = cmp$naive, log_safe = cmp$safe,
      abs_err = cmp$abs_err, rel_err = cmp$rel_err
    ))
  }
}

print(long_res)
write_csv(long_res, "04-long-branch-1overk.csv")

## ---- experiment 5: singleton subtraction with kObs in {2, k} ---------------

cat("\n=== Experiment 5: log(1 - p_const - p_singleton) ===\n")

## Use a 5-tip tree so P(constant) + P(singleton) < 1 even at k = 2 (a
## 3-tip binary site is necessarily either constant or singleton: sum to 1).
tr_sing <- read.tree(text = "(((t1:BL,t2:BL):BL,(t3:BL,t4:BL):BL):BL,t5:BL);")
n_tips_sing <- 5L
sing_res <- data.frame()
for (k in c(2L, 4L, 8L)) {
  for (bl in c(1e-12, 1e-9, 1e-6, 1e-3, 0.1, 1.0)) {
    tr_sing$edge.length <- rep(bl, nrow(tr_sing$edge))
    p_c <- mk_const_jc(
      as.integer(tr_sing$edge[, 1]),
      as.integer(tr_sing$edge[, 2]),
      as.numeric(tr_sing$edge.length),
      nTip = n_tips_sing, kStates = k,
      root_freqs = rep(1/k, k), rate_multipliers = 1)
    p_s <- mk_sing_jc(
      as.integer(tr_sing$edge[, 1]),
      as.integer(tr_sing$edge[, 2]),
      as.numeric(tr_sing$edge.length),
      nTip = n_tips_sing, kStates = k,
      root_freqs = rep(1/k, k), rate_multipliers = 1)
    p_tot <- p_c + p_s
    cmp <- log1m_compare(p_tot)
    sing_res <- rbind(sing_res, data.frame(
      n_tip = n_tips_sing, k_states = k, branch_len = bl,
      p_const = p_c, p_singleton = p_s, p_total = p_tot,
      one_minus_p = 1 - p_tot,
      log_naive = cmp$naive, log_safe = cmp$safe,
      abs_err = cmp$abs_err, rel_err = cmp$rel_err
    ))
  }
}

print(sing_res)
write_csv(sing_res, "05-singleton-subtraction.csv")

## ---- experiment 6: synthetic worst-case for log(1-p) =======================

cat("\n=== Experiment 6: log1p vs log(1-p) on synthetic p near 1 ===\n")

p_grid <- c(1 - 10^(-(1:16)))  ## p = 0.9, 0.99, ... 1 - 1e-16
syn_res <- data.frame()
for (p in p_grid) {
  cmp <- log1m_compare(p)
  syn_res <- rbind(syn_res, data.frame(
    p = p, one_minus_p = 1 - p,
    log_naive = cmp$naive, log_safe = cmp$safe,
    abs_err = cmp$abs_err, rel_err = cmp$rel_err
  ))
}

print(syn_res)
write_csv(syn_res, "06-log1p-synthetic.csv")

## ---- summary ---------------------------------------------------------------

cat("\n=== Summary ===\n")
cat(sprintf("Caterpillar max rel_err vs reference: %.3e\n",
            max(cater_res$rel_err, na.rm = TRUE)))
cat(sprintf("Balanced max rel_err vs reference   : %.3e\n",
            max(bal_res$rel_err, na.rm = TRUE)))
cat(sprintf("Ascertainment max rel_err           : %.3e (worst p = %g)\n",
            max(asc_res$rel_err, na.rm = TRUE),
            asc_res$p_const[which.max(asc_res$rel_err)]))
cat(sprintf("Singleton subtraction max rel_err   : %.3e\n",
            max(sing_res$rel_err, na.rm = TRUE)))
cat(sprintf("Synthetic log1p worst rel_err       : %.3e (at p = 1 - %g)\n",
            max(syn_res$rel_err, na.rm = TRUE),
            syn_res$one_minus_p[which.max(syn_res$rel_err)]))

cat("\nResults written to:", out_dir, "\n")
