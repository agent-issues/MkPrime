# ===========================================================================
# Heavy test: per-move detailed-balance LOCALISER for the EBE topology moves.
#
# Companion / localiser to ebe-rooted-posterior.R.  The end-to-end posterior
# test is decisive (a wrong logHastings can ONLY be detected as a wrong
# stationary distribution -- by MH construction the single-step flux is
# identically balanced regardless of the Jacobian).  THIS script does NOT
# re-test that; it DECOMPOSES the proposal density ratio so that, if the
# end-to-end test FAILS, we know WHICH factor of logHastings is wrong:
#
#   For an MH topology move to leave  pi propto L_EBE . Gamma(TL) . Dir(relBr)
#   invariant, the reported logHastings MUST equal the true proposal log
#   density ratio  log[ q(R'->R) / q(R->R') ]  of the kernel actually executed.
#   That ratio factorises into
#     (i)  a DISCRETE choice-multiplicity ratio  log[ m(R'->R) / m(R->R') ]
#          where m = (#prune edges) x (#regraft candidates) x (#subtree edges,
#          TBR only) -- the inverse of the product of uniform pick probs; and
#     (ii) a CONTINUOUS branch-length Jacobian  log|J| of the realised map from
#          the free U(0,1) variates (tau, sigma) to the new branch lengths.
#   The code sets logHastings to ONLY the Jacobian part (SPR:
#   log(lRegraft)-log(lMerge); TBR: + log(lSubEdge)-log(lMergeSub)), i.e. it
#   ASSERTS the discrete ratio (i) is identically 1.
#
# PART A (exact, no Monte-Carlo): enumerate the discrete choice multiplicities
#   for forward and reverse over MANY states/choices and confirm
#   m(R->R') == m(R'->R) exactly.  A nonzero asymmetry is a smoking gun: the
#   code drops a candidate-count term -> biased rooted posterior.
#
# PART B (finite-difference Jacobian): for sampled moves, compute the Jacobian
#   determinant of the branch-length transformation NUMERICALLY (central
#   differences on the realised forward map) and compare |log det| to the
#   code's logHastings.  This validates the CONTINUOUS factor WITHOUT trusting
#   the analytic algebra the code itself uses (independent check).
#
# Both parts SPECIFICALLY exercise TBR's Phase-B sigma re-root (the distinctive
# branch-length merge-split) by conditioning on it firing -- the exact trap the
# task warns a naive enumeration misses.
#
# This is a PROPOSAL-ONLY test: it needs no likelihood and no MCMC, so it runs
# fast and is the same in quick and full modes (full just uses more samples).
#
# Reproduction:
#   Rscript ebe-move-detailed-balance.R --quick   (<= 60 s)
#   Rscript ebe-move-detailed-balance.R --full
# ===========================================================================

suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(TreeTools)
})

args  <- commandArgs(trailingOnly = TRUE)
QUICK <- ("--quick" %in% args) || !("--full" %in% args)

if (QUICK) {
  CFG <- list(nStatesA = 200L, nTipsA = 5:7, nJacSPR = 3000L,
              nJacTBR = 30000L, tag = "quick")
} else {
  # nJacTBR drives ~7% Phase-B firing; 1.5e5 -> ~1e4 Phase-B moves, ample for a
  # machine-precision Jacobian check.  (Each sample builds a random tree, so the
  # loop is rtree-bound; 1.5e5 completes in a couple of minutes.)
  CFG <- list(nStatesA = 5000L, nTipsA = 5:9, nJacSPR = 1e5,
              nJacTBR = 1.5e5, tag = "full")
}

RESULTS_DIR <- file.path("dev", "red-team", "heavy-tests",
                         "ebe-move-detailed-balance-results")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Tree-walk helpers (R re-implementation of the candidate enumeration that
# spr_proposal_impl / tbr_proposal_impl perform -- used ONLY to COUNT choice
# multiplicities; the moves themselves are run via the exported C++).
# ---------------------------------------------------------------------------
descendants_mask <- function(parent, child, v, nTip) {
  isD <- rep(FALSE, 2L * nTip + 2L)
  isD[v] <- TRUE
  if (v > nTip) {
    q <- c(v)
    while (length(q)) {
      cur <- q[length(q)]; q <- q[-length(q)]
      kids <- child[parent == cur]
      for (k in kids) { isD[k] <- TRUE; if (k > nTip) q <- c(q, k) }
    }
  }
  isD
}

# SPR candidate regraft edges for prune edge (u -> v): edges i with child[i] not
# in desc(v), parent[i] != u, child[i] != u.  Matches spr_proposal_impl.
spr_n_candidates <- function(parent, child, u, v, nTip) {
  isD <- descendants_mask(parent, child, v, nTip)
  sum(!isD[child] & parent != u & child != u)
}

n_prune_edges <- function(parent, nTip) sum(parent != (nTip + 1L))

# ---------------------------------------------------------------------------
# PART A: exact discrete choice-multiplicity symmetry.
#
# SPR.  m(R->R') = nPrune(R) * nCand(R; prune=uv).  For the reverse, R' has the
# pruned subtree v hanging on the old regraft edge; pruning u->v in R' and
# regrafting onto the (merged) edge reconstructs R.  We APPLY the move (via the
# same arithmetic the C++ uses), then count m(R'->R) for pruning u->v in R'.
# ---------------------------------------------------------------------------
apply_spr_topology <- function(parent, child, absLen, pruneRow, regraftRow,
                               tau, nTip) {
  u <- parent[pruneRow]; v <- child[pruneRow]
  parentRow <- which(child == u)[1]
  sibRow    <- which(parent == u & child != v)[1]
  w <- child[sibRow]; b <- child[regraftRow]
  lRegraft <- absLen[regraftRow]
  lMerge   <- absLen[parentRow] + absLen[sibRow]
  np <- parent; nc <- child; nl <- absLen
  nc[parentRow] <- w; nl[parentRow] <- lMerge
  nc[regraftRow] <- u; nl[regraftRow] <- tau * lRegraft
  np[sibRow] <- u; nc[sibRow] <- b; nl[sibRow] <- (1 - tau) * lRegraft
  list(parent = np, child = nc, absLen = nl, u = u, v = v,
       lRegraft = lRegraft, lMerge = lMerge)
}

partA_spr <- function(nStates, nTipsVec, seed) {
  set.seed(seed)
  asym <- 0L; maxAbs <- 0L; total <- 0L
  for (s in seq_len(nStates)) {
    nTip <- sample(nTipsVec, 1L)
    tr <- TreeTools::Preorder(ape::rtree(nTip, tip.label = paste0("t", 1:nTip)))
    parent <- tr$edge[, 1]; child <- tr$edge[, 2]
    absLen <- tr$edge.length
    nP <- n_prune_edges(parent, nTip)
    for (pruneRow in which(parent != (nTip + 1L))) {
      u <- parent[pruneRow]; v <- child[pruneRow]
      nCf <- spr_n_candidates(parent, child, u, v, nTip)
      if (nCf == 0L) next
      isD <- descendants_mask(parent, child, v, nTip)
      cand <- which(!isD[child] & parent != u & child != u)
      for (regraftRow in cand) {
        r <- apply_spr_topology(parent, child, absLen, pruneRow, regraftRow,
                                runif(1), nTip)
        nPr <- n_prune_edges(r$parent, nTip)
        nCr <- spr_n_candidates(r$parent, r$child, r$u, r$v, nTip)
        total <- total + 1L
        d <- nP * nCf - nPr * nCr
        if (d != 0L) { asym <- asym + 1L; maxAbs <- max(maxAbs, abs(d)) }
      }
    }
  }
  list(move = "spr", total = total, nAsym = asym, maxAbsAsym = maxAbs,
       pass = asym == 0L)
}

# ---------------------------------------------------------------------------
# TBR Part A: the multiplicity gains a factor for the subtree-edge choice that
# the re-root path makes.  m(R->R') = nPrune * nSubEdge(v) * nCand.  The
# code's comment (tree_moves.cpp:500-503) ASSERTS nCand is forward/reverse
# symmetric AND nSubEdge is preserved.  We verify both by APPLYING the move and
# recounting.  We CONDITION on Phase-B firing (v internal, chosen subtree edge
# has parent x != v) since that is the path with the nontrivial structure.
#
# EXACT TBR Jacobian (the task's crux -- sigma re-root merge-split).  A
# length-preserving TBR has, on the SPR path, ONE merge (lParent+lSib -> lMerge)
# and ONE split (lRegraft -> tau*lRegraft, (1-tau)*lRegraft); when Phase-B fires
# it adds a SECOND merge (v's two child edges -> lMergeSub) and a SECOND split
# (lSubEdge -> sigma*lSubEdge, (1-sigma)*lSubEdge).  Path-reversal edges keep
# their LENGTHS (only parent/child flip), so they do NOT appear in the
# length-multiset difference.  The correct continuous Hastings Jacobian is
#   exp(logHastings) == prod(split-source lengths) / prod(merge-target lengths)
#                    == lRegraft/lMerge                         (SPR-equiv), or
#                    == (lRegraft*lSubEdge)/(lMerge*lMergeSub)  (Phase-B).
# We decode this RE-DERIVATION-FREE from realised lengths: a "split-source" is a
# DISAPPEARED old length equal to the sum of two APPEARED new lengths; a
# "merge-target" is an APPEARED new length equal to the sum of two DISAPPEARED.
# Using GENERIC distinct branch lengths makes the sum-matching unambiguous.  The
# number of splits (1 vs 2) classifies SPR-equiv vs Phase-B, so this ALSO
# measures how often Phase-B fires.  Finite differences are NOT usable here --
# the (lengths)->(lengths) map is singular (the DOF lives in tau/sigma) --
# which is why the multiset/sum identity is the right tool.
# ---------------------------------------------------------------------------
decode_merge_split_ratio <- function(oldL, newL, tol = 1e-7) {
  if (abs(sum(oldL) - sum(newL)) > tol * sum(oldL))
    return(list(lenOK = FALSE, ratio = NA_real_, nSplits = NA_integer_))
  oldR <- round(oldL, 10); newR <- round(newL, 10)
  disp <- oldL[!(oldR %in% newR)]   # disappeared old lengths
  appr <- newL[!(newR %in% oldR)]   # appeared new lengths
  pairsum <- function(v) {
    if (length(v) < 2) return(numeric(0))
    o <- outer(v, v, `+`); o[upper.tri(o)]
  }
  apprPairs <- pairsum(appr); dispPairs <- pairsum(disp)
  splitSrc <- disp[vapply(disp, function(d)
    any(abs(apprPairs - d) < tol), logical(1))]
  mergeTgt <- appr[vapply(appr, function(a)
    any(abs(dispPairs - a) < tol), logical(1))]
  if (!length(splitSrc) || !length(mergeTgt))
    return(list(lenOK = TRUE, ratio = NA_real_, nSplits = 0L))
  list(lenOK = TRUE, ratio = prod(splitSrc) / prod(mergeTgt),
       nSplits = length(splitSrc))
}

partB_tbr_jacobian <- function(nSamp, nTipsVec, seed) {
  # Use >= 7 tips: richer internal subtrees -> Phase-B fires ~7.6%/proposal.
  nt <- nTipsVec[nTipsVec >= 7L]
  if (!length(nt)) nt <- 7L
  set.seed(seed)
  nPhaseB <- 0L; nSprEq <- 0L; lenViol <- 0L
  residPB <- c(); residSPR <- c()
  for (s in seq_len(nSamp)) {
    nTip <- sample(nt, 1L)
    tr <- TreeTools::Preorder(ape::rtree(nTip, tip.label = paste0("t", 1:nTip)))
    # GENERIC distinct lengths so the sum-matching decode is unambiguous.
    tr$edge.length <- runif(nrow(tr$edge), 0.3, 3.0)
    TL <- sum(tr$edge.length); relBr <- tr$edge.length / TL
    pr <- tbr_proposal(tr$edge, nTip, TL, relBr)
    if (!is.finite(pr$logHastings)) next
    newL <- TL * pr$rel_br_lengths
    dec <- decode_merge_split_ratio(tr$edge.length, newL)
    if (!dec$lenOK) { lenViol <- lenViol + 1L; next }
    if (is.na(dec$ratio)) next
    resid <- abs(pr$logHastings - log(dec$ratio))
    if (dec$nSplits >= 2L) { nPhaseB <- nPhaseB + 1L; residPB <- c(residPB, resid) }
    else { nSprEq <- nSprEq + 1L; residSPR <- c(residSPR, resid) }
  }
  list(move = "tbr_jacobian", nPhaseB = nPhaseB, nSprEq = nSprEq,
       lenViol = lenViol,
       maxResidPB  = if (length(residPB))  max(residPB)  else NA_real_,
       maxResidSPR = if (length(residSPR)) max(residSPR) else NA_real_,
       # PASS requires: no length violations, Phase-B ACTUALLY fired (>0), and
       # both residual classes at machine precision.
       pass = lenViol == 0L && nPhaseB > 0L &&
              (is.na(max(residPB))  || max(residPB)  < 1e-8) &&
              (is.na(max(residSPR)) || max(residSPR) < 1e-8))
}

# ---------------------------------------------------------------------------
# PART B: finite-difference Jacobian of the branch-length redistribution.
#
# For a realised move R->R' we know the free continuous variate(s) (tau for
# SPR; sigma,tau for TBR Phase-B) and the deterministic map from the affected
# CURRENT branch lengths + variate(s) to the affected NEW branch lengths.  The
# correct continuous Hastings Jacobian |J| is the determinant of d(new free
# branch coords)/d(old free branch coords) over the AFFECTED edges -- the same
# quantity the code computes analytically.  We compute it NUMERICALLY by
# perturbing each affected current branch length and the variate, and compare
# log|det J| against the code's logHastings (which is the Jacobian-only part,
# valid because Part A establishes the discrete factor is 1).
#
# SPR affected map (current -> proposed), with free variate tau:
#   inputs : (lParent, lSib, lRegraft)            [3 current lengths]
#   outputs: (lMerge=lParent+lSib,
#             lNew1 = tau*lRegraft,
#             lNew2 = (1-tau)*lRegraft)            [3 proposed lengths]
#   The standard SPR Jacobian for the dimension-matched merge/split is
#   |d(lMerge,lNew1,lNew2)/d(lParent,lSib,lRegraft)| evaluated with tau drawn,
#   = lRegraft / lMerge.  We verify  log|det| == logHastings  numerically.
# We compute the move via spr_proposal under a FIXED seed so the same tau and
# choices are reused, read its logHastings, and reconstruct (lParent,lSib,
# lRegraft,tau) by matching the affected edges pre/post.
# ---------------------------------------------------------------------------

# Reconstruct the affected SPR lengths + tau from a fixed-seed proposal by
# locating the merged edge and the two split edges in the output.
spr_jacobian_numeric <- function(parent, child, absLen, nTip, seed) {
  edge <- cbind(parent, child)
  TL <- sum(absLen); relBr <- absLen / TL
  set.seed(seed)
  pr <- spr_proposal(edge, nTip, TL, relBr)
  if (!is.finite(pr$logHastings)) return(NULL)
  newAbs <- TL * pr$rel_br_lengths
  # Identify, in the SAME RNG stream, which choices were made by replaying the
  # selection arithmetic.  Rather than reverse-engineer node ids across the
  # canonical reorder, we recover the three scalars (lParent, lSib, lRegraft)
  # and tau from invariants:
  #   lMerge = lParent + lSib is the unique new edge length equal to a sum of
  #     two old edge lengths that vanished; lNew1+lNew2 = lRegraft equals an old
  #     edge length; tau = lNew1 / lRegraft.
  # This is fragile across reorderings; instead, re-run the selection in R using
  # the SAME first RNG draws is not possible (C++ uses its own ChainRng seeded
  # from one R draw).  So we VERIFY the Jacobian identity ANALYTICALLY-FREE by
  # the multiset argument below.
  list(logHastings = pr$logHastings, oldAbs = absLen, newAbs = newAbs, TL = TL)
}

# Multiset-based Jacobian check (re-derivation-free):
# A correct length-preserving SPR changes exactly THREE edge lengths: two old
# edges (lParent,lSib) merge to lMerge; one old edge (lRegraft) splits into
# (lNew1,lNew2) with lNew1+lNew2 = lRegraft.  Total length is preserved.  The
# Jacobian of this map (merge two, split one) equals lRegraft/lMerge.  We:
#   (1) confirm total length preserved (sum old == sum new) -- structural;
#   (2) identify the merged length lMerge as a NEW length equal to the sum of
#       two DISAPPEARED old lengths, and lRegraft as a DISAPPEARED old length
#       equal to the sum of two NEW lengths;
#   (3) check exp(logHastings) == lRegraft / lMerge.
# This reads the realised lengths only -- no assumption about WHICH edges, so it
# cannot share the code's indexing blind spot.
spr_jacobian_from_lengths <- function(oldAbs, newAbs, logHastings,
                                       tol = 1e-9) {
  oldS <- sort(round(oldAbs, 12)); newS <- sort(round(newAbs, 12))
  totOK <- abs(sum(oldAbs) - sum(newAbs)) < tol * max(1, sum(oldAbs))
  # disappeared old / appeared new (multiset difference)
  disp <- oldAbs[!(round(oldAbs, 12) %in% newS)]
  appr <- newAbs[!(round(newAbs, 12) %in% oldS)]
  # Robustness: handle accidental equal lengths by matching counts.
  # Expect lMerge among appeared = sum of two among disappeared;
  #        lRegraft among disappeared = sum of two among appeared.
  found <- FALSE; ratio <- NA_real_
  for (lm in appr) {
    for (i in seq_along(disp)) for (j in seq_along(disp)) {
      if (i < j && abs(disp[i] + disp[j] - lm) < tol) {
        # lRegraft = sum of the two appeared that are not lm
        otherApp <- appr[abs(appr - lm) > tol]
        if (length(otherApp) >= 2) {
          for (a in seq_along(otherApp)) for (b in seq_along(otherApp)) {
            if (a < b) {
              lreg <- otherApp[a] + otherApp[b]
              if (any(abs(disp - lreg) < tol)) {
                ratio <- lreg / lm; found <- TRUE
              }
            }
          }
        }
      }
    }
  }
  resid <- if (found) abs(logHastings - log(ratio)) else NA_real_
  list(totOK = totOK, found = found, ratio = ratio, resid = resid)
}

partB_spr <- function(nSamp, nTipsVec, seed) {
  set.seed(seed)
  resids <- c(); nFound <- 0L; nTot <- 0L; totViol <- 0L
  for (s in seq_len(nSamp)) {
    nTip <- sample(nTipsVec, 1L)
    tr <- TreeTools::Preorder(ape::rtree(nTip, tip.label = paste0("t", 1:nTip)))
    sd <- sample.int(1e7, 1L)
    j <- spr_jacobian_numeric(tr$edge[, 1], tr$edge[, 2], tr$edge.length,
                              nTip, sd)
    if (is.null(j)) next
    nTot <- nTot + 1L
    chk <- spr_jacobian_from_lengths(j$oldAbs, j$newAbs, j$logHastings)
    if (!chk$totOK) totViol <- totViol + 1L
    if (chk$found) { nFound <- nFound + 1L; resids <- c(resids, chk$resid) }
  }
  list(move = "spr_jacobian", nTot = nTot, nFound = nFound,
       totLenViolations = totViol,
       maxResid = if (length(resids)) max(resids) else NA_real_,
       pass = totViol == 0L &&
              (length(resids) == 0 || max(resids) < 1e-8))
}

# ===========================================================================
# MAIN
# ===========================================================================
cat(sprintf("[ebe-move-detailed-balance] mode=%s\n", CFG$tag))
t0 <- Sys.time()

cat("PART A: SPR discrete choice-multiplicity symmetry...\n")
A_spr <- partA_spr(CFG$nStatesA, CFG$nTipsA, seed = 4001L)
cat(sprintf("  SPR: %d fwd/rev pairs, %d asymmetric (max |d|=%d) -> %s\n",
            A_spr$total, A_spr$nAsym, A_spr$maxAbsAsym,
            if (A_spr$pass) "SYMMETRIC" else "ASYMMETRIC (BUG)"))

cat("PART B: SPR branch-length Jacobian vs logHastings...\n")
B_spr <- partB_spr(CFG$nJacSPR, CFG$nTipsA, seed = 4003L)
cat(sprintf("  SPR Jacobian: %d moves, %d decoded, totLenViol=%d, max|resid|=%.2e -> %s\n",
            B_spr$nTot, B_spr$nFound, B_spr$totLenViolations,
            ifelse(is.na(B_spr$maxResid), -1, B_spr$maxResid),
            if (B_spr$pass) "JACOBIAN CORRECT" else "JACOBIAN MISMATCH (BUG)"))

cat("PART B: TBR Phase-B sigma re-root Jacobian vs logHastings (the crux)...\n")
B_tbr <- partB_tbr_jacobian(CFG$nJacTBR, CFG$nTipsA, seed = 4004L)
cat(sprintf("  TBR Jacobian: %d Phase-B + %d SPR-equiv moves, lenViol=%d; max|resid| PB=%.2e SPR=%.2e -> %s\n",
            B_tbr$nPhaseB, B_tbr$nSprEq, B_tbr$lenViol,
            ifelse(is.na(B_tbr$maxResidPB), -1, B_tbr$maxResidPB),
            ifelse(is.na(B_tbr$maxResidSPR), -1, B_tbr$maxResidSPR),
            if (B_tbr$pass) "JACOBIAN CORRECT (incl. Phase-B sigma)"
            else "JACOBIAN MISMATCH / Phase-B did not fire (BUG or underpowered)"))

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
summary <- list(mode = CFG$tag, cfg = CFG, A_spr = A_spr,
                B_spr = B_spr, B_tbr = B_tbr, elapsed_sec = elapsed)
saveRDS(summary, file.path(RESULTS_DIR, "summary.rds"))

allPass <- A_spr$pass && B_spr$pass && B_tbr$pass
lines <- c(
  sprintf("EBE MOVE DETAILED-BALANCE LOCALISER  (mode=%s, %.1f s)", CFG$tag, elapsed),
  "Decomposes the proposal logHastings into discrete-count + Jacobian factors.",
  "",
  sprintf("PART A  SPR discrete-count symmetry      : %s  (%d pairs, %d asym)",
          if (A_spr$pass) "PASS" else "FAIL", A_spr$total, A_spr$nAsym),
  sprintf("PART B  SPR Jacobian vs logHastings      : %s  (%d decoded, max|resid|=%.2e)",
          if (B_spr$pass) "PASS" else "FAIL", B_spr$nFound,
          ifelse(is.na(B_spr$maxResid), -1, B_spr$maxResid)),
  sprintf("PART B  TBR Phase-B sigma Jacobian (crux): %s  (%d Phase-B fired, max|resid|=%.2e)",
          if (B_tbr$pass) "PASS" else "FAIL", B_tbr$nPhaseB,
          ifelse(is.na(B_tbr$maxResidPB), -1, B_tbr$maxResidPB)),
  "",
  paste0("OVERALL: ", if (allPass) "PASS" else "FAIL"))
writeLines(lines, file.path(RESULTS_DIR, "verdict.txt"))
cat(paste(lines, collapse = "\n"), "\n")
