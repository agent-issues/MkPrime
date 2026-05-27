# Repro: partLogLik partial-cache update on rate_neo / neo_joint moves
# only recomputes data->neoPartIndices, but after the fix rateNeo also
# affects TRANS partition rateScales — so trans partition log-likelihoods
# go stale, and the MH acceptance ratio uses a wrong newLogLik.
#
# This script verifies the diagnosis from source reading alone (the bug
# is internal to do_move_impl case 3 partial-CL branch, lines 4815-4831,
# and the slice sampler eval_slice_target lines 3083-3093 + 3164-3176).
#
# Setup: build a small mixed dataset, populate partition cache, then
# compare:
#   (A) ground-truth full likelihood at rate_neo' = 2 * rate_neo
#   (B) what partLogLik-cache-based partial update would give
#
# If they disagree, the MH ratio is wrong.

library("mkp"); library("TreeTools"); library("ape")

set.seed(1)
ntax <- 6L
mat <- cbind(
  matrix(sample.int(2, ntax * 6, replace = TRUE) - 1L, ntax, 6),
  matrix(sample.int(3, ntax * 4, replace = TRUE) - 1L, ntax, 4)
)
dimnames(mat) <- list(paste0("t", 1:ntax), NULL)
pd  <- MatrixToPhyDat(mat)
mkd <- MkPrimeData(pd, neomorphic = c(1L, 2L, 3L))  # 3 neo, 7 trans

tree <- Preorder(rtree(ntax, br = function(n) runif(n, 0.05, 0.3)))

# Direct evaluation: this exercises cpp_partition_log_likelihood for ALL
# partitions, neo and trans, at rate_neo and at rate_neo * 2.
rn1 <- 1.0
rn2 <- 2.0
ll1 <- MkpLogLikelihood(tree, mkd, rate_neo = rn1)
ll2 <- MkpLogLikelihood(tree, mkd, rate_neo = rn2)

# Compute the per-partition contributions explicitly to show that
# BOTH neo and trans change between rn1 and rn2.
neo_partitions   <- which(vapply(mkd$partitions, function(p) p$type, "") == "neomorphic")
trans_partitions <- setdiff(seq_along(mkd$partitions), neo_partitions)

# Compute per-partition contributions by ablating one partition at a time
# via PartialL = full - leave-one-out, or simply re-call partition LL.
# Easier: split mkd into two single-partition mkds (one per type) and eval.

cat(sprintf("rn=%.2f  full ll = %.6f\n", rn1, ll1))
cat(sprintf("rn=%.2f  full ll = %.6f\n", rn2, ll2))
cat(sprintf("delta full = %.6f\n", ll2 - ll1))

# Now simulate what the buggy partial-CL update path does:
# It would update only the neo partitions and assume trans contribution
# is unchanged from rn1.  We can emulate that by computing
#   (a) trans-only ll at rn1  (the cached, stale value)
#   (b) trans-only ll at rn2  (the true value at the new state)
#   (c) neo-only ll at rn1 and rn2
# The buggy newLogLik would be: (c@rn2) + (a)  -- skipping (b)
# The correct newLogLik is:     (c@rn2) + (b)
#
# We construct trans-only and neo-only mkds for this.

trans_idx <- which(mkd$type != "neomorphic")
neo_idx   <- which(mkd$type == "neomorphic")

# Helper: build single-partition-type mkd from the original phyDat
sub_mat_trans <- mat[, trans_idx, drop = FALSE]
sub_mat_neo   <- mat[, neo_idx,   drop = FALSE]

pd_trans <- MatrixToPhyDat(sub_mat_trans)
pd_neo   <- MatrixToPhyDat(sub_mat_neo)

mkd_trans <- MkPrimeData(pd_trans)  # all trans
mkd_neo   <- MkPrimeData(pd_neo, neomorphic = seq_along(neo_idx))  # all neo

# But ll on these splits is NOT directly comparable to the joint mkd
# because compute_partition_scales depends on JOINT nNeo/nTrans.
# So we instead measure trans ll change in the joint context by
# observing the joint partition LLs across rn1 and rn2 directly.
# Since the helper is straightforward, we'll just confirm via partition_scales:

partition_scales <- function(r, nNeo, nTrans) {
  if (nNeo == 0L || nTrans == 0L) return(c(neo = 1, trans = 1))
  denom <- 1 + r
  nTotal <- nNeo + nTrans
  c(neo = r/denom * nTotal/nNeo, trans = 1/denom * nTotal/nTrans)
}

cat("\nPartition scales:\n")
print(partition_scales(rn1, 3L, 7L))
print(partition_scales(rn2, 3L, 7L))

# transScale at rn=1: 1/2 * 10/7 = 0.714
# transScale at rn=2: 1/3 * 10/7 = 0.476
# So trans edge lengths shrink by 0.476/0.714 = 0.667 -> trans partition LL changes.
# Buggy partial-CL update keeps cached trans LL from rn=1 and adds it to
# new neo LL at rn=2 -> wrong newLogLik -> wrong MH acceptance ratio.
