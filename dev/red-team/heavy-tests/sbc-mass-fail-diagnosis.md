# SBC mass-failure diagnosis (job 17289622, HEAD a37c3a6)

**Source artefact for this diagnosis**: the *Hamilton results stored on disk*
were produced by sbc.R at commit `a37c3a6`, not by current main (`8e4e370`).
`git show a37c3a6:dev/red-team/heavy-tests/sbc.R` confirms three structural
forward!=inference prior mismatches that current main has already partially
fixed. ALL SIX arms fail (verdict.txt: MkNT_geometric is also FAIL — the task
description was stale).

## (a) Root cause per FAIL mode

### tree_length (all 6 arms, ranks pile at 0, mean rank/L ~ 0.17)

**Harness bug, now fixed in 8e4e370.** Old `.simTree` drew edge lengths from
`runif(0.05, 0.30)` giving `tl_true` median 2.46 (sd ~0.3); inference prior is
`Gamma(2, 2/50)` with mean 50, sd 35. Forward truth lies in the extreme left
tail of inference prior, so posterior tl is *always* above truth -> rank=0
dominates. (Verified: stored `tl_true[1]=2.094`; replayed under a37c3a6's
.simTree -> 2.0–3.1 range. Replayed under current 8e4e370's Gamma .simTree
-> 17.8 / matches inference prior.)

### rate_log_sd (all 6 arms, ranks pile at 0, mean rank ~17–19)

**Harness bug.** Old line 233: `rateLogSd_true <- min(rgamma(1,1,1), 3.0)`.
Forward is truncated, inference is not. Posterior > truth -> rank=0. Now fixed
in 8e4e370 (clamp removed, comment says so).

### kPrime_pooled (4 Mk' arms, bimodal at 0 and high-rank)

**Harness bug surviving in 8e4e370** — *not* a sampler bug, *not* EG-001
(the EG-001 result is just a special case of this). Inference geometric prior
is `Geo(k'-kObs)` (LogPrior:470: `u <- kPrime[transIdx] - kObs[transIdx]`),
but forward draws `kTrue = rgeom(p) + 2` (sbc.R:194, kFloor=2). The two priors
agree only when kObs==2 for every character. kObs table shows ~53% of
characters have kObs=2 and 47% have kObs in {3..8}. For kObs>2 characters the
forward prior puts mass on kTrue=2 (impossible under inference, since
inference enforces k'>=kObs). The 60% of sims with median kPrime-rank<5
(truth is below posterior mass) are kObs>2-dominated; the 8.5% with
median-rank>60 are kObs=2-dominated where forward's untruncated Geo lets
kTrue go large while posterior collapses to ~2 (per KPRIME_POSTERIOR_SHAPE.md).

### p in Mkp_geometric (rank 0), p in Mkp_empirical_geometric (rank 67)

Tied to the kPrime bug. p posterior depends on the prior counts `sum(u)` over
characters. Mismatched forward shift -> systematic posterior bias. EG has
*added* missing-Z bias (L1/EG-001) on top, pushing in the opposite direction.

## (b) Predicted post-patch (job 17295308) results

Job 17295308 uses commit 8e4e370 which fixes `.simTree` (Gamma+Dirichlet) and
the rateLogSd clamp. It does NOT fix the kPrime forward-shift.

Predicted verdict per arm:

| Arm                       | tree_length | rate_log_sd | kPrime_pooled | p           |
|---------------------------|-------------|-------------|---------------|-------------|
| MkNT_geometric            | PASS        | PASS        | —             | —           |
| MkNT_logseries            | PASS        | PASS        | —             | —           |
| Mkp_geometric             | PASS        | PASS        | **FAIL**      | FAIL (knock-on) |
| Mkp_beta_geometric        | PASS        | PASS        | **FAIL**      | —           |
| Mkp_empirical_geometric   | PASS        | PASS        | **FAIL**      | FAIL (EG-001) |
| Mkp_logseries             | PASS        | PASS        | **FAIL**      | —           |

If tree_length / rate_log_sd PASS in MkNT arms but kPrime still FAILs in Mk'
arms with the same bimodal rank shape, **diagnosis (a) is confirmed**. If
tree_length still FAILs, look for an additional Dirichlet-vs-internal-tree
mismatch (RunMkPrime may reparameterise edges).

## (c) Smallest disambiguating test

Run the existing harness with `--arm Mkp_geometric` filtered to kObs==2-only
characters (require all simulated chars to be binary in step 4; drop the sim
otherwise). If kPrime_pooled then PASSes, the forward-shift bug is confirmed;
if it still FAILs, look at the Gibbs k' sampler. Cost: ~5 min, quick mode.

## Power statement

These are interpretations of an existing FAIL, not a new test result, so no
power claim is needed. Quick mode (`Rscript sbc.R --quick`) cannot detect any
of these — N_SIM=5, N_CHAR=8 is far below the L=67, n>200 needed for AD
sensitivity. Do not let quick PASS imply anything.
