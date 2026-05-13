# Findings: slow mixing under `empirical_geometric`

**Status: COMPLETED — 2026-05-13**

Case 31 (joint p/k' MH) implemented and integrated. Local test suite
green (empirical: 38/38 pass, gibbs: 156/156 pass). Awaiting Hamilton
validation on slow-mixing tasks (t01_r01 etc.).

---

## Status at completion

Design A (joint p/k' marginal MH, case 31) **implemented and tested**:

- `src/mcmc.cpp`: `gibbs_kprime_sweep_impl` extended with optional
  `sampleAndUpdate` flag and `outLogZ` pointer returning the per-character
  Gibbs normaliser ∑_i log Z_i(p). Case 31 calls it twice — at oldP
  (no mutation) and newP (with sampling) — accepts on marginal posterior
  ratio + logit Jacobian. Snapshots k'/logLik/logPrior/partLogLik for
  rollback.
- `R/RunMkPrime.R`: `.kMoveTypes[["joint_p_kprime"]] = 31L`, registered
  in `assemble_moves` empirical_geometric branch with `weight = nTrans`
  alongside existing `mh_logit_p`, plumbed through `.DoMove`,
  `.BuildTuningMatrix`, `.AdaptTuning` targets/keys, category map, and
  scalar weight floor list.
- `R/MkPrimeMCMC.R`: default `scale_joint_p_kprime = 1.0` added.
- `tests/testthat/test-m092-adaptive-scheduler.R`: registered new name
  in the `validInMcmc` whitelist.

**Test status:** `devtools::test(filter="empirical")` → 38 pass, 0 fail.
`devtools::test(filter="gibbs")` → 156 pass, 0 fail. The adaptive
scheduler on the empirical test promotes `joint_p_kprime` to 17.3% of
move weight on its own — direct evidence the new move is doing useful
work and being rewarded by the warmup scheduler.

**Outstanding:** confirm on Hamilton that ESS-per-second improves on the
slow-mixing tasks (t01_r01 etc. from the trace analysis). The local
tests cannot reproduce the (p, k')-ridge pathology — that's a
big-problem phenomenon. Next: rerun the slow-mixing job once on Hamilton
with case 31 enabled and compare ESS to the prior tail.

## Earlier context (kept for the diagnostic story)

Last updated: 2026-05-13 (after Hamilton trace analysis + prior
formulation read).

## Root cause

**A strong negative posterior correlation between `p` and aggregate
k' creates a (p, k')-space ridge that no single-variable move can
traverse efficiently.** Measured on 6 Hamilton tail samples
(17140607/17141180, 50k thinned samples each):

```
task       cor(p, mean_k')   cor(p, max_k')   p mean   p IQR
t01_r01    −0.83             −0.62            0.53     [0.42, 0.63]
t03_r09    −0.83             −0.55            0.54     [0.43, 0.63]
t04_r09    −0.81             −0.44            0.64     [0.54, 0.74]
t05_r05    −0.83             −0.64            0.45     [0.36, 0.53]
t06_r03    −0.83             −0.62            0.52     [0.41, 0.61]
t08_r05    −0.72             −0.24            0.81     [0.72, 0.90]
```

The prior structure makes this inevitable. From src/mcmc.cpp L269–315:

```
P(k'_i = m | p) = ∑_{j=2..m}  P_emp(j) · p · (1 − p)^(m − j)
```

where m = k', j = observed-state component, (m−j) = unseen-state
component ~ Geometric(p). At p → 1 the prior concentrates k' near the
empirical body (k' ≈ kObs); at p → 0 the geometric tail dominates and
the prior tolerates k' ≫ kObs. So the data tells the chain "there are
~50 characters; some can have high k' and some can have low k'", and
the chain trades off:

- Configuration A: high p (~0.8), most k' = kObs, few unused states.
- Configuration B: low p (~0.4), some k' wandering to 30–50, many
  unused states.

Single-variable moves explore *along* the ridge slowly because
changing p requires coordinated changes in many k' values to keep the
posterior high. `gibbs_kPrime` only resamples k' *conditional on
current p*; `mh_logit_p` only resamples p *conditional on current k'*.
Neither can step *across* the ridge.

## Why this looked like "logPrior dominant + flat logL"

Per-task summary across the same 6 tasks (50k samples each):

```
metric                       median (across 6 tasks)
sd(log_likelihood)           5.5
sd(log_prior)                13.0     ← 2.4× the likelihood sd
sd(log_posterior)            14.9
cor(log_lik, log_prior)      +0.16    (positive, small)
max k' seen                  41–51
median per-char k' range     20–33
```

Variance of log_prior is dominated by `nTrans · log p + (∑u)·log(1−p)`,
which moves with the (p, k') ridge. The likelihood stays flat because
adding unused states (k' > kObs) barely changes the model fit. So the
trace looked like "k' moves through wide prior support" — true on the
surface, but the underlying constraint is the (p, k') correlation.

## Why my local toy did not reproduce it

13-tip × 50-char × n=80 toy, nChains=2, heat=0.2: the posterior
correlation is similar in structure but the chain has too few samples
and too small a problem to manifest the slow-mixing pathology. Local
acceptance-rate signals didn't separate eg from bg.

## What the Hamilton data shows for the move budget

Header of every Hamilton log preserves the post-warmup move weights:

```
Characters: mh_logit_p 55.2%, gibbs_kPrime 20.6%,
            kPrime 1.1%, block_kPrime 1.1%
Topology:   nni/spr/gibbs_spr/gibbs_subtree_swap/tbr/pspr each 1.1%
Branches:   local_dirichlet 8.3%, others 1.1% each
Rates:      slice_rate_log_sd 1.9%, rate_log_sd 1.1%, joint_tl_rls 1.1%
```

So the adaptive scheduler has rewarded `mh_logit_p` (the case-30 fix
from commit ca85725) — it gets 55% of cycles — and `gibbs_kPrime` 20%.
Each works fine on its own dimension, but neither can move along the
ridge.

Topology mixes adequately: 5137–11047 distinct topology hashes per
50k samples per task. Topology is not the bottleneck.

## Recommended fix

**Implement the Plan's Design A: joint (p, k') Metropolis-Hastings
move (case 31).**

The proposal:
1. Propose Δ in logit(p), apply to get p' (symmetric proposal on logit
   scale, Hastings factor = Jacobian of logit, same as `mh_logit_p`).
2. Resample k' from the per-character full conditional given p'
   (Gibbs draw, identical to `gibbs_kprime_sweep_impl`).
3. Accept with α = min(1, π(p') / π(p)), where
   π(p) = ∫ π(p, k') dk' is the marginal posterior on p after
   integrating out k'.

MH algebra: with q(p', k' | p, k'_old) = q_p(p' | p) · g(k' | p', data)
and q_p symmetric on logit, the k' Gibbs draw cancels via
g(k' | p', data) = π(p', k') / π(p'). The ratio reduces to the
marginal posterior ratio only.

**The marginal is cheap.** Per-character k'_i are independent given
p (under the prior) and per-character likelihoods factor, so
log π(p) = log π_p(p) + ∑_i log [∑_{k'_i} π(k'_i | p) · L_i(k'_i)].
The inner sum is exactly the normalization constant computed inside
`gibbs_kprime_sweep_impl` per character — already present, just needs
to be returned alongside the draw.

**Effort:** ~80–150 lines (modify gibbs_kprime_sweep_impl to return
log-marginal-on-p alongside the new k' values; add case 31 calling it
twice — once at current p to get π(p), once at proposed p' to get
π(p'); add R move-table entry).

**Test baseline established before any code change:**
`devtools::test(filter="empirical")` → 38 pass, 0 fail, 1 skip (slow),
2 stabilisation warnings (expected — short warmups). Green.

**What NOT to do** (my earlier recommendation that I'm retracting):
- Don't truncate the empirical_geometric prior. The ridge would still
  exist; truncation just compresses both ends slightly. Also risks
  semantic changes the package's tests may rely on (test suite has
  empirical-prior tests; need to be careful).

**Lower-priority follow-ups, once Design A is in:**
- `block_kPrime` shows ~1% scheduler weight (floor) and was 0%
  acceptance on local runs. Confirm propose-counts on a fresh run with
  acceptance dump enabled — if truly DOA, remove.
- The `kPrime` int_walk at 1.1% (floor) is also negligible weight; can
  drop to 0 without effect.

## Step-by-step status

- **Step 1 (logP split)** ✅ Hamilton tails confirm logPrior dominant,
  cor(L,Pr)=+0.16, ratio 2.44×.
- **Step 2/3 (per-move acceptance, swap rate)** ⚠️ Per-move acceptance
  not recoverable (no `.rds` summary written because tasks timed out).
  Move-weights from log header confirm plan's "mh_logit_p 55%".
  Swap_cold ≈ 1 per thinned sample (~10 iters); topology mixing OK.
- **Step 4 (topology vs scalar ESS)** — Topology proxy via hash diversity
  (5–11k distinct in 50k samples) is healthy. Topology not the
  bottleneck. Detailed tree-ESS deferred.
- **Step 5 (ESS estimator)** ✅ Sound for ρ ≥ 0.99.
- **Step 6 (wall-time profile)** ✅ R-level Rprof uninformative; move
  weights from header are the right diagnostic.
- **Step 7 (prior–data mismatch)** ✅ Mild on the body (KL median
  0.078). The geometric tail is what enables the (p, k')-ridge mixing
  pathology, not body mismatch.
- **Step 8 (model gotcha)** ✅ No strict k'-rate ridge. The actual
  ridge is in (p, k'), driven by the convolution structure of the
  prior.

## Files

- `inst/hamilton/mkp-study-retrieve.R` — fixed for tail-sample
  retrieval
- `inst/hamilton/mkp-study-inventory.R` — remote inventory helper
- `inst/hamilton/mkp-study-17140607-17141180/mkp_eg_sample/` — 6
  representative tails (~330 MB)
- `data-raw/step1_hamilton_logp_decomp.R` — Hamilton trace analysis
- `data-raw/step1_hamilton_summary.rds` — per-task summary
- `data-raw/step{5,6,6b,7}_*.R` — local toy work (superseded by
  Hamilton analysis where they conflict)
- `remote-jobs.md` — job 17140607 (timeout) + 17141180 (running)
  entries
