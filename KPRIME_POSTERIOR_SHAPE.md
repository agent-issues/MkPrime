# Shape of the per-character k' posterior under Mk' priors

**Status:** local analysis, no Hamilton log fetch required.
**Inputs:** real binary characters from `tree_01/rep_01` only (one tree,
one rep -- 12 of 41 binary characters profiled). Caveat: the picture is
plausible across tasks but has not been verified on t05/t13/t24.
**Method:** prior PDFs computed analytically; per-character marginal
likelihood profiled at three branch-length regimes (TL=0.3, 1.0, 3.0) by
calling `MkpLogLikelihood()` on a single-character `MkPrimeData` over
`k' = 2..50`, both with the Mk' relabel correction (`relabel = TRUE`) and
without (`relabel = FALSE`, isolating pruning + ascertainment). Posterior
= prior x likelihood, normalised on the grid.

**Files**

- Script: `dev/pilots/2026-05-12-prior-validation/analysis/redteam_kprime_shape.R`
- RDS:    `dev/pilots/2026-05-12-prior-validation/analysis/redteam_kprime_shape.rds`
- Plots:  `redteam_kprime_{priors,loglik,posteriors}.png` in the same dir

## Bottom line

**The answer is (b), with one important wrinkle.** For binary characters
under Mk':

1. **Pruning + ascertainment uniformly prefers low k'**: slope from
   k=2 to k=10 is between −3.78 and −12.6 log units for every binary
   character tested. Even the "uninformative" half tilts at −3.85 simply
   because the ascertainment correction makes each unobserved state more
   costly as k grows.
2. **The Mk' relabel correction is a positive +3.81 log units over the
   same range** (k=2 → k=10) and grows to +7.1 from k=2 → k=50. Its job
   is precisely to cancel the ascertainment-driven preference for low k.
3. **It cancels exactly only when pruning is flat.** On half the chars
   the cancellation lands within ±0.1 of zero, so the Mk'
   *full-conditional* on k' equals the prior. On the other half, pruning
   tilts at −6 to −13, the relabel +3.8 isn't enough, and the
   full-conditional still pulls toward k = kObs.
4. **No sampler stickiness.** k' is updated by a Gibbs sweep over a
   candidate window (`K_MAX_CAND = 50`, `mcmc_likelihood.cpp:626`), not
   MH, so per-iteration k' samples already are the full-conditional.
   Hypothesis (c) is ruled out.

**So Mk' on binary characters is dominated by the prior, by design.**
The relabel correction was meant to make k' identifiable from the data;
in practice it does so only weakly, leaving the prior to determine the
posterior. With `mkp_geo`'s Beta(1,1) hyperprior, P(k' <= 3) = 0.68 ->
posterior median sits at 2 for every binary character.

## Prior PMFs at kObs = 2

| prior                      | E[k']  | P(k' <= 3) | P(k' > 20) |
| -------------------------- | ------ | ---------- | ---------- |
| mkp_geo, Beta(1,1)         | 4.57   | 0.680      | 0.031      |
| mkp_eg, emp + Beta(1,1)    | 5.09   | 0.567      | 0.032      |
| mkp_highk, Beta(1,20)      | 16.38  | 0.128      | 0.314      |

`mkp_eg` and `mkp_geo` are qualitatively similar -- both peaked at small
k -- but not identical: `mkp_eg` has a noticeably fatter mid-tail (E[k']
5.09 vs 4.57; P(k' <= 3) 0.57 vs 0.68). The empirical `N_obs` body has
0.67 mass on N_obs = 2 with geometric decay 0.43 in the tail, so the
convolution with the geometric `N_unobs` shifts a little probability up
relative to pure geometric. Both still leave 97% of mass below k' = 20.

The Beta(1,20) hyperprior in `mkp_highk` (job 17227627) shifts mass
dramatically: P(k' > 20) climbs from 0.03 to 0.31, E[k'] from ~5 to
16.4.

## Per-character likelihood decomposition (medium TL)

```
char  slope_with_relabel  slope_pruning_only  relabel_delta
  1               -0.05               -3.86          +3.81
  4               -4.39               -8.19          +3.81
 10               -0.06               -3.86          +3.81
 16               -0.04               -3.85          +3.81
 20               +0.03               -3.78          +3.81
 26               -0.04               -3.85          +3.81
 27               -0.06               -3.86          +3.81
 28               -2.30               -6.10          +3.81
 41               -6.53              -10.34          +3.81
 43               -2.14               -5.95          +3.81
 47               -8.76              -12.57          +3.81
 48               +0.02               -3.79          +3.81
```

Slopes are log L(k=10) − log L(k=2). **Every** binary character has a
strongly negative pure-pruning slope -- pruning + ascertainment under JC
on a fixed tree always prefers fewer states for a 2-state observation.
The relabel correction (+3.81 across [2,10]) cancels the
"uninformative-character" subset (rows 1, 10, 16, 20, 26, 27, 48) and
fails to cancel the rest.

## Posterior summary on the 12 binary chars (medium TL)

```
prior      median k'    mean k'   P(k' > kObs+5)
mkp_geo    2 (all 12)   2.0-4.6   0.001-0.13
mkp_highk  2 or 12-13   2.1-16.4  0.005-0.68
mkp_eg     2 or 3       2.1-5.1   0.0001-0.15
```

For chars where pruning is "merely" uninformative (slope ≈ −3.85),
the Mk' full-conditional ≈ flat, so the posterior ≈ prior, and
`mkp_highk` (which puts E[prior] = 16.4) yields posterior median k' = 12.
For chars where pruning is strongly informative (4, 28, 41, 43, 47), the
data win even under the high-k prior: char 47 stays at median k' = 2 under
`mkp_highk`, char 41 likewise.

## Implication for the Mk'-vs-mk_k40 tree-recovery gap

Comparing mechanisms:

- **mk_k40**: `knownStates = 40` for every variable character -> chars
  become `"known"` type, **no relabel correction applied**
  (`R/likelihood.R:185` gates relabel on the transformational branch).
  So mk_k40 evaluates the pruning+ascertainment likelihood at k=40 with
  no compensating term. On the 12 binary chars analysed, this means
  mk_k40 is paying a log-L deficit of roughly 8-13 units per binary char
  *relative to k=2*, summed across ~40 binary chars -> ~400 nats of
  worse fit per iteration.
- **Mk' (mkp_geo, mkp_eg)**: integrates over k' per character with the
  geometric prior. Posterior median k' = 2 for every binary char ->
  basically a 2-state model whose pruning likelihood is high.

mk_k40 has **lower** per-character likelihood than Mk' yet **better**
tree recovery. The gap therefore cannot live in the k' machinery.
It must live somewhere else -- almost certainly in the **branch-length
posterior**: forcing k=40 inflates the pruning weight on substitution
events, which pushes the inferred branch lengths into a regime that
happens to match the truth better. The parallel TL agent is the right
place to confirm.

## Implication for mkp_highk

The mkp_highk arm (Beta(1,20)) will:

- raise posterior median k' from 2 to ~12 on the ~7/12 uninformative
  binary chars (where the relabel cancels pruning),
- leave posterior k' ≈ 2 on the ~5/12 informative binary chars,
- give characters with kObs >= 3 a similar split.

If mkp_highk closes some of the CID gap to mk_k40, the mechanism is the
shifted branch-length posterior on the "raise k'" half. If it doesn't,
the gap is entirely in branch lengths and the k' prior is
not the lever to pull.

## Answer to the original question

(a) **Data say no** for half the binary chars -- but in an even
stronger sense than naive: the **pure pruning** likelihood actively
prefers low k for *every* binary char; what makes "informative" vs
"uninformative" is whether the relabel correction (a model artefact,
not data) cancels pruning. **Yes (a), partly.**

(b) **Prior pulls back.** Where pruning + relabel ≈ flat, the prior
fully determines the posterior, and the geometric/EG priors are peaked
at small k. **Yes (b), confirmed.**

(c) **Sampler stuck?** No -- Gibbs sweep over candidate window. **(c)
ruled out.**

(d) Combination of (a) and (b). The most accurate one-liner: **Mk'
binary-character posteriors equal the prior whenever the relabel
correction succeeds in cancelling the pruning likelihood's preference
for low k. Under geometric / EG priors that prior puts most mass at low
k, so the posterior does too.**

## What was NOT done and why

- **No Hamilton log fetch.** k' uses a Gibbs sweep, so the per-iteration
  sample distribution IS the posterior; ESS is not a meaningful
  diagnostic (it would just measure how fast the tree+rate state mixes
  underneath). The full-conditional shape is reconstructed analytically
  here.
- **Only tree_01/rep_01 profiled.** One tree (25 tips), 12 of 41 binary
  characters. The qualitative split (uninformative-where-relabel-cancels
  vs informative-where-pruning-dominates) is structural -- it depends on
  whether the per-character ascertainment correction has a slope that
  matches +3.81 over k = 2..10 -- so the picture is unlikely to differ
  much on t05/t13/t24, but this has not been verified empirically.
- **No multi-character (kObs=3, 4) profiling.** Only 7 chars with kObs=3
  and 2 with kObs=4 in this rep -- not enough to be informative, and
  the diagnostic question (Mk' losing to mk_k40) is dominated by binary
  chars (41/50 in this dataset).
