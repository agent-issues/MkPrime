# SBC-TL-MIX-001 — diagnosis

**Mechanism:** *forward simulator bug in the SBC harness, not a production-code bug.*
`dev/red-team/heavy-tests/sbc.R:139` iterates the edge list in `rev(seq_len(nrow(edges)))`,
which is **postorder** on a `TreeTools::Preorder` edge list. The loop reads
`states[pa]` before that parent's state has been drawn (it is still the
`integer(...)` default `0`). Only the edge whose parent is the root sees a
correctly-drawn parent state; every other edge draws its child against
state 0. Children of internal nodes are then overwritten correctly when those
nodes are themselves reached as children later in the same reverse traversal,
but by then the *tip* draws against the wrong parent state are already locked in.

Result: tips are simulated approximately as independent Bernoulli draws
biased toward state 0 along each tip's incident edge. Phylogenetic
correlation is largely destroyed; inference (which uses correct JC pruning)
fits this by collapsing tree length downward.

## Evidence (sim 15, seed 20262541, 8 tips, 28 informative binary chars,
true topology fixed, `tl_true = 4.134`)

### Part A — likelihood profile, `rel_br` fixed at truth's normalised weights

`MkpLogLikelihood(..., kPrime = 2, nCat = 1, coding = "variable",
relabel = TRUE)` on the production data, evaluated on a TL grid:

| TL    | logL    | logPrior | logPost |
|-------|---------|----------|---------|
| 0.500 | -100.49 | -7.15    | -107.64 |
| 1.000 | -97.27  | -6.48    | -103.75 |
| **1.200** | **-97.09** | -6.30 | **-103.39** |
| 1.400 | -97.24  | -6.16    | -103.40 |
| 1.500 | -97.41  | -6.09    | -103.50 |
| 2.000 | -98.88  | -5.83    | -104.70 |
| 4.000 | -108.35 | -5.21    | -113.56 |
| **4.134 (truth)** | **-108.99** | -5.18 | -114.17 |
| 8.000 | -122.29 | -4.68    | -126.97 |
| 50.000| -135.39 | -4.53    | -139.92 |

* **Brent MLE (likelihood-only): TL = 1.199** (logL = -97.091)
* **MAP (logL + Gamma(2, 2/50) logPrior): TL = 1.200**
* logL(TL = truth) − logL(TL = MLE) = **−11.6 nats**: truth is exp(11.6) ≈ 1.1×10⁵ times less likely than the chain's mode under the inference likelihood, *given the truth's own topology and branch proportions*.
* Plain-Mk profile (`relabel = FALSE`) **also** peaks at TL = 1.2: the Mk' relabel
  correction is **not** the cause.

### Part B — sampler is healthy

20 000 sampling iterations × `autoTune = TRUE` × `fixTopology = TRUE` × start
on the true topology:

| Move           | accept |
|----------------|--------|
| tree_length    | **0.44** |
| branch_lengths | 0.60   |
| dirichlet_branch | 0.17 |
| local_dirichlet | 0.19  |
| kPrime / gibbs_kPrime / p | 1.00 |
| rate_log_sd    | 0.39   |
| joint_tl_rls   | 0.19   |

* Final adapted `scale_tree_length` = **1.61** (Bactrian scale; healthy width).
* 70/70 unique `tree_length` values in the kept window (no freezing).
* Posterior `tree_length`: mean = 1.64, sd = 0.45, 95% CI [0.88, 2.61].
  This is exactly what the analytic logPost peak predicts (≈ 1.20 with mild
  upward pull from Gamma(2, 0.04) prior whose mode is 25).

The chain is correctly sampling from the model's posterior. The posterior
just isn't centred at truth, because *the data is not what the inference
model thinks the truth generated*.

### Part C — confirmation: swap `rev(seq_len(...))` → `seq_len(...)` in
`.simJCchar`, refit MLE

Same seed, same RNG stream for char draws, same tree:

* **buggy simulator** (current `sbc.R`): n_char kept = 27/30, **MLE TL = 2.14**
* **fixed simulator** (forward iteration): n_char kept = 29/30, **MLE TL = 4.59**
  ≈ tl_true = 4.13

30-seed sweep on fresh trees:

| simulator | median MLE / tl_true | mean MLE / tl_true | sims with MLE < tl_true |
|-----------|----------------------|---------------------|-------------------------|
| buggy     | **0.80**             | 1.63                | 16/30                   |
| fixed     | 1.51                 | 2.17                | 12/30                   |

The buggy simulator's median MLE sits **20% below** truth (and below the
fixed simulator's median by a factor ~1.9). The fixed simulator's MLE has
modest positive small-sample bias (n_char ≈ 28; MLEs at 30 chars per sim are
noisy but unbiased on average is too strong a claim from 30 sims — what
matters is direction). The buggy simulator's downward shift is the direction
of the observed v5 SBC failure (`tree_length` ranks pile HIGH → posterior
sits LOW relative to truth).

## Recommended fix (do NOT apply here — surface for next campaign step)

One-line change in `dev/red-team/heavy-tests/sbc.R:139`:

```diff
-  for (e in rev(seq_len(nrow(edges)))) {
+  for (e in seq_len(nrow(edges))) {
```

The trailing comment block around lines 141-144 should be tightened to
state that forward iteration on a `TreeTools::Preorder` edge list is the
required direction for a root-to-leaves draw; the rate-convention note can
stay.

## Downstream implications

* **Every SBC submission since the harness was first run is contaminated.**
  v1-v5 results for the universal `tree_length` HIGH bias are explained by
  this single bug and do **not** indicate any defect in `RunMkPrime` /
  `MkpLogLikelihood` / `LogPrior`. SBC-WARMUP-001 (refuted) and
  SBC-WARMUP-002 (a real but independent adaptive-scheduler floor issue
  affecting `N_WARM = 20000` only) were both correctly red-herrings for
  the universal bias.
* **EG-001's empirical 0.053 borderline in `Mkp_empirical_geometric` is
  also sitting on contaminated data** — re-evaluate against a clean v6 run.
* SBC-PRIOR-CONFOUND-001 (already withdrawn 2026-05-27) and
  SBC-MASS-FAIL-001 (kPrime parameterisation, partially fixed in v3/v4)
  remain valid harness defects independent of this one.

## Update 2026-05-27 (retracted, see follow-up below) — SBC-HARNESS-005: rev() fix was necessary but not sufficient

v6 (job with the rev() fix; `sbc.R` now iterates forward) still showed
`tree_length` FAIL with mean rank ~52/67 (expected 33.5) across 4 arms
where results were collected.

**Root cause of the residual bias:** the start-tree topology, not the
simulator. With `fixTopology = TRUE`, the harness conditions MCMC on the
starting tree's *topology*. A valid SBC test of continuous parameters under
fixed topology requires the starting topology to equal the true topology.
The harness was using `AdditionTree(pd)` as the start, which is parsimony-
based. At `tl_true ~ Gamma(2, 0.04)` (mean=50), near-saturation makes data
nearly uninformative about topology:

* 10/10 sampled sims: RF(AdditionTree, true_tree) ≥ 8 (maximum = 10 for 8 taxa)
* 8/10 sims: RF = 10 (completely wrong topology)
* sbc-warmup-trace sim 1, RF=10 start: rank = 67/67 (worst case)
* sbc-warmup-trace sim 1, `--true-topo` (RF=0): rank = 61/67 (plausible under H0)

The MCMC with the wrong topology samples `p(tree_length | y, T_wrong)`;
the wrong topology forces a shorter tree to explain the data, pushing
the TL posterior systematically below truth → ranks pile HIGH.

**Claim in this document that "v1-v5 results are explained by the rev() bug"
is now understood to be incomplete.** The rev() bug (SBC-HARNESS-004)
contributed to bias in both the simulator and the likelihood surface, but
the topology conditioning was also wrong in every v1–v6 run — and was the
dominant driver of the HIGH rank signal seen in v6 after the rev() fix.

**Fix (SBC-HARNESS-005, applied 2026-05-27):** sbc.R start-tree block
replaced with:
```r
start_tree <- true_tree
start_tree$edge.length <- rep_len(0.1, nrow(true_tree$edge))
```
This makes the harness explicitly test `p(θ|y, T_true)`, which is the
correct and unambiguous SBC target for continuous parameters under a
fixed-topology MCMC. v7 is the first run with both fixes applied.

## Update 2026-05-27 (follow-up — retracts the SBC-HARNESS-005 stanza)

The SBC-HARNESS-005 stanza above is **incorrect** and is retracted. v7 also
failed for the same reason v6 failed: **SBC-HARNESS-006** — the Hamilton
SLURM submit script invoked an out-of-tree stale copy of `sbc.R` at
`/nobackup/pjjg18/mkp-study/red-team/heavy-tests/sbc.R` (last touched
2026-05-27 13:16, before SBC-HARNESS-003/004/005 were committed). Every
SLURM SBC run since that point silently used the buggy rev()-postorder
simulator and the AdditionTree start. The in-tree fixes (this document's
"Recommended fix" + SBC-HARNESS-005) were correct, but the SLURM jobs
never loaded them.

**Discriminator (decisive):** sim 1 (seed 20261527) under the buggy
simulator gives `n_var = 27`, kObs after filtering matches v7's saved
`summary.rds`; the in-tree simulator gives `n_var = 30`. The v7 forward
sim matches the buggy version exactly.

**Local re-run with in-tree sbc.R:** MkNT_geometric, 200 sims, same seeds
as v7 — `tree_length` AD p = 0.9860. Clean PASS. So:

* The rev() bug WAS the dominant driver of `tree_length` HIGH rank pile
  in every run from v1 onward (consistent with this document's original
  analysis).
* The "needs topology fix" stanza was based on v7 results that came from
  a code path that didn't include the rev() fix — so the v7 failure
  doesn't constrain whether the topology fix is needed independently.
* Whether the topology fix is *necessary* (vs nice-to-have) is now
  testable empirically: v8 will run with both fixes and the correct path.

**Lesson:** trust the smoking-gun discriminator (`n_var` per sim, kObs
pattern) over rank-pattern interpretation. The MLE flat-likelihood check
at TL=25..500 was decisive: it said "data carries no info, posterior
must be ≈ prior; v7's posterior wasn't" — pointing at the simulator/data,
not the MCMC. The path bug was upstream of both.

## Out-of-lane (ruled out)

* Sampler stuck (mechanism 1 in the brief): refuted by Part B — acceptance
  0.44, scale tuned to 1.61, 70/70 unique samples.
* Edge-Dirichlet prior pathology (mechanism 3): refuted by Part A — logPrior
  on TL is monotonic, prior pulls *up* toward 25, not down toward 1.4.
* Numerical loss in `fast_neg_exp` / `expm1` near rt → 0 (FAST-EXP-001):
  audited stable in this regime; sim 15 has TL ≈ 4 → mean rt ≈ 0.6, nowhere
  near the rt ≲ 1e-12 region.

## Reproducer

* `dev/red-team/heavy-tests/sbc-tl-mix-001-diagnose.R` — Part A + B
  (likelihood profile + instrumented MCMC). Runs in ~2 s.
* `dev/red-team/heavy-tests/sbc-tl-mix-001-confirm.R` — Part C (buggy vs
  fixed simulator). Runs in ~30 s for the 30-seed sweep.
