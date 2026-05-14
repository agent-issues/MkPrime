# REPORT: Aligned-units re-sim (corrected diagnosis)

**Date:** 2026-05-14 (run + analysis)
**Run:** `dev/pilots/2026-05-14-aligned-units-resim/`
**Seed (data):** 20260514
**Seed (MCMC):** default init
**nIter:** 100 000  thin=40  post-burnin n=1594
**Elapsed:** ~20 min

## Headline

**The "tree_length bias" framing was wrong.** The chain is not biased in its
parameter medians — it is failing to mix in topology space. With effective
ESS near zero on all continuous parameters and 1093 unique topologies in
1594 post-burnin samples, the posterior medians are essentially random
samples from a random-walk trajectory, not estimators of a stationary
posterior.

The unit-mismatch hypothesis from `dev/pilots/2026-05-13-tl-bias-investigation/REPORT.md`
is **not confirmed and not the relevant question**: aligning simulator and
model parameterisations does not improve recovery because the chain never
gets close to the true topology under either parameterisation.

## What was checked

### 1. Posterior summaries (raw, no RelabelEcology)

| Parameter | Median | IQR | Truth |
|-----------|--------|-----|-------|
| phi | 3.15 | [2.14, 4.69] | 4 |
| pi0 | 0.32 | [0.18, 0.49] | 0.75 |
| theta_1 | 0.98 | [0.86, 1.00] | 1.0 |
| tree_length | 39.1 | [25.1, 58.1] | 12.7 |
| gamma_e (derived) | 2.07 | [1.53, 3.16] | 1.75 |

### 2. ESS (batch-means)

All near zero (≤ 12 out of 1594 samples). This alone invalidates median-based
verdicts.

### 3. Ridge diagnostic (advisor-recommended)

Correlations on log-transformed posterior:

| Pair | r |
|------|---|
| log tl vs log γ_e | −0.06 |
| log tl vs log φ | −0.01 |
| log tl vs log rate_neo | −0.16 |
| log tl vs pi0 | 0.02 |
| pi0 vs theta | −0.27 |

`cor(log tl, log γ_e)` is essentially zero (well below the 0.4 plan
threshold). **There is no tl ↔ γ ridge.** The previous chain A's reported
ridge correlation of +0.34 was likely chain-A-specific noise, not a
structural identifiability problem.

### 4. Topology recovery (the real story)

Clustering-information distance to truth tree:

| Stat | Value | Interpretation |
|------|-------|---------------|
| min normalised CID | 0.491 | Chain never gets within 49% of truth |
| median normalised CID | 0.749 | Typical sample is 3/4 of the way to random |
| max normalised CID | 0.886 | Worst samples are near-random |
| Exact matches (CID=0) | 0 / 1594 | Truth topology never visited |

Bipartition support:

| Bipartition | Posterior probability |
|-------------|----------------------|
| True ecology grouping (A1-4 + C1-4) | 0.006 |
| Wrong lineage grouping (A1-4 + B1-4) | 0.000 |

Modal topology holds 0.4% of post-burnin samples; 1093 unique topologies in
1594 samples. The chain is essentially performing a random walk in tree
space, never finding a high-likelihood basin.

### 5. Comparison with chain A (default sim params)

| Metric | Chain A | Aligned |
|--------|---------|---------|
| min CID to truth | 0.332 | 0.491 |
| P(true bipartition AC) | 0.000 | 0.006 |
| Posterior tl median | 20.4 | 39.1 |

**Both chains have the same fundamental problem.** Chain A is slightly less
bad in topology proximity but neither recovers truth. The "tl bias"
difference (20.4 vs 39.1) is just where the random-walk happened to wander.

## Implications

### What the previous "unit mismatch" investigation missed

The investigation in `dev/pilots/2026-05-13-tl-bias-investigation/REPORT.md`
focused on a parameter-recovery question (does tl median equal truth?)
without first verifying that the chain has converged. With ESS near zero,
that question is undefined — there is no stationary posterior to recover
from.

The earlier framing of `cor(log tl, log γ_e) = 0.344` as a ridge was over-
interpreted from a single chain with poor mixing. The aligned re-sim shows
the correlation is −0.06; the value in chain A was likely a transient
trajectory artifact.

### What's actually wrong: the dataset has insufficient signal

**Cross-checking existing sim3 results confirms the problem is not
ecology-specific:**

| Run | Blind CID min | Aware CID min | Notes |
|-----|--------------|---------------|-------|
| sim3-mcmc-result.rds | 0.274 | 0.344 | Original v1 |
| sim3-v2-result.rds | 0.282 | 0.274 | After v2 fixes |
| sim3b-mcmc-result.rds | 0.183 | 0.189 | Negative control |

**Every blind chain on this dataset also fails to visit the true topology.**
Zero exact matches in any run, including v1 and v2, including blind and
aware, including the negative-control sim3b. The aware chains aren't
worse than blind by much — the ecology layer is roughly neutral, not
harmful, but it cannot rescue a dataset with weak phylogenetic signal.

This means the previous Sim 3 framing — "ecology-aware should recover the
convergent tree" — was asking too much of the data. The convergent-tree
Goldilocks config (16 tips × 240 chars, stemBr=0.10, rootBr=0.15) has
strong convergent character noise relative to phylogenetic signal *by
design*; that's the point of the test. But the data isn't informative
enough for *any* model to actually recover truth — the test as defined
cannot succeed.

### Likelihood is correctly implemented

For the record: I briefly mis-claimed during this analysis that the
ecology-aware likelihood was unimplemented. That was wrong — see
`src/mcmc_ecology.cpp:285` (`pruning_jc_acrv_flat_ecology`), line 507
(`pruning_mkn_acrv_flat_ecology`), line 895
(`cpp_log_likelihood_ecology`), and line 920 (gamma_e_compute call).
The likelihood code applies the per-(edge, character) rate factors
`z=0→1/γ, z=1→φ/γ, z=2→(1/φ)/γ` correctly inside the pruning loop.

### Real outcome of this pilot

This diagnoses a **fundamental MCMC-mixing problem on the Sim 3
convergent-tree dataset**, independent of any v2 prior fix. The v2 prior
fix work (asymmetric slab, theta logging, RelabelEcology) is sound and
necessary; it's just not sufficient.

## Recommendations

### Do NOT

- Dispatch Step 5 HPC multi-rep until topology mixing is understood and
  fixed.
- Update `sim3-multirep.R` simulator params to "aligned" values — that
  does not fix the underlying issue.
- Trust any parameter-recovery claims from chains where ESS is near zero.

### Should

1. **Redesign Sim 3 to have recoverable phylogenetic signal.** Either:
   - Increase data: 32+ tips, 500+ chars, or
   - Strengthen branch lengths so phylogenetic signal dominates
     convergent character noise, or
   - Reframe the success metric: instead of "recover the truth topology
     exactly", evaluate "does ecology-aware sample closer-to-truth
     topologies than blind?" via bipartition probabilities or CID
     deltas. The existing aware-vs-blind CID comparison in
     `sim3-multirep.R` already does this — that's the right metric.

2. **Truth-init topology test.** Run a chain initialised from the truth
   tree and check whether it leaves the truth basin. If it doesn't, the
   posterior is concentrated near truth but the chain can't find that
   region from random starts (mixing problem). If it does leave, the
   posterior itself doesn't concentrate on truth (likelihood landscape
   problem).

3. **Use sim3-multirep.R as designed.** The existing multi-rep driver
   already measures bipartition probabilities and CID across replicates
   — those are robust to single-chain mixing failure. The headline figure
   should be "aware better than blind in N of 20 reps" not "aware median
   tl equals truth".

4. **Document the topology-mixing limit in vignette.** Be honest that
   v2 Sim 3 cannot exactly recover the convergent-tree topology in single
   chains; report bipartition probabilities instead.

## Key files

- `dev/pilots/2026-05-14-aligned-units-resim/run.R` — this run
- `dev/pilots/2026-05-13-step3b-mode-persistence-A/` — chain A baseline
- `dev/pilots/2026-05-13-step3b-mode-persistence-B-truthInit/` — chain B
  (truth-param init, same dataset)
- `inst/simulations/ecology/sim3-pilot.R` — the Goldilocks config under
  investigation
- `inst/simulations/ecology/sim3-multirep.R` — multi-rep driver (blocked)
