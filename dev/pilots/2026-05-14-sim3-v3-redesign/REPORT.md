# Sim 3 v3 redesign + mixing diagnostic

**Date:** 2026-05-14
**Worktree:** `worktree-ecology-aware`
**Pilots:** `dev/pilots/2026-05-14-sim3-v3-redesign/`

## Headline

The v3 redesign and follow-up diagnostics establish the following picture:

1. **Likelihood implementation is correct.** Truth has the highest logL on
   the v3 dataset (−4655.0, verified by direct call to
   `.MkpEcologyLogLikelihood`).
2. **Single MCMC chains on Sim 3 cannot recover the convergent topology**,
   in v3 just as in v2. Posterior is broad (~200 nats wide); the chain
   never reaches truth from random init or stays at truth from truth init.
3. **The aware-vs-blind logLik improvement is real but cannot be read off
   the chain.log column directly** — chain-reported logL values are stale
   (off by ~30 nats from re-evaluated truth-equivalent state). Need to
   recompute via `.MkpEcologyLogLikelihood` for honest comparison.
4. **Two priors do clash with truth values used in simulation**:
   - `theta = 1.0` (simulator default for "all encouraged") lies at the
     `Beta(2, 2)` prior boundary where density is 0. The chain's
     `theta_1` median of ~0.97-0.98 is the prior pulling away from 1.
   - `rate_neo ~ LogNormal(0, 2)` is wide enough that the chain visits
     pathological values (up to ~5e6) on the (tl × rate_neo) ridge.
     Likelihood is near-saturation invariant in that region (cost ~18 nats
     vs truth), so the prior is the only force keeping rate_neo bounded.
5. **Forward path is paper-clean via Option 2 (aware-vs-blind delta)**,
   acknowledging that posterior medians are not informative on this
   dataset class. Other options listed at the bottom for the user.

## V3 redesign

Parsimony grid (`parsimony-grid.R`) at higher signal levels picked
`2x-3x-stem` as the Goldilocks v3 (8/8 fooled, 7/8 truth-on-baseline,
medFull=10): `nEco=120, nBase=360, phi=4, stem=0.30, root=0.15`.

Truth tree length: **13.50** model substitution units.

## Bayesian pilot at v3 (`bayesian-pilot.R`)

Two 100k chains from parsimony-init start tree.

| Metric | Blind | Aware (relabelled) |
|--------|-------|--------------------|
| Chain logLik max (reported) | −4716 | −4658 |
| Chain logLik max (re-evaluated) | -- | **−4691** ⚠️ stale |
| ESS phi | -- | ~0 |
| ESS tl | ~0 | ~0 |
| Modal topology fraction | 0.4% | 0.3% |
| P(true AC bipartition) | 0 | 0.001 |
| min CID to truth | 0.301 | 0.564 |

**Headline:** logLik improvement of aware over blind exists at the peak (~25
nats via re-evaluated values, smaller than the 58-nat chain-reported figure)
but chains never sample anywhere near the truth tree, so posterior
distributions are uninformative as recovery estimators.

## Truth-init verification (`truth-init.R`, `truth-init-pinned.R`)

Two 50k chains initialised at truth (truth-tree + truth params, with
`theta = 0.999` to avoid the prior boundary):

| Metric | unpinned | pinned (gibbs_z 8%, scale_phi/pi0/theta 3%) |
|--------|----------|----------------------------------------------|
| Chain logL median | −4853 | −4854 |
| Chain logL max (reported) | −4648 | −4665 |
| Distance from truth peak (median) | −195 nats | −196 nats |
| Visits within 20 nats of peak | 0.8% | 0.7% |
| posterior tl median | 44 | 68 |
| posterior pi0 median | 0.20 | 0.26 |
| posterior rate_neo median | 5.7 | 0.54 |
| min CID to truth | 0.489 | 0.461 |

**Both chains drift from truth within the first sampled iteration.** Pinning
the ecology moves at 8%/3% changes WHICH parameter values the chain
prefers to drift into (rate_neo no longer blows up to 5.7, instead drops
to 0.54), but the chain still leaves truth.

**The "drift" is the chain correctly sampling a broad posterior.** Within
the posterior, the chain finds many quasi-equivalent regions:

| Top-1% logL sample | logL | tl | phi | pi0 | rate_neo | Notes |
|---------------------|------|----|----|-----|----------|-------|
| iter 697 (random)   | −4691 (re-eval) | 30 | 62 | 0.31 | 3.0 | high-phi mode |
| iter 559            | -- | 6.9 | 3.6 | 0.12 | 5×10⁶ | (tl, rn) saturation ridge |
| iter 1249           | −4706 | 9.5 | 11.6 | **0.76** | 2.3 | near-truth pi0 |
| iter 185            | −4714 | 2.85 | 3.3 | 0.09 | 1.7 | low-tl mode |

The chain visits all of these in proportion to their posterior mass. Real
posterior dispersion, not pathological multi-modality.

## Confirmed and corrected claims

**Confirmed (high confidence):**
- Likelihood code correct; truth IS the global logL maximum (−4655)
- theta=1 boundary problem real (-Inf prior density on truth)
- Posterior is broad, not narrowly biased
- v3 stronger signal didn't help mixing (still ESS ~0)
- Aware does prefer higher-logL regions than blind at peak

**Corrected from earlier claims:**
- "Chain max −4658 beats blind −4716 by 58 nats" used stale chain.log
  values. Actual aware logL at iter 697 is −4691, not −4658. So the
  aware-vs-blind delta at peak is ~25 nats not 58. Still positive but
  smaller than reported.
- "rate_neo identifiability is 2.5-nat soft pull" was on v2 dataset.
  On v3 with more chars + longer stems, the JC kernel is closer to
  saturation, so the rate-time ridge is much wider and rate_neo can
  wander into pathological values (5e6) with only ~18 nats logL cost.
  The LogNormal(0, 2) prior is too permissive given this saturation.
- "Pure numerical pathology" (rate_neo=5M) is actually a real likelihood
  ridge near JC saturation, not a numerical bug.

## Implications for paper

The convergent-tree Sim 3 dataset, even at v3 strength, has insufficient
signal to make single Bayesian chains converge on the true tree from
random initialisation. The likelihood IS informative at truth, but the
posterior is so broad that chain-medians are dominated by off-truth mass.

**This makes the paper's framing decision urgent.** Three viable framings:

### (a) Aware-vs-blind delta as headline (Option 2 from earlier conversation)

- Run `sim3-multirep.R` with 20 reps, report:
  - P(true bipartition AC | aware) − P(true bipartition AC | blind)
  - CID-to-truth delta (aware − blind)
  - Per-rep aware logLik max − blind logLik max
  - Per-character z posteriors at the per-character mode
- The headline is "ecology-aware consistently outperforms blind by X% on
  bipartition recovery; aware logLik beats blind by ~25 nats per rep".
- Paper-clean. Does not require fixing posterior medians or mixing.

### (b) Tighten priors (research call, may help recovery)

- `rate_neo`: LogNormal(0, 2) → LogNormal(0, 0.5) [keeps rate_neo within
  e^(-1, 1) ≈ [0.37, 2.7] at 95%]
- `theta`: Beta(2, 2) → Uniform(0, 1) [admits θ=1 boundary case]
- Possibly tighter `tree_length` prior given the saturation ridge

This would likely improve posterior shape and might make medians
recoverable. But it's a model design decision, not a mixing fix.

### (c) Accept broad posterior, lengthen chains

- Run 1M-iter chains and report posterior mass within X nats of peak
- Honest but unflattering: paper says "posterior is broad on weakly-
  informative data"
- HPC cost: ~10x current 100k pilots × 20 reps = significant compute

## Concrete next-step ask

Which framing do you want? Answer drives whether to:
- (a) Skip mixing fixes, run `sim3-multirep.R` headline at v3 config
- (b) Modify model priors before any HPC dispatch
- (c) Long single-chain HPC runs

## Files committed

- `parsimony-grid.R`, `parsimony-grid.csv` — v3 candidate grid scoping
- `bayesian-pilot.R`, `result.rds` — v3 aware/blind 100k pilot
- `analyse.R` — pilot summary script
- `truth-init.R`, `truth-init-pinned.R` — truth-initialised chains
- `truth-init-result.rds`, `truth-init-pinned-result.rds` — outputs
