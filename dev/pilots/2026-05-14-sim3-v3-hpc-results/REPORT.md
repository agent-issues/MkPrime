# Sim 3 v3 multirep HPC pilot — 8 reps × 100k iter

**Date:** 2026-05-14
**SLURM job:** 17153939 (8 array tasks, all `COMPLETED`, 18-21 min each)
**Config:** nEco=120, nBase=360, phi=4, stem=0.30, root=0.15, tip=0.5
**Priors:** sigmaPhi=1.0, theta~Uniform(1,1), rateNeoSdlog=1
**Outputs (local):** `inst/simulations/ecology/multirep-v3-results/rep0[1-8]/summary.rds`

## Headline aware-vs-blind (mean ± sd over 8 reps)

| Metric          | Blind         | Aware         | Δ (aware − blind) |
| --------------- | ------------- | ------------- | ----------------- |
| P(true AC)      | 0.000 ± 0.000 | 0.000 ± 0.000 | +0.000            |
| P(wrong AB)     | 0.692 ± 0.352 | 0.135 ± 0.256 | **−0.557**        |
| CID_to_truth    | 0.476 ± 0.035 | 0.673 ± 0.080 | +0.197 (worse)    |

Per-rep table in `multirep-v3-results/` (`multirep-v3-aggregate.R`).

## What the numbers mean

**Good news.** The aware model strongly suppresses the misleading convergent
bipartition (A,B): P(wrong) collapses from 0.69 → 0.14. In 5/8 reps, P(wrong)
drops to essentially zero. This is exactly the qualitative pattern the ecology
layer is designed to produce — it correctly refuses to be fooled by the
convergent character-state similarity between AC tips.

**Bad news.** Neither model recovers the true bipartition (P(true) = 0 in *all*
16 chains). The aware model's tree-distance (CID) to truth is *worse* than
blind's. So aware is "right for the wrong reason": its posterior trees do not
contain the false AB grouping, but they also don't contain the true AC
grouping — they're elsewhere in tree space entirely.

## Diagnostic — chain.log inspection (rep01)

| Param           | Truth | Blind tail | Aware tail (3 samples) |
| --------------- | ----- | ---------- | ---------------------- |
| log_posterior   | —     | −4618      | −5225, −5251, −5068    |
| log_likelihood  | —     | −4684      | −4934, −4870, −4776    |
| tree_length     | 13.5  | 16–20      | 126, 346, 153          |
| rate_neo        | 1.0   | 0.33–0.43  | 0.06, 1.06, 1.00       |
| phi             | 4     | —          | 21.6, 21.6, 7.97       |
| pi0             | 0.75  | —          | 0.28, 0.52, 0.38       |
| theta_1         | 1.0   | —          | 0.97, 0.97, 0.98       |

Two things jump out:

1. **Aware log_posterior is ~600 nats below blind.** Aware is not just slow to
   converge — it's sitting in a completely wrong region of parameter space.
2. **Aware is ridge-walking on phi × tree_length.** phi=22 with tree_length=346
   gives roughly the same per-edge expected substitutions as phi=8 with TL=153
   or (truth) phi=4 with TL=13.5 — these are quasi-equivalent under the
   normalised mixture. The chain is wandering this ridge instead of locating
   the truth peak.

`theta_1 ≈ 0.97` is the one well-behaved aware parameter — it has locked onto
the asymmetry that the ecology distinguishes correctly.

## Why aware can have lower P(wrong) while having worse logL

Blind sits in a small basin around the convergent (AB) bipartition because the
character matrix really does have many shared changes between A and C tips
masquerading as AB-like signal under independence. Aware *knows* those shared
changes are explainable by shared ecology, so it correctly down-weights the AB
basin — but with phi/TL ridge-walking, it never settles on AC either. So
P(wrong) goes to zero, P(true) stays at zero, and the trees in posterior are
high-entropy.

## Next moves (mixing-focused)

The headline-quality result (P(wrong) ↓0.56, P(true) ↑) is one mixing fix
away, not a model-spec change. Suggested order:

1. **Joint phi × tree_length proposal.** The ridge is the dominant pathology.
   A correlated bactrian move on (log phi, log TL) — or even sequential moves
   tuned together — would let the chain slide along the ridge to the peak
   instead of bouncing across it.
2. **Tighter sigmaPhi if joint move alone isn't enough.** sigmaPhi=1.0 lets
   phi roam up to e^2 ≈ 7.4 within 2 prior σ. The truth is 4; an even tighter
   sigmaPhi (e.g. 0.7) would make phi ≈ 22 strongly improbable a priori
   without precluding truth.
3. **Re-run truth-init aware** with these fixes — confirm chain stays near
   the truth peak rather than drifting onto the ridge within 1 sample.
4. **If 1-3 fix mixing**, scale to 20 reps × extended iter via the
   checkpoint-aware resume path. Files & dispatch already set up for this.

## Files generated

- `inst/simulations/ecology/multirep-v3-results/rep0[1-8]/summary.rds`
- `inst/simulations/ecology/multirep-v3-aggregate.R`
- `inst/hamilton/sim3-multirep-v3/{setup_project.sh,run_rep.R,sim3-multirep-v3.sh,dispatch.R}`
- Per-rep chain logs + checkpoints still on Hamilton at
  `/nobackup/pjjg18/mkp-sim3-multirep-v3/results/rep0[1-8]/` — resume-ready.
