# Sim 3 v3 multirep HPC results — pre-fix vs post-fix

**Dates:**
- Pre-fix multirep: 2026-05-14 morning, job 17153939 (8 reps, all `COMPLETED`, 18-21 min each)
- Post-fix multirep: 2026-05-14 afternoon, job 17157994 (8 reps, all `COMPLETED`, 37-47 min each; ~2× slower due to periodic eco-resync + tighter convergence)

**Config:** nEco=120, nBase=360, phi=4, stem=0.30, root=0.15, tip=0.5
**Outputs (local):** `inst/simulations/ecology/multirep-v3-results/rep0[1-8]/summary.rds`

## Headline aware-vs-blind (mean ± sd over 8 reps)

| Metric          | Blind (pre)   | Aware (pre)   | Δ (pre)    | Blind (post)  | Aware (post)  | Δ (post)   |
| --------------- | ------------- | ------------- | ---------- | ------------- | ------------- | ---------- |
| P(true AC)      | 0.000 ± 0.000 | 0.000 ± 0.000 | +0.000     | 0.001 ± 0.004 | 0.002 ± 0.004 | +0.001     |
| P(wrong AB)     | 0.692 ± 0.352 | 0.135 ± 0.256 | **−0.557** | 0.551 ± 0.396 | 0.381 ± 0.341 | −0.170     |
| **CID_to_truth** | 0.476 ± 0.035 | 0.673 ± 0.080 | +0.197 (worse) | 0.484 ± 0.034 | **0.432 ± 0.049** | **−0.052** (better) |

## What the post-fix numbers say

1. **The expected sign on CID is now present.** Aware-vs-blind delta in
   CID-to-truth flipped from +0.197 (aware *worse*) to −0.052 (aware
   *better*). This is the qualitative pattern the ecology layer was
   designed to produce. The pre-fix +0.197 was the artefact of a chain
   that was destabilised by the kPrime / slice / topology-Gibbs
   accumulator leaks and ended up wandering in low-density posterior
   regions.
2. **P(wrong) suppression is muted (-0.170 vs prior -0.557) but real
   in aggregate.** High variance: 3/8 aware reps still show
   P(wrong) > 0.4 (reps 1, 5, 7). These are reps where the chain found
   the wrong bipartition despite the eco layer — i.e. the data alone
   was too consistent with AB.
3. **P(true) is still effectively zero.** Neither model recovers the
   true AC bipartition in 7/8 reps. This is a topology-mixing problem
   — the chain rarely visits the AC bipartition even though the
   ecology layer correctly discourages AB. With `gibbs_spr` and
   `gibbs_subtree_swap` gated in eco mode (S-4), topology moves rely
   on random-MH spr/nni/tbr/pspr only; weight redistribution may be
   needed.
4. **Zero `[eco-resync` drift warnings across all 8 reps × 100k iter.**
   The accumulator is rock-solid in production.

## Comparison to pre-fix multirep

The pre-fix run looked like a strong aware-vs-blind story on P(wrong)
(−0.557) but came with a catastrophically worse CID (+0.197). It was a
**false win**: the chain reported a posterior with no AB bipartition
not because it had found the truth, but because the accumulator bugs
were preventing it from settling anywhere coherent. Aware logged its
trees in a high-entropy "neither AB nor AC" region of tree space.

Post-fix, the chain reaches more focused posterior regions in both
models. The blind chain still falls for AB convergent signal 55% of
the time. The aware chain catches the AB-as-ecology illusion in most
reps but isn't immune; and although CID drops below blind, P(true)
remains low because topology mixing isn't strong enough yet.

## Bugs fixed between the two runs

See `dev/red-team/findings.md` for the full ledger. Six fixes
(commits `36e9b53`, `ad25928`, `b254b98`, `37fdbe5`):

| ID  | Severity | Title                                                              |
| --- | -------- | ------------------------------------------------------------------ |
| S-1 | HIGH     | `gibbs_kprime_sweep` wrote non-eco logLik (36% of moves)          |
| S-2 | HIGH     | `block_kprime_shift` else-branch same bug                          |
| S-3 | MED      | Pre-proposal drift diagnostic used non-eco likelihood              |
| S-4 | HIGH     | `gibbs_spr` / `gibbs_subtree_swap` used non-eco partial CL (34% of moves) |
| S-5 | HIGH     | `eval_slice_target` + `slice_scalar_impl` used non-eco likelihood  |
| L-4 | MED      | R `LogPrior` rejected theta ∈ {0, 1} + 0×log(0) NaN trap           |

Plus model-spec change: pi0 prior tightened from Beta(7, 3) (ESS=10)
to Beta(75, 25) (ESS=100). Plus production safety net: periodic
from-scratch resync in `run_mcmc_batch_cpp` every 20 iter, logs drift
> 0.5 nats to stderr.

## Next moves

Priority order:

1. **Scale to 20 reps** to tighten the sd on the P(wrong) and CID
   estimates. 8 reps is enough for the qualitative sign on CID but
   not for a paper-quality figure.
2. **Topology mixing.** The `gibbs_spr` gate makes random-MH the only
   topology channel in eco mode; weights are static at small fractions.
   Either redistribute spr/nni/tbr weight when ecologyAware, or build
   an eco-aware streaming candidate evaluator (long-term).
3. **Investigate the 3/8 reps where aware still finds AB.** Are those
   reps where the ecology-edge marginals are weak (deep stem
   confusion), or where the topology mixing simply doesn't have time
   to escape AB? Per-rep posterior tree plots would tell.
4. **Address open findings.** L-1 (wEdge root-edge spec), L-2/L-3
   (simulator/model alignment, defer until kEco > 2), L-5 (ecology
   rate hard-coded), M-1/M-2 (per_ecology mode refE handling),
   M-3 (`gibbsZEvery` unused). Each is a small unit of work.

## Files

- `inst/simulations/ecology/multirep-v3-results/rep0[1-8]/summary.rds`
- `inst/simulations/ecology/multirep-v3-aggregate.R`
- `inst/hamilton/sim3-multirep-v3/{setup_project.sh,run_rep.R,sim3-multirep-v3.sh,dispatch.R}`
- Hamilton `/nobackup/pjjg18/mkp-sim3-multirep-v3/results/rep0[1-8]/`
  retains chain logs + checkpoints for resume/extension.
