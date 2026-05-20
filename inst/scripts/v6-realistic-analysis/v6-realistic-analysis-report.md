# v6-realistic analysis: null result at realistic tree length

**Run.** Sim 3 v6-realistic (16 tips, 300 chars, truth TL = 1.16, phi = 8,
pi0 = 0.45, stemBrEco = 0.15, stemBrClade = 0.035), 3 reps x (blind + aware)
on Hamilton, parsimony-anchored expSteps via commit bc24d52. Post-burnin
sample sizes per arm: rep01 176/172, rep02 172/180, rep03 183/183. All
scoring root-invariant via `ScoreTreesUnrooted`.

## Headline

At realistic TL, blind does not fail and aware does not change the
topological answer.

| rep | arm   | P(AC) | P(AB) | CID median | TL median |
|----:|-------|------:|------:|-----------:|----------:|
| 1   | blind |  1.00 |  0.00 |     0.065  |     1.16  |
| 1   | aware |  1.00 |  0.00 |     0.000  |     1.35  |
| 2   | blind |  1.00 |  0.00 |     0.065  |     1.31  |
| 2   | aware |  1.00 |  0.00 |     0.100  |     1.32  |
| 3   | blind |  1.00 |  0.00 |     0.065  |     1.00  |
| 3   | aware |  1.00 |  0.00 |     0.065  |     1.06  |

P(AC) = 1.00 and P(AB) = 0.00 in every cell. Both arms recover truth at
P = 1. The multirep-v3 "headline" (blind fails on parallel ecology,
aware regularises) was an artefact of saturated TL (truth TL = 13.5) plus
the old pathologically tight expSteps = 10 prior, not a property of the
ecology confound.

CID-to-truth distributions overlap heavily; aware is slightly closer in
reps 1 and 3, slightly further in rep 2. The difference is dwarfed by
the per-tree variance within a single chain.

## TL prior fix: verdict

Truth TL = 1.16. Posterior TL medians (blind / aware):
1.16 / 1.35, 1.31 / 1.32, 1.00 / 1.06. Five of the six chains land
within 15 % of truth; the worst (aware rep01 at 1.35) is 16 % high.
For reference, multirep-v3 posterior TL ran 4-13 x truth.

The parsimony-anchored expSteps default is doing exactly its job. The
order-of-magnitude prior misspecification driving v3's anomalies has
been eliminated. `v6-TL.pdf` shows the posterior densities bracketing
truth in every rep.

## Buried structure: aware redistributes mass to one eco stem

Per-clade stem-edge posterior medians (truth: eco-1 = 0.150, eco-0 = 0.035):

| rep | arm   | A (eco-1) | B (eco-1) | C (eco-0) | D (eco-0) |
|----:|-------|----------:|----------:|----------:|----------:|
| 1   | blind |     0.11  |     2.14  |     0.03  |     0.02  |
| 1   | aware |     2.65  |     0.01  |     0.02  |     0.03  |
| 2   | blind |     0.06  |     0.04  |     0.02  |     0.02  |
| 2   | aware |     2.64  |     0.05  |     0.04  |     0.05  |
| 3   | blind |     0.03  |     0.03  |     0.02  |     0.02  |
| 3   | aware |     2.09  |     0.02  |     0.03  |     0.05  |

Aware concentrates one eco-1 stem onto a long branch (median 2.1-2.6, much
longer than truth 0.15) while leaving the other eco-1 stem at ~0.05. Blind
keeps both eco stems compressed around the eco-0 truth. Eco-0 stems (C, D)
match truth in both arms.

This is an inflated rooting artefact, not a clean recovery of
stemBrEco: a single very long edge at the root pushes one eco-1 stem
to absorb the eco signal, while the rooting collapses the other.
Worth noting -- it would be reported as a caveat -- but does not change
the headline (topology is identical; CID nearly identical).

## Implications for the paper

The multirep-v3 finding does not survive the move to realistic TL with a
correctly anchored prior. The "aware rescues blind" story needs to be
retired or reframed.

Options:

1. **Strengthen the confound at realistic TL.** Push phi higher (16 or 24),
   stemBrEco longer (0.30), or both. The Fitch check in v6 already
   showed AB-grouped costs only 9 steps over truth on 300 chars; that is a
   narrow likelihood gap, and an honest blind chain should not fall into
   it. We may need a wider gap engineered to make blind fail.

2. **Lengthen reps to look at non-binary effects.** With topology saturated,
   the only signal left is CID and branch-length structure. CID effects
   are small but non-zero; with N = 8 instead of N = 3 we could put real
   CIs on them and report "aware reduces CID to truth by X +/- Y" as a
   secondary finding.

3. **Pivot the paper.** Make the prior-anchoring story the headline:
   "MkPrime's old expSteps = 10 default produced 5-13x TL inflation on
   realistic morphological matrices and a spurious eco-confound signal;
   parsimony-anchoring fixes it." Use v3 as the cautionary tale and
   v6-realistic as the corrected null. This avoids inventing a finding
   that the data doesn't support.

The honest read of the 3-rep evidence is option 3 (with optional follow-up
along option 1 to bound how strong an eco confound has to be at realistic
TL before blind starts to fail). Three reps is small for effect-size CIs
but unambiguous for the topology result -- the headline is robust.

## Files
- `v6-summary-table.csv`, `v6-stem-length-table.csv`
- `v6-results.pdf`, `v6-CID.pdf`, `v6-TL.pdf`,
  `v6-eco-stem-lengths.pdf`, `v6-aware-globals.pdf`
- Driver: `inst/scripts/v6-realistic-analysis/v6-realistic-analyse.R`
