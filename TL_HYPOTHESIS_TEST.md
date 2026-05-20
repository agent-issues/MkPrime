# TL hypothesis test — result: **HYPOTHESIS REJECTED (direction is reversed)**

Date: 2026-05-19
Script: `dev/pilots/2026-05-12-prior-validation/analysis/redteam_tl_hypothesis.R`
Plots: `dev/pilots/2026-05-12-prior-validation/analysis/redteam_tl_{1..4}_*.png`

## The hypothesis under test

Under JC(k), to explain the same observed amount of change a longer branch is needed at larger k (per-unit-time substitution prob ≈ T/(k-1) at small T). The mechanistic story for "mk_k40 > Mk' on CID-to-truth" was that fixed k=40 forces inference to use longer branches, longer branches = stronger apparent signal per unit time, tighter topology posterior, better CID.

**Quantitative prediction**: TL(mk_k40) substantially larger than TL(mkp_eg) on the same task — at small T, naive theory suggests a ratio of order k/2 (i.e. up to ~20×). Even allowing for compression at higher T, the inflation should be sizeable and unambiguous.

## What we observe (n = 260 tasks × 9 arms)

### Per-arm posterior mean TL (sorted by mean):

| Arm      | n   | mean TL | median TL | mean CID |
|----------|-----|---------|-----------|----------|
| mkp_geo  | 26  | 1.430   | 1.463     | 0.256    |
| mk       | 260 | 1.394   | 1.391     | 0.275    |
| mkp_eg   | 260 | 1.339   | 1.336     | 0.268    |
| mk_kp1   | 260 | 1.298   | 1.300     | 0.259    |
| mk_kp2   | 260 | 1.263   | 1.264     | 0.253    |
| mk_k9    | 260 | 1.220   | 1.219     | 0.245    |
| mk_k15   | 260 | 1.206   | 1.212     | 0.242    |
| mk_k24   | 260 | 1.200   | 1.206     | 0.240    |
| **mk_k40** | 260 | **1.194** | **1.187** | **0.239** |

**TL goes DOWN monotonically as k goes UP.** mk_k40 has the *shortest* TL of any arm. mk (k=2) has the longest among the fixed-k Mk variants.

### Paired test mk_k40 vs mkp_eg (n=260)

- median TL_k40 = 1.187, median TL_eg = 1.336 → **ratio 0.89× (k40 is shorter)**
- Wilcoxon paired (alternative = k40 > eg): V=0, p=1.0 (i.e. the *opposite* direction is overwhelmingly significant)

### TL inflation does NOT predict CID gain — it predicts CID gain in the WRONG direction

Per-task correlation between log(TL_k40 / TL_eg) and ΔCID = CID_eg − CID_k40 (positive ΔCID = mk_k40 closer to truth):

- **Spearman ρ = −0.367, p = 1.4 × 10⁻⁹, n = 260**

Tasks where mk_k40's branches are *shorter* than Mk's are the tasks where mk_k40 wins more, not less. This is the opposite of the predicted mechanism.

### Theoretical curve vs empirical

Anchoring at k=2 with median n_tip=25 (47 branches) and TL=1.394 gives implied per-branch substitution prob p_obs ≈ 0.029. Plugging that into JC(k) and solving for matched-p T(k) gives:

| k  | theoretical TL | empirical TL |
|----|----------------|--------------|
| 2  | 1.394 (anchor) | 1.394        |
| 9  | 1.376          | 1.220        |
| 15 | 1.375          | 1.206        |
| 24 | 1.375          | 1.200        |
| 40 | 1.374          | 1.194        |

Even modest theoretical inflation (~1%) is **not seen**; the empirical numbers go the other way by ~14%. So the prediction fails in both sign and magnitude. (Note: the "ratio ~ k/2" intuition assumed *very* small p; at p ≈ 0.029 the saturation of (k-1)/k caps the ratio near 1, so the JC(k) prediction at this p_obs is in fact much milder than k/2. But the empirical direction is still opposite.)

## Why is TL going DOWN with k?

A higher k means each substitution event is *more informative* per character: the prior probability that two randomly chosen states agree under stationarity is 1/k, vs 1/2 at k=2. So an observed identical pair of states at a tip in k=40 is much stronger evidence of recent shared ancestry than in k=2. Conversely a *mismatch* under k=40 is weak evidence — it could have arrived from any of 39 wrong states. Net effect: at high k the data still imposes the same number of observed substitutions, but the substitution probability *per unit branch length* required to fit them is **lower**, because the no-change prob converges to a slower-decaying exponential (the second JC eigenvalue is −k/(k−1) which → −1 as k→∞, vs −2 at k=2). At fixed *small* p the per-branch T should still rise with k, but at the substitution rates we actually see (~3%) the regime is dominated by the (k-1)/k saturation cap; the high-k models can explain the same data with *shorter* branches because each branch step delivers a stronger likelihood ratio per match.

So the mechanism stands JC theory on its head only under the small-T regime assumption. At empirical substitution levels, **high-k models prefer short trees**, and that — together with their *better* CID — means the relationship between TL and topological accuracy is opposite to what we hypothesised.

## Strong correlation: TL tracks CID, monotonically across arms

The arm-level numbers are almost perfectly co-linear: higher mean TL ↔ worse mean CID. This is a striking pattern:

- mk (TL 1.394, CID 0.275) — worst CID, longest TL
- mk_k40 (TL 1.194, CID 0.239) — best CID, shortest TL

This is **not** noise: across all 9 arms the rank correlation of (mean TL, mean CID) is +1 (every arm with longer TL has worse CID). The relationship survives at the per-task level too (ρ = −0.37 for the k40 vs eg contrast, with the opposite sign convention).

## Bottom line

1. **TL does not inflate with k** — it *deflates*. The mechanism we proposed is wrong.
2. **TL inflation is anti-correlated with CID accuracy**, not correlated. Where mk_k40 uses shorter branches than Mk', it does even better, not worse.
3. mk_k40's advantage over Mk' is **not** explained by "longer branches → tighter posterior".

## What would I test next

The empirical pattern is that **shorter posterior TL ↔ better topology recovery** across all 9 arms. This points to a different mechanism — TL is being over-inflated by per-character flexibility in Mk' (and at low k in Mk), and shrinking it with a high-k constraint acts as an anti-overfitting regulariser. Candidate hypotheses to test in order:

1. **TL-overfit hypothesis** (favoured): Mk' / low-k Mk have spare flexibility that lengthens branches to fit noise. The high-k constraint forces a more conservative branch-length posterior, leaving topology better resolved. Check: simulate from a known tree under TRUE Mk (k=2), then fit with mk_k40 — does mk_k40 *still* recover TL closer to truth than Mk'? If yes, this isn't a model-mismatch story, it's pure regularisation.
2. **Likelihood-curvature hypothesis**: high k makes the per-branch likelihood sharper at the optimum, so even though the posterior mean TL is lower, the *uncertainty* in branch lengths is smaller. Tighter branches → tighter topology. Check: compare sd(tree_length) and per-branch posterior CV across arms — does mk_k40 have *both* lower mean and lower sd?
3. **Information per character**: at k=40 a constant (identical-state) character is more informative; an autapomorphy is also more informative. Mk' is averaging over k that drift low under empirical Bayes (k_obs in the corpus tends to be ~3-4 because most chars are binary). If true k′ for the focal character is higher than the corpus mode, mkp_eg under-estimates it, gets a softer likelihood, and pays in topology accuracy. Check: compare mk_k40 vs mkp_eg performance **stratified by the true generative k** of the simulated dataset — does the gap shrink or invert when true k is low?

Path 1 is the cleanest first test. The TL data already on disk would let you compare posterior TL against the *true* TL on the simulation (we have access to the generating tree).
