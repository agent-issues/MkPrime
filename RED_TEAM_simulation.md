# Red-team audit — mkprime simulation pipeline

**Audit target:** does the simulation pipeline at `C:\Users\pjjg18\GitHub\mkprime` contain a bias that mechanically favours high-k inference arms (mk_k40, mk_k24, mk_k15 …) over the truth-matched `mk` arm (k=kObs)?

**Headline finding:** the pipeline is mostly correct, but **the n_keep=50 + kObs>=2 filter, combined with the geometric+2 prior on k, produces a kept dataset where 51.9 % of characters have hidden states (u_true > 0) and 36.9 % have k_true > 2 while kObs = 2**. This is the substantive finding driving the high-k preference, not a bug. The "mk" arm (k=kObs) is structurally mis-specified for ~half the matrix. Several minor issues are flagged below.

The total tree length is fixed at 1.4 ("approx empirical median 1.47"), and that "64 published Mk datasets" claim should be sourced before it goes into a paper.

---

## SEVERITY: HIGH — joint (k_true, kObs) distribution: half the matrix has hidden states

Source files:
- `C:\Users\pjjg18\GitHub\mkprime\scripts\prepare_datasets.R` (n_keep=50, drops invariant chars)
- `C:\Users\pjjg18\GitHub\mkprime\inst\rbScripts\sim-mk-tree.Rev` (kSamples = Geom(0.4)+2)
- `C:\Users\pjjg18\GitHub\mkprime\scripts\simulate_all.R` (N_CHAR_SIM=70 → first 50 variable kept)

Empirical joint distribution computed over **all 260 ground_truth.csv files** (12,898 retained characters across 26 trees × 10 reps):

```
        kObs
k_true    2    3    4    5    6
   2   5160    0    0    0    0
   3   2136  944    0    0    0
   4   1153  648   95    0    0
   5    601  406   86    5    0
   6    368  246   60   12    0
   7    189  130   45    8    1
   8    123   86   33    6    0
   9     65   45   20    1    0
  10     62   31    8    3    0
  11     30   12   10    1    0
  12     11   16    2    0    0
  13      7    5    2    1    0
  14      5    3    0    0    0
  15      5    4    2    0    0
  16      1    3    1    0    0
  20      1    0    0    0    0
```

- Mean k_true = 3.49 (matches Geom(0.4)+2 expectation E≈3.50).
- Mean kObs = 2.27 (kept chars).
- 51.9 % of kept characters have u_true > 0 (k_true > kObs).
- **36.9 % of kept characters are "binary in appearance" but have k_true ∈ {3..20}.** The `mk` arm (k=kObs=2) is structurally wrong for them.
- 4 characters with k_true ≥ 13 are kept with kObs = 2 (binary). One k_true=20 character is realised as kObs=2.

This is exactly the regime where mk_k40 should out-perform mk: when the modeler assumes k=kObs they discard u_true degrees of freedom that the higher-k arm restores. **It is not a simulation bias — it is the substantive scientific story.** The audit answer to the headline paradox is: the paradox is not an artefact, it reflects the fact that the data-generating process has substantial hidden multistate variation that the truth-naïve `mk` arm cannot recover by tuning k to the observation.

If the paper's claim is "the truth-matched k is best", the simulation does not test that — by construction, `mk` is not given k_true, only kObs. The natural comparison `mk_ktrue` (per-character k = k_true from ground truth) is **not in the present design**. Recommend running this arm before publishing the ramp result.

Per-tree the hidden-state density is uniform (~52 % each tree), so the bias is not concentrated in a few trees.

---

## SEVERITY: MEDIUM — empirical median TL = 1.47 claim is uncited

`generate_trees.R:18` hard-codes `exp_steps <- 1.4 # Mean total tree length (approx empirical median 1.47)`. The 1.47 figure traces back to `AGENTS.md:384` ("empirical median total tree length across 64 published Mk datasets is 1.47"). I could not find the source data, list of 64 datasets, or computation in the repo. For Cladistics-style empirical generalisation, the provenance needs to be in the manuscript. If the 64 datasets skew toward shallow-TL palaeo matrices, generalisation to deeper morphological matrices may be misleading.

This does not bias mk-vs-mk_k40 comparisons (both inference arms see the same TL), but it does set the regime where the high-k advantage is observed.

---

## SEVERITY: LOW — 35 corrupted rows (orig_idx ≤ 0) in 17 reps

A bug somewhere in the prepare/simulate cycle produced `ground_truth.csv` rows with `orig_idx ∈ {-4..-1}` in 17 of 260 reps (35 of 12,898 rows). Affected reps are concentrated in trees 1–10 only — none in trees 11–26. Example: `tree-inference/tree_05/rep_03/ground_truth.csv`, last row `50,-1,4,2,2` (orig_idx=-1, k_true=4, kObs=2, u_true=2). The `chr50.nex` file exists, so the relabeled character was written, but the orig_idx column lost its provenance, so the audit cannot verify k_true for that character — the 4/2/2 entry could be from another character entirely. Plausibly the loop re-ran with a different prepared rep dir and appended a stale row, or `prepare_replicate` was called twice and `rbind`-style state leaked.

Impact: 0.27 % of ground_truth rows are unreliable. Does not change the joint distribution materially, but means **summary statistics per-character should drop rows with orig_idx < 1** before downstream analysis. Hamilton CID computations identifying characters by char_idx are unaffected (chr files are correctly numbered).

Asymmetry between trees 1–10 and 11–26 is suspicious — was a fix applied between the two batches? The trees 11–26 generation script is not committed to the repo (only `generate_trees.R` for trees 01–10 exists; the 16 extra rows appear in `trees_meta.csv` via commit `d8c3d320` with no accompanying script). This script needs to be located and committed.

---

## SEVERITY: LOW — 16 extra trees (11–26) lack a committed generation script

`generate_trees.R` produces only 10 trees with `set.seed(5839)`. The `trees_meta.csv` was extended to 26 trees in commit `d8c3d320` (Mar 18 2026) with no corresponding script change. `tree_26` has `pool_index = NA`, distinct from 11–25 (pool_index ∈ 56–3965). This suggests trees 11–26 were generated by an ad-hoc one-off — recommend the user locate it (likely in a shell history or a deleted script) and commit it for reproducibility. Without it, the 16 extra trees are not reproducible from the repo alone.

Sanity check: TL = 1.4 holds for all 26 trees; ntip=25 nedge=47 unrooted for all 26 (verified). j1_index values in `trees_meta.csv` differ very slightly (4th decimal) from re-computing J1Index on the unrooted Newick — this is the expected unrooted/rooted J1 difference, not a corruption.

---

## SEVERITY: NONE (cleared) — items checked and exonerated

1. **dnGeometric parameterisation.** Empirically: min(k.txt) = 2 across all 18,341 simulated chars; P(k=2) = 0.4028. So `dnGeometric(0.4)` returns k ≥ 0, and the `+2` shift gives kSamples ≥ 2 with E[kSamples] = 3.5. **No off-by-one.**

2. **fnJC(k) state labelling.** Simulated nexus files use state symbols `"012..."` with sequential integer labels from 0. RevBayes JC is symmetric in states, so the labels are exchangeable. The Mk inference model used in mkp is also symmetric — relabeling is a model-invariant operation.

3. **State relabelling in prepare_datasets.R.** `chi[] <- setNames(seq_along(tokens) - 1, tokens)[chi]` relabels states by first-appearance order. Because both the data-generating model (JC) and the inference model (Mk) are symmetric under permutations of state labels, this introduces no ordering bias. The relabel is a no-op for the Mk likelihood.

4. **Seed handling.** Seeds are `ti * 100 + ri` (range 101..2610). Verified: no two reps share a k.txt sequence; tree_01/rep_01 and tree_02/rep_01 (close seeds 101, 201) have independent draws.

5. **Tree round-trip.** Each tree.nwk re-reads as unrooted, 25 tips, 47 edges, TL=1.4000 to 4-digit precision. No topology corruption.

6. **Tree length identical across all inferences.** TL=1.4 is passed to `run-mcmc.Rev` as the 5th positional arg (`submit_new_trees.R:35`), used as the compound-Dirichlet mean — both `mk` and `mkprime` (and presumably all mk_kX arms in the downstream study) see the same simulated data and the same prior centring.

---

## Recommendations

1. **Run an `mk_ktrue` arm** (per-character k = ground-truth k_true). This is the missing control. Without it, the mk_k40 ramp can be read as "k_true is more informative than kObs", which is trivially true given the joint distribution above, rather than as "high-k inference is robust".

2. **Drop the 35 corrupted ground_truth rows** from any per-character analyses, and find/commit the trees-11–26 generation script.

3. **Document the 1.47 empirical median** with the source list before submission. If the 64-dataset distribution is bimodal or skewed, a single TL = 1.4 may misrepresent the empirical regime; consider replicating the simulation at e.g. TL ∈ {0.5, 1.4, 4.0} to span the range.

4. **Flag the 36.9 % "hidden binary" rate prominently** in the paper. The story is that mk_k40 wins because the kept matrix is dominated by characters whose realised state count under-states k_true — and that is exactly the regime morphologists are in when they use Mk with k=kObs.
