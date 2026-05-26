# Heavy tests — TreeESS root-dependence + subtree-swap detailed balance

**Lane:** D2 (mcmc-diagnostician)
**Scripts:** `tree-ess-rooting.R`, `subtree-swap-db.R`
**Author:** mcmc-diagnostician (agent run 2026-05-26)
**Targets two OPEN findings from Wave 1 (math-prover L3):**

- SWAP-001 (`src/mcmc.cpp:1874`) — `GibbsSubtreeSwap` commits without MH
- SWAP-002 (`src/mcmc.cpp:2980-2983`) — `WeightedSubtreeSwap` `logHR`
  omits partner-set ratio

and a parallel root-dependence audit of `src/tree_ess.cpp`.

---

## (a) `tree-ess-rooting.R` — TreeESS root-dependence audit

### What this tests
`MkPrime::TreeESS` (R/treeESS.R) uses
`TreeDist::RobinsonFoulds` as its default distance function. RF distance is
defined on unrooted bipartitions, so per-tree re-rooting must leave the
full distance matrix bit-identical. The downstream Frechet correlation ESS
and median pseudo-ESS are deterministic functions of that matrix (plus a
seeded subsample for `maxRows < n`), so they should be identical too.

Project memory `project_scoring_bug` reminds us that `ape::prop.part` is
root-dependent — a sibling tool that produced exactly this kind of
silent bias elsewhere. The audit here is the corresponding test for the
ESS code path.

### Theoretical basis
RF distance is an unrooted invariant: `RF(T, T') = |splits(T) ⊕ splits(T')|`
where `splits` are root-independent bipartitions
(Robinson & Foulds 1981). Re-rooting changes the parent pointer of the
root edge and reorders internal node IDs but leaves the bipartition set
identical. Therefore `dmat_unrooted == dmat_rerooted`.

### Pass criterion
1. `identical(dmat_orig, dmat_rerooted) == TRUE`, or
   `max |dmat_orig - dmat_rerooted| < 1e-12` if floating-point equality
   leaks via the C++ pipeline.
2. `|ESS_orig - ESS_rerooted| < 1e-8` for both `medianPseudoESS` and
   `frechetCorrelationESS`.

A FAIL on (1) is a real bug in `src/tree_ess.cpp` (the routine fails to
canonicalise / unroot before computing distances). A PASS on (1) but
FAIL on (2) suggests an undocumented non-determinism (a random subsample
seed, an unguarded NaN path).

### How to run
- Small-N: `Rscript dev/red-team/heavy-tests/tree-ess-rooting.R --quick`
  (n_tip=6, n_trees=200, ~30 s on a laptop)
- Full-scale: `Rscript dev/red-team/heavy-tests/tree-ess-rooting.R`
  (n_tip=8, n_trees=10000, ~3-5 min on Hamilton login node, no SLURM
  required — RAM and CPU usage modest; a SLURM script would be overkill)

### Output interpretation
Inspect `dev/red-team/heavy-tests/tree-ess-rooting-results/verdict.txt`.
The first non-blank line is `VERDICT: PASS|FAIL`. `summary.rds` contains
the full ESS pair, the dmat absolute diff, and the configuration.

### What a failure would mean
- If `dmat_identical = FALSE` → `TreeDist::RobinsonFoulds` (or the way
  `TreeESS` invokes it) is implicitly root-sensitive, e.g. via implicit
  rooting choice when reading trees. This would invalidate any ESS-based
  early-stopping criterion that compares across rooted vs unrooted
  posteriors — including the `minTreeEss` stopping path in `RunMkPrime`
  if checkpoint/resume re-roots trees.
- If `dmat_identical = TRUE` but ESS differs → the ESS routine has an
  uncontrolled RNG path; this would be a milder bug (deterministic given
  the same seed, but cross-machine differences possible). File as a new
  finding.

---

## (b) `subtree-swap-db.R` — β=0 detailed-balance test for two subtree-swap moves

### What this tests
At β = 0 the Gibbs and weighted subtree-swap moves both reduce to
candidate-uniform sampling from `{self} ∪ partners(A)`. A correct chain
under a flat topology prior must visit each of the 15 labelled binary
unrooted topologies on 5 tips equally often. The test is the same
chi-squared uniformity check that
`tests/testthat/test-tbr-detailed-balance.R` applies to NNI, SPR, TBR.

The two arms of interest:

- **GibbsSubtreeSwap (`src/mcmc.cpp:1654-1878`).** Per L3, this commits
  the chosen partner without an MH step. If `|partners(A)|` and
  `|partners_T'(B)|` differ — which they generally do, because the
  partner set depends on the local ancestor/descendant structure that a
  swap changes — the chain's stationary distribution is biased away
  from uniform.

- **WeightedSubtreeSwap (`src/mcmc.cpp:2825-3009`).** Has an MH step,
  but `logHR` at lines 2980-2983 contains only the branch-fraction
  Beta-density component. The missing
  `log(|partners_T'(B)| / |partners_T(A)|)` term means MH accept ratios
  are biased.

### Theoretical basis
Detailed balance under flat target requires
`q(T → T') = q(T' → T)`. For the Gibbs draw at β = 0:

```
q(T → T') = (1/nEdge) · 1/(|partners_T(A)| + 1)
q(T' → T) = (1/nEdge) · 1/(|partners_T'(B)| + 1)
```

These are equal iff `|partners_T(A)| = |partners_T'(B)|`. The math-prover
proof (`dev/red-team/proofs/hastings-tree-moves.md` §6) gives the
algebra; this harness gives the empirical confirmation.

### How the harness ensures only one move is applied
The driver is a pure-R MH loop that invokes one proposal function per
arm. The proposal functions for the two subtree-swap arms call the
exported helpers `MkPrime:::get_valid_swap_partners_cpp` and
`MkPrime:::swap_subtrees_cpp` — the same C++ routines invoked by
`gibbs_subtree_swap_impl` and `weighted_subtree_swap_impl`. No
`RunMkPrime` is constructed; there is no `moveWeights` to mis-configure;
no other move can fire. The control arm uses `MkPrime:::ProposeNni`
verbatim (matching the existing TBR detailed-balance test).

### β = 0 reduction (key derivation)
**Gibbs:** at β = 0, `exp(β · log L) = 1` for every candidate; the
sampler at `mcmc.cpp:1820-1841` becomes a uniform draw over
`{self, partners[0], ..., partners[K-1]}`. The chosen partner is
committed unconditionally (`mcmc.cpp:1874`; no `log α` test). The harness
emulates this by drawing uniformly and returning `logHastings = Inf`
(always-accept).

**Weighted:** at β = 0, `candW[pi][b] = 1` and `wOrig = 1`, so
`mCand[pi] = nBins` for all pi. In the C++ code (`mcmc.cpp:2917-2933`),
`sumM = 1 + K·nBins`, so `P(self) = 1/(1 + K·nBins)` and
`P(partner_i) = nBins/(1 + K·nBins)`. **Conditional on a non-self draw**,
the partner choice is uniform: `P(partner_i | move) = 1/K`. The harness
draws topology uniformly over `{self} ∪ partners` for simplicity (i.e.
`P(self) = P(partner_i) = 1/(K+1)`). The two schemes differ only in the
self-stay probability — a mixing-rate difference, not a stationary
difference: the topology-conditional partner choice (the only thing the
detailed-balance bug can affect) is `1/K` in both. The MH `logHR` from
`mcmc.cpp:2978-2983` reduces (at β = 0) to
`dbeta(f_default; α_old, β_old, log) - dbeta(f_new; α_new, β_new, log)`.
The harness reproduces this exact expression — no partner-count
correction — and submits it to a normal MH accept/reject.

A consequence: each harness iteration corresponds to ~`nBins · K / (K+1)`
"effective" C++ iterations (in the sense of real proposals), so the
harness's nominal n_iter is conservative relative to the C++ chain.

### Pass criterion
- 15 distinct topologies visited.
- Chi-squared p > 0.005 against uniform over 15 buckets.
- `max(freq) / min(freq) < 2`.

Threshold 0.005 matches the existing
`test-tbr-detailed-balance.R::.expect_uniform` helper, which has
empirically passed for NNI/SPR/TBR. Bonferroni-adjusted across the three
arms here, the strict family-wise threshold would be 0.005 / 3 ≈ 0.00167,
but any deviation strong enough to reject at the 0.005 level on 5000
samples (full-scale) is a real signal.

### How to run
- Small-N: `Rscript dev/red-team/heavy-tests/subtree-swap-db.R --quick`
  (n_iter = 50 000, thin = 25 → 2000 samples per arm; ~1-2 min total).
- Full-scale: `Rscript dev/red-team/heavy-tests/subtree-swap-db.R`
  (n_iter = 1 000 000, thin = 200 → 5000 samples per arm; ~15-30 min
  on a laptop, ~5-10 min on a Hamilton compute node).

**Hamilton submission (optional, per arm):** the test is fast enough that
SLURM is unnecessary in normal use — running it directly on a login node
is the recommended pathway. If a SLURM submission is required for batch
campaigns, a 1-CPU, 4 GB-RAM, 1 h walltime job suffices; no special
resource ask. Resumability is not relevant at this scale.

### Output interpretation
Inspect `dev/red-team/heavy-tests/subtree-swap-db-results/verdict.txt`:

```
VERDICT: PASS | FAIL | WARN

  nni_control:           ... PASS expected
  gibbs_subtree_swap:    ... FAIL ⇒ confirms SWAP-001
  weighted_subtree_swap: ... FAIL ⇒ confirms SWAP-002
```

`summary.rds` contains the full topology frequency tables, chi-squared
statistics, and timing.

### What a failure would mean

- **NNI control FAIL.** Harness bug, not an MCMC finding. The TBR
  detailed-balance test in `tests/testthat/` already passes NNI/SPR/TBR;
  if the analogue here fails the topology-key routine or the chain
  driver has drifted. Reported as `VERDICT: WARN`.

- **Gibbs swap FAIL while control passes.** Confirms SWAP-001: the chain
  is not detail-balanced because Gibbs commits without correcting for
  partner-set count asymmetry. Filing recommendation: upgrade SWAP-001
  status from OPEN to CONFIRMED in `dev/red-team/findings.md`. Fix
  candidates documented in the existing finding.

- **Weighted swap FAIL while control passes.** Confirms SWAP-002:
  `logHR` is missing the `log(|partners_T'(B)| / |partners_T(A)|)` term.
  Same upgrade recommendation. Fix outline: compute
  `|partners_T'(B)|` after the swap, add the log-ratio to `logHR`
  before the MH ratio at `mcmc.cpp:2980-2983`.

- **Both swap arms pass at full-scale.** Either (i) partner-set sizes
  happen to coincide often enough on 5 tips that the bias is undetectable,
  or (ii) the moves are detailed-balanced through a mechanism the
  proof did not see. In case (i), the bug remains real but undetectable
  at n=5 — the harness should be re-run at n=6 or n=7 with longer
  chains. In case (ii), L3's proof needs revisiting. Report
  `VERDICT: PASS` and note the open ambiguity in the finding.

---

## Joint constraints

- **Reproducibility.** All RNG seeds are hardcoded constants
  (4711, 9931, 2025, 3071, 4099). The two harnesses generate independent
  output directories so reruns are non-destructive.
- **No regular test-suite pollution.** Both scripts live under
  `dev/red-team/heavy-tests/`; nothing is added to `tests/testthat/`.
- **No commits.** Per the orchestrator brief, results are local artefacts
  only.
- **Worktree state.** Branch base at `2bc4da7` (current main HEAD as of
  the campaign). Wave 1+2 artefacts read via `git show main:...`.
