# M-157: Fused ascertainment correction

## Problem

`constant_site_prob_jc` accounts for 6.6% of total CPU (VTune S-PROF4).
It runs a separate full tree traversal that recomputes the exact same
`exp()` per edge per rate category as the main pruning pass. The MkN
and Het ascertainment paths have the same redundancy.

## Findings during investigation

**Bug in `het_constant_site_prob`:** Uses a per-edge P_ss(t) product
shortcut that is only correct for star trees. For general binary trees,
the off-diagonal CL elements at internal nodes contribute to the
likelihood at ancestor nodes, so a proper Felsenstein pruning traversal
is required. The non-Het functions (`constant_site_prob_jc`,
`constant_site_prob_mkn`) correctly do full Felsenstein traversal. No
tests exist for `het_constant_site_prob`. This will be fixed as part of
the fused implementation (proper CL propagation inherently correct).

## Design

### Core idea

Fuse the constant-site CL propagation into the existing flat pruning
functions as a side-channel. In each edge loop iteration, after computing
the per-edge transition parameters (which require `exp()`), also
propagate a small set of constant-site "pseudo-characters" using the
same parameters. This eliminates the redundant `exp()` calls and loop
overhead of the standalone ascertainment functions.

### JC optimisation: 1 pseudo-character

Under JC symmetry, all k constant-site patterns yield the same site
likelihood. Instead of propagating k pseudo-characters (as
`constant_site_prob_jc` does), propagate **one** pseudo-character (all
tips in state 0) and multiply the result by k at the end.

Side buffer per node: k doubles (the k CL states for the single
pseudo-character). Total: `(maxNode + 1) × k` doubles. For the typical
case (k=2, 54 tips): 107 × 2 = 214 doubles = 1.7 KB.

### MkN: 2 pseudo-characters

Rates are asymmetric, so the all-0 and all-1 constant-site patterns
have different probabilities. Two pseudo-characters with 2 states each.
Side buffer: `(maxNode + 1) × 4` doubles.

### Het (F81): k pseudo-characters

Each constant-state pattern under F81 gives a different likelihood
because the frequency vector is non-uniform. Use k pseudo-characters
with k states each. Buffer: `(maxNode + 1) × k²` doubles. (For k=2
with nRot=1, this is `(maxNode + 1) × 4` doubles; for k=5, 25 per
node.)

The fused implementation replaces the buggy analytical shortcut with
proper Felsenstein pruning, fixing the correctness issue.

### Buffer management

Allocate the side buffer as local `std::vector<double>` and
`std::vector<uint8_t>` inside each flat pruning function, gated on
`outConstProb != nullptr`. The buffers are small (< 10 KB typical) and
allocated once per partition per likelihood eval. The current standalone
functions already allocate similar-sized vectors, so this is no worse.

## Call sites affected

### Flat pruning functions to modify (mcmc_likelihood.cpp)

All are `static` (file-scope), so signature changes have no ABI impact.

| Function | Pseudo-chars | New parameter |
|----------|-------------|---------------|
| `pruning_jc_flat` | 1 (JC sym) | `double* outConstProb = nullptr` |
| `pruning_jc_acrv_flat` | 1 (JC sym) | `double* outConstProb = nullptr` |
| `pruning_mkn_flat` | 2 | `double* outConstProb = nullptr` |
| `pruning_mkn_acrv_flat` | 2 | `double* outConstProb = nullptr` |
| `pruning_f81_het_acrv_flat` | k | `double* outConstProb = nullptr` |

When `outConstProb` is nullptr, the function behaves exactly as before
(zero overhead — the side-buffer allocation and pseudo-char loops are
entirely skipped).

### `cpp_partition_log_likelihood` (mcmc_likelihood.cpp)

For each partition type, when `coding != 0`:
1. Pass `&constSiteProb` to the flat pruning call.
2. Use the returned value for the Lewis correction.
3. Remove the separate `constant_site_prob_jc` / `constant_site_prob_mkn`
   / `het_constant_site_prob` calls.

The fallback non-workspace path (heap-allocated temp buffer) also gets
the fused parameter.

### Partial CL path (node_cl_cache.h)

The partial CL path in `apply_dirty_set_and_score` calls
`constant_site_prob_jc` / `constant_site_prob_mkn` separately (lines
609, 615, 628). These are NOT fused by this task because the partial CL
path uses a completely different evaluation strategy (dirty-set
incremental update rather than full traversal). The standalone
ascertainment functions remain available for this path.

**Future work (M-161):** per-partition cache validity could eventually
cache the constant-site probability alongside the partition likelihood,
eliminating redundant ascertainment calls in the partial CL path too.

### Gibbs kPrime sweep (gibbs_partial_cl.h, mcmc.cpp)

The Gibbs paths already have their own `evaluate_const_prob` and
`constProbAccum` mechanism. **Not touched** (Agent D territory).

### Standalone functions (ascertainment.cpp)

`constant_site_prob_jc`, `singleton_site_prob_jc`,
`constant_site_prob_mkn`, `singleton_site_prob_mkn` remain as-is. They
are:
- Exported to R (`// [[Rcpp::export]]`) for testing
- Used by the partial CL path
- Used by `const_site_prob_for_k()` for the Gibbs sweep
- Needed for future singleton correction (Phase 7)

`het_constant_site_prob` and `het_singleton_site_prob` (static helpers in
mcmc_likelihood.cpp) will be kept but marked deprecated once the fused
path replaces their call sites. A correctness test will be added that
validates the fused result against the proper Felsenstein traversal.

## Implementation steps

1. **Add fused side-channel to `pruning_jc_acrv_flat`** (the highest-
   traffic JC path). Test against `constant_site_prob_jc`.

2. **Add fused side-channel to `pruning_jc_flat`** (non-ACRV JC). Test.

3. **Add fused side-channel to `pruning_mkn_acrv_flat`** and
   **`pruning_mkn_flat`**. Test against `constant_site_prob_mkn`.

4. **Add fused side-channel to `pruning_f81_het_acrv_flat`**. Test
   against a reference value computed from a proper Felsenstein traversal
   on the constant-site pseudo-characters (NOT against the buggy
   `het_constant_site_prob`). This also fixes the Het bug.

5. **Wire up `cpp_partition_log_likelihood`** to use fused ascertainment
   for all partition types. Remove standalone ascertainment calls from
   the full-evaluation path.

6. **Full test suite** — verify zero numeric drift.

7. **Benchmark** — microbench before/after on Sun2018.

## What NOT to change

- `_persite` functions (Agent D territory)
- `gibbs_partial_cl.h` (Agent D territory)
- Standalone ascertainment functions in `ascertainment.cpp`
- Singleton correction paths (Phase 7)
- Workspace allocation in `mcmc.cpp` (no stride changes needed)
- Partial CL path in `node_cl_cache.h`

## Expected savings

- `constant_site_prob_jc` was 6.6% of CPU; ~83% of that is redundant
  `exp()` calls. Fusing eliminates those → **~5.5% total wall-time
  reduction** for the full evaluation path.
- Additional savings from eliminating the CL propagation loop overhead
  and cache thrashing of the second traversal → estimated **~6% total**.
- The partial CL path (NNI, beta_simplex, Dirichlet) still calls the
  standalone functions, so the savings apply only to full-evaluation
  moves (SPR, TBR, pSPR, parameter moves).

## Risk assessment

**Low risk.** The fused code path mirrors the existing standalone
functions exactly — same exp_term, same CL propagation arithmetic, same
root accumulation. The `outConstProb == nullptr` default preserves
existing behaviour. The existing test suite (608+ tests covering
likelihood, ascertainment, rate matrix paths) validates correctness.
