# F81-Het collapse: design stub

**Branch:** `wip/f81-het-collapse-design` (off PR #2 tip — assumes JC
collapse commits b99b0ae + 8d4638e + a966d8b are merged).
**Status:** Design stub. No implementation yet.
**Why deferred:** Need a math pass before code; not a mechanical port of
the JC collapse.

## Problem

The MCMC hot path's F81-Het kernel (`pruning_f81_het_acrv_flat` in
`src/mcmc_likelihood.cpp`, dispatched at the three sites in
`cpp_partition_log_likelihood` and at the Gibbs k′ sweep persite call in
`src/mcmc.cpp:3537/3560`) does not enjoy the JC lumpability collapse
that stages 1 and 2 of this PR landed for JC. When Het is active and
`kFull > kObsMax + 1`, the kernel still pays full `nChar × kFull` stride
per node.

Per-call cost ratio if F81-Het lumpability were available:
`kFull / kEff` where `kEff = kObsMax + 1`. For kObs=2, kFull=12 this is
a potential 4× saving on every accept/reject in the Het arm.

## Why JC's argument doesn't carry over

JC is *strongly lumpable under any state partition* because Q is
`(1/k − I)·μ` — all off-diagonals equal and all rows sum to 0. Any
coarsening (including the kObs observed states + 1 lumped) produces a
new valid CTMC with the same form, and `kFull` enters only through the
analytic `p_same`/`p_diff` — so the math identity in
`src/likelihood.cpp::pruning_jc_collapsed` works without conditions on
the data.

F81-Het uses a discretised Dirichlet over per-character stationary
frequencies. Strong lumpability for an F81 chain with stationary
frequencies `π_0, …, π_{kFull−1}` requires the lumped class `U` to
satisfy `π_s = π_t` for all `s, t ∈ U` (uniform within U). In F81-Het:

- The MCMC marginalises over `nBetaCat` mixture components per character.
- Each component's frequencies are bin midpoints of `Beta(α, (k−1)α)`,
  rotated to enforce labelling symmetry.
- Within one component, the frequency assigned to a never-seen state is
  determined by the rotation index — generally distinct across the
  states in U.

So *plain* lumping of U into one column fails: the column would carry
different effective π depending on which member of U you're notionally
tracking, and the F81 transition matrix on the coarsened state space
isn't well-defined.

## Candidate approaches (sorted by expected ROI)

### A. Component-conditional lumpability

For a fixed Het component `b`, the rotation defines a permutation of
frequencies over kFull states. The unseen states in U receive some
subset of those frequencies. **If** the rotation logic for kObs < kFull
were redefined so that all states in U always receive the *same*
frequency in every rotation, then F81-Het becomes JC-like on U and
lumping is valid per-component.

Whether this is consistent with the symmetry argument that motivates
the rotation in the first place needs careful thought — the rotations
exist to avoid label-dependence of the likelihood. A "U gets the
average frequency" coarsening might be defensible because labels within
U are unobserved anyway, but this needs a proof, not a hand-wave.

**Pros:** Largest realised speedup; reuses the JC-collapse dispatch
plumbing directly.
**Cons:** Requires modifying `compute_het_bins` semantics; needs a
correctness proof against the un-coarsened mixture posterior; downstream
implications for `het_constant_site_prob` and `het_singleton_site_prob`.

### B. Block-diagonal coarsening (multiple lumped classes)

If the kFull states split into J ≤ kEff groups where within-group
frequencies are equal in every Het component (e.g., by binning the
rotation outputs), the chain is lumpable on those J classes. For the
existing rotation scheme this likely gives `J ≈ kFull / nBetaCat`,
which is a partial saving but cleaner than (A).

**Pros:** No change to the F81-Het mixture math; lumpability proof is
direct.
**Cons:** Realised speedup is far smaller than (A) — only when
`kFull / nBetaCat` is meaningfully less than kFull. For typical
nBetaCat=4 and kFull=12, you'd lump to ~3 classes per group ⇒ kEff
might be ~kObs + 3 not kObs + 1.

### C. Skip — Het arm has a separate optimisation backlog

The Het arm hasn't been profiled in the post-collapse world yet. It's
plausible that other costs (computing per-character bins; integrating
over mixture components in singleton_site_prob) dominate the pruning,
in which case state-space collapse won't move the wall clock much.

**Action:** Before designing (A) or (B), profile a Het-enabled MCMC run
(e.g., one of the mkp-eg-het arm tasks) and measure what fraction of
runtime is in `pruning_f81_het_acrv_flat` vs the surrounding F81-Het
machinery. If pruning is <30% of Het wall time, this whole project is
low-value.

## Files / dispatch sites that would change

If A or B were implemented:

- New kernel: `pruning_f81_het_acrv_flat_collapsed` mirroring the
  existing `pruning_f81_het_acrv_flat` (~600 lines including ACRV and
  the persite variant).
- New ascertainment: `het_constant_site_prob_collapsed`,
  `het_singleton_site_prob_collapsed`.
- Dispatch in `cpp_partition_log_likelihood` type==2, transformational
  allSame, transformational heterogeneous sub-batch — currently all
  branch on `useHet` BEFORE checking `useCollapse`. Add `useCollapse`
  inside the `useHet` branch too.
- Dispatch in the Gibbs sweep at `src/mcmc.cpp:3537/3560` — currently
  the collapsed path only fires when `!useHet`.
- Workspace check: collapsed asc stride `ascNeed = kEff` (or `kEff*kEff`
  for fused-asc Het) instead of `kStates` (or `kStates*kStates`) —
  already automatically satisfied since `kEff ≤ kStates`.

## What to do first

1. **Profile.** Run a Het MCMC under the existing setup, instrument
   pruning vs total. If pruning is the bottleneck, proceed; otherwise
   close this branch as "investigated, not worth it".
2. **Math.** If proceeding, sit down and prove (A) is valid — does
   constraining U to share a common frequency-per-rotation give the
   same posterior as the current scheme? Write the proof in this file
   before writing code.
3. **Implement** mirroring the JC collapse pattern; the dispatch
   infrastructure is already in place.

## Cross-references

- JC collapse commits on PR #2 branch: b99b0ae (Gibbs sweep), 8d4638e
  (cpp_partition_log_likelihood), a966d8b (dead-code removal).
- PR #2 original ("kObs < kFull JC collapse, non-flat path"): 380a356.
- F81-Het kernel implementation:
  `src/mcmc_likelihood.cpp::pruning_f81_het_acrv_flat` (~line 1290 onwards).
- Het ascertainment: `het_constant_site_prob` (~line 1684),
  `het_singleton_site_prob` (~line 1800).
- Rotation / bin computation: `compute_het_bins` (~line 1673).
