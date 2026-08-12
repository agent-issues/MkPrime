# gibbs-spr-db — exact π-invariance magnitude measurement and regression gate for `gibbs_spr`

## What this tests

That `gibbs_spr` is not π-invariant is **already confirmed analytically**, twice
and independently of this harness (`dev/red-team/proofs/esjd-allocation.md` §Q4,
plus a separate verifier that rebuilt the argument from source and ruled out
five defences). The decisive argument needs no simulation: under the flat
Dirichlet(1,…,1) branch prior (`src/mcmc.cpp:323–324`) the set **E** of trees
carrying two exactly-equal incident edges is π-null; every accepted `gibbs_spr`
lands in **E**; therefore `(πK)(E) > 0 = π(E)`.

So this harness is deliberately **not** evidence-of-existence. It exists for two
jobs:

* **(a) Magnitude.** How far is the sampled distribution from the exact target,
  in units someone can act on — what fraction of posterior samples sit on the
  π-null manifold, and how badly is the branch-length distribution distorted?
* **(b) Regression gate.** Fail clearly on current behaviour, pass on a
  corrected kernel. The `spr_fixed_surrogate` arm is what a fixed `gibbs_spr`
  must reduce to (uniform τ + Jacobian + MH); it passes the same gate, so a
  green gate after the fix is meaningful rather than vacuous.

It also settles the subsidiary dispute in `hastings-tree-moves.md` §5, which
called the move "watertight with caveats" on the grounds that the τ = 1/2
restriction "is compensated by separate branch-length moves" — that is the
harness's **Q2**.

### Which code path is driven

The `0.5 * lReg` halving is present in **all three** commit paths:

| Path | Lines | Driven here? |
|---|---|---|
| main partial-CL | `src/mcmc.cpp:1463–1466` | **yes** |
| Q-heterogeneity | `src/mcmc.cpp:1777–1780` | no |
| full-evaluation fallback (the lines §Q4 cites) | `src/mcmc.cpp:1921–1924` | no |

The fixture uses no Q-heterogeneity and `coding = "variable"`, so
`gibbs_spr_impl` (`src/mcmc.cpp:1192–1205`) routes to the **main partial-CL
path** — the one a default run actually executes. A harness that only exercised
the fallback would prove less than it looks like it does. Because all three
paths share the same halving and the same "commit without MH", a fix must be
applied to all three and re-gated here.

## Theoretical basis

### The exact target

Set **β = 0**. `do_move_cpp(..., beta = 0)` passes β straight to the move
implementations, so every Gibbs weight becomes `exp(0 · (LL − max)) = 1` and
every MH log-ratio loses its likelihood term. The target reduces to the prior,
and in `cpp_log_prior` the branch-fraction prior is

```cpp
// Dirichlet(1,...,1) = log((n-1)!) = lgamma(n)
lp += std::lgamma(static_cast<double>(relBrLengths.size()));
```

i.e. **flat on the simplex**, with no topology term. None of the kernels tested
alters `treeLength`, so we condition on `treeLength = 1`. The exact target is

> π = Uniform{labelled unrooted binary topologies on n tips} ⊗ Dirichlet(1,…,1)
> on the 2n−3 edge fractions, the two factors independent.

For n = 6: 105 topologies × the uniform 8-simplex. Crucially **π is exactly
i.i.d.-sampleable**, so the reference is *truth*, not another sampler: a uniform
labelled unrooted topology comes from sequential random-edge insertion (attaching
taxon k+1 at a uniformly chosen edge of a uniform tree on k taxa is uniform over
the (2k−3)!! trees on k+1 taxa), and Dirichlet(1,…,1) from normalised Exp(1)
draws. The script *verifies* the generator in-run — a 105-bucket chi-squared on
the reference draw plus a check that it produces zero exact ties — rather than
asserting it.

A side effect worth recording: because the branch-fraction prior is *constant*,
`gibbs_spr`'s failure to update `state->logPrior` after reallocating branch
lengths is harmless. That is one candidate mechanism ruled out, not a defence of
the move.

### The test: πKᴮ = π from an exact stationary start

A kernel K is π-invariant iff πKᴮ = π for every B. So: draw `x_r ~ π`
**exactly**, apply B sweeps of K, record `y_r`. The `y_r` are **i.i.d. draws
from πKᴮ** — independent across replicates, no burn-in, no thinning, and no
ergodicity or mixing assumption anywhere. Comparing `{T(y_r)}` against `{T(z_j)}`
for an independent exact π sample is then an ordinary two-sample test with valid
p-values. This is strictly stronger than a long thinned chain, whose nominal
p-values are invalid under autocorrelation, and it makes every arm directly
comparable at equal cost.

Bias grows and then saturates in B (πKᴮ → the kernel's own stationary law), so B
is kept moderate and the budget spent on replicates.

### Why β = 0 is decisive for a FAIL, and only for a FAIL

The branch-length map — merge `parentRow` + `sibRow` into their sum, split the
regraft edge into two halves — is **β-independent**; β touches only the candidate
weights. A kernel that is not π-invariant at β = 0 is not π-invariant, full stop.
The converse does not hold, and the magnitude numbers inherit that asymmetry:
they are distortions of the β = 0 target and need not equal the β = 1 posterior
distortion. The **π-null mass transfers most directly**, since it is set by the
acceptance rate and the branch-move dose rather than by the likelihood. A β = 1
magnitude estimate is obtainable by prior importance sampling (draw from π
exactly, weight by `exp(logLik)`) and is left as an extension; it cannot change
the existence verdict.

### Statistics — and why the topology marginal is the wrong test on its own

The defect is in the *continuous* part, and the measured results confirm the
earlier proof's concession that the topology marginal is fine. The panel is
therefore dominated by branch-fraction statistics, each with an exact reference:

| Statistic | Exact law under π |
|---|---|
| `topo` — topology bucket (bipartition key) | Uniform over (2n−5)!! buckets |
| `int_frac` — total fraction on internal edges | Beta(n−3, n) |
| `ord_min`, `ord_med`, `ord_max` — order statistics of the fraction vector | uniform spacings |
| `simpson` — Σ xᵢ² | exactly simulable |
| `bal_min` — min over *adjacent* edge pairs of \|xᵢ/(xᵢ+xⱼ) − 0.5\| | exactly simulable; > 0 a.s. |
| `n_tie` — count of adjacent edge pairs with **bit-exactly equal** lengths | **exactly 0 with probability 1** |

`n_tie` is the sharpest diagnostic and is measure-theoretic rather than
statistical: under any continuous π, P(two specified edges exactly equal) = 0.
Both a bit-exact count and a relative-difference (≤ 1e-12) count are reported so
that a last-bit artefact from `preorder_weighted_impl` cannot be mistaken for
the signal, or vice versa. (In practice the two agree exactly.)

Topology keys come from descendant-tip bitmasks over the `parent`/`child`
arrays — root-invariant by construction. `ape::prop.part` is **never** used: it
is root-dependent and has caused a real scoring bug in this project
(`project_scoring_bug`). Adjacency is likewise computed from the arrays
(rooting-invariant: two edge rows are adjacent iff they share a node).

Total-variation distance on the topology marginal is reported against a
**matched H0 baseline** — the TV of an exact π subsample of the same size —
because TV is noise-dominated at small n/nTopos and would otherwise look like a
large deviation for a perfectly correct kernel.

## Arms

Every C++ arm drives the **real compiled kernel** through
`do_move_cpp(dataPtr, statePtr, code, 0L, 0.5, 0.5, 1L, beta = 0)`. Nothing
about those moves is re-implemented in R.

| Arm | Kernel | Role | Gated? |
|---|---|---|---|
| `gibbs_spr` | code 10 alone | the accused (Q1) | yes |
| `gibbs_spr+br_1to1` | code 10, code 4 alternating | composite (Q2) | yes |
| `gibbs_spr+br_1to10` | code 10 then 10 × code 4 | composite, defence's better shot | yes |
| `gibbs_spr+br_1to30` | code 10 then 30 × code 4 | branch-mixing dose-response (full mode only) | no |
| `spr` | code 6 → `spr_proposal_impl` | **primary positive control** | yes |
| `branch_lengths` | code 4 alone | **positive control** for the branch-fraction statistics | yes |
| `tbr` | code 17 → `tbr_proposal_impl` | second correct reference, corroborating | no |
| `spr_fixed_surrogate` | R-level MH on `spr_proposal`, accept at `log U < logHastings` | **gate validity**: the corrected kernel | yes |
| `weighted_spr` | code 13 | second accused (GSPR-003) | no |
| `Rspr_c0.95 … c0.00` | R-level MH accepting at `log U < c · logHastings` | **power ladder** | no |

**Why `spr` is the control.** `spr_proposal_impl` (`src/proposals.cpp:110–166`)
draws `tau = unif_rand()` (`:121`) and returns
`logHastings = log(lRegraft) − log(lMerge)` (`:161`). The symmetric fraction
density cancels, leaving the Jacobian as the only required term — exactly
correct. `tbr_proposal_impl` (`src/tree_moves.cpp:492–494`) is a second correct
reference, kept *corroborating rather than gating* only because it is a compound
move with more surface area (path reversal, the degenerate `x == v` branch that
consumes an RNG draw and sets `lMergeSub = lSubEdge`), so a deviation there
would be more ambiguous to attribute.

**`weighted_spr` is not a control.** It is independently confirmed defective as
**GSPR-003** — it omits the required Jacobian `log(lRegraft) − log(lMerge)`.
It is included as an independent corroboration arm and **never** enters the
`gibbs_spr` verdict. (A separate observation, not adjudicated here: its Hastings
ratio at `src/mcmc.cpp:3116–3119` also substitutes the *sampled* bin's Beta
density for the bin *mixture*, and `concentration = 2 · nBins`
(`src/mcmc_state.h:46`) leaves the components heavily overlapping, so there may
be a second independent error in the same expression.)

**Why the power ladder matters.** `c = 1` is exactly correct; `c < 1` is a
*smooth, non-atomic* bias of tunable size. It measures the detection floor for
biases that do **not** leave the atomic tie fingerprint, which is what makes a
null result on the composite arms interpretable instead of merely quiet.

## Pass criterion

Stated before the code was written.

* **Continuous statistics** (`int_frac`, `ord_min`, `ord_med`, `ord_max`,
  `simpson`, `bal_min`): two-sample Kolmogorov–Smirnov against the exact
  reference. An arm fails a statistic if **p < 2.08e-4**.
* **Topology**: chi-squared against uniform over the (2n−5)!! buckets; fails if
  **p < 2.08e-4**. Full bucket occupancy is additionally required *only* when
  the expected count per bucket is ≥ 20 — below that, empty buckets are ordinary
  Poisson noise under H0 and requiring all of them spuriously fails correct
  kernels (the same trap recorded as SBC-HARNESS-001 in `subtree-swap-db.R`).
* **Ties (`n_tie`)**: fails if the fraction of replicates with at least one
  bit-exactly-equal adjacent edge pair exceeds **1e-3**, against an exact
  reference that must show **exactly 0**. This is a floating-point-coincidence
  guard, not a statistical threshold: two independent doubles coincide with
  probability ~2⁻⁵².

**Threshold rationale.** 6 gated arms × 8 statistics = 48 tests; a family-wise
α of 0.01 gives Bonferroni 0.01/48 = 2.08e-4. The same threshold applies to
accused and control arms alike — which makes the controls *easier* to pass
(correct behaviour) and detection of the accused *harder*, i.e. conservative in
the direction that matters. No tolerance was widened anywhere to obtain a
verdict; the one criterion that was changed after first execution was the
topology-occupancy rule above, which was **loosening a spurious FAIL on correct
kernels**, not loosening a FAIL on the accused (`gibbs_spr` fails on `n_tie`,
`int_frac`, `ord_min` and `bal_min` regardless).

**Verdict logic** — Q1 and Q2 are kept separate, since conflating them is how
this dispute started:

1. Generator self-check fails, **or** either primary control (`spr`,
   `branch_lengths`) fails, **or** `spr_fixed_surrogate` fails → **WARN**: the
   harness is suspect (or the gate cannot certify a fix) and nothing else is
   interpretable.
2. Otherwise `gibbs_spr` failing any gated statistic → **Q1 = FAIL**.
3. Otherwise, PASS is reported **only** if the power ladder established a
   detection floor; without one the verdict is **INCONCLUSIVE**. PASS is never
   reported without a power statement.
4. `gibbs_spr+br_*` failing → **Q2 = FAIL**. Their passing is always reported
   against the ladder's detection floor, never as a bare PASS.

Exit status is 1 on FAIL, so the script can be used directly as a CI-style gate.

## How to run

The driver is `do_move_cpp`; **`RunMkPrime()` is never called**, so the mandatory
triple-guard (`maxTime` + `setTimeLimit()` + bash `timeout`) is **vacuous** here.
There is no engine invocation, no streaming output, and no convergence criterion
that can run away — the cost is the fixed reps × sweeps product, printed at
startup. Exactness comes from tiny trees, not long chains.

* **Small-N (execution check):**
  `Rscript dev/red-team/heavy-tests/gibbs-spr-db.R --quick`
  n = 6, 600 replicates, 40 sweeps, 4 800 exact reference draws.
  **Measured: 36 s** single-threaded (i7-10700), and already decisive on Q1.
* **Headline configuration (what the committed results were produced with):**
  `Rscript dev/red-team/heavy-tests/gibbs-spr-db.R --reps 5000 --sweeps 150`
  **Estimated ~19 min** by linear scaling from the measured quick run
  (state construction 1.6 ms/replicate; 0.016–0.036 ms per C++ move;
  0.28 ms per `weighted_spr` move).
* **Full default:** `Rscript dev/red-team/heavy-tests/gibbs-spr-db.R`
  (20 000 replicates, 150 sweeps) — **estimated ~75 min** local, single-threaded.
  Adds nothing to the Q1 verdict; it buys a lower power floor and tighter
  composite magnitudes.
* Overrides: `--reps N --sweeps B --n NTIP --ref-mult M`.
* **Hamilton is not recommended and deliberately not scripted.** Nothing here
  needs it: Q1 is settled at a few hundred replicates, and Q2's power is
  replicate-limited rather than wall-clock-limited. If a detection floor an
  order of magnitude below the local one were ever wanted, the ask would be a
  16-task array (one seed per task, `--reps 200000 --sweeps 150`, ~2 h/task),
  each task writing only its per-replicate statistic matrix (~20 MB), so the
  `feedback_no_oversample` /nobackup budget is untouched.

## Output interpretation

`dev/red-team/heavy-tests/gibbs-spr-db-results/`:

* `verdict.txt` — headline `VERDICT:`, then blocks in this order: reference
  generator self-check, controls, **gate validity**, Q1, Q2, **magnitude**,
  power ladder, second accused, limitation, interpretation.
* `summary.rds` — config (including which commit path was driven), per-arm seeds
  and acceptance rates, the per-replicate statistic matrices, topology tables,
  the power ladder and every test result.

Read in this order: (1) the generator self-check and the control block — if
either is red, stop, the harness is broken; (2) `gate-on-correct-kernel` — if
that is red the gate cannot certify a fix; (3) the `n_tie` line of the
`gibbs_spr` row, which settles Q1 by itself; (4) the magnitude block; (5) the
composite rows against the power floor.

## What a failure would mean

* **`gibbs_spr` fails, controls and the corrected surrogate pass** — the gate is
  live and red, as expected on current code. The mechanism fingerprint is a
  `n_tie` rate near `nCand/(nCand+1)`: the deterministic `0.5 * lReg` split at
  `src/mcmc.cpp:1463–1466` committed with no MH step and no Jacobian. Ruled-out
  alternatives: the stale `state->logPrior` (harmless, the prior is flat in the
  fractions); the two-edge merge alone (`spr` performs the same merge and
  passes).
* **Composite arms fail** — the "compensated by separate branch-length moves"
  defence does not hold, and default runs (`gibbsSpr = TRUE`,
  `R/MkPrimeMCMC.R:356`) sample something other than the posterior.
* **A control fails** — the harness is the first suspect: non-uniform topology
  generator (checked in-run), a mis-set β, or a state not rebuilt from an exact
  π draw. A genuine defect in `spr` or `beta_simplex` would be a much larger
  finding needing its own harness.
* **`spr_fixed_surrogate` fails** — the gate is unusable and must be repaired
  before it can certify any fix; this is why it is gated rather than merely
  reported.
* **After the fix, `gibbs_spr` still fails** — check that all three commit paths
  were patched, not just the one the reviewer cited.
