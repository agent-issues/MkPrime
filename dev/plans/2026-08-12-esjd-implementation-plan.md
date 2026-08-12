# ESJD/s proposal scheduling — implementation plan

**Companion to:** `2026-08-12-esjd-proposal-scheduling.md` (the assessment).
**Status:** plan. Nothing implemented.

---

## 0. The architecture, in one sentence

**ESJD/s replaces the *within-block* score; the min-ESS/s bandit becomes the
*across-block* allocator.**

The assessment diagnosed two broken mechanisms — `.AdaptMoveWeights`'s
acceptance-based score and the tuning bandit's min-ESS/s objective. They are
broken *as currently composed*, both trying to be global optimisers. Separated
by role they are complementary:

| level | mechanism | why it belongs there |
|---|---|---|
| within block | ESJD/s | this is where the real contests are — 6 topology moves, 5 branch-length moves — and where "which kernel moves this parameter furthest per second" is the right question |
| across blocks | min-ESS/s (repaired) | MkPrime's convergence criteria (`minEss`, `minTreeEss`, `maxRhat`) are **bottleneck** objectives; sum-ESJD is a **throughput** objective and would pour budget into the 100-dimensional branch block while starving `rate_log_sd` |

This split also disposes of the commensurability problem for free: RF² is never
compared to whitened σ², because tree moves only ever compete with each other.

**v1 does not attempt optimal global allocation.** That is a research problem.
Bottleneck-across-blocks + ESJD-within-block is defensible, matches the
package's own stopping rules, and maps onto the existing `pinnedWeights`
machinery.

### Block partition

| block | moves | metric |
|---|---|---|
| topology | `nni`, `spr`, `pspr`, `tbr`, `gibbs_spr`, `gibbs_subtree_swap`, `weighted_spr`, `weighted_subtree_swap` | RF distance² (§3) |
| branch lengths | `branch_lengths`, `dirichlet_branch`, `local_dirichlet`, `weighted_branch_lengths`, `block_gibbs_branch`, `tree_length`, `slice_tree_length` | whitened log-scale jump², summed over edges |
| k′ | `kPrime`, `block_kPrime`, `gibbs_kPrime` | integer jump², summed over trans characters |
| p / hyper | `mh_logit_p`, `gibbs_p_marginal`, `slice_kprime_s`, `slice_kprime_r` | whitened jump² |
| rates | `rate_loss`, `rate_neo`, `rate_log_sd`, `beta_scale`, `slice_*`, `joint_tl_*` | whitened log-scale jump² |

Joint moves (`joint_tl_rls` etc.) touch two blocks. Assign each move to the
block of its *primary* target and let its jump be measured over that block
only — a joint move earns its keep by moving its primary parameter more
efficiently than the marginal move does.

### The auto-pin is in scope, and is a decision for S3

`R/RunMkPrime.R:914–928` auto-pins the truly always-accept types (`gibbs_p`,
`slice`, `gibbs_kprime_sweep`, `slice_kprime_hyper`, `gibbs_p_marginal`) at
their initial weights, because `accept × dim / cost` would otherwise hand them
astronomical scores. They are therefore **outside the current criterion
entirely**, and `gibbs_kPrime` needs a bespoke warmup throttle
(`.WarmupGibbsCap`, `M-171`, ~200 ms/call) and a matching restore step on top.

Under a measured criterion these could be **unpinned and scheduled**, which is
the most tangible prize in this work: it would retire a hand-tuned cap and an
assumption ("one draw per cycle is already optimal") in favour of a
measurement. But it also widens the blast radius well beyond a scoring change.

**Decision: keep the auto-pin through S1/S2, and treat unpinning as a separate
flag in S3, evaluated on its own.** S1 must nonetheless *instrument* the pinned
moves — measuring `gibbs_kPrime`'s actual ESJD/s is the evidence that decides
whether unpinning is worth attempting, and it costs nothing to collect while
the moves stay frozen.

Note the auto-pin does **not** cover `gibbs_spr`, `gibbs_subtree_swap`,
`weighted_*` or `block_gibbs_branch`; those are scored today and sit in the
topology and branch blocks above.

## 1. `dim` comes out of the score

Non-negotiable and easy to miss in review. ESJD sums over coordinates, so a
`block_gibbs_branch` move touching `nEdge` coordinates is credited by `nEdge`
terms **automatically**. Retaining `dimAdj` (`R/RunMkPrime.R:4464`) would
double-count dimension. The new within-block score is

```
rate_m = (Σ over interval of squared whitened jump) / (moveTimeNs_m / 1e9)
```

with a rejected proposal contributing exactly **0** to the numerator while its
time still counts in the denominator.

## 2. Instrumentation — snapshot the target block, not the state

`src/mcmc.cpp:6163–6197` is a single clean point: one weighted move per
iteration per chain, `t0`/`t1` and the `accepted` flag already in hand.

Do **not** snapshot whole state. Each move already declares a `target`
(`R/RunMkPrime.R:3513` onward), so a static move→block table lets us snapshot
only the affected coordinates. Costs:

| block | snapshot cost |
|---|---|
| rates / p / hyper | 1–3 doubles |
| k′ | `nTrans` ints (225 on Sun2018) |
| branch lengths | `nEdge` doubles (~105) |
| topology | bipartition set or parent/child arrays |

All negligible against a likelihood evaluation, but the table keeps cheap
scalar moves genuinely cheap. **Care needed:** topology moves also change
branch lengths (SPR reattachment), `tree_length` scales all edges, and the
`gibbs_*` moves change topology *and* branches — the block table must reflect
what each move actually mutates, not what its `target` field nominally says.

New accumulators, mirroring `moveTimeNs` exactly:
`jumpSq(ch, moveIdx)` and `jumpSqSq(ch, moveIdx)` (the sum of squares is what
gives the rate a standard error, which the noise gate needs), returned
alongside `move_time_ns`.

**Whitening scales need no new architecture.** `chain_rhos` is already
estimated in R from the warmup snapshot buffer (`.EstimateJointRhos`,
`.AccumulateRhoSnapshot`) and pushed into C++ per batch as `jointRhos`. Mirror
that mechanism for per-coordinate scales. A coordinate with no usable scale
contributes 0 — reads as zero jump, which closes the gate rather than
inventing a rate.

## 3. The tree metric: RF, decided

**Robinson–Foulds.** Settled on MS's `treess`/`TreeDist` experience: for tree
ESS, RF outperforms clustering-information distance because its behaviour
*corresponds to the nature of tree autocorrelation* — the very property an
ESJD numerator needs. CID's better-behaved geometry is an advantage for
summarising tree-to-tree *dissimilarity*, not for measuring how far a Markov
chain has travelled in tree space. `TreeESS()` already defaults to
`RobinsonFoulds` for the same reason.

So the earlier saturation worry is deprioritised, but note what it implies
about signal, which S1 should confirm rather than assume:

- For `nni`, an accepted move changes exactly one bipartition, so the RF jump
  is deterministic (= 2) and carries no information beyond acceptance.
- For `spr`/`tbr`/`pspr`, the jump scales with reattachment distance — **this
  is the signal worth measuring**, and precisely what `accept × dim / cost`
  throws away.
- `gibbs_spr` can "accept" the current attachment point, i.e. move nowhere.
  ESJD scores that honestly at 0; acceptance cannot see it at all.

**S1 records both metrics anyway** — it is one extra accumulator over a run
that is being done regardless, and it converts "RF suffices here too" from an
expectation into a measurement. RF is the default and the tie-break; CID is
recorded only as a cross-check, and would only be revisited if RF's per-move
distribution turns out degenerate for `tbr` specifically.

## 4. Reproducibility — offer an ESJD-per-iteration mode

Wall clock is already in the denominator, so ESJD adds no *new*
irreproducibility (`.AdaptMoveWeights` has traded bit-reproducibility for
wall-clock cost normalisation since inception — a deliberate call, not to be
"fixed"). But scoring on **ESJD per iteration** rather than per second is
fully deterministic, and that gives the test suite something it currently
lacks: schedule assertions that can be pinned exactly. Ship it as an internal
flag used by tests and SBC, not as a user-facing option.

## 5. Correctness argument — and why SBC is not needed

Every registered move is a complete π-invariant kernel, so for any fixed
weights `π(Σ wᵢKᵢ) = Σ wᵢ(π Kᵢ) = π`. Changing the mixture cannot break
invariance. MkPrime freezes weights at the warmup boundary, so the sampling
phase remains a **fixed-kernel** MCMC and not even the diminishing-adaptation
argument is required.

The correctness risk is therefore confined to two things, both testable
cheaply and neither needing SBC:

1. **Instrumentation must not perturb the chain.** Same seed, jump measurement
   on vs off → **bit-identical chains**. This is the single most valuable test
   in the plan. It catches an accidental RNG draw inside the jump computation,
   and it catches the class of bug MkPrime has already been bitten by twice
   (marginal_k Bug B: an evaluator reading `state->parent/child` instead of its
   passed arguments). Jump measurement reads state at exactly the hazardous
   moment — before and after a mutation — so this is the highest-risk diff in
   the whole plan.
2. **The freeze boundary must not move.** Assert it.

Run the existing fast SBC as cheap regression insurance (weights are frozen
post-warmup, so `p` AD ≈ 0.86 should be unchanged); if it moves, something in
(1) is wrong. A full Hamilton SBC campaign is **not** warranted for a weight
rule change.

## 6. Staging

Harness before instrumentation. This is the non-obvious ordering and it is the
right one: S3 cannot be evaluated without a trustworthy instrument, and
building the harness first forces the multimodal dataset to exist.

### S0 — measuring instrument (no package changes)

`dev/benchmarks/proposal-schedule-ab.R`.

- **Multimodal target, constructed not found.** Nothing in
  `.AGENTS/memory/validation-datasets.md` is annotated as topology-multimodal.
  Simulate ~60 characters on tree A and ~60 on tree B, 16–20 taxa, A and B
  separated by a distant SPR. Gives a genuinely bimodal topology posterior
  *with ground truth on the island count* — better than a hard empirical
  matrix, and reusable. Add it to `validation-datasets.md`.
- **Equal wall clock, not equal iterations.** One shared `nIter` sized from the
  fastest arm's pilot, every run truncated to the shortest run's elapsed time.
  Adaptation pinned in **absolute** iterations so a cheaper arm does not also
  get a longer warmup.
- **Metrics, in order:** tree pseudo-ESS/s and min-ESS/s computed from
  **pooled independent runs** (never per-chain — a chain confined to one island
  reports healthy ESS, tree-ESS included); island coverage and switch rate;
  rank-normalised split-R̂ across independent runs at several prefixes.
  Tempering rungs are one run's ladder, not replicates.
- **A/A pass first.** #767's A/A found 5/5 seed pairs favoured whichever arm ran
  *second*, mean +9.2%, on a configuration identical to itself. So: pair by
  seed **and alternate arm order across pairs**, or the harness inherits that
  bias. Report the noise floor before any A/B number is quoted.

**Gate:** A/A returns a null with a quantified noise floor.

### S1 — C++ accumulators, measurement only

Jump accumulators + block table + whitening-scale push. **No weight rule
changes.** ESJD/s logged next to accept/propose/time.

**Gate:** bit-identical chains with instrumentation on vs off; 100% covr line
coverage on new R and C++ (house rule — budget tests for every accumulator
path, including the always-accept ones); full testthat FAIL = 0.

### S2 — the empirical question

On Sun2018 and the S0 target: does ESJD/s rank moves differently from
`accept × dim / cost`? Confirm RF's per-move distribution is non-degenerate
for `tbr` (§3).

**Gate — this is the real decision point.** If the rankings agree, stop:
the answer is "no change", and S1's diagnostics are still worth keeping. Only
disagreement licenses S3.

### S3 — new scoring rule, flag-gated, default off

Within-block ESJD/s; `dim` removed; `.DecayLowAcceptMoves` retired under the
same flag; fixed-share floors (`wMin` 0.01 / `wMinScalar` 0.02) replaced by
draw-count floors `min(kMinDraws / intervalDraws, 1/nValid)`; noise gate added
to both the within-block score and the bandit. Softmax-vs-`sqrt` resolved in
favour of the existing annealed softmax — do not stack both diversifiers.
Bandit window enlarged past the current `nrow >= 10L` admission.

### S4 — A/B and flip

Equal-wall-clock A/B on the S0 harness, paired seeds with alternating order.
Flip the default only on a signed win outside S0's noise floor.

**Note:** 5 seed pairs × 2 arms at equal wall clock is likely Hamilton
material rather than local, per the long-runs rule. Local is for smokes and
the S1/S2 diagnostics.

## 7. Fable-tier review checkpoints

Three specialists already exist for exactly these hazards. Named checkpoints,
each vetoable:

| when | agent | question |
|---|---|---|
| after this plan, before S1 | `math-prover` | Is bottleneck-across-blocks + ESJD-within-block coherent? Is the fixed-weight invariance argument airtight given the freeze boundary, and does the joint-move primary-block assignment introduce a bias? |
| after S1 | `numerical-auditor` | The C++ instrumentation diff. State-read ordering around mutation, cache-coherence interaction, whitening underflow when a posterior SD → 0, the M-159 cache-boost interaction. **Highest-value review in the plan** — this is the marginal_k Bug B failure mode. |
| at S0 | `mcmc-diagnostician` | Harness design: is the constructed bimodal target genuinely bimodal, are the pooled-run metrics sound, does the A/A design actually neutralise the order effect? |

If any run with `isolation: worktree`, check `git merge-base` when they
return — they inherit main HEAD, not this worktree.

## 8. Open items

- ~~Tree metric~~ — closed: RF (§3). CID recorded as a cross-check only.
- The M-159 cache-aware weight boost (`src/mcmc.cpp:6136`) multiplies weights
  *outside* the score. Its interaction with a measured criterion needs stating:
  the boost changes which moves are drawn, hence what gets measured.
- Joint-move block assignment (§0) is a judgement call flagged for
  `math-prover`.
