# ESJD/s as the proposal-scheduling optimand — assessment

**Status:** analysis + recommendation. Nothing implemented.
**Prompted by:** StratoBayes/StratoBayes-R-package#767, which replaces
"distance from 0.234 acceptance" with measured whitened ESJD/s as the
criterion for **mixture weights over heterogeneous kernels** (leaving 0.234
in place as the target for each proposal's **scale**).

---

## 1. The bug #767 fixes does not exist here

`.AdaptMoveWeights()` (`R/RunMkPrime.R:4427`) scores

```
logScore = log(accept_rate) + log(dim) - log(mean_cost_s)
```

softmaxed at an annealed temperature (`tStart = 2.0` → `tEnd = 0.5`).

`accept_rate` sits in the **numerator**, so a move that never accepts scores
→ 0. StratoBayes's pathology — an out-of-support draw returning early is the
*cheapest* path, and `1/sqrt(cost)` then up-weights it for failing — has no
analogue. There is nothing to fix by direct transfer.

## 2. But acceptance is uninformative here for a different reason

MkPrime also runs `.AdaptTuning()` (`R/RunMkPrime.R:4787`), which drives each
MH move **that has a step-size knob** to a fixed per-move acceptance target:

| target | tuning key | moves |
|---|---|---|
| 0.35 | `scale_*` / `int_walk_window` | `tree_length`, `kPrime`, `rate_loss`, `rate_neo`, `rate_log_sd`, `beta_scale`, `mh_logit_p` |
| 0.23–0.234 | `beta_simplex`, `*dirichlet_alpha`, `int_walk_window` | `branch_lengths`, `dirichlet_branch`, `local_dirichlet`, `block_kPrime` |
| 0.25 | `scale_joint_*` | `joint_tl_rls`, `joint_tl_rl`, `joint_tl_rn` |

Note what is **absent**: `nni`, `spr` and `pspr` appear in the `targets`
vector but map to `NA_character_` in `tuningKeys`, and `tbr` is in neither.
Topology MH moves have no step size to adapt, so their acceptance is *not*
driven anywhere — the entries are expected rates, not enforced ones.

For everything in the table above, once scale adaptation has converged
`accept_rate` is pinned at its target **by construction**, independent of how
well the move mixes. The weight score then collapses to

```
score ≈ target_p × dim / cost
```

i.e. a pure dimension-over-cost rule with no mixing signal at all. The two
adaptation loops work against each other: `.AdaptTuning` removes exactly the
variance `.AdaptMoveWeights` is trying to read.

This is #767's own diagnosis in a different costume — one number doing two
jobs. Acceptance is the right target for **scale** (Roberts, Gelman & Gilks
1997), and using it again for **weights** buys nothing once scale has
converged.

### The three move classes fail three different ways

Taking §2, §3 and §5 together, no class of MkPrime move is well served by an
acceptance-based weight criterion, but for three unrelated reasons:

| class | acceptance signal | why the weight score fails |
|---|---|---|
| tunable scalar MH (`scale`, `beta_simplex`, `dirichlet_simplex`, `int_walk`, joint) | driven to a fixed target by `.AdaptTuning` | signal *destroyed* — score → `target × dim / cost` |
| topology MH (`nni`, `spr`, `pspr`, `tbr`) | free, not adapted | signal alive but measures *frequency* of movement, never *magnitude* |
| always-accept (Gibbs, weighted, slice) | accept by construction | signal never existed — score is exactly `dim / cost` |

All three point the same way: the quantity the criterion needs is displacement.

## 3. Always-accept kernels are structurally unscoreable

The `.AdaptTuning` target table lists `NA_real_` for every Gibbs, weighted,
block and slice move — because they accept by construction:

`gibbs_kPrime`, `gibbs_p_marginal`, `gibbs_spr`, `gibbs_subtree_swap`,
`block_gibbs_branch`, `weighted_branch_lengths`, `weighted_spr`,
`weighted_subtree_swap`, `slice_rate_loss`, `slice_rate_neo`,
`slice_rate_log_sd`, `slice_tree_length`, `slice_beta_scale`,
`slice_kprime_s`, `slice_kprime_r`

These accept by construction, so `accept_rate` carries no information and the
score is effectively `dim / cost` for all of them. A
`gibbs_kPrime` sweep that leaves every k'ᵢ unchanged scores **identically** to
one that resamples the whole vector. A slice sampler that never leaves its
bracket scores identically to one that traverses the marginal. That is a large
and growing fraction of the registry, and it is precisely the set of moves the
`marginal_k` audit has been chasing.

This is a **stronger** argument for ESJD in MkPrime than anything in #767's own
context: ESJD measures displacement, which is defined for always-accept
kernels, and acceptance is not.

## 4. `.DecayLowAcceptMoves` is #767's error with the sign flipped

`R/RunMkPrime.R:4556`: any move whose batch acceptance falls below
`accept_floor = 0.02` is multiplied by `decay = 0.7` per batch, down to
`0.1 × initialWeight`.

Where #767 rewarded cheap failure, this punishes expensive success. A `tbr`
move on a large tree can genuinely sit below 2% while being the only kernel
that crosses topology islands; per acceptance it displaces far more than `nni`.
And because topology moves have no step size to adapt (§2), nothing pulls their
acceptance back up when a dataset drives it down — the registry's own expected
rate for `spr`/`pspr` is 0.10, only 5× the decay floor, so the rule is one hard
dataset away from firing on a move that is behaving exactly as designed.

Currently latent rather than active, but it is the same category error:
acceptance rate used as a mixing proxy across heterogeneous kernels.

## 5. The tuning bandit's min-ESS/s objective — the most direct transfer

`.MinEssPerSec()` (`R/RunMkPrime.R:4983`) is the objective for the Tuning-phase
bandit, which perturbs one free weight three ways and keeps the argmax. #767's
two objections to per-chain ESS both apply, and the second harder than in
StratoBayes:

1. **Multimodality.** ESS rewards a chain that stays put. Topology islands are
   *the* canonical multimodality in phylogenetics. Partly mitigated already —
   `TreeESS(..., RobinsonFoulds)` is folded into the minimum when ≥ 20 trees
   are buffered (`tuneWithTreeEss`) — but a chain confined to one island still
   reports a healthy tree-ESS within it.
2. **Noise.** The window admits from `nrow >= 10L`, and the scaled budget aims
   at only ~100 samples per candidate. ESS at n = 10–100 is mostly noise, so
   the bandit is partly optimising seed luck. Unlike #767's replacement, it has
   **no noise gate at all**: it always adopts the argmax, however unresolved.

ESJD/s has an honest standard error at that window length (a sum and a sum of
squares over draws), which is what makes #767's empirical-Bayes shrinkage gate
possible.

## 6. The blocker: topology has no free ESJD

**Do not adopt naively.** #767's whitening works because every StratoBayes
parameter is continuous. A whitened ESJD over MkPrime's continuous parameters
scores `nni`, `spr`, `tbr`, `pspr`, `gibbs_spr`, `weighted_spr`,
`gibbs_subtree_swap` at *exactly zero* — they would all be starved to their
floors, which is strictly worse than the status quo.

Two workable designs:

**(a) Single commensurate metric.** Accumulate a tree-space jump (Robinson–
Foulds or SPR distance) alongside the whitened parameter jump, paid only on
acceptance, normalised by a running mean pairwise RF to make it dimensionless.
`src/tree_ess.cpp` + `TreeDist` already provide the machinery. Risk: RF² and
whitened σ² are being added in units nothing justifies.

**(b) Two-block competition — recommended.** Topology moves compete among
themselves on tree-ESJD/s; parameter moves compete among themselves on
whitened ESJD/s; the topology/parameter budget split is pinned, or adapted on
its own schedule. This never compares RF² to whitened σ², and it composes with
the existing `pinnedWeights` machinery rather than fighting it.

Whitening details for MkPrime specifically:
- `kPrime` is integer-valued with `dim = nTrans` — squared jump in k' units is
  natural and needs no whitening beyond a per-character scale.
- `rel_br_lengths` is a simplex; whiten on the log scale or use a
  compositional metric, not raw Euclidean.
- `tree_length`, `rate_loss`, `rate_neo`, `rate_log_sd`, `beta_scale` whiten on
  the log scale against posterior SD, as #767 does.

## 7. What does *not* transfer

- **#767 item 2 (post-adapt weight updating).** MkPrime freezes weights at the
  warmup boundary deliberately: the sampling phase is a fixed-kernel MCMC.
  Keep it frozen. Weight adaptation is already the only adaptation whose
  wall-clock dependence costs bit-reproducibility; do not extend its reach.
- **#767's `∝ sqrt(rate)` diversification.** MkPrime's annealed softmax
  temperature already does that job, and more flexibly. Adopting ESJD means
  choosing one, not stacking both.
- **The M-159 cache-aware weight boost** (`src/mcmc.cpp:6136`) multiplies
  weights outside the score entirely. Any change to the criterion has to state
  how it interacts with that boost.

## 8. What transfers cleanly and cheaply

**Floors as draw counts, not budget shares.** `wMin = 0.01` /
`wMinScalar = 0.02` are fixed shares of the budget, exactly the design #767
replaced with `min(kMinRateDraws / intervalDraws, 1 / nValid)`. The point is
that the *cost* of keeping a move measurable should fall as the interval
lengthens. Same criticism, same fix, no new instrumentation needed.

**A noise gate on the bandit.** Independent of ESJD: the bandit should refuse
to adopt a candidate whose advantage is inside its own standard error.

---

## Recommendation — instrument before deciding

1. **Add per-move squared-jump accumulators as measurement only.**
   `src/mcmc.cpp:6163–6197` is a single clean instrumentation point: one move
   per iteration per chain, selected by weight, with `t0`/`t1` and the
   `accepted` flag already in hand. Add `jumpSq(ch, moveIdx)` (and its sum of
   squares) beside `moveTimeNs`, surfaced through the same return list as
   `move_time_ns`. Use design (b) — separate continuous and tree blocks.
   **Change no weights.** Log ESJD/s next to accept/propose/time.

2. **Compare the two rankings on real runs.** Sun2018 (54 taxa, 225 chars) for
   the parameter block; a dataset with genuine topology islands for the tree
   block — none of `.AGENTS/memory/validation-datasets.md` is annotated as
   multimodal, so one needs identifying or constructing first. If
   `accept × dim / cost` and ESJD/s rank the moves the same way, there is
   nothing to buy and the answer is "no change."

3. **Only if they disagree:** swap the score's numerator, retire
   `.DecayLowAcceptMoves`, replace the fixed-share floors with draw-count
   floors, and decide the softmax-vs-sqrt question.

4. **Decision mechanism:** equal-wall-clock A/B, paired by seed, on a
   multimodal target — never an ESS comparison, for the reasons in §5.

Step 1 is cheap, reversible, and answers the question with data instead of
argument. Steps 2–4 are only worth budget if step 1 shows disagreement.
