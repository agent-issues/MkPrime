# GSPR-001 fix design — and why "targeted fix vs unification" is a false choice

**Status:** design settled; `weighted_spr` Jacobian applied, `gibbs_spr` not yet.
**Gate:** `Rscript dev/red-team/heavy-tests/gibbs-spr-db.R --quick` (23–36 s, exit 1 on FAIL).

---

## 1. What the correct kernel must contain

`gibbs_spr` merges `(l_parent, l_sib) → lMerge` and splits `lReg → (f·lReg,
(1−f)·lReg)`. Any correct version needs three things the current one lacks:

1. **A drawn split fraction.** A deterministic `f = 0.5` maps a positive-measure
   set into a π-null set. No Hastings term can repair that; the fraction must
   have a density.
2. **The SPR Jacobian** `log(lReg) − log(lMerge)`. Derived from the bijection
   `(l_parent, l_sib, lReg, f_new) ↔ (lMerge, a, b, f_old)` with
   `f_old = l_parent / lMerge`; block-diagonal with
   `|∂(lMerge, f_old)/∂(l_parent, l_sib)| = 1/lMerge` and
   `|∂(a, b)/∂(lReg, f_new)| = lReg`. Present in `spr_proposal_impl`
   (`src/proposals.cpp:161`), `tbr_proposal_impl` (`src/tree_moves.cpp:492`) and
   pspr (`src/mcmc.cpp:3871`); absent from both `gibbs_spr` and `weighted_spr`.
3. **A real accept/reject**, with the proposal-density ratio for the
   likelihood-weighted candidate choice.

## 2. The reverse-proposal probability is the hard part — and it is already solved

Computing `q(x|y)` naively means re-evaluating the candidate weights from the
proposed state: a second full set of likelihood evaluations, doubling the cost.

`weighted_spr_impl` avoids that with an argument worth stating explicitly,
because it is what makes a correct `gibbs_spr` affordable:

- Both directions prune the **same** subtree at the same edge, so the residual
  backbone and the candidate edge set are identical (Lemma S1 of
  `hastings-tree-moves.md`).
- "Self" from `y`'s perspective **is** the chosen candidate; the original
  position is one of `y`'s candidates.
- The merge/split is local, so every other edge length is untouched.

So the *set of enumerated configurations* is the same from `x` and from `y`.
The shared normaliser `sumM` therefore cancels, and the reverse weight of `x`'s
own configuration is exactly the `selfW[oldBin]` already computed in the forward
pass — which is why `weighted_spr` evaluates self at bin midpoints of `lMerge`
(`src/mcmc.cpp:2962–2963`) rather than only at its current fraction.

## 3. Base the fix on `spr_proposal_impl`, not on `weighted_spr`

**An earlier draft of this section decided the opposite — to route `gibbs_spr`
through `weighted_spr`'s machinery — on the argument that assembling §1's three
requirements using §2's cancellation reproduces `weighted_spr` anyway, so
unification is the smaller diff. Measurement refuted that.**

Adding the Jacobian to `weighted_spr` improved **every** statistic in the gate
harness, by 2× to 19×:

| statistic | relBias before | relBias after |
|---|---|---|
| `int_frac` | −0.293 | **−0.015** |
| `ord_min`  | −0.744 | −0.379 |
| `ord_med`  | −0.346 | −0.093 |
| `ord_max`  | +0.354 | +0.092 |
| `simpson`  | +0.460 | +0.118 |
| `bal_min`  | +0.634 | +0.238 |

(acceptance 0.649 → 0.399, consistent with a correctly tightened ratio; `n` also
differs, 5000 full vs 600 quick, so compare `relBias`, not KS p-values.)

That settles the Jacobian: necessary, correctly signed, and worth keeping. It
also shows it is **not sufficient** — residual biases of 9–38 % mean
`weighted_spr` carries at least one further defect, most plausibly in the
`selfW`/`sumM` cancellation of §2, since "self" is enumerated over bins of
`lMerge` while candidates are enumerated over bins of `lReg`, and the claim that
both directions enumerate the same configuration set with the same weights has
not actually been verified.

**So `weighted_spr` is not a safe base.** Building `gibbs_spr` on it would
inherit an unidentified defect into a default-on kernel — trading a known bug
for an unknown one.

**Revised decision: build the `gibbs_spr` fix on `spr_proposal_impl`**
(`src/proposals.cpp:110–166`), which this same harness certifies clean (KS p
0.093–0.962 across all six statistics, zero π-null mass), with the
likelihood-weighted candidate selection layered on top and its reverse selection
probability computed honestly rather than assumed to cancel.

Cost: that may require a second evaluation pass to obtain `q(x|y)`, where
`weighted_spr` assumed the cancellation for free. Accept it. Correctness first,
and the ESJD instrumentation folded into this pass is exactly what will measure
whether the extra cost matters — which is a better position than inheriting a
cancellation nobody has checked.

Consequence for `weighted_spr`: it stays default-off and GSPR-003 stays open
with its severity unchanged, now recording a *second* defect beyond the
Jacobian. Verifying the §2 cancellation is the next step there, but it is not on
this fix's critical path.

## 3a. The construction — and why it needs no second evaluation pass

§3 accepted that computing `q(x|y)` honestly might cost a second pass. It does
not. The reason is a structural fact neither existing routine exploits
explicitly:

> **`x` and `y` have a common pruned intermediate.** Pruning subtree `S` at the
> same edge from `x` and from `y` yields the *identical* tree `R` — same
> topology, same edge lengths, with the vacated pair merged. Regrafting is the
> inverse of pruning, so both endpoints prune to one shared `R`.

Hence the candidate edge set is the edges of `R`, identical from both
directions, and `x`'s own position is simply one of those edges (the `lMerge`
one). So a selection distribution that is **a function of `R` alone** has the
same normaliser in both directions, and it cancels exactly — no assumption
required, unlike §2's unverified appeal to bin-grid symmetry.

The one condition: the candidate weights must not depend on the drawn fraction.
Write `w_e` for the likelihood of regrafting `S` at edge `e` of `R` **at a fixed
reference fraction of ½**. Then with `τ ~ U(0,1)` drawn independently:

```
q(y|x) = p(τ_y) · w_{e_y} / Σ_e w_e
q(x|y) = p(τ_x) · w_{e_x} / Σ_e w_e
```

`Σ_e w_e` cancels because the weights are shared; `p(τ_y) = p(τ_x) = 1` because
`τ` is uniform. What survives is

```
logHR = log w_{e_x} − log w_{e_y} + log(lReg) − log(lMerge)
```

plus the usual `β·(logLik_y − logLik_x)` and prior ratio in the MH step.

**Had the weights been evaluated at the drawn `τ` instead of a fixed ½, the
normalisers would be `Σ_e w_e(τ_y)` versus `Σ_e w_e(τ_x)` and would not
cancel** — that is what would have forced a second pass. Fixing the reference
fraction is what buys the cancellation, and it costs nothing in proposal quality
because `τ` is drawn afterwards.

### Where the current code goes wrong, precisely

The existing weights are already `w_e` at ½ for the candidates — those
evaluations are reusable as-is. The defect is the **self** weight:

```cpp
double wOrig = std::exp(beta * (llOrig - maxLL));   // llOrig = state->logLik
```

`state->logLik` is `x` evaluated at **its own** branch fractions, while every
candidate is evaluated at ½. Self is scored on a different footing from the
alternatives, so the selection distribution is not a function of `R` alone and
the normaliser does not cancel. The free `state->logLik` is exactly the
shortcut that breaks it — the same shape of error as `weighted_spr`'s
`lMerge`-bins-versus-`lReg`-bins asymmetry (§3).

### Two symmetry conditions, checked against the code

**Prune-edge count: symmetric, no correction needed.** `eligible` is every edge
with `parent != root` (`src/mcmc.cpp:1211–1215`). The root is always
trifurcating in this representation, so exactly 3 edges have `parent == root`
and `|eligible| = nEdge − 3` **independently of topology**. The prune edge
`(u → v)` has `parent = u ≠ root` in both `x` and `y`, so it is eligible in
both, and the subtree `S` ↔ edge correspondence is 1–1. Probability `1/|eligible|`
cancels exactly. (TBR argues the same thing for itself at
`src/tree_moves.cpp:245–247`.)

**Candidate set: NOT symmetric as written — this is a second defect.**
`src/mcmc.cpp:1255–1259` excludes the three edges incident to `u`:

```cpp
if (isDesc[state->child[i]]) continue;                            // S's edges
if (state->parent[i] == u || state->child[i] == u) continue;      // incident to u
```

In `R` those three collapse to the single merged edge, so from `x` the candidate
set is `E(R) \ {e_x}`. Apply the same filter from `y`, where `u` sits on `e_y`,
and it is `E(R) \ {e_y}`. **Each direction excludes its own current position.**
The two sets have equal size but different membership, so

```
Σ_{e ≠ e_x} w_e   ≠   Σ_{e ≠ e_y} w_e
```

and the normaliser does not cancel. The §3a cancellation fails under the filter
as written.

**Fix: include the merged edge as a candidate**, making the set exactly `E(R)`
from both directions. In unpruned indexing the merged edge is the
`(parentRow, sibRow)` pair, so it needs a small special case evaluated at ½ of
`lMerge` — which is the *same* evaluation §3a already requires for the self
weight. So including self is not a tidiness choice: **it is what makes the
normaliser cancel**, and it and the `state->logLik` correction are one change,
not two.

### The gate cannot validate the selection term — a real limitation

`spr_fixed_surrogate` is "R-level MH on `spr_proposal` accepting at
`log U < logHastings`", i.e. uniform τ + Jacobian + MH, with **no** likelihood
weighting. At the gate's β = 0 every weight is `exp(0) = 1`, so:

- the selection ratio `log w_{e_x} − log w_{e_y}` is identically 0, and
- the set asymmetry above is invisible, because both sets have size
  `|E(R)| − 1` and uniform weights.

At β = 0 the fixed kernel therefore *reduces exactly to*
`spr_fixed_surrogate` — reassuring, since it means a correct implementation is
guaranteed to pass and any red gate indicates an implementation bug rather than
a design flaw. But it also means **the gate is necessary and not sufficient**:
it certifies the Jacobian, τ and MH machinery while being structurally blind to
the two things §3a and this section actually add.

Those need their own test, and it can be deterministic rather than statistical:
assert that the enumerated candidate set from `x` and from `y` are the **same
set of edges of `R`**, and that `w_{e_x}` is computed by the identical code path
as every `w_{e_y}`. A fixture that prunes both endpoints and compares the two
enumerations settles it without any sampling.

### Net cost

| | evaluations |
|---|---|
| now | `nCand` at ½ (self is free but wrong) |
| fixed | `nCand` at ½, **+1** for self at ½, **+1** for the chosen candidate at the drawn `τ` |

Two extra likelihood evaluations per move. Cheap, and far cheaper than the
`2·nCand` a naive reverse pass would cost.

### What changes semantically

Drawing the original position is **no longer a guaranteed no-op**: `τ` will
generally differ from `x`'s current fraction, so that outcome is a legitimate
branch-fraction move at unchanged topology. The current `if (rnd < wOrig) return
false` early-out must go.

## 4. Cost and behaviour changes to expect

- **Acceptance falls from ~0.833 to below 1.** The measured 0.833 ≈
  `nCand/(nCand+1)` was never an acceptance rate — it was the probability of not
  drawing "self". A real accept step will reject some of those.
- **RNG streams change**, so bit-reproducibility against pre-fix runs is gone.
  Unavoidable for any correctness fix here; no checkpoint format changes.
- **Mixing may drop slightly** on datasets where the τ=½ point value happened to
  be a good proposal. That is the price of targeting π, and the ESJD
  instrumentation folded into this pass is what will measure it rather than
  guess.

## 5. All three paths

The halving appears in **all three** commit paths and every one must be patched
and re-gated:

| path | commit site |
|---|---|
| main partial-CL | `src/mcmc.cpp:1463–1466` |
| Q-heterogeneity | `:1777–1780` |
| full fallback | `:1921–1924` |

The gate fixture drives the **main partial-CL** path — what a default run
executes — so a green gate is meaningful for the common case, but paths 2 and 3
need their own coverage before this is called done.

## 6. Order of work

1. ~~`weighted_spr` Jacobian~~ — applied and **empirically validated**: every
   statistic improved 2–19×, `int_frac` from −29.3 % to −1.5 % (§3). Kept.
   Residual deviation means a second defect remains; GSPR-003 stays open.
2. Fix `gibbs_spr` on the `spr_proposal_impl` pattern (§3), all three paths.
3. Re-gate: `--quick` must flip FAIL → PASS, and the `spr` /
   `branch_lengths` / `tbr` controls must stay clean (they are untouched, so any
   movement there means the harness or the build is wrong, not the fix).
4. Full `testthat`; expect topology-move test churn from the RNG change.
5. ESJD jump accumulators at the same commit sites, per
   `2026-08-12-esjd-implementation-plan.md` §2.

## 7. Not in this fix

- **GSPR-002** (stale `logPrior`, likelihood-only candidate weights): the
  harness confirmed these are *harmless under the flat Dirichlet prior* — one
  candidate mechanism ruled out rather than a defence. Routing through
  `weighted_spr`, which computes `compute_log_prior_at` and updates
  `state->logPrior`, retires GSPR-002 as a side effect. Verify, don't assume.
- `gibbs_subtree_swap`: separate defect (missing neighbourhood normaliser), not
  a deterministic-split problem. Default-on; needs its own fix and its own gate
  arm.
