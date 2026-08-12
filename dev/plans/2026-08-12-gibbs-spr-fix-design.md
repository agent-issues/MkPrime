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
