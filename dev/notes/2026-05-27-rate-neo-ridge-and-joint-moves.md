# Rate_neo / tree_length ridge & joint moves after audit Issue 1 fix

**Status:**
- **A1 (`joint_tl_rn`) — DONE.** Implemented in `src/mcmc.cpp` (case 33) and fully wired in R (commit `e08d328`, main, 2026-05-27). See task checklist below.
- **A2 (neo_joint fate) — OPEN.** Needs a 5–10k-iter warmup on a mixed dataset to measure `cor(log rate_loss, log rate_neo)`, `cor(log T, log rate_neo)`, `cor(log T, log rate_loss)`. Harvest from the A3 Hamilton smoke run.
- **A3 (SBC on mixed data) — OPEN.** SBC harness with mixed-partition simulations has not yet been run against the patched chain. Planned for Hamilton.
- **B (deterministic trans-invariant move) — DEFERRED.** Profile first.

**Context:** Branch `fix/partition-rate-normalisation` (commits `c739d0c`,
`9e68cad`, `4a38326`), now merged, landed RB-style partition-rate normalisation. Audit
at `dev/rb-equivalence/notes/partition-rate-and-acrv-audit.md` Issue 1.
NEWS.md has the user-facing summary. Quick recap:

    neoScale   = r / (1+r) · (n_neo + n_trans) / n_neo
    transScale = 1 / (1+r) · (n_neo + n_trans) / n_trans

The fix makes `rate_neo` (= r) control the **split** between neo and trans
effective rates, with `tree_length` (= T) carrying the overall scale and
the joint-weighted mean held at 1. Changing r now moves both partitions'
effective branches in opposite directions.

That has two follow-up consequences that this note unpacks.

---

## Item A (urgent): the (T, r) posterior ridge has rotated

### What changed

**Pre-fix**, `rate_neo` was a one-sided multiplier on neo edges:
`neo_eff = T · b · r`, `trans_eff = T · b`. The ridge between T and r in
the posterior was diagonal in the neo direction only — doubling r and
halving T kept neo evolution unchanged, while trans evolution was unaffected
by r at all. So **for mixed datasets the (T, r) ridge was essentially the
neo-only ridge**, and the trans data acted as an independent anchor that
strongly constrained T directly.

**Post-fix**, `rate_neo` enters both partitions. Trans no longer anchors T
directly; T has to balance against the joint `n_neo·neoScale + n_trans·transScale`
constraint. The new ridge geometry, computed from the formula:

- A *deterministic* invariance line exists for *trans* evolution:
  `T·transScale = const` ⟺ `T/(1+r) = const`, i.e. `T ∝ (1+r)`.
- A *deterministic* invariance line for *neo* evolution:
  `T·neoScale = const` ⟺ `T·r/(1+r) = const`, i.e.
  `T ∝ (1+r)/r = 1/r + 1`.

The posterior ridge is the data-weighted compromise between these two lines.
For typical mixed datasets (Casali-style with n_trans ≫ n_neo) the
trans-invariance line dominates: the ridge is approximately `T ∝ (1+r)`,
which is **a different ridge from the pre-fix one**.

### What that breaks for existing mixing aids

There are four joint/ridge-aware moves currently in the codebase:

| Move (moveType) | Pre-fix purpose | Status post-fix |
|---|---|---|
| `neo_joint` (18) — scale rate_loss & rate_neo by same factor | Empirical: rate_loss and rate_neo posteriors observed to be correlated | **Validity unclear.** Under the new parameterisation rate_loss only enters stationary frequencies (not rate scale), while rate_neo controls the rate split; their posterior correlation is no longer guaranteed to be on a 1:1 log-diagonal. Probably still useful but probably mistuned. |
| `joint_tl_rls` (21) — 2D Bactrian on (log T, log rate_log_sd), ρ adapted | T × ACRV variance trade-off | **Unaffected** (rate_log_sd doesn't enter the partition-rate formula) |
| `joint_tl_rl` (22) — 2D Bactrian on (log T, log rate_loss), ρ adapted | T × rate_loss trade-off when neo present | **Unaffected** for the same reason |
| *(no move)* — (T, rate_neo) joint | n/a — pre-fix r didn't constrain T much | **NEW RIDGE.** No move currently rides it. |

The headline problem: the strongest posterior correlation introduced by the
fix is on (T, r), and **nothing currently mixes along it**. Independent
Bactrian/slice moves on T and r will see acceptance drops and ESS drops
proportional to how tight the ridge is.

### What to do (urgent)

**A1. Add `joint_tl_rn` (implemented as moveType 33 — case 25 was already
taken by `gibbs_kprime_sweep`): 2D Bactrian on (log T, log r) with adaptive ρ.** Mechanics mirror `joint_tl_rls` exactly; the
infrastructure already exists in `bactrian_2d_perturbation` and the
`jointRhos` matrix plumbing from M-120 (`dev/plans/2026-03-29-1716-m-120-...`).
Concretely:

- New case in `do_move_impl` at `src/mcmc.cpp:~4350` (next to cases 21, 22):
  ```cpp
  case 33: { // joint_tl_rn: 2D Bactrian on (log T, log r)
    double z1, z2;
    bactrian_2d_perturbation(jointRho, z1, z2);
    double multT = std::exp(scaleTuning * z1);
    double multR = std::exp(scaleTuning * z2);
    state->treeLength = oldTL * multT;
    state->rateNeo    = oldRN * multR;
    logHastings = std::log(multT) + std::log(multR);
    break;
  }
  ```
- R-side wiring: add `joint_tl_rn` to the move-name allowlist in
  `R/MkPrimeMCMC.R::.ValidateMoveWeights` and to the move construction in
  `R/RunMkPrime.R` (look for the existing `joint_tl_rls` definition and
  copy/adapt). Gate on `hasNeo`.
- ρ adaptation: piggyback on the existing warmup correlation estimator
  used by `joint_tl_rls` / `joint_tl_rl`. Track `cor(log T, log r)` in
  the chain matrix; cap at ±0.95.
- PLC: falls through to `default` (full recompute). Both T and r affect
  every partition; partial-cache update is impossible for this 2D move.
- nodeCL cache: invalidate all. Same reason.

**A2. Revalidate `neo_joint` and decide its fate.**

Re-derive what posterior correlation `neo_joint` was riding. Under the new
parameterisation:

- rate_loss enters the neo Q-matrix via stationary frequencies
  (π₀ = rl/(1+rl), π₁ = 1/(1+rl)), but NOT the rate (the eigenvalue of Q
  is fixed at 2 regardless of rl — see `src/node_cl_cache.h:225-241`).
  So scaling rl shifts equilibrium, not speed.
- rate_neo enters the partition rate scales as above.

The pre-fix `neo_joint` (scale both by the same factor) made geometric
sense when rate_neo was a pure neo-rate scalar: doubling rate_neo halved
the effective T-budget for neo, which doubling rate_loss compensated for
in some sense. Post-fix this rationale is murkier.

**Concrete action:** run a 5–10k-iter warmup on a representative mixed
dataset (the smoke `out/mkprime_635_by_nt_9v.rds` is a reasonable target;
Casali is too small on the neo side to give a strong neo signal). Estimate

    cor(log rate_loss, log rate_neo)
    cor(log T,        log rate_neo)
    cor(log T,        log rate_loss)

from the sampled posterior. Three plausible outcomes:

1. `(rate_loss, rate_neo)` correlation is high. `neo_joint` is still
   earning its keep; possibly retune it as a 2D Bactrian instead of a 1D
   diagonal (i.e. replace case 18 with a 2D-ρ-adapted version).
2. `(rate_loss, rate_neo)` correlation is weak but `(T, rate_neo)` is
   strong. Drop `neo_joint`, rely on the new `joint_tl_rn` from A1.
3. All three correlations are weak. Drop `neo_joint`, no replacement
   needed beyond A1.

Until A2 is done, leave `neo_joint` in place — it can't hurt, it's just
possibly inefficient.

**A3. SBC the patched chain on a mixed dataset.**

The audit predicts the cross-sampler (T, r) rhat between MkPrime and
RevBayes should approach the within-sampler rhat once both fixes are in.
Run the SBC harness with the patched MkPrime against either (a) the patched
RB once the matching ACRV fix lands, or (b) a known-correct toy posterior
derivable analytically on a 3-tip 2-character mixed dataset.

If SBC fails *post-A1*, that's a hint the ridge geometry I sketched above
is wrong — please dig.

---

## Item B (microopt, defer until measured): the deterministic trans-invariant joint move

### The idea

For a single-parameter move that touches `rate_neo`, the partial cache for
**every partition** becomes stale (the F1 finding from the external-reviewer
review on 2026-05-27). So the cheapest available move type for rate_neo
currently does an `O(n_partitions × n_branches × k²)` full recompute.

But trans evolution is invariant under the joint transformation

    (r, T)  →  (r', T'),    T' = T · (1+r') / (1+r)

because `T'·transScale(r') = T·transScale(r)` identically. Under this
transformation:

- Trans `partLogLik` entries are exactly unchanged → don't recompute.
- nodeCL `unit.rateScale` and `topo.edgeLen` change individually but
  their product (which is what `apply_transition` consumes) is invariant →
  stored CLs are still correct in value.
- Neo `partLogLik` and nodeCL entries change non-trivially (by factor
  r'/r in effective branch) → must invalidate.

So a deterministic "trans-invariant rate_neo" move could refresh ~half the
work of a full recompute. Microopt savings: roughly the trans share of
total CL work, weighted by the rate_neo move's share of total iterations
(~6% under default weights when hasNeo).

### Maths

Hastings: Bactrian draw on log r gives r' = r · exp(σε). Push T → T'
deterministically as above. The Jacobian of the bijection (r, T) → (r', T')
is

    |J| = ∂(r', T')/∂(r, T) = | (r'/r)·0       0           |
                              | ∂T'/∂r         (1+r')/(1+r) |
        = (r'/r) · (1+r')/(1+r)
                                        ── but the (r,T) coords each carry their own scale →
                                           log|J| = log(r'/r) + log((1+r')/(1+r))

So log-Hastings:

    log H = log(r'/r)              // Bactrian on r
          + log((1+r')/(1+r))       // Jacobian from the deterministic T push

Prior delta: both r (LogNormal) and T (Gamma) priors enter:

    Δ log prior = log π_r(r') - log π_r(r)
                + log π_T(T') - log π_T(T)

Reverse move: same Bactrian draw with sign flipped; the bijection is its
own inverse modulo sign. Detailed balance holds because the joint move is
a deterministic bijection on (r, T) composed with the Bactrian kernel on r,
and the Jacobian is correctly accounted.

### Cache-update plumbing (the unpleasant part)

The current nodeCL cache invalidation API has two granularities:
`invalidate_neo_cls` and `invalidate_all_cls`. To exploit the deterministic
trans invariance, we need a third mode: "trans units' rateScale and
topo.edgeLen change in lockstep, the product is invariant, do NOT mark
clValid = false."

Concretely, after the joint move accepts:

1. For trans units (isMkN == false): update `unit.rateScale = new_transScale`,
   update `topo.edgeLen[e] = T' · relBr[e]`. Leave `unit.clValid = true`.
2. For neo units (isMkN == true): update `unit.rateScale = new_neoScale`,
   update `topo.edgeLen[e] = T' · relBr[e]`. Set `unit.clValid = false`
   (regular invalidation).
3. For partLogLik: refresh only neo partitions; trans entries stay.

The risk: the "do NOT mark stale even though edge lengths changed" invariant
is fragile. If a future code change updates trans `rateScale` *without*
updating `topo.edgeLen` in lockstep (or vice versa), trans CLs become
silently wrong. That's a worse failure mode than the F1 partial-cache bug
because it can't be caught by a "cache vs direct" parity test on a single
state — it only diverges across moves.

### When to do it (decision criterion)

Profile first. Specifically, on a representative chain post-A1:

1. What % of wall-clock is spent in rate_neo / neo_joint / slice_rate_neo
   moves?
2. Within those, what's the breakdown between neo vs trans CL work in
   the full recompute?

If `rate_neo moves > 5%` of wall-clock AND `trans share > 50%` of CL work,
the deterministic move saves real time. Otherwise leave it.

The standalone `joint_tl_rn` from A1 doesn't get this optimisation — both
T and r change *and* their changes don't cancel for trans, so full recompute
is genuinely required for the 2D Bactrian. The deterministic move is a
*different* move that would coexist with `joint_tl_rn`, used adaptively
(e.g., when the chain is in a regime where the deterministic constraint
is close to the local posterior gradient direction).

### Risks / open questions

1. **Mixing.** The deterministic move follows a 1D line in (r, T) space.
   If the actual posterior ridge isn't aligned with that line, the move
   doesn't help mixing — it just explores along the wrong direction at
   higher acceptance. The 2D Bactrian `joint_tl_rn` is strictly better for
   mixing because ρ adapts to the actual posterior shape. The deterministic
   move is *only* worth implementing if (a) the cache-saving optimisation
   is measurably worth it, and (b) the deterministic line happens to be
   close to the posterior ridge under typical priors (i.e. trans-dominated
   datasets, where the trans-invariance line nearly is the posterior ridge).
2. **Cache fragility** as noted above.
3. **Detailed balance verification.** The math is straightforward but the
   project's tradition (see `test-tbr-detailed-balance.R`) is to
   numerically verify each new move on a flat-posterior toy. Budget for
   that.

---

## Suggested execution order

1. **Day 1.** Implement A1 (`joint_tl_rn`). Same shape as `joint_tl_rls`.
   Test that it runs, mixes, and doesn't break existing tests. **Ship.**
2. **Day 1.** Run the smoke warmup, dump the posterior correlations
   needed for A2. **Decide neo_joint's fate.** Either patch it (case 18 →
   2D Bactrian) or remove it. **Ship.**
3. **Week 1.** SBC the patched chain on mixed data (A3). If it passes,
   the partition-rate fix is fully validated. If it fails, the bug is
   either in the patch or in the ridge framing above — both worth catching.
4. **Deferred.** Item B only if profiling shows rate_neo moves are a
   measured bottleneck. The cache-fragility risk makes this a poor
   speculative investment.

## Files & symbols to know

| Concern | Location |
|---|---|
| Partition-rate helper | `src/mcmc_state.h::compute_partition_scales` |
| The new RB-style scaling applied | `src/mcmc_likelihood.cpp::cpp_partition_log_likelihood` (scaledEdge at top) |
| Existing 2D-joint move code (template for A1) | `src/mcmc.cpp` case 21 / case 22, around line ~4370; `bactrian_2d_perturbation` in same file |
| Existing 1D `neo_joint` (to revalidate) | `src/mcmc.cpp` case 18, line ~4364 |
| Joint-move R wiring (template for A1) | `R/RunMkPrime.R`, search for `joint_tl_rls` |
| ρ adaptation | `R/RunMkPrime.R`, search for `jointRhos` |
| nodeCL cache invalidation surface | `src/node_cl_cache.h::NodeCLCache::invalidate_*` |
| Partition-LL partial-cache switch | `src/mcmc.cpp::do_move_impl`, search for `newPC = state->partLogLik` |
| F1 regression test (cached vs direct after rate_neo moves) | `tests/testthat/test-partition-rate.R::partLogLik stays in sync...` |
| The audit document | `dev/rb-equivalence/notes/partition-rate-and-acrv-audit.md` |
| The original M-120 design note (joint Bactrian infrastructure) | `dev/plans/2026-03-29-1716-m-120-*.md` |

## Task status

### A1 — DONE (2026-05-27)

`joint_tl_rn` implemented and wired end-to-end:

| Concern | What was done |
|---|---|
| C++ move | `src/mcmc.cpp` case 33: 2D Bactrian on (log T, log r), `logHastings = log(mult1) + log(mult2)` |
| Cache invalidation | Falls to `default: invalidate_all` — correct (both T and r changed) |
| PLC (partLogLik) | Falls to `default: full recompute` — correct (rate_neo invalidates all partitions per Issue 1 fix) |
| Move type map | `joint_tl_rn = 33L` in `.moveTypes` |
| R BuildMoves | Added gated on `hasNeo && joint2d` with `weight = 1, dim = 2L` |
| rho estimation | `rho_tl_rn` tracked in `.EstimateJointRhos`; uses `cor(log T, log rate_neo)` |
| rho initialisation | `chainRhos` initialised with `rho_tl_rn = 0.0` |
| rho pass-through | `.BuildChainRhoMatrix` passes `rho_tl_rn` for `joint_tl_rn` moves |
| Allowlist | Added to `.ValidateMoveWeights` valid-name vector |
| Default scale | `scale_joint_tl_rn = 0.5` in `MkPrimeMCMC` defaults |
| Default MH tuning | `joint_tl_rn = 0.25` in `.defaultMhTargets` |
| Tuning-name map | `joint_tl_rn = "scale_joint_tl_rn"` in `.scaleTuningParam` |
| Move category | Added to `.moveCategories` under "Rates" |
| ResumeMkPrime | `scale_joint_tl_rn` tuning restored on resume |

Smoke-test ESS comparison vs independent Bactrian: **pending A3 Hamilton run**.

### A2 — OPEN: neo_joint fate

Run a 5–10k-iter warmup on a representative mixed dataset (pid 635 or pid
3832 from the rb-equivalence matrix set). Extract and report:

    cor(log rate_loss, log rate_neo)
    cor(log T,         log rate_neo)
    cor(log T,         log rate_loss)

Use the `rhoSampleBuf` dump from `RunMkPrime` or a short standalone warmup.
Three outcomes govern neo_joint's fate (see Item A above). Record verdict as
`dev/notes/YYYY-MM-DD-neo-joint-fate.md`. **Do not change case 18 without
committing to one of the three outcomes.**

Harvest these correlations from the A3 Hamilton smoke warmup — they come
for free from the same run.

### A3 — OPEN: SBC on mixed-partition data

Run the existing SBC harness (`dev/red-team/heavy-tests/sbc.R`) with a
mixed-partition dataset (neo + trans chars). Confirm:
- All scalar parameters pass rank-histogram uniformity (p > 0.05 after
  Bonferroni or FDR correction at the campaign's standard 5% level).
- `rate_neo` and `tree_length` are not correlated failures.

If A3 passes, the partition-rate fix is fully validated end-to-end for
mixed datasets. Write verdict to `dev/red-team/sbc-results-mixed/verdict.txt`.

**Note:** the rb-equivalence smoke run (Hamilton, pid 950/635 × by_nt_9v)
is a *cross-sampler* equivalence check — not SBC. Both are needed;
they answer different questions.

### B — DEFERRED

Profile first. Only implement the deterministic trans-invariant move if
`rate_neo` moves exceed 5% of wall-clock and the trans share exceeds 50%
of CL work in that budget. See Item B above for the full decision criterion.
