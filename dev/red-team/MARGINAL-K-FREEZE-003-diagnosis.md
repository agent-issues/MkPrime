# MARGINAL-K-FREEZE-003 — free-topology freeze, REPRODUCED + LOCALIZED + FIXED (2026-06-02)

> Repro driver: `dev/red-team/numerical/marginal-k-freeze-repro.R`
> Results: `dev/red-team/numerical/marginal-k-freeze-results/` (freeze-repro.log,
> sweep-r02.rds, dataset-dependence.rds, partA-runmkprime.rds).
> Fix verified 2026-06-02 — see RESOLVED below. (§forks retained for history.)

## RESOLVED (2026-06-02) — fix implemented + verified

Both bugs fixed on `feat/marginal-k`:

- **Bug B fix** (`src/mcmc.cpp`, `src/mcmc_state.h`): `compute_per_kprime_log_lik`
  now takes `parent`/`child` parameters, and `cpp_log_likelihood_marginal`
  threads its own `parent`/`child` args through — so a *proposed* topology is
  evaluated correctly rather than `state`'s pre-commit (old) tree. The 7 internal
  `state->parent/child` reads in the helper now use the passed topology; the
  Gibbs-kPrime-sweep caller (case 25) passes `state->parent/child` (unchanged
  behaviour). **Verified:** the move-type gap-sweep now shows **nni/spr/pspr
  baseline_gap = 0.000** (was 3.305 / 1.302 / 4.843), and the tracked
  `test-marginal-k-cache-option-a.R` NNI warm≠cold test now **PASSES**.
- **Bug A mitigation** (`R/RunMkPrime.R` `.BuildMoves`): `gibbs_spr`,
  `gibbs_subtree_swap`, `weighted_spr`, `weighted_subtree_swap` are disabled
  under `marginal_k` (they write a fixed-kPrime `state->logLik` and never call
  the marginal evaluator, so the Bug-B threading fix does not reach them — the
  gap-sweep confirmed residuals: gibbs ~10.6, weighted 0.17/1.5 nat). Topology is
  still searched by the now-correct `nni`/`spr`/`pspr`/`tbr`. A marginal-aware
  Gibbs candidate eval, and the weighted-move residual root-cause, are deferred.
- **Guards:** opt-in C++ coherence assertion (`-DMKPRIME_CHECK_MARGINAL_COHERENCE`)
  + always-on `tests/testthat/test-marginal-k-free-topology.R`.

**End-to-end** (RunMkPrime, free topology, the formerly-frozen `n16_c16_r02`):
marginal_k **no longer freezes** — sd(tl/rls/p) = 47.4 / 3.79 / 0.059 (was
0/0/0); all moves accept (nni 0.91, spr 0.45, mh_logit_p 0.84, slice 1.0); Gibbs
moves absent from the schedule. marginal/loglik/cache testthat: 153 pass / 0 fail.

**Follow-ups:** (1) RB-equivalence VALIDATED deterministically and (2) the
weighted-move residual ROOT-CAUSED, with the audit extended to every move — see the
FOLLOW-UP section immediately below. (3) re-enabling the disabled moves efficiently
(force-cold-per-eval / cache fingerprint) remains deferred ("efficiency tomorrow").

---

## FOLLOW-UP (2026-06-02) — RB-equivalence validation + exhaustive cache-coherence audit

Two correctness follow-ups closed ("validate; explore the gap").

### (A) Free-topology RB-equivalence — VALIDATED deterministically

The freeze fix changed *which* topology the marginal evaluator reads, so the decisive
correctness check is the Rao-Blackwell likelihood identity generalised across topology
space — not a noisy MCMC overlap:

    lse_{k' in [kObs, K]} [ logPrior_S(k') + logLik_S(k') ] == logPrior_M + logLik_M

`tests/testthat/test-marginal-k-rb-free-topology.R` checks this to 1e-7 on 8 random
topologies x p in {0.1, 0.4} x rate_log_sd in {0, 0.6} (nCat = 4, exercising BOTH the
homogeneous and the ACRV pruning paths the fix rerouted). **64 assertions PASS.** This
proves the marginal evaluator is RB-correct on *arbitrary* trees: target-equivalence
with sampled_k is established deterministically, with no MCMC noise. The complementary
"an accepted topology MH move commits the correct marginal LL" facet is guarded by
`test-marginal-k-free-topology.R` (committed `state->logLik` == cold recompute after
nni/spr/tbr/pspr). A signal-bearing MCMC posterior-overlap run is *confirmatory* (the
existing T-OVL uses signal-free data + an over-powered KS bar); given the deterministic
identity it is not load-bearing, and remains a local-smoke + Hamilton-prep item.

### (B) The weighted-move residual — a marginal-cache-coherence bug

`compute_full_loglik_at` does NOT force the marginal `charLLCache` cold: under
marginal_k it calls `cpp_log_likelihood_marginal`, which uses the fast-path when
`state.charLLCacheReady` is true and sets it true after filling. The `useCache` gate
(mcmc.cpp:4306) has **no topology/branch fingerprint** — the per-(char,k') raw LLs it
stores depend on topology, branches, and rate, so the cache is valid ONLY across a
pure-`p` change. A move that evaluates multiple topology/branch configs per call
without forcing the cache cold between them reads the FIRST config's cache for every
subsequent eval (state==warm!=cold).

Every `compute_full_loglik_at` multi-eval site (mcmc.cpp):

| move (type)                | sites              | gated under marginal_k? |
|----------------------------|--------------------|-------------------------|
| gibbs_spr (10)             | 1816               | yes (Bug A)             |
| gibbs_subtree_swap (11)    | 2447               | yes (Bug A)             |
| weighted_branch_scale (12) | 2568               | **NOW (this audit)**    |
| block_gibbs_branch (15)    | 2705, 2761         | **NOW (this audit)**    |
| weighted_spr (13)          | 2891, 2934, 3028   | yes                     |
| weighted_subtree_swap (14) | 3149, 3228         | yes                     |

The exhaustive coherence sweep (`marginal-k-freeze-repro.R`, SWEEP_ORDER extended to
every move type) on the frozen cell confirms the committed-gap signatures:

    nni/spr/pspr/tbr/mh_p ......  0.000    single eval / cache-by-design  -> SAFE
    weighted_branch_scale (12) . -1.261    state==warm!=cold  (cache)  -> was UNGATED
    block_gibbs_branch (15) .... -1.604    state==warm!=cold  (cache)  -> was UNGATED
    gibbs_spr/subtree (10/11) .. +10.5     state!=warm==cold  (fixed-k')
    weighted_spr/subtree (13/14) -0.9/-2.3 state==warm!=cold  (cache)
    gibbs_kprime_sweep (25) .... +6.848    state!=warm==cold  (fixed-k')

**Action:** `weighted_branch_scale` (12) and `block_gibbs_branch` (15) were UNGATED
under marginal_k (default-off, so latent — but a user enabling either would silently
corrupt the chain). Both now gated in `.BuildMoves`. The complete set of six
incoherent moves is now disabled under marginal_k.

**k'-sampling moves (kPrime 7, gibbs_kprime_sweep 25, block_kprime_shift 26):** the
sweep fires case 25 via the low-level `do_move_cpp` and it shows the +6.85 fixed-k'
gap — but `.BuildMoves` ALREADY excludes every k'-sampling move under marginal_k (k'
is integrated out; only `mh_logit_p` is added; RunMkPrime.R:3603). The sweep result
confirms that exclusion is *load-bearing* (those moves would corrupt if scheduled),
not redundant; `test-marginal-k-free-topology.R` Test 2 now also asserts it.

**Enabled marginal_k schedule is coherent.** After gating 12 + 15, every move that can
appear under marginal_k either evaluates the proposal exactly once (the MH topology /
branch / rate / joint moves, after the `do_move_impl` entry invalidation at
mcmc.cpp:4706), forces the cache cold per eval (the slice sampler,
`eval_slice_target`:3347), or uses the cache by design (`mh_logit_p`, valid across a
pure-`p` change). The coherence sweep reads 0.000 for every such move.

**Deferred (efficiency, "tomorrow"):** re-enable the six moves under marginal_k by
forcing the `charLLCache` cold per intermediate eval (NOT a blanket force-cold in
`compute_full_loglik_at` — `mh_logit_p` needs the warm fast-path), and/or add a
topology/branch fingerprint to the `useCache` gate.

---

(Original diagnosis follows — retained for the record.)

## TL;DR

Under `likelihoodMode="marginal_k"` with **free topology**, every
topology-changing move writes an **incorrect `state->logLik`**. That wrong value
becomes the slice `logY0` / MH baseline for the continuous parameters
(`tree_length`, `rate_log_sd`, `p`), so their proposals — evaluated at the true
(lower) marginal LL — can never beat the inflated baseline ⇒ 0% acceptance ⇒
freeze. The Gibbs topology moves (which always commit) re-inflate the baseline
every iteration; the chain locks. `mh_logit_p` (case 30) is the **only** move
that leaves `state->logLik` coherent.

This is the **v1 marginal-k gating bug**. The Stage-2 `sampled_k` SBC PASSED but
ran `fixTopology=TRUE`, which never exercises topology moves — so it
**structurally cannot** catch this. Free topology is the production regime
(marginal_k is the intended ≥100-tip default).

## Reproduction (RunMkPrime, free topology, 16-tip × 16-char, autoTune=FALSE)

| cell | mode | sd(tree_length) | sd(rate_log_sd) | sd(p) | verdict |
|------|------|----------------:|----------------:|------:|---------|
| n16_c16_r02 | **marginal_k** | **0.0000** | **0.0000** | **0.0000** | **FROZEN** |
| n16_c16_r02 | sampled_k | 23.06 | 1.675 | 0.046 | mixes |
| n16_c16_r03 | marginal_k | 2.495 | 1.207 | 0.068 | mixes |
| n16_c16_r03 | sampled_k | 0.898 | 1.507 | 0.027 | mixes |

**Frozen-chain (r02) per-move acceptance — the fingerprint:**

```
       tree_length 0.0000   branch_lengths 0.0000   nni 0.0000   spr 0.0000
         gibbs_spr 0.6678   gibbs_subtree_swap 0.5595   tbr 0.0000   pspr 0.0000
  dirichlet_branch 0.0000   local_dirichlet 0.0000   mh_logit_p 0.0000
       rate_log_sd 0.0000   slice_rate_log_sd 0.0000   joint_tl_rls 0.0000
```

Only `gibbs_spr` and `gibbs_subtree_swap` accept (Gibbs moves always commit a
sample); **everything else is pinned at 0%**. Healthy chains (r03, or any
sampled_k) accept 0.4–0.9 across the board.

## Localization — discriminating move-type coherence-gap sweep

For each move type: build a fresh marginal_k state (INIT-001 makes the init
`state->logLik` correct), fire that move until it accepts, then measure

```
baseline_gap = state->logLik (what the move LEFT) − cold marginal LL
warm_gap     = warm eval (cache fast-path) − cold marginal LL
```

`probe()` is validated by two controls reading ≈0: the post-init state, and
`mh_logit_p` (below).

| moveType | name | baseline_gap (nat) | warm_gap (nat) | reading |
|---------:|------|-------------------:|---------------:|---------|
| 30 | mh_p (MH) | 4e-11 | 4e-11 | **clean** (control) |
| 5  | nni (MH) | **3.305** | 3.305 | state==warm ≠ cold (wrong) |
| 6  | spr (MH) | **1.302** | 1.302 | state==warm ≠ cold (wrong) |
| 20 | pspr (MH) | **4.843** | 4.843 | state==warm ≠ cold (wrong) |
| 13 | weighted_spr | **0.168** | 0.168 | state==warm ≠ cold (wrong) |
| 14 | weighted_subtree_swap | **2.862** | 2.862 | state==warm ≠ cold (wrong) |
| 10 | **gibbs_spr** | **+10.65** | **0.000** | fixed-kPrime baseline; cache invalid→warm=cold |
| 11 | **gibbs_subtree_swap** | **+10.59** | **0.000** | fixed-kPrime baseline; cache invalid→warm=cold |

(`warm_gap=0` for the Gibbs rows means the cache is left **invalid**, so the warm
eval recomputes cold — i.e. the *cache* is fine; the bug is purely the wrong
`state->logLik` the move left behind, not a cache that holds a wrong value.)

The `warm_gap` column cleanly separates **two distinct bugs**:

### Bug A — Gibbs topology moves write a fixed-kPrime `state->logLik`
`gibbs_spr_impl` / `gibbs_subtree_swap_impl` (cases 10/11, **default ON** via
`gibbsSpr`/`gibbsSubtreeSwap=TRUE`) group transformational characters at the
**single** current `state->kPrime` (`src/mcmc.cpp:1254-1282`) through the
sampled-k partial-CL machinery, producing a **fixed-kPrime** likelihood — no
marginal-over-k′ logsumexp, no geometric `P(u|p)` weight. They then write
`state->logLik = candLL[chosen]` (`mcmc.cpp:1410`, `:2132`) and **return early**
from the dispatch (`:4976`, `:4979`), before any marginal re-eval. `do_move_impl`
invalidates the cache at its top (`:4694`) — so the cache is left **invalid**
and the next warm eval simply recomputes cold (`warm_gap=0`; the cache never
holds a wrong value) — but nothing restores the correct **marginal** baseline.
Result: `state->logLik` inflated (the same INIT-001 inflation, now recurring on
every accepted Gibbs move). This is the dominant in-situ driver: the frozen r02
chain accepts ONLY gibbs moves, so it re-inflates every iteration.

**The inflation magnitude is the missing geometric weight — confirmed
quantitatively.** The omitted term is `−Σ_char log P(u_char | p)`; with all
`u=0` (k′=kObs) that is `−n_char·log(p)`. Sweeping p (PART B.4) the gibbs gap
tracks this prediction with correlation **1.0000**:

| p | gibbs_gap (nat) | −n_char·log(p) | mh_p accept after one gibbs_spr |
|---:|---------------:|---------------:|--------------------------------:|
| 0.30 | 18.71 | 19.26 | 0.000 |
| 0.50 | 10.67 | 11.09 | 0.000 |
| 0.70 | 5.45 | 5.71 | 0.073 |
| 0.90 | 1.60 | 1.69 | **0.950** |
| 0.97 | 0.46 | 0.49 | **0.950** |

(n_char=16; measured gap sits just below the all-`u=0` bound because realized
u>0 characters contribute slightly less than `−log p` each.) This is independent
confirmation that Bug A is exactly "fixed-k likelihood with the geometric `P(u|p)`
normalization dropped," and it explains why the B.2 gap is ~uniform across cells
(shared n_char and shared init p=0.5).

### Bug B — MH/weighted topology moves leave the accepted-move likelihood wrong
For `nni`(5), `spr`(6), `pspr`(20), `weighted_spr`(13), `weighted_subtree_swap`(14),
`baseline_gap == warm_gap` ≠ 0 (1.3–4.8 nat here): after an accepted move,
`state->logLik` equals the **warm** eval, and **both disagree with a cold rebuild**
of the committed tree. So the likelihood written on accept is simply **wrong** for
the resulting tree. **Measured facts (solid, and new):** warm ≠ cold after an
accepted MH/weighted topology move, and the failure spans SPR / pSPR / weighted —
not just the already-tracked NNI case
(`tests/testthat/test-marginal-k-cache-option-a.R` Test 3, warm−cold 0.164 nat).
**Root NOT pinned:** `do_move_impl` invalidates *both* cache tiers at its top for
non-p moves (`:4694`), so the proposal eval *should* be fully cold — which is in
tension with a pure "stale-cache" story. Two cold recomputes of the same committed
tree disagreeing by ~3 nat points instead to an **eval-time-vs-commit-time tree /
edge-length bookkeeping** mismatch (the proposal is evaluated against a
(parent,child,edgeLen) that differs from what gets committed). Could also be a
cache-invalidation gap. Localizing this is fix-stage work. (In a chain already
corrupted by a Gibbs move these MH moves also reject — r02 fingerprint — so Bug B
compounds the freeze rather than rescuing it.)

### Causal close
On the frozen dataset: `mh_logit_p` accepts **95.3%** from a clean baseline;
after **one** `gibbs_spr` (baseline_gap = 10.57 nat) it accepts **0.0%**. One
fixed-kPrime topology move freezes the continuous sampler. (`probe()` invalidates
the cache and does not touch `state->logLik`, so the only changed quantity is the
baseline.)

## Why ~7/16 freeze and not 16/16 (p-gated bistability)

The freeze is **gated by p**, via the Bug-A gap `≈ −n_char·log(p)` (PART B.4
above). A "self-healing" move — `mh_logit_p` on accept, or a slice on completion
— rewrites `state->logLik` to the correct marginal (e.g. `slice_scalar_impl`,
`src/mcmc.cpp:3441`); its accept probability is `~exp(−gap)`. So self-healing
fires readily when the gap is small (high p: at p=0.90 mh_p heals and accepts
0.95) and essentially never when the gap is large (low/mid p: p≤0.5 → 0%). The
chain is therefore **bistable**: if p sits in the large-gap (low/mid) region when
a Gibbs move corrupts the baseline, no self-healing move can fire and the chain
locks; if p has drifted into the small-gap (high-p) region, a self-healing move
breaks the inflation and the chain mixes. r02 froze with p stuck at init
(sd(p)=0); r03 mixed with p moving (sd(p)=0.068). Which basin a chain lands in is
seed/data-dependent (the gap is ~uniform ~10.5 at the shared init p=0.5, so it is
**not** a gap-magnitude threshold across cells) → a ~7/16 incidence, not all-or-
none. ⚠ The exact escape dynamics (how a chain reaches high p before locking,
given that p-moves are themselves gap-suppressed) are **not** fully pinned —
plausibly an escape during the brief uncorrupted window before the first accepted
Gibbs move. Stated as: tested **p-gated bistability**; precise tipping dynamics a
follow-up.

## Why the SBC missed it
`fixTopology=TRUE` ⇒ no topology moves ⇒ `state->logLik` stays correct after the
INIT-001 fix ⇒ calibration SBC passes. The bug lives entirely on the
topology-move path, which only the free-topology (production) regime exercises.
**Lesson:** pair a calibration SBC (which often fixes topology) with a
free-topology mixing check.

## state->kPrime under marginal_k
`get_mcmc_state()$kPrime` returns all 2's on this binary dataset (kObs=2). k′ is
analytically integrated out in marginal_k, so `state->kPrime` is a vestigial
field. This matters for the fix: the Gibbs moves grouping by `state->kPrime`
(`:1258`) are not just leaving a wrong baseline — they also **sampled the new
topology using fixed-k weights**, i.e. against the wrong target. So a
"recompute the baseline after the Gibbs move" patch fixes the freeze but NOT the
topology-proposal correctness.

## Fix forks (NOT chosen here — needs user/advisor call)
1. **Gate Gibbs/weighted topology moves OFF under marginal_k** (route topology
   through the marginal-correct path). Simplest; loses Gibbs efficiency. But the
   MH topology path (Bug B) must be fixed too, or it is no rescue — see below.
2. **Make the Gibbs candidate eval marginal-aware** (marginalize k′ per
   candidate). Correct target + correct baseline; larger change, perf cost.
3. **Recompute `state->logLik = compute_full_loglik` (marginal) after every
   topology move.** Fixes the freeze (Bug A baseline) cheaply, but leaves the
   Gibbs topology *proposal* using fixed-k weights (wrong target) and does not
   touch Bug B (the wrong accepted-move LL on the MH/weighted topology path).
4. **Bug B must be fixed regardless** (the cache-option-a warm≠cold root on the
   MH/weighted topology path): a cold recompute after the move, or fixing the
   proposal-eval/commit edge-length bookkeeping so warm == cold.

A defensive **assertion** (`|state->logLik − compute_full_loglik| < tol` after
each accepted move under marginal_k, in a debug build) would have caught both
bugs immediately and is cheap insurance for any fix.
