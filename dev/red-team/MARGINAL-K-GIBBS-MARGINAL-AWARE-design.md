# Design note — marginal-aware Gibbs candidate eval (re-enable gibbs_spr / gibbs_subtree_swap under marginal_k)

Status: DESIGN / OPEN. The deferred half of FREEZE-003 Phase 2. **Consult advisor
BEFORE coding** (per the overnight plan). Written 2026-06-02 eve.

## Why these two moves are NOT fixed by the Phase-2 scratch-eval

`gibbs_spr`(10) / `gibbs_subtree_swap`(11) do NOT call the marginal evaluator. They
group the transformational characters **by their current `state->kPrime`**
(`mcmc.cpp:1261-1289`), allocate per-`(kPrime, partition)` CLGroups, run ONE caching
downpass (`caching_downpass`), compute residual CLs (detach the pruned subtree), then
score each candidate reattachment by an O(depth) partial-CL update — a **fixed-kPrime
likelihood**: `LL(y | tree, k = state->kPrime_i)`, with NO marginal Σ over
k' ∈ [kObs_i, K] and NO geometric `P(u_i | p)` weight. They then write
`state->logLik = candLL[chosen]` and return.

Under marginal_k, `state->kPrime` is not a sampled quantity (k' is integrated out;
the chain pins it at a placeholder). So the candidate weights are the WRONG TARGET:
the move samples topology ∝ fixed-k LL instead of ∝ marginal LL. A cache-coherence
fix (the scratch-eval flag) cannot repair a wrong proposal *target* — hence these two
stay disabled while 12/13/14/15 were re-enabled.

The coherence gap-sweep confirms the mechanism: gibbs_spr/gibbs_subtree_swap leave
`baseline_gap ≈ +10.5` (= the omitted `-n_char·log p` geometric weight), and
`warm_gap = 0` (the warm recompute is correct; only the committed fixed-k value is
wrong) — the `state!=warm==cold` signature of Bug A, distinct from the
`state==warm!=cold` cache signature the scratch-eval fixed.

## What "marginal-aware" requires

Each candidate topology T_c must be scored by its MARGINAL log-likelihood
  L_marg(T_c) = Σ_i logSumExp_{k' ∈ [kObs_i, K]} [ LL_i(T_c, k') + log P(k'-kObs_i | p) ]
(plus the Model-A/B truncation normaliser, exactly as `cpp_log_likelihood_marginal`
does). The Gibbs/weighted candidate weight is then `exp(beta · L_marg(T_c))`.

## Options

### Option A — route candidates through `compute_full_loglik_at` (scratch)
For each candidate, build (parent, child, edgeLen) and call
`compute_full_loglik_at(..., fillCharLLCache=false)` — the marginal evaluator.
- Correct and tiny: reuses the exact machinery the weighted moves now use.
- Cost: O(nCand) FULL marginal evals (each a full pruning over all k'-slots),
  losing the partial-CL O(depth) reuse that is gibbs_spr's whole reason to exist.
- **Observation: this makes gibbs_spr essentially weighted_spr WITHOUT the
  branch-fraction integration.** weighted_spr (already re-enabled, marginal-correct)
  is strictly richer. So Option A delivers little NEW capability — it is a redundant,
  slower cousin of a move we already ship. Likely NOT worth landing.

### Option B — wire FU-3b: per-k' partial-CL slots (the real marginal-aware Gibbs)
Run the caching downpass / residual / candidate-update **per k'-slot**
(k' ∈ [kObs, K]) using the dormant Tier-2 cache `perKpClSlots[ko]` (the
`McmcState::perKpCl*` fields landed structurally in FU-3, never wired). For each
candidate, combine the per-k' candidate LLs by the geometric-weighted logSumExp.
- Preserves the O(depth) partial-CL reuse PER k'-slot → keeps the efficiency niche.
- Cost: ~(K-kObs+1)×, i.e. up to `kMaxKprimeCand` (~10-30) × the fixed-k path, but
  still far cheaper than nCand full re-prunings on large trees.
- Large change: per-k' residual bookkeeping, the Tier-2 cache contract (allocation,
  invalidation rhythm — note the scratch-eval invariant must extend to Tier-2 once
  it is no longer dormant!), ascertainment pseudo-groups per k'.
- This is the only version that justifies a distinct marginal gibbs move.

## Recommendation (for the advisor / Martin)

Marginal topology search is ALREADY correct and available under marginal_k via
nni / spr / pspr / tbr (MH) + the four re-enabled weighted/block moves. The open
question is purely **mixing efficiency on large trees**, where gibbs_spr's
partial-CL reuse would help — but only Option B realises that, at the cost of wiring
Tier-2 (and extending the scratch/coherence contract to it). Option A is redundant
with weighted_spr.

So the decision is binary:
1. **Leave 10/11 deferred indefinitely** — rely on weighted_*/MH for marginal topology;
   document that gibbs_spr/gibbs_subtree_swap are sampled_k-only by design. (Cheapest;
   no correctness gap since the moves are gated, not silently wrong.)
2. **Invest in Option B (FU-3b)** if large-tree marginal_k mixing proves inadequate
   in practice — gated behind an empirical mixing need, not done speculatively.

Suggest NOT coding Option A. Pick (1) now and revisit (2) only if/when a large-tree
marginal_k run shows a topology-mixing deficit that the weighted moves don't cover.

## If Option B is chosen — verification gates (same as Phase 2)
- gap-sweep: cases 10/11 `baseline_gap → 0` (currently +10.5).
- A Test-3 analogue: fire 10/11 many times under marginal_k; committed==cold AND
  warm==cold on accept+reject; and the Tier-2 ready-flag coherence (now live).
- A deterministic per-candidate check: the move's candidate weights == the marginal
  LL of each candidate from `compute_full_loglik_at` (scratch), to ~1e-8.
- moves-on overlap with gibbsSpr/gibbsSubtreeSwap ON, marginal vs sampled.
- full testthat FAIL=0; sampled_k path bit-identical (the per-k' grouping must reduce
  to the current fixed-k path when only one k'-slot is occupied).
