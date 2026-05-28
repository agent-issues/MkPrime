# 2026-05-28 — Implementation plan: marginalise k'_i out of the Mk' likelihood

Status: **PLAN — no code yet.** Scoped for a future implementer (probably
me on a later run) to pick up cold. Companion to
`dev/notes/2026-05-28-kprime-viability.md` (the diagnostic that justifies
this work).

## TL;DR

Per-character k'_i is *not identified* by the data at empirical tree
lengths — finding **EG-003** in `dev/red-team/findings.md` confirms
Spearman(u_post, u_true) ≈ −0.1 under every k'-prior arm. The MCMC's
per-character u_post numbers are dressed-up prior summaries. At
150 tips and beyond we cannot afford to (a) keep paying the per-cycle
cost of the case-25 Gibbs sweep for non-information, or (b) report those
numbers to users.

The fix is **Rao-Blackwellisation**: drop k'_i from the sampled state
and replace the per-character likelihood with

    L_marg(y_i | tree, μ, p) = Σ_{u=0}^{u_max} L(y_i | tree, μ, kObs_i + u) · P(u | p)

The posterior is then over (tree, μ, σ, p) only — the same marginal-on-θ
target as the joint sampler in expectation, with the slow discrete
coordinate removed.

The plan below is a **mode switch**, not a replacement
(`likelihoodMode = c("sampled_k", "marginal_k")`), starting with the
**geometric** arm only. The other three Mk' arms ride the same
evaluator and land in follow-ups once geometric mode is validated.

## §1 — Why "geometric only" for v1

Five arms exist (`R/MkPrimeModel.R:152–155`): `geometric`,
`beta_geometric`, `empirical_geometric`, `logseries`, plus the degenerate
MkNT baseline (k = kObs, no hyperparameter — out of scope, nothing to
marginalise).

`geometric` wins for v1 because:

- **Cleanest math.** `P(u | p) = p (1-p)^u` is one-parameter and exactly
  normalised on `u ∈ {0, 1, 2, …}`. No `Z_i(p)` to worry about (Model B
  conditional / Model A unconditional both reduce to a flat-on-u sum
  weighted by the same geometric pmf in the marginal evaluator — the
  Model A vs B distinction reappears in the **truncation rule on the
  prior**, not the per-character marginal arithmetic).
- **Model A SBC works.** `dev/red-team/proofs/kprime-priors.md` §1 and
  `dev/notes/2026-05-28-kprime-viability.md` T3/T6 establish that under
  unconditional `k' ~ 2 + Geo(p)`, p is identifiable from
  N_TIP = 8 × N_CHAR = 30, MLE p̂ unbiased, log BF over flat = 14–259 nats.
  Mini-SBC at N_SIM=50 passed both `p` and `kPrime_pooled` rank tests
  (`T6-modelA-recovery-demos`, AD p = 0.038 and 0.011). So the marginal-k
  geometric mode comes pre-equipped with a validated SBC harness.
- **Each other arm carries an open arithmetic or sampler question.**
  `empirical_geometric` just landed EG-001 (per-character Z_i(p) fix,
  2026-05-28) and its Model A SBC is pending Hamilton submission.
  `beta_geometric` just landed its (s, r) reparameterisation 2026-05-28
  and its T9 SBC is in-flight (Hamilton job 17306210).
  `logseries` carries the latent LS-001 (missing Z_i(c)) — currently
  hidden because c is fixed. Mixing the marginal-k refactor with any of
  those open questions makes the validation lane muddy. Better: get
  geometric marginal-k validated against geometric sampled-k first, then
  port to the other arms once their own SBC is settled.

Defending it in one line: **validate the mechanism against the arm with
the cleanest math and a coherent generative joint, then port.**

The plan covers `geometric` end-to-end. For each other arm a §11.x
follow-up sub-section describes the diff against the geometric
implementation.

## §2 — Mode switch, not replacement

Add `likelihoodMode = c("sampled_k", "marginal_k")` to `MkPrimeModel()`,
parallel to `priorVariant = c("conditional", "unconditional")` (which
landed 2026-05-28 — see `R/MkPrimeModel.R:149` and `:168–175` for the
precedent shape).

Rationale:

- **Comparability.** Posterior-overlap on `tree_length` / `p` between
  sampled-k and marginal-k on the same dataset is a load-bearing
  validation step (§7). That requires both modes coexist in one build.
- **Bit-reproducibility of historical runs.** Existing checkpoints,
  pre-EG-001 production samples, and the SBC v8/v10 results are all on
  sampled-k. Keeping sampled-k as a callable mode means we don't need to
  pin a git tag to reproduce them.
- **Bisection on regressions.** If a marginal-k run shows surprising
  posterior behaviour at scale, switching to sampled-k on the same model
  / data / chain config is the natural A/B.
- **Maintenance burden is small.** The marginal-k path piggybacks on
  the existing case-25 batching infrastructure (see §3 + §4) — there's
  no parallel pruning kernel to maintain.

The state vector schema does **not** shrink under marginal mode —
`state->kPrime` remains as a no-op field initialised to `kObs`. The
operational win is "no slow discrete coordinate to mix," not schematic.
This is worth flagging in the user-facing release note so the
motivation ("drop k'_i from the MCMC state vector") is honoured in
spirit if not strictly in struct layout.

## §3 — Architectural insight: case 25 IS the marginal evaluator (minus one step)

The critical reading from `src/mcmc.cpp:3694–4082`
(`gibbs_kprime_sweep_impl`) is that **phase 1 of the Gibbs sweep already
produces exactly the per-(char, k') log-weights the marginal evaluator
needs.** Phase 2 (lines 4084–4122) is the only step that differs:

- **Phase 1 (lines 3906–4082).** Loops `ko = 0 .. K_MAX_CAND-1`. For
  each `ko`, computes site-pattern likelihoods for all transformational
  characters at `k = kObs + ko`, applies ascertainment correction
  (`logAscCorr`), applies Mk' relabelling
  (`mk_prime_relabel_log(k, kObs)`), and stores
  `charLogW[ti * K_MAX_CAND + ko] = β · ll + logPrior_k`. Includes:
  - JC-lumpability collapse for `ko ≥ 2` via
    `pruning_jc_acrv_persite_collapsed` (line 3987) — the per-k'
    evaluations get the same speedup that landed in PR #2;
  - Pattern dedup via `partAct[pi].patTrans` — characters with identical
    tip columns share one pruning call;
  - Prior-ceiling pre-filter (M-164) at `LOG_CUTOFF = -25.0` —
    characters whose `(β · LL_best + logPrior_k)` falls below the running
    max cutoff are terminated, so per-char `nCand` is typically ≪
    K_MAX_CAND;
  - ACRV via `cpp_acrv_rates`;
  - Het via `pruning_f81_het_acrv_persite`.

- **Phase 2 (lines 4084–4122).** Categorical sample from the `charLogW`
  weights, write `state->kPrime[gi] = kObs_i + chosen`.

**The marginal-k evaluator replaces phase 2 with**

    for each transformational char ti:
      log_L_marg[ti] = logSumExp(charLL[ti, 0..nCand-1])
    return sum(log_L_marg)

where `charLL` is `charLogW` without the prior or `β` factor — see §3.1.

### §3.1 — Refactor shape (no code yet, just signatures)

Extract phase 1 into a helper that returns the **raw** per-(char, ko)
log-likelihoods, not the weights `β · LL + logPrior_k`. Suggested
signature (illustrative only):

```cpp
// New helper in src/mcmc_likelihood.cpp.
// Populates charLL[ti * K_MAX_CAND + ko] for ko in 0..maxKo-1 (or up to
// early termination), and charNCand[ti] = number of evaluated k' values.
// Same batched / collapsed / dedup / cutoff infrastructure as case 25
// phase 1, but no prior or β baked in.
void compute_per_kprime_log_lik(
    McmcData* data, McmcState* state,
    const NumericVector& edgeLen,
    const NumericVector& acrvRates,
    double cutoffMaxLL,        // for prior-ceiling pre-filter
    int    K_MAX_CAND,
    std::vector<double>& charLL,    // length nTrans * K_MAX_CAND
    std::vector<int>&    charNCand  // length nTrans
);
```

Case 25 (sampled-k) calls this, then runs phase 2.

The new marginal-k evaluator calls this, then runs

```cpp
double log_L_marg = 0.0;
for (int ti = 0; ti < nTrans; ++ti) {
  double mx = -INF;
  for (int c = 0; c < charNCand[ti]; ++c)
    if (charLL[ti*K_MAX_CAND + c] > mx) mx = charLL[ti*K_MAX_CAND + c];
  double s = 0.0;
  int ko = 0;
  for (int c = 0; c < charNCand[ti]; ++c, ++ko)
    s += std::exp(charLL[ti*K_MAX_CAND + c] - mx)
         * std::exp(logPriorByU[ko]);
  log_L_marg += mx + std::log(s);
}
```

For neomorphic chars and any chars with `kPrime` known/fixed, fall
through to the existing `cpp_partition_log_likelihood` path — they are
not marginalised.

### §3.2 — Where the chain currently calls `cpp_log_likelihood`

Every non-case-25 move that needs a likelihood ratio. Grep for
`cpp_log_likelihood` / `cpp_log_likelihood_partitioned` in
`src/mcmc.cpp` (≈ lines 3243, 3325, 4132, 4206, 4218 — and the proposal
sites in `mcmc_likelihood.cpp:2786`). All of these currently pass
`state->kPrime` through. Under marginal mode they must instead dispatch
to a marginal evaluator.

Suggested API: a single `cpp_log_likelihood_marginal(...)` orchestrator
in `mcmc_likelihood.cpp` that loops partitions, calls
`compute_per_kprime_log_lik` on the transformational partitions, and
falls back to `cpp_partition_log_likelihood` (with `kPrime = kObs`,
since u = 0 always) on neomorphic. Then in `mcmc.cpp` introduce a thin
dispatcher:

```cpp
inline double total_log_lik(...) {
  return data->marginalK
    ? cpp_log_likelihood_marginal(...)
    : (data->usePartitioned
        ? cpp_log_likelihood_partitioned(...)
        : cpp_log_likelihood(...));
}
```

Every existing call-site changes from the explicit dispatcher to
`total_log_lik(...)`. About 6 sites.

## §4 — Critical files to touch

### C++ (`src/`)
- **`src/mcmc_state.h`**.
  - `McmcData`: new `bool marginalK = false;` flag.
  - `McmcState`: per-(partition, node, k_offset) cache slot for partial CLs
    (see §5 for cache strategy). New field
    `PerKpNodeCl perKpCl;` or extend `nodeCL` with a k-axis.
- **`src/mcmc_likelihood.cpp`**.
  - `prepare_mcmc_data(...)`: add `bool marginalK` parameter, plumb to
    `data->marginalK`.
  - `cpp_log_prior(...)` (lines ~217–460): under `marginalK == true`,
    drop the per-character `k'_i` prior term entirely; only the
    hyperprior on `p` (Beta) is kept. The k'_i prior is consumed by the
    marginal-likelihood sum, not the prior.
  - `cpp_log_likelihood(...)`, `cpp_log_likelihood_partitioned(...)`:
    not changed directly; a new sibling `cpp_log_likelihood_marginal`
    lives alongside them and dispatches per partition.
  - **New: `cpp_log_likelihood_marginal(...)`.** Per partition,
    transformational → call `compute_per_kprime_log_lik` then logSumExp;
    neomorphic → forward to `cpp_partition_log_likelihood` unchanged.
- **`src/mcmc.cpp`**.
  - Extract phase 1 of `gibbs_kprime_sweep_impl` (lines 3779–4082) into
    the new `compute_per_kprime_log_lik` helper (§3.1).
  - Case 25 (`gibbs_kprime_sweep`, line 4630) calls the helper, then
    runs the old phase-2 categorical sampling. Behaviour unchanged
    under sampled-k mode.
  - **Move-weight dispatch** at the .BuildMoves layer (R side, see
    `R/RunMkPrime.R` move-weight assembly): under marginal mode,
    register case 25 and case 26 (block_kprime_shift, line 4633) at
    weight 0. They become no-ops. R-side also drops case 7 (single-char
    kPrime MH) — already legacy, but check.
  - Replace every `cpp_log_likelihood(*data, ..., state->kPrime, ...)`
    /  `cpp_log_likelihood_partitioned(...)` call with a dispatcher
    that picks marginal vs sampled by `data->marginalK`. About 6 sites.
  - **Case 30 (`mh_logit_p`, line 4643).** Under marginal mode, a p-move
    changes every transformational character's likelihood (because
    `P(u|p)` shifts globally). Cache invalidation: any p-move must
    re-run the marginal evaluator from scratch (or re-weight a cached
    per-(node, k') CL — see §5). Under sampled-k mode this move only
    touches the prior, since kPrime is the relevant state variable.
- **`src/RcppExports.cpp`** + **`src/mcmc_likelihood.cpp`** export tables:
  one new R-exported entry point if we want
  `cpp_log_likelihood_marginal_xptr` for unit testing.

### R (`R/`)
- **`R/MkPrimeModel.R`**.
  - `MkPrimeModel()` constructor: add `likelihoodMode = c("sampled_k",
    "marginal_k")` argument, default `"sampled_k"` for v1.
    `match.arg`, store in returned list. Validation: error or warn if
    `likelihoodMode == "marginal_k"` is combined with
    `qHeterogeneity = TRUE` or `kPrimePrior != "geometric"` (v1 scope).
  - `LogPrior()` (lines 518–745): under `marginalK == TRUE`, skip the
    geometric-arm `u` term at lines 595–597; keep only the
    `dbeta(p, ...)` hyperprior.
  - `print.MkPrimeModel()`: print mode if non-default.
- **`R/RunMkPrime.R`**.
  - Plumb `model$likelihoodMode` into the C++ `prepare_mcmc_data` call
    (the new `marginalK` boolean).
  - `.BuildMoves` / move-weight assembly: under marginal mode, set
    weight 0 on move types 25 and 26 (and 7, if still registered);
    redistribute that weight onto case 30 (`mh_logit_p`) since `p` now
    dominates the chain's exploration of u_max.
  - Trace column construction (line 325 comment lists the format):
    under marginal mode, omit the `kPrime_i...` columns entirely. The
    trace becomes `iter, swap_cold, topo_hash, tree_length, ...`
    without the per-char kPrime columns. Downstream
    `R/Convergence.R:55–57, :309–351` filter on `kPrime_` prefix — those
    code paths still work; they just yield an empty kPrime group under
    marginal mode. Add a one-line "marginal mode — kPrime not sampled"
    note in `Convergence` output.
- **`R/likelihood.R`** (lines 44–215).
  - `MkpLogLikelihood()` user-facing API: add `likelihoodMode` argument,
    or detect it from a passed `model` if the model is supplied. Adds a
    code path that returns the marginal likelihood for direct user
    queries.
- **`R/MkPosterior.R`** + **`R/Convergence.R`**.
  - Filter `kPrime_*` columns out of `scalarParams` under marginal
    mode (they don't exist in the trace).
- **`R/Checkpoint*.R`** — see §9.
- **Tests under `tests/testthat/`**:
  - New `test-marginal-k-geometric.R` (the bit-identity-strong unit
    test of §7 + targeted regression tests).

## §5 — Cache invalidation strategy

This is the make-or-break decision. Two options:

### Option A — per-(partition, node, k_offset) partial CL cache  (RECOMMENDED)

Extend `state->nodeCL` (currently per-(partition, node)) with a k-axis,
so that for each transformational partition we cache the Felsenstein
partial likelihood at every internal node for every candidate `ko ∈
0 .. K_MAX_CAND-1`. Size grows by a factor of `K_MAX_CAND` (≤ 50, often
≪ 50 thanks to early termination).

Cache invalidation rules:

| Move | Action |
|---|---|
| NNI / SPR / TBR (topology) | Invalidate affected subtree across **all** k-slices for the affected partitions. Existing partial-CL invalidation logic just iterates over the new k-axis. |
| Branch-length moves (cases 4, 12, 23, 24, etc.) | Invalidate affected subtree across all k-slices for affected partitions. |
| `tree_length` scale (case 2) | Invalidate everything. (Same as today.) |
| `rate_log_sd` (cases 31, 34) | Invalidate everything. (Same as today.) |
| `rate_neo` (case 33) | Invalidate neomorphic partitions only (transformational unaffected). |
| **Case 30 (`mh_logit_p`)** | **Cache is fully valid.** The per-k' partial CLs do not depend on `p`. A `p`-move only changes the `logPriorByU[ko]` weights consumed by the per-character `logSumExp` step. Re-running just the logSumExp (over the cached `charLL` if we cache those too — see below) is O(nTrans × u_max), no pruning. |

The win is that `p`-moves — which under marginal-k will dominate
mixing of `u_max(p)`-effective parameters — are free. The cost is
K_MAX_CAND× larger CL cache.

**Sizing.** Current `nodeCL` for a 150-tip dataset: ≈ 150 internal
nodes × ≈ 300 chars × (k = 2..) double-precision arrays per partition.
With K_MAX_CAND ≈ 50 (worst case; early termination cuts it to ~5–15
for typical p), the multiplier is bounded above by 50× and typically
~10×. For a 150-tip × 300-char dataset that's ~ 600 MB in worst case;
we should size the workspace lazily so only actually-evaluated k' slices
allocate.

**Sub-option A'.** Cache `charLL[ti, ko]` (the per-character per-k'
log-likelihoods after ascertainment + relabelling) instead of the
per-node partial CLs. Smaller (`nTrans × u_max` doubles per chain) but
needs full re-pruning whenever any branch / topology / rate changes —
the cache is invalidated by almost every move except `p`. So it's only
useful as a `p`-move accelerator. **Use both:** node-level per-k' CLs
for the main path (Option A), plus `charLL[ti, ko]` derived from them.
A `p`-move re-runs only the per-character logSumExp from the cached
`charLL`.

### Option B — full re-evaluation (no per-k' cache)

Marginal evaluator is called from scratch on every move. Acceptable for
small datasets (≤ 25 tips). At 150 tips this is ~5× slower than Option
A on tree/branch moves and ~10–50× slower on `p`-moves. Not viable for
the user's headline 150-tip target.

**Plan: implement Option A.** Sub-option A' (with the `charLL` cache for
`p`-move acceleration) is a cheap bolt-on that should land in the same
PR.

## §6 — Adaptive u_max(p) per chain step

Geometric `P(U ≥ u | p) = (1-p)^u`. Want the truncation tail mass below
ε = 1e-6 (or 1e-8 if we want to be paranoid). Solve for u_max:

    u_max(p) = ceil( log(ε) / log(1-p) )

| p      | u_max(p) at ε = 1e-6 |
|--------|----------------------|
| 0.1    | 132                  |
| 0.3    | 39                   |
| 0.5    | 20                   |
| 0.7    | 12                   |
| 0.9    | 6                    |
| 0.99   | 3                    |

The case-25 hard cap is `K_MAX_CAND = 50` plus a `LOG_CUTOFF = -25.0`
prior-ceiling pre-filter, which together adaptively prune unused
candidates. **Reuse this directly**: at small p the tail extends past
50 but the prior-ceiling cutoff (which combines per-character LL with
the prior) will terminate most characters well before u = 50 because
the LL of "this character has 30 more hidden states than were observed"
is far below the optimum.

**Per-character adaptivity.** kObs varies per character (range 2–20).
For chars with kObs already large, the marginal LL is sharp around
u = 0 — likelihood pushes `P(u > 0 | y_i, tree)` toward 0 very fast.
The M-164 prior-ceiling pre-filter handles this by terminating
characters as soon as their candidate weight falls below the
running-max cutoff. Confirmed by inspection of lines 3927–3961.

**Concrete plan.** Keep `K_MAX_CAND = 50` as the hard cap (matches the
sampled-k path so the validation lane is symmetric). Keep
`LOG_CUTOFF = -25.0`. If the cap is binding for some character at
small p, raise K_MAX_CAND with a build-time `#define` — this is
exposed in the same header for both paths so they stay in sync.

For the analytic tail, the marginal evaluator should add a per-character
**tail correction** to compensate for the truncated mass:

    log L_marg(y_i) ≈ logSumExp(LL_evaluated + logPrior_evaluated) + log_tail_corr_i

where `log_tail_corr_i ≈ log(1 + tail_mass · max_LL_in_tail / sum_evaluated)`.
For the conservative cutoff this is < 1e-6 nats per char — negligible
for headline LL but worth adding for SBC bit-identity. Defer to v1.1
unless the bit-identity test (§7) fails.

## §7 — Validation: detailed-balance / correctness

Rao-Blackwellisation preserves the marginal posterior on (θ_marginal) =
(tree, μ, σ, p) by construction. Bit-identity to sampled-k is **not**
achievable — different state space. The validation strategy uses three
independent rungs, in increasing cost order. All live under
`dev/red-team/heavy-tests/marginal-k/` (new subdir).

### §7.1 — Rung 1: bit-identical sum-of-products test (run first)

Direct test of the marginal evaluator's arithmetic, independent of the
sampler. For fixed (tree, μ, p, kObs vector):

    L_marg_native = cpp_log_likelihood_marginal(tree, μ, p)
    L_marg_naive  = logSumExp_{u_vec} [
                      cpp_log_likelihood(tree, μ, kObs + u_vec)
                      + Σ_i log P(u_i | p)
                    ]

where `u_vec` ranges over a product grid for small `nTrans` (say
nTrans ≤ 4, u ∈ {0, 1, …, 6}; 7^4 = 2401 evaluations). These must
agree to numerical tolerance (1e-9 in log-LL on a 4-char × 8-tip
synthetic dataset). Discriminates between

- correct evaluator arithmetic;
- correct ascertainment/relabelling-correction placement (inside vs
  outside the marginal sum — they must be **inside**, per character per
  k');
- correct tail-correction (§6) under tight cutoffs.

Targeted file: `tests/testthat/test-marginal-k-geometric.R`. Move into
heavy-tests if the brute-force loop is slow.

### §7.2 — Rung 2: posterior overlap, sampled-k vs marginal-k

For each of several synthetic datasets (kObs mix, n_tip ∈ {8, 16, 25},
n_char ∈ {16, 32, 50}), run both modes with matched seeds:

- N_iter = 12000, N_warm = 4000;
- KS test on marginal posteriors of `tree_length`, `rate_log_sd`, `p`;
- KS test on marginal RF distance of the topology;
- pass bar: KS p > 0.01 on every parameter in every cell.

Heavy-test driver: `dev/red-team/heavy-tests/marginal-k/T-OVL-sampled-vs-marginal.R`.

This is the primary correctness check.

### §7.3 — Rung 3: SBC under Model A on (tree_length, rate_log_sd, p)

Reuse the existing geometric Model A SBC harness — it already passes
under sampled-k (see `dev/notes/2026-05-28-kprime-viability.md` T6 and
the `experiment/model-a-prior` worktree mini-SBC). Flip the
`likelihoodMode = "marginal_k"` switch and re-submit.

Expected outcome: AD p > 0.4 on `tree_length`, `rate_log_sd`, `p` at
N_SIM = 200, N_iter = 12000, N_tip = 8, N_char = 30 — matching or
exceeding the sampled-k Model A result. (The `kPrime_pooled` rank test
is irrelevant under marginal-k since kPrime is no longer sampled.)

Hamilton driver: `dev/red-team/heavy-tests/marginal-k/T-SBC-marginal-geometric.R`
plus `dev/red-team/heavy-tests/submit-marginal-k-sbc.sh` (single 8h
task, mirrors the existing geometric Model A submit script). Per the
`feedback_pkgload_prebuild` memory, pre-build the package on the login
node before sbatch (and use `${SRC}/dev/red-team/...` paths per
`feedback_slurm_inscript_path`).

### §7.4 — Cross-cutting smoke tests

- `cpp_log_prior` under marginal mode returns the same value as
  sampled-k with `kPrime = kObs` — that is, the per-character u-term
  drops out and only the Beta(p) hyperprior remains.
- `log_likelihood + log_prior` at any (tree, μ, p) under marginal mode
  equals `logSumExp over u of [LL(tree, μ, kObs+u) + log P(u | p)]` +
  hyperprior — same as §7.1 but checks the joint posterior.
- Trace columns under marginal mode contain no `kPrime_*` and the
  filter logic in `R/Convergence.R` produces the expected empty group
  without error.

## §8 — Performance estimate (rough)

Per-cycle cost (rough, order-of-magnitude):

| Component | Sampled-k | Marginal-k (Option A) |
|---|---|---|
| Per-non-25-move likelihood eval | 1× pruning per char | ~5–15× pruning per char (cached after first eval) |
| Case 25 (kPrime sweep) | full case-25 cost | **removed** |
| Case 30 (p-move) | scalar prior eval | scalar weighted re-sum over cached per-k' LLs (~ nTrans × u_max ops, no pruning) |
| Net per cycle | A | A × 5–15 minus the (substantial) case-25 cost |

Naive accounting suggests marginal-k is 2–10× slower **per cycle**.
But ESS/sec is the right metric, and the predicted gain comes from:

- Eliminating an O(nTrans) discrete-coordinate sweep whose ESS on
  `tree_length` is bounded by the worst-mixing per-character kPrime
  (per Rao-Blackwell intuition: the marginalised chain mixes at least
  as well on `tree_length` / `p` as the joint chain projected to those
  coordinates);
- Removing the auto-correlation that the case-25 sweep induces between
  `tree_length` and `p` (every kPrime move shifts the implied per-char
  rate, dragging tree_length on the next cycle);
- Eliminating the slow-coordinate trap when one character's kPrime
  walks off into the tail and the chain spends thousands of cycles
  shifting it back.

**Quantitative expectations are hedged.** The expected wall-clock-to-
converged-tree on 150 tips is "comparable or better than sampled-k"
with a non-trivial probability of being substantially better
(2–10× ESS/sec gain on `tree_length`, `p`). Cannot put a hard number
without running it.

**Benchmark plan.** After §7.2 passes, run a head-to-head on the
Casali / Allain2012 / AllainAquesbi2008 cells that
`project_hyperprior_ess_bench` already covers — same datasets, same
N_iter, both modes, measure ESS/sec on (tree_length, p, rate_log_sd).
Drop into `dev/red-team/heavy-tests/marginal-k/T-BENCH-ess-vs-sampled.R`.

## §9 — Checkpoint migration

`R/MkPrimeMCMC.R:78–217` describes the checkpoint format. The current
checkpoint state vector includes `kPrime` (an `IntegerVector` per
transformational character).

Under marginal mode:
- **Write side.** Checkpoint serialisation can omit `kPrime` (it's
  uninteresting — always `= kObs`). Add a `likelihoodMode` field to the
  checkpoint header so resume can detect mismatch.
- **Read side.** `ResumeMkPrime()` must check
  `checkpoint$likelihoodMode == model$likelihoodMode` and abort with
  a clear error if they disagree. Cross-mode resume is **not**
  supported — different state spaces, no honest projection.
- **Backwards-compat.** Pre-marginal checkpoints have no
  `likelihoodMode` field; treat absent as `"sampled_k"`. Old
  sampled-k checkpoints resume into sampled-k mode as before, no
  migration needed.
- **Mode-mismatched resume.** If the user calls `RunMkPrime(...,
  likelihoodMode = "marginal_k")` and finds a sampled-k checkpoint, the
  error message should suggest either (a) deleting the checkpoint
  (`overwrite = TRUE`), or (b) keeping `likelihoodMode = "sampled_k"`.

## §10 — `dev/STATUS.md` update sketch

A new top-level §marginal-k section (between `## Reference: MkNT` and
`## Maintenance`) once landed:

```
## marginal-k (likelihood-mode flag)

**Implements.** `MkPrimeModel(likelihoodMode = "marginal_k")` swaps the
sampled-k k'_i state for an analytic sum over u_i ∈ {0, .., u_max(p)},
using the case-25 batched per-(char, k') likelihood machinery wrapped
with logSumExp instead of categorical sampling. v1 supports the
`geometric` arm only; other arms accept the flag and error out with a
"not implemented in v1" message.

**Motivation.** EG-003: per-character u_post is prior-dominated at
empirical tree lengths (Spearman(u_post, u_true) ≈ -0.1). The k'_i
state is a slow discrete coordinate that carries no information at
scale (≥ 50 chars, mean edge ≤ 0.05). Removing it via
Rao-Blackwellisation preserves the marginal posterior on (tree, μ, σ,
p), removes the case-25 sweep cost, and (per Rao-Blackwell intuition)
should improve ESS/sec on tree_length / p.

**Proof.** Rao-Blackwell preserves the marginal on θ \ k by
construction. Bit-identity test at fixed (tree, μ, p) confirms the
evaluator arithmetic against a brute-force sum-of-products at
nTrans ≤ 4 (tol 1e-9). Heavy-tests at `dev/red-team/heavy-tests/
marginal-k/`.

**Testing.**
- Bit-identity unit test (§7.1): PASS / FAIL / hash.
- Posterior overlap vs sampled-k (§7.2): KS p > 0.01 on
  (tree_length, p, rate_log_sd) on 9 cells (n_tip × n_char grid).
- Model A SBC (§7.3): AD p > X on (tree_length, p, rate_log_sd) at
  N_SIM = 200, N_iter = 12000, N_tip = 8, N_char = 30.

**Empirical performance.** Hyperparameter level: matches Model A
sampled-k (T6 reference). Per-character level: per-character u not
sampled — explicit non-feature. If users want per-character u, use
sampled-k mode.

**Evaluation.** Default mode for the geometric arm at scale (≥ 100
tips). Sampled-k retained as a callable mode for comparability and
backwards-compat.
```

Each existing arm's "Empirical performance" / "Evaluation" subsection
should add a one-line cross-reference to §marginal-k once the
respective arm gets ported. For geometric: "v1 marginal-k mode lives
at §marginal-k; default is sampled-k until that mode's SBC suite
lands."

## §11 — Follow-ups for the other arms

Once §1–§10 land, each other arm is a small diff against the geometric
implementation. The work per arm:

### §11.a — empirical_geometric

The marginal sum becomes

    L_marg = Σ_{m = max(2, kObs_i)}^{m_max} L(y_i | tree, μ, m) · P_trunc(m | p, kObs_i)

where `P_trunc` is the per-character truncated EG pmf (computed by
`.LogPriorEmpiricalGeometric` and `cpp_log_prior`'s EG branch lines
299–407 — the `egLogPriorByK` precomputation already used by case 25 at
3819–3861 is exactly what we need). The only new arithmetic is that
the per-(char, m) prior weights become character-specific (truncation
floor differs), where geometric's are character-independent in `u`.
Reuse case-25's `egLogPriorByK` table verbatim. Honours `priorVariant`.

### §11.b — beta_geometric

`P(u | α, β)` is BetaGeo, table already precomputed at lines 3767–3777
(`bgLogPrior[]`). Marginal sum is `Σ_u L(y_i | kObs+u) · BetaGeo(u; α, β)`.
α and β are still sampled (the marginal-k mode does not eliminate
those — they have a hyperprior of their own). The `(s, r)` slice on
(α, β) lives unchanged. The case 30 `p`-move is replaced by the
existing α / β moves; their cache-invalidation behaviour mirrors
case 30 (all per-(char, k') prior weights change globally, but the
cached per-(char, k') LLs are unchanged).

### §11.c — logseries

`P(k' | c)` is Logseries with the same truncation pattern as EG. `c`
is currently fixed, so under marginal-k the prior weights are
chain-constant — Option A's `charLL` cache is enough; no analogue of
case 30 needs to fire. Landing logseries marginal-k is a strict
performance win at any scale. **Must close LS-001 first** (per the
arm-level §logseries entry in `STATUS.md`) since the marginal sum
needs the correctly normalised per-character truncation; but only if
we also want to enable a sampled `c` later.

## §12 — Risks and unknowns

- **Cache memory at 150 tips × 300 chars × K_MAX_CAND = 50.** Worst-case
  600 MB. Need lazy allocation of k-slices, and fail-loud with a clear
  message if the workspace exceeds (say) 4 GB. Probably need a profvis
  / VTune pass once a 150-tip benchmark exists.
- **Het + marginal-k.** Defer. v1 errors out if
  `qHeterogeneity = TRUE && likelihoodMode = "marginal_k"`. Het
  multiplies the per-evaluation cost by ~nBetaCat ≈ 4; the case 25
  Het path exists (line 3978) so the marginal evaluator can adopt it,
  but it 10×s the cache memory in the worst case. Out of v1 scope.
- **Partition-API + marginal-k.** Defer. The `cpp_log_likelihood_
  partitioned` path is what mixed-character-type datasets need. The
  marginal evaluator must respect the per-partition rate (classRate)
  scaling. Adding marginal-k to the partition API is mechanically
  straightforward (just thread `classRate[]` through the per-(char,
  k') pruning calls — case 25 already does this), but the per-class
  prior on σ_c interacts with cache invalidation in ways that need
  re-thinking. Defer to v1.1 after geometric marginal-k validates.
- **Numerical conditioning of the per-character logSumExp.** Should
  be fine; charLL values for a single character span ≤ 25 nats by
  construction (the M-164 cutoff). Subtract max-per-char before
  exp/sum/log. Add an explicit numerical-auditor pass under
  `dev/red-team/numerical/` after v1 lands (per the red-team agent
  rotation memory).
- **Tail correction (§6).** v1 truncates at K_MAX_CAND = 50 + the
  M-164 cutoff. If the bit-identity test (§7.1) shows the truncation
  bites, add an analytic tail correction. Quantify in the §7.1 driver.

## §13 — Out of scope for v1

- Other Mk' arms (beta_geometric, empirical_geometric, logseries) —
  see §11.
- Q-matrix heterogeneity (`qHeterogeneity = TRUE`) — see §12.
- Partition-API (`usePartitioned = TRUE`) — see §12.
- Per-character u_post posterior output — explicit non-feature.
  Sampled-k mode is the answer when the user needs per-char u_post.
- Cross-mode checkpoint resume — see §9.

## §14 — Suggested PR sequencing

One refactor PR + one mode-add PR keeps the diff reviewable:

1. **PR-A (refactor only, bit-identical).** Extract `compute_per_kprime
   _log_lik` from `gibbs_kprime_sweep_impl` phase 1. Case 25 calls it
   then runs phase 2. Behaviour, RNG draws, and bit-output unchanged.
   Existing test suite must pass unchanged. **No new mode flag yet.**
2. **PR-B (marginal-k geometric arm).** Add `likelihoodMode` flag to
   `MkPrimeModel()` and plumb to C++. Implement
   `cpp_log_likelihood_marginal` (geometric only). Implement Option A
   cache. Add §7.1 bit-identity test. Add §7.2 posterior-overlap
   heavy-test driver. Update STATUS.md per §10. Disable case 25 / 26 /
   7 in marginal mode.
3. **PR-C (Hamilton SBC validation).** Submit §7.3 SBC, land results
   into STATUS.md. Per `feedback_pkgload_prebuild` and
   `feedback_slurm_inscript_path` — pre-build on the login node, use
   `${SRC}/dev/red-team/...` paths.

Per project memory `feedback_no_prs.md`, "PR" here means **branch +
direct commit on a worktree branch**; Mk-prime/r is private with GHA
disabled and we do not open PRs. Per `feedback_no_branch_switch_main_
checkout.md`, do the implementation in a worktree (e.g.
`C:\Users\pjjg18\GitHub\mkp.worktrees\marginal-k`), not on the shared
main checkout.

## §15 — One-paragraph elevator pitch for the implementer

The case-25 `gibbs_kprime_sweep_impl` in `src/mcmc.cpp:3700–4122`
already computes every per-(char, k') log-likelihood the marginal
evaluator needs — it just samples categorically at the end instead of
doing a logSumExp. Refactor phase 1 of that function into a helper,
add a `likelihoodMode = "marginal_k"` flag on `MkPrimeModel()`,
implement the new `cpp_log_likelihood_marginal` (geometric arm only
for v1) on top of the helper, extend the partial-CL cache with a k-axis
so `p`-moves stay cheap, disable case 25 / 26 / 7 in marginal mode,
and validate against a bit-identical sum-of-products test, a
sampled-vs-marginal posterior overlap test, and the existing geometric
Model A SBC harness flipped to marginal mode. About 6 cpp_log_likelihood
call-sites in `src/mcmc.cpp` need rerouting through a small dispatcher.
Other arms (`beta_geometric`, `empirical_geometric`, `logseries`) port
in follow-ups by swapping the `P(u | hyperparams)` lookup table; the
machinery is the same.
