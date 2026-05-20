# Red-team: state->logLik / state->logPrior drift in ecology-aware MCMC

Scope: state synchronisation in `do_move_impl` and helpers. Worktree
`C:\Users\pjjg18\GitHub\mkp\.claude\worktrees\ecology-aware`.

## Summary

There is a **non-ecology log-likelihood being written into `state->logLik`** by
two Gibbs/MH paths that are dispatched **even when `data->ecologyAware == true`**:

- `gibbs_kprime_sweep_impl` (moveType 25) — `src/mcmc.cpp:3811-3823`
- `block_kprime_shift_impl` (moveType 26) — `src/mcmc.cpp:3895-3927`

Both call `cpp_partition_log_likelihood` (the non-ecology pruning), sum the
parts and overwrite `state->logLik`. The ecology rate-modifier terms
(`phi`, `z`, `pi0`, `gamma_e`) are dropped from the stored value. From that
point on the running `state->logLik` is a non-ecology surrogate, while every
subsequent ecology MH move (cases 34/35/37) computes `newLL` via
`cpp_log_likelihood_ecology` and compares against the corrupted baseline.

Magnitude consistent with smoking-gun: the ecology correction over 480 cells
(84 z=0, 393 z=1, 3 z=2) with phi≈26 contributes hundreds of nats. The
reported +307-nat drift at sample 1 is exactly the order produced by removing
the ecology multiplier from the partition sum.

A second, less-likely-but-related bug: the Gibbs sweep in case 25 also draws
the new k' values from a **non-ecology** full conditional (the per-character
log-weights at `mcmc.cpp:3406-3808` are built from `cpp_partition_log_likelihood`
intermediates, not the ecology pruning). The k' samples are themselves biased,
but the visible accountancy drift comes from the final `state->logLik = newLL`.

## Bug 1 — gibbs_kprime_sweep writes non-ecology logLik

File: `src/mcmc.cpp`
Lines: 3810-3823

```
ClWorkspace* wsPtr = state->clWs.ready() ? &state->clWs : nullptr;
int nParts = (int)data->parts.size();
std::vector<double> newPLC(nParts);
double newLL = 0.0;
for (int pi = 0; pi < nParts; ++pi) {
  newPLC[pi] = cpp_partition_log_likelihood(    // <-- NON-ecology
    *data, pi, state->parent, state->child, edgeLen,
    state->kPrime, state->rateLoss, state->rateLogSd, state->rateNeo,
    state->betaScale, wsPtr);
  newLL += newPLC[pi];
}
state->logLik = newLL;                          // <-- baseline corrupted
state->partLogLik = std::move(newPLC);
```

Move case 25 is dispatched unconditionally (`mcmc.cpp:4478-4480`); there is
no `if (data->ecologyAware)` guard, and there is no parallel
`gibbs_kprime_sweep_eco_impl`. When ecology mode is on, this runs and
silently replaces the ecology log-lik with the bare Mk-prime sum.

**Fix.** In ecology mode, after the k' updates, call
`cpp_log_likelihood_ecology` to set `state->logLik` (matching the gibbs_z
sweep pattern at `mcmc.cpp:4052-4055`); leave `state->partLogLik` empty (the
ecology path does not use it — see `fill_partition_cache` early return at
line 497). The k'-conditional sampling itself should also be reworked to use
the ecology pruning, but the **accountancy fix is to recompute logLik via
the ecology call before returning true**.

**Verification.** Add one printf before line 3822:
```
double freshEco = cpp_log_likelihood_ecology(*data, state->parent, state->child,
  edgeLen, state->kPrime, state->rateLoss, state->rateLogSd, state->rateNeo,
  state->phi, state->zMatrix, state->pi0, state->theta);
REprintf("kprime_sweep: partition=%.2f eco=%.2f diff=%.2f\n",
         newLL, freshEco, freshEco - newLL);
```
Expect a few-hundred-nat gap whenever phi differs meaningfully from 1.

## Bug 2 — block_kprime_shift mirrors the same mistake

File: `src/mcmc.cpp`
Lines: 3889-3927 (path active when `hasPLC` true; ecology mode actually has
`partLogLik` empty so `hasPLC` is false, BUT the else-branch at lines 3908-3918
still calls `cpp_partition_log_likelihood`).

In ecology mode `state->partLogLik` is empty (set by `fill_partition_cache`
returning at line 497), so execution falls into the `else` block at
`mcmc.cpp:3908`:

```
} else {
  newPC.resize(nParts);
  newLogLik = 0.0;
  for (int pi = 0; pi < nParts; ++pi) {
    newPC[pi] = cpp_partition_log_likelihood(*data, pi, ...);   // <-- WRONG
    newLogLik += newPC[pi];
  }
}
```

The MH acceptance ratio at line 3921 then compares this non-ecology newLogLik
against the (currently still-correct) `state->logLik` ecology value, producing
a bogus ratio; on acceptance the corrupted value is written to
`state->logLik` at line 3925.

**Fix.** Branch on `data->ecologyAware`: if true, call
`cpp_log_likelihood_ecology` for `newLogLik`, and do not populate
`state->partLogLik` on acceptance.

## Bug 3 — pre-proposal drift diagnostic is also non-ecology

File: `src/mcmc.cpp`
Lines: 4122-4126

```
double freshLL = cpp_log_likelihood(*data, state->parent, state->child, ...);
double drift = std::abs(state->logLik - freshLL);
```

In ecology mode this compares the (correct-or-corrupted) ecology logLik to a
non-ecology baseline, so `diagDriftCount` always reports drift in eco mode
and is useless for catching the real bug. Switch on `data->ecologyAware` and
use `cpp_log_likelihood_ecology` when true.

## Looks correct

- Cases 34 (scale phi), 35 (scale pi0), 37 (scale theta): full MH with
  `cpp_log_likelihood_ecology`; on accept assigns `state->logLik = newLL`,
  `state->logPrior = newLP` (mcmc.cpp:4596-4601, 4645-4651). Symmetric MH
  with explicit Hastings — fine.
- Case 36 (`gibbs_z_sweep_impl`): full ecology recompute after the sweep
  (mcmc.cpp:4051-4061). Correct pattern.
- General move dispatcher `!hasPLC` branch (mcmc.cpp:4874-4889): routes to
  `cpp_log_likelihood_ecology` when `eco` is true. Topology + branch-length
  + rate scalar moves go through this path in ecology mode and produce a
  consistent `newLogLik`.
- `init_mcmc_state` (line 444-469): stores the R-computed `log_lik`,
  `log_prior` from `.InitState`. `fill_partition_cache` is a no-op in eco
  mode (line 497) so the R value is preserved. The R initialiser uses
  `.MkpEcologyLogLikelihood` — same kernel as `cpp_log_likelihood_ecology`,
  so an init mismatch of more than a few nats from numeric drift alone is
  unlikely (and the 7000-iter post-init growth points to per-move
  accumulation, not init).
- `wEdgeDirty` is declared at `mcmc.cpp:211` and set once at line 483 but
  never read — `wEdge` is recomputed inside every
  `cpp_log_likelihood_ecology` and `gibbs_z_sweep_impl` call, so there is no
  stale-cache exposure on topology / branch-length moves.

## Could-not-verify

- The exact frequency of moveTypes 25 / 26 in the truth-init schedule —
  if they run only every K iterations, a single sweep contaminates
  everything downstream. The `moveWeights` table in `RunMkPrime.R` and the
  schedule that drives `do_move_cpp` would tell you how many times bug 1
  fires during the 7020-iter warmup.
- Whether `block_kprime_shift` was actually invoked (it requires
  `transIdx > 0` and a non-zero `intWalkWindow`); but bug 1 alone is
  sufficient to produce a one-shot ~300-nat shift.
- No periodic resync logic exists. The only consistency check
  (`diagDriftCount` at 4126) is broken for ecology mode (bug 3) and does
  not re-sync even when triggered.

## Recommendation

Add an `if (data->ecologyAware)` branch to both kPrime helpers calling
`cpp_log_likelihood_ecology` for the final `state->logLik` update. Add a
periodic resync block inside `do_move_impl` (every ~500 iter) that
recomputes the full ecology logLik / logPrior from scratch and overwrites
`state->logLik` / `state->logPrior` — at minimum during warmup. This both
hides any residual accumulation bug and turns the existing
`diagDriftCount` into a useful guardrail.
