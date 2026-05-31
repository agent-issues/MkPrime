# MkPrime arm status

Per-arm implementation / proof / validation status for the Mk' family.
Updated as arms progress through implementation, proofs, SBC, and sampler
tuning.

Note: this page is created in the `feat/marginal-k` worktree to cover the
v1 marginal-k landing. Cross-merge with `main`'s existing STATUS.md (when
it lands) by adding the §marginal-k section below; pre-existing per-arm
sections (`geometric`, `empirical_geometric`, `beta_geometric`,
`logseries`) are not modified here.

## marginal-k (likelihood-mode flag)

**Implements.** `MkPrimeModel(likelihoodMode = "marginal_k")` swaps the
sampled-k k'_i state for an analytic sum over u_i ∈ {0, .., u_max(p)},
using the case-25 batched per-(char, k') likelihood machinery
(`compute_per_kprime_log_lik`, extracted in PR-A) wrapped with
`logSumExp` instead of categorical sampling. v1 supports the `geometric`
arm only; combinations with `empirical_geometric` / `beta_geometric` /
`logseries`, `qHeterogeneity = TRUE`, or partition-API
(`usePartitioned = TRUE`) abort at the `MkPrimeModel()` / `RunMkPrime()`
boundary with a "not implemented in v1" error pointing to plan §11/§13
in `dev/notes/2026-05-28-marginal-k-plan.md`.

**Motivation.** EG-003: per-character u_post is prior-dominated at
empirical tree lengths (Spearman(u_post, u_true) ≈ -0.1, see
`dev/red-team/findings.md`). The k'_i state is a slow discrete
coordinate that carries no information at scale (≥ 50 chars, mean edge
≤ 0.05). Removing it via Rao-Blackwellisation preserves the marginal
posterior on (tree, μ, σ, p), removes the case-25 sweep cost, and (per
Rao-Blackwell intuition) should improve ESS/sec on tree_length / p.

**Proof.** Rao-Blackwell preserves the marginal on θ \ k by construction.
Bit-identity test at fixed (tree, μ, p) confirms the evaluator arithmetic
against a brute-force sum-of-products at nTrans = 4, kObs = 2, p = 0.7
to tolerance 1e-6 nats (truncation-bias-limited; the marginal evaluator
runs u ∈ {0..49} via kMaxKprimeCand = 50 while the brute force in
`tests/testthat/test-marginal-k-geometric.R` runs u ∈ {0..11} — tail
mass beyond the brute-force cap is < 5.3e-7 per char, bounding the
discrepancy).

**Cache strategy.** Plan §5 sub-option A' (charLL cache) is implemented:
the marginal evaluator stores per-(char, k') raw LL on `state` after
each rebuild, and `mh_logit_p` (case 30) reuses it (no re-pruning,
O(nTrans × u_max) ops). The full Option A per-(node, k') CL cache is
deferred to v1.x — not in v1.

**Testing.**
- Bit-identity unit test (§7.1, `tests/testthat/test-marginal-k-geometric.R`):
  PASS as of feat/marginal-k landing.
- Posterior overlap vs sampled-k (§7.2): driver at
  `dev/red-team/heavy-tests/marginal-k/T-OVL-sampled-vs-marginal.R`.
  Not yet run; companion sbatch script is Stage 2 work.
- Model A SBC on p (§7.3): **VALIDATED CLEAN (2026-05-31).** Two bugs were
  found and fixed en route: MARGINAL-K-INIT-001 (init-logLik freeze, 8a26872)
  and MARGINAL-K-TRUNC-001 (geometric truncation normaliser, e5ebf90). Post-fix
  local N=400 SBC: pooled Anderson-Darling rank-uniformity p = 0.862 (p),
  0.711 (tree_length), 0.517 (rate_log_sd); low-p (<0.10) gross frac rank<=2
  0.92 -> 0.05. The residual low-p *conditional* rank deflation is a
  CORRECT-posterior shrinkage artifact (conditional sub-binning of SBC ranks is
  not a valid calibration test — SBC guarantees only marginal uniformity), NOT a
  forward/inference mismatch; established by independent review + the C-i check
  below.
- C-i truncation guard (`tests/testthat/test-marginal-k-truncation.R`,
  "marginal-k LL == full uncapped/uncutoff reference at low p"): **PASS.** The
  package marginal-k summation (candidate cap + log-cutoff + analytic Z_A)
  equals the explicit full logSumExp over k' in [kObs,K] to <= 7.3e-12 across
  kObs {2..6} x p {0.005..0.5}; the -25 log-cutoff drops <= 1.9e-11 mass. So
  p < 0.05 is a weakly-identified regime, not a calibration defect. (Caveat:
  this is a forward/inference *consistency* check; a shared-wrong truncation
  model would need a large-data recovery check.)

**Known limit / Stage 1b watch-item.** The C++ candidate cap
`kMaxKprimeCand = 50` (mcmc_state.h) currently exceeds the SBC truncation
`K = 30`, so it does not bite. Raising the model's `K` (e.g. the intended
K = 200 real-data default) ABOVE `kMaxKprimeCand` would silently cap the
marginal numerator while `Z_A` normalises the full [2,K] — reintroducing
MARGINAL-K-TRUNC-001 across a wide p range. Stage 1b must wire `kprimeTruncK`
from the model AND couple `kMaxKprimeCand >= K` (plus a loud `K >= max(kObs)`
guard). The C-i guard test above would catch a regression of this kind.

**Empirical performance.** Hyperparameter-level identity vs sampled-k:
heavy-test PENDING (T-OVL). Per-character u not sampled — explicit
non-feature; sampled-k mode remains for users who need per-char u_post.

**Evaluation.** Default mode for the geometric arm at scale (≥ 100 tips)
once T-OVL and SBC validate. Sampled-k retained as a callable mode for
comparability, backwards-compat, and the per-char-u use case.

**Disabled moves under marginal-k.** Case 25 (`gibbs_kprime_sweep`),
case 26 (`block_kprime_shift`), and case 7 (`int_walk` on `kPrime`) are
dropped from the `.BuildMoves` schedule. Case 30 (`mh_logit_p`) gets the
redistributed weight since `p` now dominates `u_max(p)`.

**Open follow-ups.**
- PR-C: T-OVL + SBC validation on Hamilton.
- v1.x: full per-(node, k') CL cache (plan Option A) once T-OVL passes
  and a 150-tip benchmark identifies the cache as the bottleneck.
- v1.x: marginal-k port to `empirical_geometric` /
  `beta_geometric` / `logseries` arms (plan §11). Each is a small
  diff against the geometric evaluator (swap the `P(u | hyperparams)`
  lookup table).
- v1.x: checkpoint cross-mode resume migration (plan §9). Currently
  abort is the behaviour — a sampled-k checkpoint cannot be resumed
  into marginal-k mode and vice versa; the checkpoint header carries
  no `likelihoodMode` field yet so the check fires on the first call
  to `RunMkPrime(..., likelihoodMode = "marginal_k")` against a stale
  checkpoint.
