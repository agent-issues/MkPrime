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

**Stage 1b — truncation cap K wired from the model (DONE 2026-06-01).**
`MkPrimeModel(kprimeTruncK = K)` (default **200** real-data; SBC pins **30**)
now flows through `.InitMcmcData` → `set_kprime_trunc_k()` (a setter mirroring
`set_branch_bins`, so the 30-arg `prepare_mcmc_data` signature is untouched)
into `McmcData.kprimeTruncK`. The C++ candidate cap `kMaxKprimeCand` was raised
**50 → 256** so the marginal numerator can sum the full support `[2, K]` at the
K = 200 default; `set_kprime_trunc_k` enforces `2 <= K <= 256` and
`.InitMcmcData` enforces `K >= max(kObs)` (friendly abort, else a character has
empty truncated support → −Inf). Cost: `charLogW`/`charLLCache` grow to
`nTrans·256` doubles (~0.6 MB at 300 chars); the heavy per-(node,k′)
`PerKpClCache` is dormant in v1 (forward-decl only), so this does NOT scale that
allocation. Any `K <= 30` result is bit-identical to the old cap (the per-char
cap `nEff = min(nCand, K − kObs + 1)` is unchanged), so SBC at K = 30 is
unaffected. Verified end-to-end:
- New guard `test-marginal-k-truncation.R` ("numerator reaches the full support
  at K=200"): the package marginal at K = 200 equals the explicit logSumExp over
  the full `[2, 200]` (diff 0; a cap-50 numerator would land on the `[2,51]`
  reference, off by ~1.3e-3 ≫ 1e-7 tol → FAIL), with G1/G2 anti-vacuity guards
  and a live-K contrast (K = 200 vs 30 moves the LL ~0.80 nats).
- Finding: the numerator's `k′ > 51` tail is **likelihood-suppressed**
  (`P(data|k′)` decays ~1 nat / 10 states), so the cap's practical impact is
  small for low-kObs characters even though the *prior* mass there is large — a
  useful bound. K matters chiefly via `Z_A(p)` at small p (real-data headroom).
- `RunMkPrime` end-to-end smoke at K = 200 (fixed-topology, marginal_k): PASS.

**Stage 2 — sampled_k geometric prior truncated to match marginal_k (DONE 2026-06-01).**
`likelihoodMode = "sampled_k"` now targets the SAME posterior as `"marginal_k"`
(Rao-Blackwell consistency): the sampled-k geometric **prior** is truncated at K
and renormalised by Z(p), exactly matching the marginal-k normaliser (Model A:
−logZ_A = −log(1−(1−p)^(K−1)), shared across chars; Model B: per-character
−logZ_B,i = −log(1−(1−p)^(K−kObs_i+1))). Changes:
- **Prior** (`src/mcmc.cpp cpp_log_prior` plain-geometric branch + R
  `MkPrimeModel.R LogPrior`): add the truncated renormalised per-character mass
  under both Model A/B (was untruncated, Model-B-only); support guard rejects
  k′ > K under sampled_k. The R and C++ priors mirror each other bit-for-bit.
- **Gibbs conjugacy** (`src/mcmc.cpp` case 9 `gibbs_p`): the truncated geometric
  is no longer Beta-conjugate (Z(p) is p-dependent), so case 9 is guarded to
  return false for the plain geometric; p is sampled via case 30 `mh_logit_p`
  instead (the `.BuildMoves` scheduler switched gibbs_p → mh_logit_p, weight 3,
  mirroring empirical_geometric). The legacy conjugate-Beta gibbs_p tests were
  retired in favour of a non-conjugacy guard test.
- **Gibbs sweep cap** (`src/mcmc.cpp` case 25 `gibbs_kprime_sweep`): the
  always-accept sweep caps its candidate range at `nEff = min(nCand, K−kObs+1)`
  (recomputing maxW over the retained range) so it never draws k′ > K. Cases
  7/26 (MH) reject k′ > K via the −Inf prior — no explicit cap needed.
- **K wiring** (`R/RunMkPrime.R .InitMcmcData`, `R/MkPrimeModel.R`):
  `set_kprime_trunc_k` now fires for the geometric arm under BOTH modes (was
  marginal_k only), with the `K >= max(kObs)` guard and `[2,256]` validation
  broadened to the geometric arm. The C++ struct default `kprimeTruncK` was
  aligned **30 → 200** to match the `MkPrimeModel` default, so a dataPtr built
  directly via `prepare_mcmc_data` (bypassing `.InitMcmcData`) still agrees with
  the R prior; SBC and the truncation tests pin K = 30 explicitly via the model.

Verified:
- **Proof** (`dev/red-team/proofs/marginal-k-sampled-rb-consistency.md`,
  math-prover, rated "Watertight"): summing the truncated sampled-k joint over
  k′ ∈ [kObs,K] equals the marginal-k value per character, both models; every
  formula checked line-by-line against source (no discrepancies); algebra
  confirmed to machine precision in R. Caveat: the multi-character lift is by the
  product-of-sums factorisation (now also checked empirically, below).
- **Deterministic bit-check** (`test-marginal-k-truncation.R`, Stage-2 RB +
  multi-char RB): logSumExp over k′ of the full sampled-k JOINT
  (`eval_log_prior_cpp` + `eval_full_loglik_cpp`) == marginal-k joint to 1e-7,
  single-char (Model A/B × kObs {2,6} × p {0.02,0.08}) AND a 2-character (k′_1,
  k′_2) grid; with G1/G2 anti-vacuity + stale-binary (K-sensitivity) +
  support-cap (k′ = K+1 → −Inf) guards.
- **R/C++ prior parity** (`test-partition-prior.R`, `test-partition-hyperprior.R`):
  R `LogPrior` == C++ `eval_log_prior_cpp` to 1e-10 on the geometric arm.
- **Scheduling / guard** (`test-gibbs.R`, `test-beta-geometric-prior.R`,
  `test-bg-hyperparameter-moves.R`): geometric schedules mh_logit_p not gibbs_p;
  C++ refuses a misrouted gibbs_p under the truncated geometric.
- Full `testthat` suite green except the pre-existing out-of-scope red test below.

**SBC on sampled_k (Stage 2, PENDING Hamilton).** Driver
`dev/red-team/heavy-tests/marginal-k/T-SBC-sampled-geometric.R` (mirrors the
marginal-k SBC with `likelihoodMode = "sampled_k"`; forward truncates at
K_MAX_PRIOR = 30 by rejection-redraw, inference pins kprimeTruncK = 30; pass bar
AD > 0.4 on tree_length / rate_log_sd / p — k′_pooled excluded as a Talts
boundary artefact, project memory `project_sbc_kprime_structural`). Full run is a
Hamilton array job (pre-build on the login node first per `feedback_pkgload_prebuild`);
the posterior-overlap heavy test (T-OVL) is the same class of follow-up.

**Out-of-scope red test (pre-existing).** `test-marginal-k-cache-option-a.R`
"NNI (case 5) refreshes the cache" FAILS (warm − cold = 0.164 nats). Proven
**not** a Stage 1b regression: the failure is bit-for-bit IDENTICAL at K = 30 and
K = 200, so the K-wiring did not cause it. This is the tracked Tier-2 partial-CL
warm≠cold bug (a topology move does not refresh the per-(node,k′) cache);
SBC-irrelevant because SBC runs `fixTopology = TRUE`. Tracked separately.

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
