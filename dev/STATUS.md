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

**SBC on sampled_k + T-OVL (Stage 2) — SBC PASS (fixed-topo); T-OVL EXPOSES a marginal_k FREEZE under free topology (Hamilton 2026-06-02, f18a091).**
Driver `dev/red-team/heavy-tests/marginal-k/T-SBC-sampled-geometric.R`
(`likelihoodMode = "sampled_k"`; forward truncates at K_MAX_PRIOR = 30 by
rejection-redraw, inference pins kprimeTruncK = 30; k′_pooled excluded as a Talts
boundary artefact, `project_sbc_kprime_structural`). Jobs (login-node pre-build
done; stale-`.so` guard active):
- `17333100_[0-39]` sampled SBC **batch 1** (seedBase 20260528) → `sbc-results-sampled-b1/`
- `17333101` batch-1 aggregator (afterany)
- `17333102_[0-39]` sampled SBC **batch 2** (seedBase 20260901) → `sbc-results-sampled-b2/`
- `17333103` batch-2 aggregator (afterany)
- `17333104` T-OVL posterior-overlap grid (§7.2) → `marginal-k/T-OVL-verdict.txt`

**Pass bar = the PRE-REGISTERED 2-batch criterion, NOT the per-run all-3-AD>0.4
strict gate** (which rejects a perfect sampler ~78% = 0.6³, so it would chase
tl/rls noise). Apply `marginal-k/pool-sampled-batches.R` once both aggregators
finish: (a) pooled (N=400) AD>0.05 on each of tl/rls/p AND no batch <0.01;
(b) low-p (p_true<0.10) frac p-rank≤2 ≤0.15 & mean normRank ∈ [0.30,0.70];
(c) tl/rls freeze extreme(0|L) ≤~0.03. Given the RB proof ("Watertight") + the
both-variant prior bit-check, an SBC/T-OVL **miss ⇒ mixing, not target** — do not
re-litigate the prior.

**RESULT 2026-06-02.** SBC: pre-registered 2-batch verdict **PASS** (N=400,
200/200 good each) — pooled AD tl=0.37 rls=0.21 p=0.76 (all >0.05, no batch
<0.01); low-p (p_true<0.10, n=39) frac p-rank≤2=0.103, mean normRank=0.366 (the
truncation corner is calibrated under sampled_k); freeze ≤0.030. The per-run
strict gates FAIL on lone tl/rls/p dips (the 0.6³ noise) — pooled is clean.
**⇒ Stage 2 sampled_k *calibration* VALIDATED — but only for the FIXED-topology
regime the SBC exercises (`fixTopology = TRUE`, driver line 262).**

T-OVL (FREE topology): verdict FAIL (42/48 param-cells KS≤0.01), and the cause is
a **real bug, not test noise** (*correcting an earlier note here that had it
backwards*). **marginal_k whole-chain FREEZES in 7/16 cells** — `sd = 0` on tl,
rls AND p, stuck at init — while **sampled_k freezes in 0/16**. In the frozen
marginal chains EVERY Metropolis + slice move has **0% acceptance** (only the
Gibbs topology moves "accept"); non-frozen cells accept ~0.6–0.9. This is the
documented marginal-k `state→logLik`-inconsistency freeze (resume doc: "stuck at
init 0.5"; "sampled_k NEVER freezes"), firing under FREE topology — which the
`fixTopology` SBC structurally cannot reach. **REPRODUCED + LOCALIZED 2026-06-02**
(`dev/red-team/numerical/marginal-k-freeze-repro.R`; full write-up in
`dev/red-team/MARGINAL-K-FREEZE-003-diagnosis.md`): under marginal_k, **every
topology-changing move writes an incorrect `state->logLik`; `mh_logit_p` is the
only clean move.** A move-type coherence-gap sweep (validated by init + mh_p
controls reading ≈0) separates **two bugs** via warm-vs-cold: **(A)** the
**default-ON Gibbs moves** `gibbs_spr`/`gibbs_subtree_swap` write a **fixed-kPrime**
`state->logLik` (grouping at the single `state->kPrime`, `mcmc.cpp:1254-1282/1410`;
no marginal-over-k′ sum), cache left invalid (warm recomputes cold). The inflation
**= the omitted geometric weight `−n_char·log p`**: the gibbs gap tracks `−16·log p`
across p with **corr 1.0** (+10.6 nat at init p=0.5), nailing the mechanism.
**(B)** the MH/weighted topology moves (`nni`/`spr`/`pspr`/`weighted_*`) write a
**wrong accepted-move likelihood** (+1.3–4.8 nat; warm==state≠cold) = the
cache-option-a warm≠cold bug, now shown **broader than NNI**; root unpinned
(cache-invalidation gap vs proposal/commit edge bookkeeping — fix-stage). Causal
close: one `gibbs_spr` drops `mh_p` acceptance **95.3% → 0.0%**. In-situ r02
fingerprint: only `gibbs_spr` (0.67) + `gibbs_subtree_swap` (0.56) accept; all 12
other moves at 0.0000. The **7/16 incidence is p-gated bistability** (gap is
p-gated, ~uniform ~10.5 at the shared init p=0.5 — NOT a cross-cell magnitude
threshold): a self-healing move (mh_p-on-accept / slice, `mcmc.cpp:3441`) fires
only when the gap is small (high p; at p=0.9 mh_p heals & accepts 0.95) and never
at low/mid p (≤0.5 → 0%), so a chain locks iff p sits in the large-gap region when
a Gibbs move corrupts. Precise escape dynamics a follow-up.
**RESOLVED 2026-06-02 (commits `7d34892` fix + `24ad4f2` tbr guard).**
- **Bug B fixed.** `compute_per_kprime_log_lik` / `cpp_log_likelihood_marginal`
  now read topology from the passed `parent`/`child` (signature + 7 internal
  reads rerouted), so a topology-changing MH proposal is evaluated on the
  PROPOSED tree, not `state`'s old one. `gibbs_kprime_sweep_impl` (case 25)
  passes `state`'s own tree -> bit-identical (sampled_k path unchanged).
- **Bug A + cache-coherence audit (6 moves gated; 4 since RE-ENABLED — see Phase 2).**
  Six moves write an incoherent
  committed `state->logLik` under marginal_k and are disabled in `.BuildMoves`: the
  gibbs pair (`gibbs_spr`/`gibbs_subtree_swap`, fixed-k') plus four multi-eval moves
  (`weighted_spr`/`weighted_subtree_swap`/`weighted_branch_scale`/`block_gibbs_branch`)
  that read a STALE marginal `charLLCache` — `compute_full_loglik_at` does not force
  it cold between topology/branch evals and the `useCache` gate has no fingerprint
  (valid only across a pure-`p` change). `weighted_branch_scale`+`block_gibbs_branch`
  were UNGATED before this audit (default-off, so latent — a user enabling either
  would have silently corrupted the chain). k'-sampling moves (kPrime/gibbs_kPrime/
  block_kPrime) are likewise excluded under marginal_k (k' is integrated out).
  nni/spr/pspr/tbr still search topology.
- **Opt-in guard.** Compile `-DMKPRIME_CHECK_MARGINAL_COHERENCE` to warn when a
  committed marginal `logLik` != a cold recompute; the always-on CI guard is
  `tests/testthat/test-marginal-k-free-topology.R`.
- **Per-move coherence VERIFIED (exhaustive sweep).** Coherence gap-sweep over
  EVERY move type: nni/spr/pspr/tbr/mh_p -> **0.000**; the six gated moves show the
  corruption (weighted_branch_scale -1.26, block_gibbs_branch -1.60, gibbs +10.5,
  weighted_spr/subtree -0.9/-2.3; gibbs_kprime_sweep +6.85 if force-fired). Every
  move that CAN appear under marginal_k reads 0.000. End-to-end r02/r03 **unfrozen**
  (sd > 0; tbr accepts 0.43/0.32). Full testthat **FAIL=0 PASS=5858**.
- **RB-equivalence VALIDATED (deterministic, free topology).** The decisive check is
  the Rao-Blackwell likelihood identity `lse_k'[joint_S] == joint_M`, now generalised
  across topology space: `tests/testthat/test-marginal-k-rb-free-topology.R` holds it
  to 1e-7 on 8 random topologies x p x rate-heterogeneity (64 assertions PASS). This
  establishes marginal_k == sampled_k TARGET-equivalence on arbitrary trees
  deterministically (no MCMC noise). A signal-bearing MCMC posterior-overlap run is
  now **confirmatory, not load-bearing** (the legacy T-OVL's signal-free `sample(0:1)`
  data + over-powered per-cell KS are why it false-FAILed); local-smoke + Hamilton-prep
  item. **marginal_k free topology: freeze fixed, schedule coherent, RB-equivalence
  validated at the likelihood level.**
- **Phase 2 — 4 weighted/block moves RE-ENABLED (2026-06-02 eve).** Candidate-weight
  classification (advisor-driven): `weighted_branch_scale`(12)/`weighted_spr`(13)/
  `weighted_subtree_swap`(14)/`block_gibbs_branch`(15) all select candidates via the
  MARGINAL evaluator `compute_full_loglik_at` + a full MH accept (prior+Hastings) →
  pure cache-coherence (Bug B); the gibbs pair (10/11) select via the fixed-kPrime
  partial-CL path (wrong target) → stay deferred. Fix = a **scratch-eval** flag
  `fillCharLLCache` (default true) on `cpp_log_likelihood_marginal`/
  `compute_full_loglik_at`: when false the per-(char,k') `charLLCache` is neither READ
  (forced cold) nor WRITTEN. Threaded `false` into all 8 intra-move evals; tier-1-only
  suffices (tier-2 `perKpCl` cache is dormant — grep-proven never read). Case-12 commits
  via the generic MH path (coherent ready cache on accept); 13/14/15 self-accept leaving
  it cold; the next mh_logit_p reads a coherent cache either way (tested as warm==cold).
  `.BuildMoves` drops `&& !marginalK` for the four (all default-OFF/opt-in → production
  schedule unchanged); cli_inform now lists only the 2 deferred gibbs moves. VERIFIED:
  gap-sweep 12/13/14/15 baseline_gap & warm_gap → **0.000** (were -1.26/-0.90/-2.26/-1.60);
  `test-marginal-k-free-topology.R` Test 3 (300 accept/reject fires/move: committed==cold
  & warm==cold) + Test 4 (mh_p-after-weighted committed==cold) PASS. **DEFERRED:**
  marginal-aware gibbs candidate eval (10/11). PROPOSAL-SELECTION confirmed
  deterministically: Test 5 (`eval_preorder_paths_cpp`) asserts the candidate-eval
  (`preorder_into`) and commit (`preorder_weighted_impl`) paths give bit-identical marginal
  LL (1e-9) over 90 topologies + arbitrary edge orders -> selection LL == landed-tree LL.
  Full suite FAIL=0. **Overnight overlap (confirmatory):** structurally UNDERPOWERED on
  marginal_k's skewed/slow params -- marginal-p-ESS fell to 45-72 at 200k (no Gibbs-p;
  mh_logit_p mixes p poorly), p d_sd swung 1.25-8.62 on ESS noise; one moves-on rate_log_sd
  FAIL (d_sd=6.86 @ ESS~300) is the same moderate-ESS d_sd-test artifact (settled by Test 5),
  not a re-enabled-move bug; where ESS is adequate all params agree. Load-bearing evidence =
  RB identity + deterministic checks; overlap is supporting. Follow-ups (non-blocking): raise
  the harness sd-test ESS floor; a better marginal-p move.

**Deploy note (Hamilton).** GitHub auth is dead on the cluster (SSH publickey
denied; HTTPS creds empty) and `/nobackup` had purged the stale `mkp-source`
worktree — so deploy is **auth-free**: `git archive <sha> | scp | extract` into
`${SRC}`, then a login-node `pkgload::load_all` to build the `.so`. NOT `git pull`.

**Tracked cache bug — NOW IN SCOPE (pre-existing) = freeze Bug B.**
`test-marginal-k-cache-option-a.R` "NNI (case 5) refreshes the cache" FAILS
(warm − cold = 0.164 nats). Proven **not** a Stage 1b regression: the failure is
bit-for-bit IDENTICAL at K = 30 and K = 200, so the K-wiring did not cause it.
This is the partial-CL warm≠cold bug (a topology move does not refresh the
per-(node,k′) cache). ⚠ Previously filed "SBC-irrelevant because SBC runs
`fixTopology = TRUE`" — **that dismissal is now disproven** (FREEZE-003 repro,
above). It is **Bug B** of the freeze, and the FREEZE-003 sweep shows it is
**broader than NNI** (also `spr`/`pspr`/`weighted_*`). NB the *dominant in-situ*
freeze driver is **Bug A** (the default-ON Gibbs moves' fixed-kPrime
`state->logLik`, +10.6 nat) — frozen chains accept only Gibbs moves — so a
cache-option-a fix alone will **not** unfreeze marginal_k; Bug A must be fixed
too. **No longer Tier-2; together these block the marginal_k production default
and are the gating bugs for v1 marginal-k.**

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
- PR-C: T-OVL + sampled_k SBC validation — SUBMITTED to Hamilton 2026-06-02
  (jobs `17333100`–`17333104`; see the Stage-2 block above). Collect with
  `marginal-k/pool-sampled-batches.R` (pre-registered 2-batch criterion) and read
  `marginal-k/T-OVL-verdict.txt`; flip this arm to the v1 default if both pass.
- Stage 2 residual: the R-fallback `.DoMove` `gibbs_p` (`R/RunMkPrime.R`) still
  does an UNtruncated `Beta(a+nTrans, b+sumU)` draw — wrong for the truncated
  geometric, but currently UNREACHABLE (the `.BuildMoves` scheduler emits
  `mh_logit_p` for the geometric, and the C++ engine case 9 is guarded). Mirror
  the C++ non-conjugacy guard in the R-fallback if `gibbs_p` is ever
  re-scheduled, so a misrouted call can't silently sample the wrong p.
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
