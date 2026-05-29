# EBE design-lock spec (Phase 0 contract)

This is the shared contract for the Phase-1 build agents. All math, the single shared
P-helper, the entry-point list, and the test contracts are fixed here. Do not deviate
without updating this file.

## 1. State / data conventions (existing, unchanged)
- Neomorphic binary character: **state 0 = absent, state 1 = present**.
- Base rates (from `rate_loss`, the loss:gain ratio): `r01_base = 2/(1+rate_loss)`
  (gain 0→1), `r10_base = 2·rate_loss/(1+rate_loss)` (loss 1→0). Total `λ_base = 2`.
- Base equilibrium: `π1_base = 1/(1+rate_loss)` (P present), `π0_base = rate_loss/(1+rate_loss)`.
  This is also the existing neomorphic root distribution (`mcmc_ecology.cpp:1250-1251`) — **reused unchanged**.
- `z_{c,e} ∈ {0,1,2}` per (neomorphic character, non-reference ecology); reference
  ecology has no column (treated as z=0). `φ` = shared magnitude (global mode: `phi[0]`;
  per_ecology: `phi[s]`). `pi0`/`theta` priors on `z` — unchanged.

## 2. The equilibrium tilt (the ONLY new math)
For neomorphic character `c` on an edge of ecology `s`, define the **tilted equilibrium**:

    s_z = +1 if z=1 (toward present), −1 if z=2 (toward absent), 0 if z=0 or s=refEcology
    logit(π1) = logit(π1_base) + s_z · log φ_s
    π1_{z,s} = expit( logit(π1_base) + s_z·log φ_s ),   π0_{z,s} = 1 − π1_{z,s}

Then build the 2×2 transition matrix with the **decay rate FIXED at the base
λ = 2** (i.e. tilt the *attractor*, not the *tempo*):

    ex   = exp(−2 · tEff)            // tEff = edge_length · rate_neo · ACRV_rate  (as today, via neoEl)
    P00  = π0 + π1·ex      P01 = π1·(1−ex)
    P10  = π0·(1−ex)       P11 = π1 + π0·ex

This is the EXISTING MkN P formula (`pruning_mkn_acrv_flat_ecology_raw` ~L672-681) with
`(pi0_, pi1_) = (π0_{z,s}, π1_{z,s})` and `lam = 2` hard-fixed (NOT recomputed from
rate-scaled rates). Two consequences, both required:
- **Exact graceful degradation:** at z=0, `π1 = π1_base` ⇒ `r01=2π1_base=r01_base`,
  `r10=r10_base`, `lam=2` ⇒ byte-identical to baseline MkN. (Test 1.)
- **No global rate–TL ridge** (NOT full rate-orthogonality): the tilt changes only
  `(π0,π1)`, `λ` stays 2, so a `φ`/`z` move applies no global multiplicative `t·φ`
  rescale. BUT the stationary substitution rate `2λπ0π1 = 4π0π1` falls with the tilt
  (1.0 at π1=0.5 → 0.36 at π1=0.9), so a near-constant character trades off *locally*
  against low `rate_loss`/`rate_neo`. Pooled `φ` is identified by the ensemble of
  ecology-correlated characters; per-character interpretation needs ecology dispersion.
  **The recovery test MUST confirm `φ` and `rate_loss`/`rate_neo` are recovered
  separately, not traded off.**

Mixture over edge ecology (unchanged structure): `P_mix(e) = Σ_s wEdge(e,s)·P_{z(c,s),s}`.
Root: dot product `Σ_s rf[s]·clRoot[s]` with `rf = (π0_base, π1_base)` — unchanged.

`φ` symmetry note: the full benign 2-fold map is `(φ→1/φ, flip all z=1↔2, theta→1−theta)`
— it preserves BOTH likelihood and prior (the z-prior is asymmetric unless theta=0.5).
Equivalent modes for the tree; the conditional Gibbs cannot make the global flip by
chance, so no trap/bias (unlike the rate model's *inequivalent* modes). v1 keeps `φ` free
under the existing `LogNormal(0,sigmaPhi)`; **canonicalise to φ≥1 in reporting, flipping
z labels and theta accordingly**. (Folding the prior to φ≥1 is a v2 option.)

## 3. The single shared helper (R1 mitigation — ONE source of truth)
Add to `src/mcmc_ecology.cpp`, called by BOTH the main kernel and const-site:

    static inline void ebe_mkn_P(
        int z, int ecoState, int refEcology,
        double pi1_base, const double* phi, int magnitudeMode,
        double tEff, double* P /*len 4: P00,P01,P10,P11*/);

It computes `π1` per §2 (using `phi[0]` if magnitudeMode==0 else `phi[ecoState]`;
`s_z=0` when `ecoState==refEcology`), `π0=1−π1`, `ex=exp(−2·tEff)`, fills `P`.
`pi1_base` is passed in (= `1/(1+rate_loss)`), computed once per call. NEVER inline this
math anywhere else.

## 4. `gammaE` → identity
`gammaE` is a rate-mean normaliser with no meaning under equilibrium semantics. Force it
to 1.0 at both sites: `compute_gamma_e_ecology` (return all-1) and the inline Gibbs block
(`gibbs_z_sweep_impl` ~L4686-4694). Keep the plumbing/signatures; neutralise the value.

## 5. Entry points that MUST route through `ebe_mkn_P`
**Build order (walking skeleton):** land `ebe_mkn_P` + the main kernel (entry 1) and pass
T1 + T2 on a tiny tree FIRST; only then propagate to entries 2–8. Catches a tilt-math
error once rather than replicated across 8 sites.
1. `pruning_mkn_acrv_flat_ecology_raw` (~L610) — main neomorphic kernel.
2. `const_site_prob_mkn_eco_single_raw` (~L1019) — two passes (all-0, all-1), same helper, same rf.
3. `per_char_log_lik_ecology` (~L1930, type==0) — Gibbs path.
4. `cpp_partition_log_likelihood_ecology_raw` (~L1181, type==0) + Rcpp wrapper (~L1406).
5. `compute_eco_work_item_ll` (~L1482) — OpenMP work-item core.
6. `cpp_log_likelihood_ecology` (~L1672) — orchestrator (build work list; gammaE=1).
7. `CppLogLikelihoodEcologyCached` (~L2389) → `populate_eco_cache_full` + total.
8. `CppPartialEvalEcologyNNI` (~L2437) → `eco_cache_partial_eval_nni` + restore. **R3 audit.**
JC ecology pruner (~L306): transformational/known carry NO ecology effect — confirm it
applies no tilt in EBE (z ignored). **R8.**

## 6. Move / cache obligations
- `scale_phi` (case 34, mcmc.cpp ~L5291): the accept path MUST `invalidate_all_cls()`
  — `φ` now changes the kernel directly. **R2 audit.**
- `scale_pi0`/`scale_theta` (35/37): now affect only the prior (fed `gammaE`, now 1);
  keep their `invalidate_all_cls()` as harmless insurance in v1.
- `gibbs_z_sweep` (36): unchanged; still `invalidate_all_cls()` after the sweep.
- Priors (`LogPrior` / `cpp_log_prior`): **no formula change in v1**; keep R↔C++ in sync.

## 7. Rooting  (REVISED 2026-05-29 — infer the root by default)
- Root *distribution*: unchanged (base equilibrium). No new machinery.
- Root *placement*: with z≠0 the likelihood depends on it (T3). EBE is non-stationary,
  so the **rooted** tree is the inferential object.
- **DECISION (Martin's steer): the root is INFERRED by default, not fixed.** Rationale:
  (a) the tree-move Hastings ratios are *structural* (branch-length Jacobians + candidate
  counts; e.g. SPR `log(lRegraft/lMerge)`, TBR adds `log(lSubEdge/lMergeSub)`) —
  **likelihood-independent**, so MH already proposes over rooted trees and applies the EBE
  likelihood ratio on top; (b) the **prior is root-invariant** — branch term is the
  constant `lfactorial(nEdge-1)` (flat Dirichlet) and tree-length is Gamma on the preserved
  total; *proven empirically* in `dev/ecology/prior-reroot-check.R` (maxdiff=0 over
  redistribute + reroot; sensitive control). So nothing but the *data* (directional signal)
  drives the inferred root. The "root floats" behaviour (`dev/ecology/rootprobe.R`: 5
  distinct root bipartitions) is the **desired inference**, not a bug.
- **OPTION — fix the root** via a supplied outgroup taxon `g` (for known outgroups /
  tighter mixing when signal is weak): re-root start tree on `g`; constrain moves to keep
  the `{g | ingroup}` split — exclude root-adjacent NNI (`u==root`), root-children subtree
  swaps, and regraft onto `g`'s pendant edge. Local per-move guards, ~irreducible over
  ingroup topologies by the n-taxon↔(n−1)-rooted bijection. **Not the default.**
- **VERIFICATION (gates "infer-root is trustworthy"):** each topology move active in
  ecology mode must have a *rooted* Hastings ratio that is a correct proposal-density ratio
  over rooted trees (untested under the old root-invariant model). Active eco moves: NNI,
  SPR, TBR, pSPR (gibbs_spr/gibbs_subtree_swap self-disable at `mcmc.cpp:1024,1791`).
  - **Analytic leg — DONE (2026-05-29, math-prover → `dev/red-team/proofs/ebe-rooted-moves.md`).**
    All four default-active moves PROVEN to return `log[q(R'→R)/q(R→R')]` over *rooted*
    trees incl. branch Jacobians. TBR (prime suspect) confirmed: `nSubEdge=2m−2` and `nCand`
    are reroot-invariant (m preserved; 804/804 + 510/510 numeric), two disjoint Jacobian
    factors sum correctly; the "candidate-count ratio = 0" comment holds *as a rooted
    statement*. pSPR correct (residual tree identical f/r, Fitch root-invariant). Load-bearing
    **Assumption 3** — `preorder_weighted_impl` is a root-*preserving* relabelling, identity
    Jacobian — explicitly checked (root node + degree-2 invariant across thousands of moves).
    **No move needs fixing or disabling. Infer-root valid as-is on the proposal side.**
    `weighted_spr`/`weighted_subtree_swap` are OFF by default (no eco self-disable but zero
    weight); their guided-move `Z=Z'` cancellation is UNPROVEN → must verify IF ever enabled
    in eco mode.
  - **Empirical leg — DONE (mcmc-diagnostician → `dev/red-team/heavy-tests/`;
    autocorrelation-aware grader: TVD vs null band at ESS, every PASS requires its positive
    control to fire). Reconciled per-move:**
    - **NNI — CONFIRMED:** proof + empirical end-to-end posterior match, powered control.
    - **SPR — CONFIRMED:** proof + empirical end-to-end + exact localiser (count symmetry
      433,744 fwd/rev pairs / 0 asym; Jacobian log(lRegraft/lMerge) to 1.8e-15 over 100k moves).
    - **TBR — CONFIRMED (strong); end-to-end `--full` confirmatory, pending:** proof (count +
      Jacobian + reachability; 804/804 + 510/510 numeric) + empirical Jacobian incl. the
      Phase-B σ re-root (the crux) to 8.88e-16 over 13,877 fired moves. Count symmetry rests on
      the analytic proof (diagnostician declined to re-measure — a recount probe risks a false
      FAIL). `ebe-rooted-posterior-hamilton.sh` (~6-10h) gives the belt-and-suspenders
      end-to-end verdict; CONFIRMATORY, not load-bearing.
    - **pSPR — analytic-only:** proof (Fitch root-invariant, residual tree identical f/r,
      candidate symmetry from SPR); empirically UNVERIFIED. Default-on (`pSpr=TRUE`).
    - Methodology note: a first raw-χ² grader FALSELY flagged correct SPR (MCMC autocorrelation
      inflates χ²); rebuilt ESS-aware. Reassuring re: grader soundness.
  - **GATE STATUS (advisor 2026-05-29): infer-root trustworthy for Phase 2b — which is
    MOVE-INDEPENDENT (prior/Gibbs only), so it cannot depend on the TBR/pSPR residuals at all
    (decoupled, not merely "closed enough").** NNI/SPR confirmed both legs; TBR count (proof) +
    Jacobian (proof + machine-precision empirical) correct. What the `--full` run ACTUALLY
    gates is the **production infer-root claim (rodent, Phase 4)** — NOT Phase 2b: TBR's
    *factors* are verified but its *integrated end-to-end stationary* behaviour is not
    (factor-correctness ≠ assembled-move-correctness), and **pSPR (default-on, data-dependent
    density, ZERO empirical coverage) is the real residual** — not TBR. PLAN: (1) proceed Phase
    2b now; (2) fold pSPR into the SAME `--full` run (a wiring change to the existing
    end-to-end harness, NOT a new one — SendMessage the diagnostician) so one HPC run closes
    TBR end-to-end + pSPR together; (3) **do NOT claim/ship infer-root rodent results until
    `--full` closes.**
- **R3/R4 (still relevant):** `eco_cache_partial_eval_nni` root-invariance — 1A verified
  partial==full at z≠0 (max|diff| 0–7e-15); SPR/TBR use full recompute in eco mode.
- Caveat: weak directional signal → diffuse root posterior, slow mixing over rootings →
  that is exactly when the fix-root option earns its keep.

## 8. Test contracts (the oracle and the kernel agent build to these)
- **T1 graceful (exact):** all z=0 ⇒ `ebe` logLik == baseline non-ecology MkN logLik,
  tol 1e-10, for codingType 0 and 1.
- **T2 correctness:** mixed z, φ≠1 ⇒ `ebe` C++ == independent R oracle, tol 1e-10,
  ±ACRV, **for codingType 0 AND 1** — coding=1 is the REAL ascertainment-correctness
  gate (the oracle independently computes the `−log(1−pConst)` correction).
- **T3 non-stationarity:** re-root same topology/branch-lengths ⇒ logLik CHANGES when
  z≠0, INVARIANT (tol 1e-10) when z=0.
- **T4 const-site PLUMBING only:** Σ over all 2^nTip tip patterns of `P(pattern)` == 1 on
  a ≤6-tip tree. This checks only normalisation (holds for ANY stochastic P); it is NOT a
  tilt-correctness test — that is T2-with-coding=1.
- **T-rec separation:** in the recovery sim, simulate with known `φ` AND known
  `rate_loss`/`rate_neo`; confirm all are recovered separately (no φ↔rate trade-off).
- **T5 reconciliation:** partial-eval logLik == full recompute after an accepted `phi`
  move (R2) and after an NNI (R3).
- **T6 prior-invariance:** under a flat likelihood, `scale_phi/pi0/theta/gibbs_z` leave
  the posterior == prior (no dropped normaliser).

## 9. Oracle (Phase 1C) contract
Independent R function `ebe_loglik_R(tree, tipStates, rate_loss, phi, z, wEdge or eco,
refEcology, rate_neo, rateLogSd/nCat, coding)` that: reconstructs ecology marginals →
`wEdge` (may reuse `.EcologyEdgeWeights`), builds per-edge mixed P from §2 in plain R,
prunes (postorder), applies root `(π0_base,π1_base)`, applies the const-site correction
by the same two-pass logic, averages over ACRV. Derived from THIS spec, not from the C++.
**Enforced: agent 1C must NOT read `src/mcmc_ecology.cpp` or any C++ kernel source** —
otherwise T2 becomes tautological. The oracle must implement the const-site
`−log(1−pConst)` correction itself so T2-with-coding=1 is a genuine check.
