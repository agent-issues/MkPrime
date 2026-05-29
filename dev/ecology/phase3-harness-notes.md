# Phase 3 validation-harness notes (EBE)

Phase 3 = end-to-end validation harnesses + adversarial review. Decisive tests
are spec'd in `dev/ecology/ebe-spec.md §8` and the plan's "Verification" block.
Harnesses live in `inst/ecology/simulations/`. Run via `devtools::load_all`
in a subprocess (installed pkg is stale). **All harnesses are downstream of
Phase 2b + GATE 2 — the model must be the wired EBE before they mean anything.**

## sim-ebe-null.R — RAN 2026-05-29, null degradation CONFIRMED
Spec tests 1 + 6 (graceful degradation + null-on-structure). z=0 everywhere,
`ecoEquilMode=TRUE`; runs aware-EBE AND blind MkNT on the same data.
**Result (16 tips, 60 neo + 60 trans, 100k iter each):** DECISIVE checks PASS —
pi0=0.751 (right at prior mean 0.75 → no spurious sparsity collapse / no invented
ecology); aware≈blind: logLik −550.86 vs −550.57, tree_length [0.97,2.06] vs
[0.74,2.13], tree-recovery CID 0.515 vs 0.521 (all overlap → ecology layer does
NOT distort the phylogeny = the "honest regulariser" property). Validates the
Phase 2b neo-z fix at full-MCMC level (trans chars not inflating pi0).
**phi caveat (calibrated + paper-relevant):** under the null phi is WEAKLY
IDENTIFIED (pi0~0.75 ⇒ few slab cells), so it reverts toward its LogNormal(0,1.5)
prior — mean 3.72, **CI [0.08, 40.5] INCLUDES 1**. NOT a manufactured effect (wide
CI, includes 1); the point estimate alone looks like an effect, so **phi must be
read WITH its CI**, and the recovery sim must show phi becomes TIGHTLY identified
(CI excludes 1) only under real dispersed signal — that contrast is the null↔effect
discriminator (spec §2 "per-character interpretation needs ecology dispersion").
Harness phi criterion recalibrated: "phi 95% CI includes 1" (was a wrong
`|log phi|<0.20` "phi→1"). tree/logLik/pi0 thresholds still provisional.

## sim-ebe-recovery.R — DESIGN ONLY (do NOT fabricate the signed-z check blind)
Model on `sim2-recovery.R` BUT three EBE-critical changes:
1. **Ecology-driven chars must be NEOMORPHIC**, not transformational. sim2 put
   the z≠0 signal on transformational chars (right for the OLD rate model, WRONG
   for EBE — EBE tilts only neomorphic equilibria). This is the single biggest
   adaptation.
2. **Dispersed ecology** (spec §2: φ is identified by the ensemble of
   ecology-correlated neomorphic chars across independent origins). Random tip
   ecology on a coalescent tree; require ≥4 origins per state (count edge-ecology
   transitions via `.AssignEdgeEcology`; resample if too clumped). Matches the
   rodent structure measured in `dev/eco-origins-check.R`.
3. **Signed z**: include both z=1 (toward present in eco 1) and z=2 (toward
   absent) neomorphic chars, plus z=0 neutral neomorphic AND z=0 transformational
   (the trans block confirms the neomorphic-only-z fix: they must NOT be flagged).

Known simulating params: `phi` (e.g. 4), `rateLoss` (e.g. 1.5 — asymmetric so
pi1_base≠0.5, exercising the reflect-about-base subtlety, spec R5), `baseRate`
(rate_neo). nTip ≈ 24 (dispersion + signal). nIter ≥ 100k (aware has more params;
may need longer — calibrate).

### Decisive checks
- **φ recovered** (confident: read `phi` col): canonicalise to φ≥1 (spec §2:
  flip z labels + theta when reporting), posterior covers true φ, |log φ| away
  from 0.
- **φ↔rate SEPARATION** (spec §2 / T-rec, confident: read `rate_neo`/`rate_loss`
  cols if present): each covers its true value; posterior cor(phi, rate_neo) and
  cor(phi, rate_loss) NOT near ±1 (no trade-off ridge). This is the headline
  "no global rate–TL ridge" claim — must be shown empirically.
- **pi0** below prior mean (reflects the fraction of neo chars with z≠0).
- **Tree accuracy ≥ blind** (run blind too; dispersed/random ecology → aware ≈
  blind on tree is acceptable; the win is φ + direction, not topology here).
- **SIGNED direction recovered** — ⚠ **OPEN API QUESTION, resolve post-GATE-2
  before writing this check.** Need the per-character POSTERIOR over z (P(z=1) vs
  P(z=2) vs P(z=0) per neomorphic char), then confirm z=1-simulated chars favour
  "toward present" and z=2 favour "toward absent". `s$zMatrix` (RunMkPrime.R:1584)
  is the per-state (final/current) z, NOT a posterior. Determine whether the chain
  streams z states (aggregate to marginals) or whether a z-trace/z-marginal output
  must be enabled. Do not guess the interface — confirm it on a live run, since
  signed recovery is the EBE headline ("directional, sign-identifiable") and a
  wrong accessor would silently mis-score it.

## Adversarial review (Phase 3, concurrent — AFTER the diagnostician's current
## root run finishes; do not double-book the compile lock)
- **numerical-auditor**: F81 conditioning, log-sum-exp at extreme tilts
  (π1→0/1 when φ large), const-site underflow (spec R6), and the math-prover's
  tagged caveat (4) zero-length-edge robustness under τ,σ→0. → `dev/red-team/numerical/`.
- **math-prover** already cleared the rooted moves; the const-site / non-stationary
  likelihood formula correctness is covered by GATE 1's T2-with-coding=1 oracle.
