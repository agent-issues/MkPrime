# EBE rooted-topology posterior correctness (end-to-end)

## What this tests
The Ecology-Biased-Equilibrium (EBE) neomorphic likelihood is **non-stationary**:
ecology tilts the F81 equilibrium, so the likelihood depends on ROOT PLACEMENT
(spec `dev/ecology/ebe-spec.md` §7, contract T3). MkPrime infers the root by
default: the topology moves produce rooted trees and the Metropolis-Hastings
step applies the EBE likelihood ratio. The correctness claim under challenge is
that the topology moves active in ecology mode (**NNI, SPR, TBR**, and pSPR)
sample the **correct posterior over ROOTED trees** — i.e. each move's *rooted*
Hastings ratio is a valid proposal-density ratio. These ratios were built and
tested only under the OLD root-invariant model, so they are unverified for
rooted inference. A wrong Hastings ratio yields a biased root posterior that
*still runs and still looks plausible* — exactly what a unit test cannot catch.

## Theoretical basis
The prior over a rooted-tree state `(topology, TL, relBr)` factorises
(`R/MkPrimeModel.R:613-622`) as

    pi(R)  propto  L_EBE(R) . Gamma(TL; shape, rate) . Dirichlet(relBr; 1,...,1)

with the relative-branch term a pure constant `lfactorial(nEdge-1)` and **no
topology term**. The prior over rooted topology is therefore FLAT (proven
empirically root-invariant, maxdiff 0, in `dev/ecology/prior-reroot-check.R`).
Hence the EXACT marginal posterior over a rooted topology `t` is the
prior-predictive evidence per topology — a KNOWN TRUTH computable by forward
likelihood evaluation alone:

    P(t)  propto  INT INT  L_EBE(t, TL, relBr) . Gamma(TL) . Dir(relBr)  d relBr d TL
          propto  E_{relBr ~ Dir(1..1)} [ L_EBE(t, TL0, relBr) ]   (at fixed TL0)

We CONDITION on a fixed total tree length `TL0`: the SPR and TBR Hastings ratios
are TL-invariant to machine epsilon (verified — TL cancels in
`log(lRegraft) - log(lMerge)`; same for TBR's `log(lSubEdge) - log(lMergeSub)`),
so fixing TL loses ZERO power on the branch-length Jacobian under test while
making the target a tractable 1-D-per-edge integral. The chain is a hand-rolled
single-site MH using the **actual exported proposal kernels** (`spr_proposal`,
`tbr_proposal`, `nni_proposal` — byte-identical `_impl` code to ecology-mode
MCMC, `src/mcmc.cpp` cases 6/17 fall through to `*_proposal_impl`) and the exact
acceptance form `logAlpha = dLogPi + logHastings` (`src/mcmc.cpp:5876`).

**The crux (a naive test misses it):** a topology-only chain FREEZES the
branch-length degrees of freedom the target integrates over — NNI permutes a
fixed length multiset (`logHastings = 0`, never redraws lengths), and SPR/TBR
preserve total TL. So every chain ALSO runs a `dirichlet_simplex` relBr move so
the branch-length VALUES vary, exercising exactly the merge/split that the
SPR/TBR Jacobians act on. TBR additionally needs an internal pruned subtree for
its distinctive Phase-B sigma re-root to fire (~3 %/proposal at 5 tips,
~6 % at 6 tips), so TBR uses 6 tips.

Reference: standard reversible-jump / fixed-dimension MH detailed balance
(Green 1995; Yang & Rannala 1997 for tree-space MCMC). The TVD-vs-ESS scaling
of an unbiased sampler is the Glivenko–Cantelli / multinomial-CLT rate.

## Pass criterion
The grading is **autocorrelation-aware** because a naive chi-squared on MCMC
samples is INVALID (this was empirically confirmed to false-FAIL a correct SPR
chain — see "What a failure would mean"). Two criteria, both required, every
replicate:

1. **TVD null-band (primary).** Under the correct target, the expected
   total-variation fluctuation of a multinomial at effective size `ESS` is
   `nullMean = (1/2) sum_i sqrt(2 P_i(1-P_i)/(pi*ESS))` with
   `nullSd = (1/2) sqrt(sum_i P_i(1-P_i)/ESS)`. PASS iff
   `TVD <= nullMean + 4*nullSd`. An UNBIASED move's TVD decays as `1/sqrt(ESS)`
   and sits in the band; a BIASED move's TVD PLATEAUS at a nonzero bias floor.
   The `4*nullSd` (≈ p > 6e-5 one-sided Gaussian) is a deliberately wide band:
   the discriminator is the *plateau vs decay* behaviour, and the band only has
   to separate the bias floor from sampling noise (the injected-bias control
   sits at ~3–4x its band).

2. **ESS-corrected G-test (secondary).** The likelihood-ratio statistic `G`
   over cells with expected count >= 5 (Cochran), rescaled by `ESS/nStored` so
   the chi-squared reference is valid. PASS iff `p > 1e-3`. The 1e-3 threshold
   (vs 0.05) is conservative for the handful of moves x replicates compared
   (a loose Bonferroni over ~8 grades keeps family-wise error < 1 %).

`ESS` is estimated per high-probability topology (mass >= 1 %) via
`coda::effectiveSize` on the 0/1 indicator series, taking the conservative
MINIMUM (the worst-mixing coordinate bounds the multinomial ESS).

**Positive controls (calibration + power statement).** An SPR chain and a TBR
chain whose Jacobian is forcibly zeroed (`logHastings := 0`, dropping the
branch-length term) MUST be flagged biased. Each control POWERS its move
family: the SPR control powers the NNI and SPR verdicts (NNI's `logHastings` is
identically 0, so a dropped-Jacobian failure is exactly what the SPR control
injects); the TBR control powers the TBR verdict (and exercises the Phase-B
sigma term). **A move is reported CORRECT only if its TVD/G verdict passes AND
its powering control fired.** A move that runs clean but whose control did not
fire is reported `EXECUTED-UNDERPOWERED`, never CORRECT — PASS-without-power is
worse than INCONCLUSIVE. In quick mode the TBR control rarely fires (Phase-B is
~3–6 %/proposal on a short tiny-tree chain), so TBR is expected
`EXECUTED-UNDERPOWERED` there; its Jacobian (incl. Phase-B sigma) is separately
pinned to machine precision by `ebe-move-detailed-balance.R`, and `--full`
provides the powered end-to-end TBR verdict.

## How to run
- **Small-N (execution gate, ~25 s):** `Rscript ebe-rooted-posterior.R --quick`
  4-tip tree (15 rooted topos) for NNI/SPR/control on a 5-tip tree (105 topos)
  for TBR/control; 800 target draws; 1 chain x 1.5e4 iter. TBR runs here too so
  the gate covers EVERY move + both controls + the verdict logic (a reference
  error must surface here, not hours into SLURM) — but TBR is
  **EXECUTED-UNDERPOWERED** in quick mode by design (see Pass criterion).
- **Full-scale:** `Rscript ebe-rooted-posterior.R --full`
  or `sbatch ebe-rooted-posterior-hamilton.sh`. 5-tip (NNI/SPR/control) and
  6-tip (TBR/control) trees; 4e5 target draws; 4 chains x 3e7 iter each, thinned
  to 30k stored topologies (feedback_no_oversample). At 6 tips Phase-B accrues
  ~1.8M events over the run, firing the TBR control. Resources: 1 node, 8 cores,
  ~6-10 h walltime. Resumable: each (move, replicate) chain is independent and
  seeded. Memory: < 2 GB (tiny trees).

## Output interpretation
`ebe-rooted-posterior-results/verdict.txt`: per-move line shows
`TVD=<obs>(thr <band>) ESS=<eff> corrP=<G p-value> acc=<rate>` and a verdict
(`SAMPLES CORRECT ROOTED POSTERIOR` / `BIASED`). The `POSITIVE CONTROL ...
detected bias: YES/NO` line is the power statement. `OVERALL:` is PASS / FAIL /
INCONCLUSIVE / PASS(partial). `summary.rds` holds every replicate's full
diagnostics (TVD, threshold, ESS, raw and corrected G, p, accept rate) and the
exact target `P` per topology.

To confirm a PASS is not a power artefact, check that `TVD` is well below `thr`
AND that for the control `TVD > thr`. A genuine FAIL shows `TVD` exceeding `thr`
with the deviation persisting (not shrinking) as `--full` lengthens the chain.

## What a failure would mean
If a move's `TVD` exceeds its band AND the ESS-corrected `corrP < 1e-3` AND the
deviation does NOT shrink as the chain lengthens, the most likely cause is a
**wrong rooted Hastings ratio for that move** (a dropped candidate-count term or
an incorrect branch-length Jacobian), biasing the inferred root posterior.
Localise WHICH factor with the companion `ebe-move-detailed-balance.R`.

Plausible alternatives to a true move bug, to rule out first:
- **Insufficient ESS / poor mixing** (the deviation shrinks with longer chains;
  TVD decays rather than plateaus) — NOT a move bug. This is the dominant
  false-positive mode: a naive raw-count chi-squared was empirically observed to
  drive a CORRECT SPR chain's p-value from 3.6e-2 to 1.2e-17 as nStored grew
  7.5k -> 300k, while its TVD DECAYED 0.023 -> 0.008 (the unbiased signature)
  and the companion detailed-balance localiser showed the SPR Jacobian exact to
  1e-15. The ESS-corrected grading exists specifically to defeat this.
- **Target MC error** (too few `targetDraws`) — inflates apparent deviation
  uniformly; mitigated by 4e5 draws in full mode and a common draw set across
  topologies.
- **Fixture degeneracy** (z all-zero, or an ecology vector giving a near-flat
  posterior) — guarded: `make_z` forces z=1 and z=2; the target entropy ratio
  is printed.
