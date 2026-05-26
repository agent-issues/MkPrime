# Lane D1 — Simulation-Based Calibration (SBC) for MkNT and Mk'

## What this tests

The MkPrime sampler's joint prior × likelihood × MCMC pipeline. SBC
challenges the correctness claim: *if the inference algorithm is correct
and the prior used at inference time matches the prior used to generate
the data, then for any one-dimensional posterior summary the rank of the
true value within `L` posterior samples is Uniform{0, …, L}* (Cook,
Gelman & Rubin 2006; Talts, Betancourt, Simpson, Vehtari & Gelman 2018).

The harness is **parameterised per prior** so it can demonstrate
calibration under matched priors *and* explicit miscalibration when the
sampler's prior is missing its truncation normaliser — i.e. it functions
as a regression demonstrator for **EG-001** (already filed; see
`dev/red-team/proofs/kprime-priors.md` §4.3).

## Theoretical basis

For a generative model with parameter `θ`, observation `y` with
likelihood `p(y | θ)`, and prior `π(θ)`:

1. Draw `θ̃ ~ π`.
2. Draw `ỹ ~ p(· | θ̃)`.
3. Run inference to obtain `L` posterior samples `θ_1, …, θ_L ~ p(· | ỹ)`.
4. Compute rank `r = #{i : θ_i < θ̃}`.

Then `r ~ Uniform{0, 1, …, L}` whenever the inference is *exact* (Talts
et al. 2018, theorem in §3). Repeating across `N_sim` simulations and
testing the pooled ranks against uniformity gives a frequentist
goodness-of-fit test for the entire inference pipeline.

**Why EG-001 shows up here.** The L6 proof
(`dev/red-team/proofs/kprime-priors.md`) shows the empirical-geometric
prior in `R/MkPrimeModel.R::.LogPriorEmpiricalGeometric` is missing the
per-character truncation normaliser
`Z_i(p) = Σ_{k ≥ kObs_i} P(k | p)`.  Under the correct (truncated) prior
the MCMC posterior on `p` would absorb the data correctly; under the
implemented (untruncated) prior the posterior is biased by
`∑_i log Z_i(p)`, which depends on `p`. The forward simulator in this
harness draws `(p, k'_i, character)` from the truncation-aware joint that
the L6 proof identifies as the *intended* model — so the rank histogram
of `p` and pooled `k'` will skew non-uniformly under the buggy sampler.
The geometric, beta-geometric (`Z_i ≡ 1`) and logseries-with-fixed-`c`
(latent LS-001) arms should remain calibrated.

## Pass criterion

Anderson–Darling test of uniformity on each parameter's pooled rank
distribution (Talts et al. 2018 §4.1 recommend a rank-uniformity test;
AD is more sensitive than KS in the tails, which is where SBC failures
typically manifest).

For each arm, let `K_arm` be the number of monitored parameters.
**Arm verdict** `PASS` iff every monitored parameter has

`p_AD > 0.001 / K_arm`     (Bonferroni-corrected across parameters)

The base threshold `α = 0.001` is conservative — six arms × ~5
parameters gives ~30 tests, and at `α = 0.001` the family-wise false
positive rate is still < 5%. Per-arm correction prevents one noisy
parameter (e.g. topology summary) from dragging the whole arm down.

`INSUFFICIENT` (fewer than 4 ranks pooled) does not count toward FAIL.

## Arms

Six arms, with predicted PASS/FAIL under the **current** (unpatched)
code:

| arm | model | prior | K | expected | reason |
|-----|-------|-------|---|----------|--------|
| `MkNT_geometric` | MkNT | geometric | 2 | PASS | k' pinned by knownStates; geometric is no-op |
| `Mkp_geometric` | Mk' | geometric | 4 | PASS | `Z_i ≡ 1`; hierarchical geometric |
| `Mkp_beta_geometric` | Mk' | beta_geometric | 3 | PASS | `Z_i ≡ 1`; BetaGeo identity |
| `Mkp_empirical_geometric` | Mk' | empirical_geometric | 4 | **FAIL** (EG-001) | missing `log Z_i(p)`; bias scales with kObs spread |
| `Mkp_logseries` | Mk' | logseries | 3 | PASS (latent LS-001) | `c` is fixed → bias is per-char additive constant, cancels in MH ratios |
| `MkNT_logseries` | MkNT | logseries | 2 | PASS | k pinned; prior reduces to support indicator |

`MkNT × empirical_geometric` is *not* run — MkNT pins `k = kObs` via
`knownStates`, so the empirical-geometric prior has no degrees of
freedom under MkNT. Skipping this combination is explicit in
`sbc.R::ALL_ARMS`.

Per-arm `K`:

- `tree_length`, `rate_log_sd` always monitored.
- `kPrime_pooled` (per-character `k'_i` ranks pooled across characters
  under the hierarchical-prior exchangeability assumption — see L6 §2.1)
  monitored for Mk' arms.
- `p` monitored for arms whose prior actually samples `p` (geometric,
  empirical_geometric); not monitored for beta_geometric (no scalar `p`
  to sample — alpha, beta are sampled instead and not exposed in the
  samples matrix) or logseries (`c` is fixed).

**Topology rank is intentionally excluded** from the rank histogram
test. Talts et al.'s rank-uniformity is for continuous parameters;
bucketing RF distance would smuggle in an arbitrary discretisation
threshold. A topology recovery diagnostic can be added in a follow-up
test under a separate pass criterion (e.g. expected RF ≤ some bound) but
is out of scope here — the SBC test is for the *parameters*.

## How to run

### Quick mode (laptop, ≤ 60 s, execution smoke test only)

```bash
Rscript dev/red-team/heavy-tests/sbc.R --quick
```

Confirms the harness executes for every arm. Per-arm verdicts at
quick scale are reported as `EXEC_OK` not `PASS`/`FAIL` because
`N_sim = 5` and `L = 34` (after warmup) give the AD test insufficient
power. The reported AD p-values at quick scale should be inspected for
*qualitative* skew direction (e.g. `Mkp_empirical_geometric` should
already show smaller `k'` p-values than `Mkp_geometric` — visible at
quick scale even if not statistically significant).

**Measured execution time** (this worktree, win-x64, `--quick`, all 6
arms):

```
real  0m6.146s
user  0m0.015s
sys   0m0.030s
```

Output at `dev/red-team/heavy-tests/sbc-results/verdict.txt`. All six
arms `EXEC_OK (good=5)`. The full quick run completes in well under the
60 s budget.

### Full scale (Hamilton, 8 h walltime per arm)

```bash
# On Hamilton, after copying sbc.R to /nobackup/$USER/mkp-study/
sbatch dev/red-team/heavy-tests/sbc-hamilton.sh
```

Resources per arm (one SLURM array task each):

| resource | value | justification |
|----------|-------|---------------|
| CPUs | 1 | RunMkPrime is single-threaded at fixTopology=TRUE |
| memory | 8 GB | conservative; MCMC state ~10s of MB |
| time | 8 h | 200 sims × ~45 s/sim = 2.5 h nominal; 8 h covers warmup variance + Hamilton's slow-start tax |
| scratch | 4 GB | per-task checkpoints |
| partition | shared | standard analysis class |

Full scale: `N_sim = 200`, `N_iter = 6000`, `thin = 60` so
`L = 100` posterior samples per simulation. This matches Talts et al.'s
recommended L ≈ 100.

To run a single arm at full scale locally (e.g. to debug):

```bash
Rscript dev/red-team/heavy-tests/sbc.R --full --arm Mkp_empirical_geometric
```

## Output interpretation

Output root: `dev/red-team/heavy-tests/sbc-results/`

```
sbc-results/
├── verdict.txt                          # top-level summary
├── MkNT_geometric/
│   ├── verdict.txt                      # per-arm AD p-values + decision
│   ├── summary.rds                      # full sim outputs
│   └── rank-matrix.csv                  # pooled ranks, one col per param
├── Mkp_geometric/...
├── Mkp_beta_geometric/...
├── Mkp_empirical_geometric/...
├── Mkp_logseries/...
└── MkNT_logseries/...
```

To plot rank histograms (per-arm post-hoc):

```r
r <- read.csv("dev/red-team/heavy-tests/sbc-results/Mkp_empirical_geometric/rank-matrix.csv")
op <- par(mfrow = c(2, 2))
for (nm in names(r)) hist(r[[nm]], breaks = 20, main = nm)
par(op)
```

A skewed histogram (e.g. systematic mass at low or high ranks) is the
qualitative signature of miscalibration; the AD p-value formalises it.

## What a failure would mean

**Predicted FAIL on `Mkp_empirical_geometric` `k'` and/or `p` ranks:**
EG-001 is confirmed by SBC — the missing `Z_i(p)` truncation normaliser
is producing a real, measurable bias in the joint posterior. The
expected fix is to add the `log Z_i(p)` term to
`R/MkPrimeModel.R::.LogPriorEmpiricalGeometric` and
`src/mcmc.cpp::cpp_log_prior`. Re-running the harness with the patch
applied should flip the arm to PASS.

**Unexpected FAIL on `Mkp_geometric` or `Mkp_beta_geometric`:** these
arms have `Z_i ≡ 1` (L6 §4.1 and §4.2) — so a FAIL here points to a bug
*elsewhere* in the inference pipeline. Most likely causes, in order:
(a) the MCMC mixing is poor enough that fixed-topology runs aren't
converging in `N_iter = 6000`; (b) the relabelling correction
(`src/corrections.cpp::mk_prime_relabel_log`) is double-counted or
missing somewhere; (c) a likelihood arithmetic error orthogonal to
either of the above.

**Unexpected PASS on `Mkp_empirical_geometric`:** the harness's forward
simulator is conditioning on `kObs` somewhere it shouldn't (defeating
the EG-001 demonstration). The most likely place is the kPrime draw —
verify `.drawKPrime()` does *not* truncate at `kObs` before simulation.
A reconcile call to the math-prover lane is warranted in that case.

**FAIL on `tree_length` or `rate_log_sd` in any arm:** these are
shared infrastructure parameters. A FAIL here is likely a mixing
problem rather than a prior bug, since neither parameter is implicated
in EG-001 or LS-001. Inspect ESS for that parameter in
`summary.rds$sims[[i]]` to confirm.

## Constraints honoured

- `feedback_no_oversample`: thinning chosen so `L ≈ 100` post-warmup,
  well below the 30k cap.
- `feedback_start_tree`: starting tree uses `TreeSearch::AdditionTree`
  (NJ fallback only on `TreeSearch` unavailable).
- `feedback_prior_modelling`: the forward simulator draws `k'_i` from
  the untruncated prior; `kObs_i` is observed *after* simulation. The
  prior at inference time depends on `kObs_i` only as a support bound
  (per L6 §2 assumption 2) — *not* as a data-conditioned truncation.
- `feedback_model_scope`: MkNT arms pin `k = kObs` via
  `knownStates = setNames(kObs, names)`. Mk' arms leave `k'` free.
- `feedback_no_prs`: no PR will be opened; this is a heavy-test
  artefact only.
- Role-spec constraints: standalone Rscript under
  `dev/red-team/heavy-tests/`, no edits to `R/`, `src/`, or
  `tests/testthat/`. Quick-mode wall recorded above. Pass criterion
  stated upfront (top of this file). Outputs under
  `dev/red-team/heavy-tests/sbc-results/`. `verdict.txt` written per
  arm + top-level.

## Trivial-fix policy

No patch file is produced. EG-001 is a code fix outside the red-team
scope (`R/` is read-only per the role spec). The fix is described above
("Predicted FAIL" section) and is already documented in
`dev/red-team/proofs/kprime-priors.md` §4.3. The orchestrator may route
the actual fix to a separate lane that has `R/` write access.

## References

- Cook, S. R., Gelman, A., & Rubin, D. B. (2006). Validation of
  software for Bayesian models using posterior quantiles. *Journal of
  Computational and Graphical Statistics*, 15(3), 675–692.
- Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A.
  (2018). Validating Bayesian inference algorithms with
  simulation-based calibration. *arXiv:1804.06788*.
- L4 proof: `dev/red-team/proofs/relabelling-correction.md` (Mk'
  relabelling correction lives in the likelihood, not the prior).
- L5 proof: `dev/red-team/proofs/ascertainment.md` (closed-form site
  probabilities under JC; commutation of relabelling × ascertainment).
- L6 proof: `dev/red-team/proofs/kprime-priors.md` (per-prior
  normalisers; EG-001 and LS-001 derivations).
