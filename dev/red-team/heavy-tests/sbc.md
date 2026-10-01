# Lane D1 — Simulation-Based Calibration (SBC) for MkNT and Mk'

## What this tests

The MkPrime sampler's joint prior × likelihood × MCMC pipeline. SBC
challenges the correctness claim: *if the inference algorithm is correct
and the prior used at inference time matches the prior used to generate
the data, then for any one-dimensional posterior summary the rank of the
true value within `L` posterior samples is Uniform{0, …, L}* (Cook,
Gelman & Rubin 2006; Talts, Betancourt, Simpson, Vehtari & Gelman 2018).

The harness is **parameterised per prior**. Every arm draws each
character's `k` from the arm's prior and simulates the character at that
`k`, so the Mk' arms test the `k'` prior they are named after.

## Forward model

Per simulation, matching the inference model term by term:

1. `p ~ Beta(1, 1)` (geometric and empirical-geometric arms; logseries fixes
   `c = 0.5`), `rate_log_sd ~ Gamma(1, 1)`, tree length
   `~ Gamma(2, 2 / expSteps)` with `expSteps = 50` spread over edges by a
   flat Dirichlet, on a random topology that inference then holds fixed.
2. `k_j` per character from the arm's prior, before any data exist:
   - Mk' geometric: `2 + Geometric(p)` truncated to `[2, 30]` by rejection,
     with `kprimeTruncK = 30` on the inference side (the prior is truncated and
     renormalised; a `pmin()` clamp would not match, cf. #57).
   - Mk' empirical-geometric: `N_obs + Geometric(p)`, `N_obs` from the packaged
     `empiricalNObs`, untruncated like the inference convolution.
   - Mk' logseries: `P(k) ∝ c^k / k` on `k ≥ 2`.
   - MkNT: `k` is pinned per character through `knownStates = kTrue`, so it is
     data rather than a parameter; the arm's prior merely supplies it.
   Mk' arms pin `priorVariant = "unconditional"` (Model A), the prior this
   forward draws from.
3. Each character is simulated under JC(`k_j`) at a rate multiplier from one
   of `nCat = 4` equiprobable discretised-lognormal categories
   (`DiscreteLognormalRates(rate_log_sd, 4)`), and `(category, data)` are
   redrawn together, with `k_j` held fixed, until the character is variable.
   That matches `coding = "variable"`: given `k`, inference divides the
   category-averaged likelihood by the category-averaged probability of
   variability (a ratio of sums over categories), and sums over `k'` outside
   that ratio.

The earlier design simulated every character at `k = 2` and dropped
invariant characters, so `kObs ≡ 2`, Models A and B coincided and `p` and
`k'` were not ranked; at `nCat = 1` `rate_log_sd` never reached the
likelihood. Those arms could not fail on the prior they were named after.

Seeds are `seedBase + ArmSeedOffset(arm) + i`, keyed on the arm's name, so
adding or removing an arm does not re-seed the others.

## Pass criterion

Anderson–Darling test of uniformity on each parameter's pooled rank
distribution (Talts et al. 2018 §4.1 recommend a rank-uniformity test;
AD is more sensitive than KS in the tails, which is where SBC failures
typically manifest).

For each arm, let `K_arm` be the number of monitored parameters.
**Arm verdict** `PASS` iff every monitored parameter has

`p_AD > 0.001 / K_arm`     (Bonferroni-corrected across parameters)

The base threshold `α = 0.001` is conservative — five arms × ~5
parameters gives ~25 tests, and at `α = 0.001` the family-wise false
positive rate is still < 5%. Per-arm correction prevents one noisy
parameter (e.g. topology summary) from dragging the whole arm down.

`INSUFFICIENT` (fewer than 4 ranks pooled) does not count toward FAIL.

## Arms

| arm | model | prior | monitored |
|-----|-------|-------|-----------|
| `MkNT_geometric` | MkNT | geometric | `tree_length`, `rate_log_sd` |
| `Mkp_geometric` | Mk' | geometric (truncated at 30) | + `p`, `kPrime_sum` |
| `Mkp_empirical_geometric` | Mk' | empirical_geometric | + `p`, `kPrime_sum` |
| `Mkp_logseries` | Mk' | logseries (`c` fixed) | + `kPrime_sum` |
| `MkNT_logseries` | MkNT | logseries | `tree_length`, `rate_log_sd` |

All five are expected to PASS under a correct sampler. There is no
beta-geometric arm: that prior is on `k' - kObs_i`, so it conditions on the
data and is not SBC-calibratable.

`kPrime_sum` is `Σ_j k'_j`, one scalar per simulation; as a function of the
parameters its rank is uniform too. Being discrete, its ties are broken
uniformly at random (as are all ranks), which keeps the rank
Uniform{0, …, L} under correct inference.

**Residual mismatch.** The `k'` sampler never proposes beyond
`kObs + 255` (`kMaxKprimeCand`), while the empirical-geometric and logseries
priors are untruncated. The forward cannot match a cap that depends on
`kObs`; the lost mass `P(k' > kObs + 255)` is negligible except for
empirical-geometric simulations with `p` below about 0.01.

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
power; quick-scale AD p-values mean nothing.
`test-heavy-test-gates.R` runs two quick arms and checks that characters
are simulated at their drawn `k` and that `p` and `k'` are ranked.

### Full scale (Hamilton, 8 h walltime per arm)

```bash
BUILD=$(sbatch --parsable dev/red-team/heavy-tests/submit-build.sh)
ARR=$(sbatch --parsable --dependency=afterok:$BUILD dev/red-team/heavy-tests/submit-sbc.sh)
```

then aggregate with `--run-id $ARR` as the footer of `submit-sbc.sh` shows.
Each task stamps `run_id:` into its `verdict-<arm>.txt`;
`aggregate-verdicts.R` ignores files from any other run, so an arm whose task
died reports INCOMPLETE rather than last run's verdict.

Resources per arm (one SLURM array task each):

| resource | value | justification |
|----------|-------|---------------|
| CPUs | 1 | RunMkPrime is single-threaded at fixTopology=TRUE |
| memory | 8 GB | conservative; MCMC state ~10s of MB |
| time | 8 h | set for the old `nCat = 1` design; `nCat = 4` costs more per sim, and the new per-sim wall has not been measured on Hamilton |
| scratch | 4 GB | per-task checkpoints |
| partition | shared | standard analysis class |

Full scale: `N_sim = 200`, `N_iter = 6000`, `thin = 60` so
`L = 100` posterior samples per simulation. This matches Talts et al.'s
recommended L ≈ 100. The 8-tip, `expSteps = 50` regime is near saturation
(`T-SBC-marginal-geometric.R` moved off it for that reason); poor mixing
there can pile ranks at the extremes even with a correct forward model.

To run a single arm at full scale locally (e.g. to debug):

```bash
Rscript dev/red-team/heavy-tests/sbc.R --full --arm Mkp_empirical_geometric
```

## Output interpretation

Output root: `dev/red-team/heavy-tests/sbc-results/`

```
sbc-results/
├── verdict.txt                          # all arms: one process, or aggregate-verdicts.R
├── verdict-<arm>.txt                    # one array task's arm, stamped with run_id
├── MkNT_geometric/
│   ├── verdict.txt                      # per-arm AD p-values + decision
│   ├── summary.rds                      # full sim outputs
│   └── rank-matrix.csv                  # pooled ranks, one col per param
├── Mkp_geometric/...
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

**FAIL on `p` or `kPrime_sum` in an Mk' arm:** the `k'` prior, its
hyperprior update, or the `k'` moves target the wrong distribution for that
prior. Check first that the forward in `.drawKPrime()` still matches the
inference prior (support, truncation, `priorVariant`); this harness family
has had four forward/inference mismatches.

**FAIL on `tree_length` or `rate_log_sd` in any arm:** shared
infrastructure parameters. Either a likelihood or ACRV error, a forward
mismatch in the rate or ascertainment model (step 3 above), or poor mixing.
Inspect ESS for that parameter in `summary.rds$sims[[i]]` before blaming
the sampler.

## Constraints honoured

- `feedback_no_oversample`: thinning chosen so `L ≈ 100` post-warmup,
  well below the 30k cap.
- The starting tree is the true topology with uniform edge lengths
  (SBC-HARNESS-005); inference holds the topology fixed.
- `feedback_prior_modelling`: the forward draws `k'_j` from the
  unconditional prior; `kObs_j` is observed *after* simulation and enters
  inference only as the likelihood's support bound.
- `feedback_model_scope`: MkNT arms pin `k` via
  `knownStates = setNames(kTrue, names)`. Mk' arms leave `k'` free.
- Standalone Rscript under `dev/red-team/heavy-tests/`, no edits to `R/`,
  `src/`, or `tests/testthat/`.

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
