# Red-team findings — open + closed bugs from rotation rounds

Trivial fixes happen inline during a round and don't land here. This
file collects non-trivial issues (bug or perf) flagged by an Opus
reviewer that need a dedicated fix commit, design decision, or follow-up.

Severity scale: HIGH (correctness; affects results), MED (correctness;
narrow conditions or off-by-ε), LOW (perf, hygiene, latent), DOC (docs
mismatch with code, no immediate impact).

| ID  | Severity | Status   | Kind  | Title                                                                 | Reference                                                                                       |
|-----|----------|----------|-------|-----------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| L-1 | MED      | OPEN     | Bug   | wEdge root-edge departs from `1/K` spec                               | `src/mcmc_ecology.cpp:104-106, 234-241, 936-942` — uses `0.5*(margParent+margChild)`. Vignette says `1/K`. Suspected ~50-100 nat impact on Sim 3 deep stem branch. |
| L-2 | MED      | OPEN     | Bug   | Simulator scalar `theta` vs model per-ecology `theta_e`               | `inst/simulations/ecology/sim3-simulate.R:180-184` — harmless at kEco=2, biased at kEco>2.       |
| L-3 | LOW      | OPEN     | Hygiene| z column-count contract: simulator has kEco, model has kEco-1        | `inst/simulations/ecology/sim3-simulate.R:32, 211`; `src/mcmc_ecology.cpp:692, 762, 947, 1154`. Add assertion or reshape simulator to match model. |
| L-4 | MED      | CLOSED   | Bug   | R `LogPrior` rejected theta ∈ {0, 1} strict + 0×log(0) NaN trap        | Fixed in commit `b254b98`. R now uses closed [0, 1] with zero-count guard, mirroring C++.       |
| L-5 | MED      | OPEN     | Bug   | Ecology rate hard-coded at 1.0 in `compute_ecology_node_marginals`    | `src/mcmc_ecology.cpp:205, 933, 1244`. Smears root marginal under long stem branches. Design decision: expose as parameter or stay fixed? |
| M-1 | LOW      | OPEN     | Bug   | `scale_phi` perturbs `phi[refEcology]` in per_ecology mode             | `src/mcmc.cpp:4544-4554`. Defer until per_ecology mode is actually used.                          |
| M-2 | LOW      | OPEN     | Bug   | `cpp_log_prior` sums phi prior over all kEco entries (includes refE)  | `src/mcmc.cpp:410-412`. Same defer as M-1.                                                       |
| M-3 | DOC      | OPEN     | Bug   | `gibbsZEvery` field plumbed but never read by `mcmc.cpp`              | Docs claim "Gibbs every N gens" gating; reality is weighted-draw only. Either honour or document. |
| S-1 | HIGH     | CLOSED   | Bug   | `gibbs_kprime_sweep_impl` wrote non-eco logLik into `state->logLik`    | Fixed in commit `36e9b53`. Now branches on `data->ecologyAware` and uses `cpp_log_likelihood_ecology`. Fires at ~36% of moves so was the dominant smoking gun. |
| S-2 | HIGH     | CLOSED   | Bug   | `block_kprime_shift_impl` else branch same bug as S-1                  | Fixed in commit `36e9b53`. Identical pattern.                                                    |
| S-3 | MED      | CLOSED   | Bug   | Pre-proposal drift diagnostic used non-eco `cpp_log_likelihood`        | Fixed in commit `36e9b53`. Now uses `cpp_log_likelihood_ecology` when `ecologyAware`.            |
| S-4 | HIGH     | CLOSED   | Bug   | `gibbs_spr_impl` / `gibbs_subtree_swap_impl` used non-eco partial CL  | Gated in commit `ad25928` — both now no-op when `data->ecologyAware`. Topology mixing falls back on spr/nni/tbr/pspr. Long-term: build eco-aware streaming candidate evaluator. |
| S-5 | HIGH     | CLOSED   | Bug   | `eval_slice_target` + `slice_scalar_impl` accept used non-eco logLik   | Fixed in commit `ad25928`. Both now route to `cpp_log_likelihood_ecology` when `ecologyAware`.   |

## Closing protocol

When a bug is fixed:
1. Update `Status: CLOSED`.
2. Reference the fix commit in the row.
3. Leave the row in the table (so future reviewers see what was found and where).

When a bug is deferred:
1. Update `Status: DEFERRED` and add a short reason.
2. Reference any tracking issue or upstream work.

When a finding turns out to be a false alarm:
1. Update `Status: DISMISSED` and add a one-line reason.
