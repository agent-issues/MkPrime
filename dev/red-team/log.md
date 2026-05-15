# Red-team rotation log

One entry per round. Bottom of file holds `last_focus:` pointer (single
integer = the area number just reviewed). The skill picks
`(last_focus mod N) + 1` next.

Format per round:

```
## Round R — area #N <name> (YYYY-MM-DD)
- Files reviewed: ...
- Scenarios traced: ...
- Findings filed: <count> (refs to dev/red-team/findings.md rows)
- Trivial fixes applied: ...
- Notes for next reviewer of this area: ...
```

---

## Round 1 — area #1 ecology likelihood + prior math (2026-05-14)

Files reviewed:
- `src/mcmc_ecology.cpp` (1100+ lines)
- `src/mcmc.cpp` ecology block in `cpp_log_prior` (lines 388-433)
- `R/MkPrimeModel.R` `LogPrior` (lines 504-686)
- `inst/simulations/ecology/sim3-simulate.R`
- `vignettes/ecology-details.qmd`

Scenarios traced:
- `gamma_e` formula across simulator / R / C++.
- Per-edge, per-character mixture in pruning loops.
- Reference-ecology skip indexing.
- Boundary cases: theta ∈ {0, 1}, pi0 ∈ {0, 1}.
- wEdge marginal root-edge handling.
- Ecology rate hard-coded at 1.0 in marginals computation.

Findings filed: 5
- B1: wEdge root-edge departs from `1/K` spec, uses `0.5 * (margRoot + margChild)`.
- B2: Simulator scalar theta vs model per-ecology theta_e (kEco > 2 only).
- B3: z column-count contract — simulator has kEco cols, model has kEco - 1.
- B4: R `LogPrior` rejects theta ∈ {0, 1} strict AND has 0 * log(0) = NaN.
- B5: Ecology rate hard-coded at 1.0 in `compute_ecology_node_marginals`.

Trivial fixes applied during round: none (round 1 was research-only).
Fixes since (separate commits): B4 fixed in `b254b98` (theta boundary +
zero-count guard); B1, B2, B3, B5 still open.

Notes for next reviewer of this area:
- After fixing B1 + B5, re-trace `wEdge` on Sim 3's deep stem branch
  vs the simulator's hard-edge labels. Suspect ~50-100 nat impact.
- B2/B3 dormant at kEco = 2; revisit before any kEco > 2 simulation.

## Round 2 — area #2 MCMC proposals + MH ratios (2026-05-14)

Files reviewed:
- `src/mcmc.cpp` move dispatch (cases 30, 34-37), Bactrian kernel (line 78),
  scale_phi/scale_pi0/scale_theta accept blocks
- `R/RunMkPrime.R` `.kMoveTypes` (line 3127), `.BuildMoves` (line 2820+)
- gibbs_z `_impl` at `src/mcmc.cpp:4011-4116`
- log-sum-exp at `src/mcmc.cpp:4023-4045`

Scenarios traced:
- Jacobian correctness for log-phi and logit-pi0/theta moves.
- Boundary reflection (logit always stays in (0,1)).
- gibbs_z behaviour at theta = 1 (zero probability for z = 2).
- scale_phi sampling refE in per_ecology mode.
- `gibbsZEvery` plumbing.
- Move ID renumbering 30 → 34-37 after main merge.

Findings filed: 4
- M1: `scale_phi` perturbs `phi[refEcology]` in per_ecology mode (breaks
  identifiability). Defer until per_ecology mode is actually used.
- M2: `cpp_log_prior` sums LogNormal phi prior over all kEco entries
  (compounds M1). Same defer.
- M3: `gibbsZEvery` plumbed from R → `McmcData` but never read in
  `mcmc.cpp`. Documented behaviour doesn't exist. Move fires purely on
  weighted draw.
- S1 (suspect): `scalarTypes` omits `scale_theta` from 2% weight floor.
  Fixed in `b254b98`.

Trivial fixes applied during round: none.

Notes for next reviewer of this area:
- Adaptive scale tuning during warmup not exhaustively audited — see
  area #7.
- Warmup-end anomaly (chain visible sampling at iter ~7020 with
  `maxWarmup=1000`) explained by Tuning-phase budget; not a moves bug.
  Worth area #7 follow-up.
- M1/M2 dormant until per_ecology mode is used.

## Round 3 — area #3 state synchronisation (2026-05-14)

Files reviewed:
- All `state->logLik =` and `state->logPrior =` assignments in `src/mcmc.cpp`
- `fill_partition_cache` (line 491)
- `gibbs_kprime_sweep_impl`, `block_kprime_shift_impl`, gibbs_z accept,
  scale_phi/pi0/theta accept, slice_scalar_impl
- gibbs_spr / gibbs_subtree_swap and their _het / _full variants
- Pre-proposal drift diagnostic (line 4114)

Scenarios traced:
- Which moves write non-eco partition sums into `state->logLik` in eco mode.
- `wEdgeDirty` lifecycle.
- Per-character Gibbs z cell-by-cell drift.
- Init mismatch between R `.MkpEcologyLogLikelihood` and C++
  `cpp_log_likelihood_ecology`.

Findings filed: 4 (all critical, all fixed)
- S1: `gibbs_kprime_sweep_impl` (case 25, 36% of moves) used
  `cpp_partition_log_likelihood` (non-eco), corrupting `state->logLik`.
  **FIXED `36e9b53`.**
- S2: `block_kprime_shift_impl` (case 26) same bug in `else` branch.
  **FIXED `36e9b53`.**
- S3: Pre-proposal drift diagnostic used `cpp_log_likelihood` in eco
  mode. **FIXED `36e9b53`.**
- S4: `gibbs_spr_impl` + `gibbs_subtree_swap_impl` (cases 10/11, ~34%
  of moves) used partial-CL machinery with no eco awareness. **GATED
  `ad25928`** — these moves now no-op in eco mode.
- S5 (added round 4-ish): `eval_slice_target` + `slice_scalar_impl`
  accept block used non-eco `cpp_log_likelihood`. **FIXED `ad25928`**.

Trivial fixes applied during round: none (large enough to warrant a dedicated commit).

Notes for next reviewer of this area:
- Periodic from-scratch resync now lives in `run_mcmc_batch_cpp` main
  loop (every 20 iter). Drift > 0.5 nats logged to stderr.
- Residual stale-log at sample 1 was still ~113 nats after S1-S5; the
  cause was likely **init mismatch** rather than mid-chain accumulator
  drift — but that diagnostic lives in area #4 (state init).
- If new accept blocks are added in future, force them through
  `compute_full_loglik_at` which already branches on `data.ecologyAware`.

## Round 4 — area #4 R-side state init (2026-05-15)

Files reviewed:
- `R/RunMkPrime.R`: `RunMkPrime` (lines 84-340), `.InitRun` (535-601),
  `.RunMkPrimeSingleRun` initialisation block (618-647), `.InitState`
  (2734-2841), `.InitMcmcChain` (3232-3253), `.SaveCheckpoint`
  (2184-2241), `ResumeMkPrime` (2288-2436), `.PerturbStart` (2516-2532)
- `R/likelihood.R`: `.MkpEcologyLogLikelihood` (140-256), `.MkpLogLikelihood`
- `R/MkPrimeModel.R`: `LogPrior` (510-708) — focused on ecology block
- `R/MkPrimeData.R`: `.ExtractEcology` (335-386)
- Cross-checked: `src/mcmc.cpp` `init_mcmc_state` (444-487),
  `fill_partition_cache` (491-517), `cpp_log_prior` ecology block
  (388-433), `get_mcmc_state` (598-640)

Scenarios traced:
- Does R `.MkpEcologyLogLikelihood` use `state$tree_length` /
  `state$rel_br_lengths` or `tree$edge.length`?
- Path of `initOverrides$tree`, `tree_length`, `rel_br_lengths` through
  `.InitState` → likelihood call.
- What chain-state fields `.InitRun` flattens vs what `.InitMcmcChain`
  expects.
- Checkpoint round-trip for `kprime_alpha`, `kprime_beta`, `beta_scale`.
- C++ vs R boundary semantics for theta ∈ {0, 1}, pi0 ∈ {0, 1}, phi=0.
- Preorder invariant under `initOverrides$tree` and `.PerturbStart`.
- Whether `fill_partition_cache` provides a safety-net resync in eco mode.

Findings filed: 9 (R4-1 .. R4-9) — 3 HIGH, 4 MED, 2 LOW.

Most exposed file: `R/RunMkPrime.R` (`.InitState`, `.InitRun`,
`.SaveCheckpoint`). R-side init now suspected to be the dominant source
of the ~113-nat gap at sample 1 in dev pilot logs.

Trivial fixes applied during round: none. Findings warrant a dedicated
init-correctness commit.

Notes for next reviewer of this area:
- R4-1 + R4-7 together explain the smoking gun: R log_lik computed from
  original `tree$edge.length`, C++ rebuilds from overridden
  `treeLength * relBrLengths`, and `fill_partition_cache` eco branch
  returns early so the mismatch survives until the first 20-iter resync.
- R4-3 contradicts the L-4 "fix mirrors C++" claim; revisit L-4 status.
- R4-4/R4-5 dormant until `beta_geometric` or qHeterogeneity is used in
  production but easy to slip through.
- Suggest single-place fix: in `.InitState`, after applying overrides,
  rebuild `tree$edge.length <- state$tree_length * state$rel_br_lengths`
  AND `tree <- TreeTools::Preorder(tree)` (re-applying relBr after
  reorder). This addresses R4-1 and R4-2 in one pass.

## Round 4 (C++) — area #4 C++ state init (2026-05-15)

Files reviewed:
- `src/mcmc.cpp`: `McmcState` struct (180-232), `init_mcmc_state` (444-487),
  `fill_partition_cache` (491-517), `allocate_cl_workspace` (528-594),
  `get_mcmc_state` (601-640), `compute_full_loglik_at` (759-783),
  periodic eco-resync block in `run_mcmc_batch_cpp` (5311-5350)
- `src/mcmc_state.h`: full file — `EcologyState`, `EcologyInfo`, `McmcData`
  ecology fields, `ClWorkspace`
- `src/mcmc_ecology.cpp`: `cpp_log_likelihood_ecology` entry (~900-1073),
  `recompute_w_edge` (1230-1257), `compute_ecology_node_marginals` call
  sites at 933, 1244, `CppLogLikelihoodEcologyPerChar` (1264-1320)
- Cross-checked R: `.InitMcmcChain` (R/RunMkPrime.R:3232), `.InitRun`
  chain reconstruction (633-647)

Scenarios traced:
- Whether `state->logLik` at the end of `fill_partition_cache` is a
  fixed point of `cpp_log_likelihood_ecology` in eco mode (it is NOT;
  the function returns early).
- All read sites for `state->wEdge` and `state->wEdgeDirty` (zero reads).
- Earliest iter at which the periodic resync fires (iter 20 for
  `startIter = 1`).
- Whether `init_mcmc_state` validates eco state shape against
  `data->ecology.kEcology` (it cannot — no dataPtr).
- Whether `cpp_log_likelihood_ecology` reuses `state->wEdge` or allocates
  fresh (always fresh).
- Storage-mode / cloning hygiene for phi, theta, zMatrix in
  `init_mcmc_state`.

Findings filed: 6 (R4C-1 .. R4C-6) — 1 HIGH, 2 MED, 3 LOW.

Most exposed file: `src/mcmc.cpp` (`fill_partition_cache`). The eco
early-return on line 497 is the single highest-leverage C++ defect and
the direct counterpart to R-side R4-7: it materially differs from the
non-eco contract (line 516 unconditionally overwrites `state->logLik`)
and is the reason any R-side init bug survives invisibly until iter 20.

Trivial fixes applied during round: none. Findings warrant a dedicated
init-correctness commit (likely shared with the R-side round-4 fix).

Notes for next reviewer of this area:
- R4C-1 + R4-1/R4-7 are the same bug seen from two sides. The cleanest
  fix touches both: rebuild `tree$edge.length` in R AND have
  `fill_partition_cache` overwrite `state->logLik` with a fresh C++ eco
  recompute (or at least warn loudly on mismatch — silent overwrite
  hides upstream bugs).
- After R4C-1 lands, the 20-iter resync block (line 5320) should
  briefly continue logging — if it goes quiet for 100+ iter on a
  truth-init replay, the smoking gun is fully resolved.
- `state->wEdge` cache (R4C-3) is genuine perf headroom; revisit after
  correctness is solid.
- R4C-4 should be folded into the same fix commit since
  `fill_partition_cache` is being edited anyway.

---

last_focus: 4
