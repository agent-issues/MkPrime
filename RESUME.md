# mkp — hand-off 2026-05-18 (late evening, post-collect)

## What this repo is

MkPrime — R package for Bayesian phylogenetic inference under the Mk′ model
for discrete morphological characters. The current focus is a 2026-05-12
prior-validation pilot comparing arms on 26 simulated trees × 10 reps to find
which k-specification best recovers the true tree (CID-to-truth). A parallel
real-data pilot (in neotrans; see below) tests whether the winning arm from
simulation transfers to six empirical matrices compared against well-corroborated
reference trees.

## Where we left off

This session was a hand-off-resume + collect cycle:

- **f6a91b7** — `dev/m9-pilot/process_pilot.R` extended to loop over all three
  models on 07203. CID-vs-asher results:

  | model    |   n | CID_mean | CID_sd |
  |----------|----:|---------:|-------:|
  | by_nt_9v |  52 |   0.5703 | 0.0157 |
  | t_9v     |  56 |   0.5737 | 0.0179 |
  | t_kv     | 242 |   0.5929 | 0.0226 |

  Real-data echo of the sim finding: fixed k=9 (both Bayesian by_nt_9v and t_9v)
  beats the t_kv baseline by ~0.02 CID against asher WCT.

- **mk_k15 / mk_k24 pilot reps completed** (jobs 17216080, 17216081 — both
  COMPLETED exit 0 on t01_r01). Timing:
  - mk_k9:  ~11.88M iter in 8h → ~118k thinned samples (thin=100)
  - mk_k15: ~10.14M iter in 8h → ~101k thinned samples
  - mk_k24:   ~235k iter in 8h → ~2,350 thinned samples (43× slower than k15)

- **mk_k15 full array submitted** as job **17218354** (`--array=0-259`, 8h each).
  Resume capability means we can extend walltime later if results warrant.

- **Preliminary 7-way CID (mk_k15 partial: 153 COMPLETED + 107 RUNNING summarised
  from partial logs).** Job 17222171 ran the summariser; results in
  `dev/pilots/2026-05-12-prior-validation/analysis/cid_seven_prior.R`:

  | Arm | n | mean CID |
  |---|---:|---:|
  | mk | 260 | 0.2754 |
  | mkp_eg | 260 | 0.2677 |
  | mk_kp1 | 260 | 0.2595 |
  | mkp_geo | 26 | 0.2562 |
  | mk_kp2 | 260 | 0.2535 |
  | mk_k9 | 260 | 0.2451 |
  | **mk_k15** | 260 | **0.2418** |

  Paired mk_k9 vs mk_k15 (260 common tasks): k15 wins 206/260, mean diff
  −0.00335, binomial p = 4.3e-22. **The flexibility ramp extends past k=9**,
  monotonically (mk_kp1 → mk_kp2 → mk_k9 → mk_k15) with diminishing
  returns (Δk_p2→k9 = 0.0084; Δk9→k15 = 0.0033). Numbers will tighten when the
  remaining 107 reps complete; ordering is very unlikely to flip.

## Pending jobs

### Hamilton HPC — mkp sim (8h walltime, resumable)

| Job ID | Arm | Status | Note |
|--------|-----|--------|------|
| 17218354 | mk_k15 (--array=0-259) | RUNNING | Full array; ~3.5 days end-to-end given shared partition |
| 17222516 | mk_k24 (--array=0-259) | QUEUED | Submitted 2026-05-19; 43× slower → ~2.3k samples/8h run; resumable to 3d |

mk_k24 queued now that k15 results confirm ramp extends past k=9 (k15 beats k9, p=4.3e-22). Both jobs resumable; extend walltime if needed after first pass.

### Hamilton HPC (long-form M9 real-data, 3-day walltime)

| Job ID | Matrix | Model | Status | Elapsed | Mem |
|--------|--------|-------|--------|---------|-----|
| 17217093 | syab07200 | by_nt_9v | RUNNING | ~6.5h | 32G |
| 17217094 | syab07200 | t_kv | RUNNING | ~6.5h | 32G |
| 17217095 | syab07200 | t_9v | RUNNING | ~6.5h | 32G |
| 17217096 | syab07202 | by_nt_9v | RUNNING | ~6.5h | 32G |
| 17217097 | syab07202 | t_kv | RUNNING | ~6.5h | 32G |
| 17217098 | syab07202 | t_9v | RUNNING | ~6.5h | 32G |
| 17217100 | syab07204 | t_kv | RUNNING | ~6.5h | 32G |
| 17217103 | syab07205 | t_kv | RUNNING | ~6.5h | 32G |
| 17217105 | syab07206 | by_nt_9v | RUNNING | ~6.5h | 32G |
| 17217106 | syab07206 | t_kv | RUNNING | ~6.5h | 32G |
| 17217107 | syab07206 | t_9v | RUNNING | ~6.5h | 32G |
| 17217660 | syab07204 | t_9v | RUNNING | ~4h | 64G |
| 17217661 | syab07205 | by_nt_9v | RUNNING | ~3.5h | 64G |
| 17217662 | syab07205 | t_9v | PENDING | — | 64G |

All complete in ~2.5d (counting from job start). On completion: `scp` the
`.trees` files from `/nobackup/pjjg18/m9-long/<matrix>/` to a per-matrix dir
under `dev/m9-pilot/`, then run an extended `process_pilot.R` looping over
matrix × model.

### Known problem — needs decision

**17217659 (syab07204 by_nt_9v at 64G) OOM'd again** after 1m43s. 64G was not
enough for this matrix+model combination. Options: resubmit at 128G, or skip
syab07204 by_nt_9v from the 6-matrix comparison and note as infeasible.

## Open items / next steps

1. **When mk_k15 array (17218354) completes**: pull summaries, extend
   `cid_six_prior.R` to a 7-way table. Then decide on k24 — submit at 8h to
   see if even ~2k samples suffices, or skip if k15 already settles the
   "does ramp continue past k=9?" question.
2. **Decide 17217659 retry** at 128G vs skip.
3. **Wait for long-form M9 to accumulate ≥50 trees/run** (earliest meaningful
   data ~09:00 BST 2026-05-19) then extend `process_pilot.R` to loop over
   matrix × model and produce the full CID table.
4. **Pull neotrans `by_nt_kv` baseline CID** for the 6 matrices, so we can do
   the (i) and (ii) comparisons from the original plan.
5. **Pre-pub `_9i` (informative coding) variants**: deferred, after pilots
   confirm M9 direction.

## Technical pointers

- **SSH**: PowerShell ssh to `hamilton8.dur.ac.uk` works; Bash-tool ssh does NOT
  (different config). Always route Hamilton commands through PowerShell.
- **PowerShell quoting**: use single-quoted strings for ssh commands containing
  `$` to avoid PowerShell variable interpolation. PowerShell 5.1 lacks `&&` /
  `||`; chain with `;` only when failure is acceptable.
- **Background SSH on Windows can hang** if the ssh session is initiated as a
  Bash background task; route through PowerShell instead. Even PowerShell ssh
  can be slow to first-handshake (~30s) when the cluster is under load.
- **Hamilton paths (mkp sim study)**:
  - Run scripts: `/nobackup/pjjg18/mkp-study/{run_one.R, summarize_streamed.R, *_array.slurm}`
  - Raw streamed logs: `/nobackup/pjjg18/mkp-study/results/t{NN}_r{MM}/{arm}_run_1.log`
  - Tree files: `/nobackup/pjjg18/mkp-study/results/t{NN}_r{MM}/{arm}_trees.nwk`
  - Summaries: `/nobackup/pjjg18/mkp-study/summary/{arm}_t{NN}_r{MM}.rds`
  - Data root: `/nobackup/pjjg18/mkprime-files/tree-inference/tree_{NN}/rep_{MM}/chr*.nex`
- **Hamilton paths (M9 real-data pilot)**:
  - Pilot (4h): `/nobackup/pjjg18/m9-pilot/syab07203/{model}_run_{1,2}.trees`
  - Long runs (3d): `/nobackup/pjjg18/m9-long/<matrix>/{model}_run_{1,2}.trees`
  - Logs: `/nobackup/pjjg18/m9-long/logs/m9L-<matrix>-<model>_<jobid>.{out,err}`
  - Source nex files: `/nobackup/pjjg18/07200_by_t_kv/`, `/nobackup/pjjg18/07202_ns_n_kv/`, etc.
- **Local tree files are gitignored**: `dev/m9-pilot/syab*/` is in `.gitignore`
  (intentional — big outputs); CID summaries are too. Record numerical results
  in commits / RESUME.md, not files.
- **RB .trees format**: NOT NEXUS. Tab-separated: `Iteration Posterior Likelihood Prior phylogeny`. Read with `read.delim(f)$phylogeny` then `ape::read.tree(text = s)` per row.
- **WCT pruning**: `KeepTip(wcTree, intersect(sampleTips, wcTree$tip.label))` before CID. Matrix→WCT mapping: `paste0("0720", 0:6)` → `wellCorroboratedTrees.nwk` lines 1–7. syab07203→asher (outgroup: Didelphis, Macropus).
- **OOM threshold (Mk′ M9 real-data)**: fnJC(9) trans-only model needs ≥64G for
  matrices with ≥214 trans.nex lines (syab07204/07205). 64G itself is **not
  always enough** — see 17217659 (syab07204 by_nt_9v at 64G OOM'd at 1m43s).
  by_nt models may need 128G+ for these larger matrices.
- **Per-arm MCMC speed (mkp sim)**: rough iter/sec on t01_r01:
  - mk_k9:  ~1,500 iter/sec
  - mk_k15: ~1,300 iter/sec
  - mk_k24:    ~30 iter/sec (43× slower — matrix exp scales like k³)
- **kObs alignment**: `chr*.nex` filenames are unpadded. `run_one.R` uses plain `sort()` → lex order (chr1, chr10, chr11, …). Analysis must match.
- **EG arm result (corrected, 2026-05-17)**: median u_post = 1.01 (binary chars, n=9917); EG posterior on k′ is calibrated near kObs+1. mk_kp1 collapses this to a point estimate, beating EG on CID.
- **No oversampling**: thin streamed MCMC to ~30k samples/run (thin_iters ≥ 16 for 8h job). See `feedback_no_oversample.md`.

## Things ruled out

- **u_post anchor hypothesis (first pass, pre-bug-fix)**: thought mk_kp1 ≈ EG
  posterior mean. Initial audit (with broken alignment) showed u_post ≈ 0.1,
  rejecting. Post-fix, hypothesis **revived** — EG sits near kObs+1.
- **Sampler bug in EG prior**: investigated negative u_post. Confirmed constraint
  correct (`R/MkPrimeModel.R:406`). The bug was a script-level alignment error.
- **Updating the simulation to match the new prior**: ruled out as circular.
- **Per-pattern vs per-character k′ confusion**: phyDat compresses to unique
  patterns; `mkd$kObs` and `state$kPrime` are per-character. MkPrime expands back.
- **kObs+1 sweet-spot hypothesis**: mk_k9 > mk_kp2 > mk_kp1 monotonically in
  6-way CID pilot. Confirmed extended in 7-way (2026-05-19, this session):
  mk_k15 > mk_k9 (paired p=4.3e-22). Ramp is monotonic up to k=15 with
  diminishing returns — open question whether k24 plateaus or continues.

## Suggested first action

```powershell
# Survey long-form M9 progress — how many trees per matrix×model so far?
ssh pjjg18@hamilton8.dur.ac.uk 'for d in /nobackup/pjjg18/m9-long/syab*/; do echo "=== $(basename $d) ==="; for f in $d/*_run_1.trees; do echo "  $(basename $f) $(wc -l < $f) lines"; done; done'
```

Then decide on mk_k15/mk_k24 full-array submission and 17217659 128G resubmit
(see "Open items" 1–2).
