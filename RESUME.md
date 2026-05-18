# mkp — hand-off 2026-05-18 (evening)

## What this repo is

MkPrime — R package for Bayesian phylogenetic inference under the Mk′ model
for discrete morphological characters. The current focus is a 2026-05-12
prior-validation pilot comparing arms on 26 simulated trees × 10 reps to find
which k-specification best recovers the true tree (CID-to-truth). A parallel
real-data pilot (in neotrans; see below) tests whether the winning arm from
simulation transfers to six empirical matrices compared against well-corroborated
reference trees.

## Where we left off

Six-way CID pilot **completed and analysed** this session:

- **8a1df6d** — 6-way CID script `dev/pilots/2026-05-12-prior-validation/analysis/cid_six_prior.R`
  + mkp_geo u_post audit (`u_post_geo.R`). **Headline result**: mk_k9 (fixed k=9
  for all variable chars) won outright — mean CID 0.231 vs mk_kp2 0.238,
  mk_kp1 0.245, mkp_eg 0.250, mkp_geo 0.257, mk 0.260. Sign-test mk_k9 vs
  mkp_eg: 23/26 wins, p ≈ 4×10⁻⁵.

- **8ba40a8** — Extended `data-raw/hamilton/run_one.R` and slurm scripts with
  `mk_k15` and `mk_k24` arms, to test whether flexibility ramp plateaus beyond k=9.
  Pilot array jobs 17216080 (k15) and 17216081 (k24) submitted; still running as
  of hand-off.

- **8985a9c** (this session) — `dev/m9-pilot/process_pilot.R` — validated
  end-to-end pipeline: RB .trees TSV → ape → KeepTip pruning to WCT taxa →
  ClusteringInfoDistance. `.gitignore` entry added for pulled tree files.

Real-data M9 pilot lives in **neotrans**, branch `m9-pilot` (see neotrans
RESUME.md), but analysis and collection script live here.

## Pending jobs

### Hamilton HPC

#### Sim k=15/k=24 extension (pilot arrays — 1 task each)
| Type | Job ID | Model | Status | ETA | On completion |
|------|--------|-------|--------|-----|---------------|
| HPC | 17216080_0 | mk_k15 | RUNNING | ~2h39m | Check `/nobackup/pjjg18/mkp-study/results/` for k15 RDS; if present, run `dev/m9-pilot/process_pilot.R`-style sanity check. If 1 rep fit in 8h, submit full `--array=0-259`. |
| HPC | 17216081_0 | mk_k24 | RUNNING | ~2h39m | Same. If OOM or too slow (>8h for 1 rep), k=24 is not feasible; note in `Things ruled out`. |

#### syab07203 M9 pilots (4h budget)
| Type | Job ID | Model | Status | ETA | On completion |
|------|--------|-------|--------|-----|---------------|
| HPC | 17216901 | by_nt_9v on 07203 | RUNNING | ~1h | Pull trees, run `dev/m9-pilot/process_pilot.R`, extend to include all three models for CID table. |
| HPC | 17216902 | t_9v on 07203 | RUNNING | ~1h | Same. |

#### Long-form M9 real-data runs (3-day walltime, 71h srMaxTime)
All produce output in `/nobackup/pjjg18/m9-long/<matrix>/`. Script to collect: `dev/m9-pilot/process_pilot.R` (needs extension for per-matrix + per-model loop). CID target: well-corroborated trees in `../neotrans/inst/wct/wellCorroboratedTrees.nwk` (matrix→WCT mapping in `../neotrans/R/AsherSmithSetup.R`).

| Type | Job ID | Matrix | Model | Status | ETA | Mem |
|------|--------|--------|-------|--------|-----|-----|
| HPC | 17217093 | syab07200 | by_nt_9v | RUNNING | ~2d22h | 32G |
| HPC | 17217094 | syab07200 | t_kv | RUNNING | ~2d22h | 32G |
| HPC | 17217095 | syab07200 | t_9v | RUNNING | ~2d22h | 32G |
| HPC | 17217096 | syab07202 | by_nt_9v | RUNNING | ~2d22h | 32G |
| HPC | 17217097 | syab07202 | t_kv | RUNNING | ~2d22h | 32G |
| HPC | 17217098 | syab07202 | t_9v | RUNNING | ~2d22h | 32G |
| HPC | 17217099 | syab07204 | by_nt_9v | OOM→resubmitted | — | — |
| HPC | 17217100 | syab07204 | t_kv | RUNNING | ~2d22h | 32G |
| HPC | 17217101 | syab07204 | t_9v | OOM→resubmitted | — | — |
| HPC | 17217659 | syab07204 | by_nt_9v | PENDING | ~3d | 64G |
| HPC | 17217660 | syab07204 | t_9v | PENDING | ~3d | 64G |
| HPC | 17217102 | syab07205 | by_nt_9v | OOM→resubmitted | — | — |
| HPC | 17217103 | syab07205 | t_kv | RUNNING | ~2d22h | 32G |
| HPC | 17217104 | syab07205 | t_9v | OOM→resubmitted | — | — |
| HPC | 17217661 | syab07205 | by_nt_9v | PENDING | ~3d | 64G |
| HPC | 17217662 | syab07205 | t_9v | PENDING | ~3d | 64G |
| HPC | 17217105 | syab07206 | by_nt_9v | RUNNING | ~2d22h | 32G |
| HPC | 17217106 | syab07206 | t_kv | RUNNING | ~2d22h | 32G |
| HPC | 17217107 | syab07206 | t_9v | RUNNING | ~2d22h | 32G |

Note: syab07201 excluded (no WCT entry). syab07203 pilot done separately (see above).

### GitHub Actions
| Run ID | Workflow | Branch | Status |
|--------|----------|--------|--------|
| 26022493849 | mem-check | main | in_progress |

## Open items / next steps

1. **Collect syab07203 pilot trees** once jobs 17216901/17216902 complete (~16:30
   BST): run `dev/m9-pilot/process_pilot.R` for all three models, produce a 3-row
   CID table vs asher WCT.

2. **Decide on k=15/k=24 feasibility** once 17216080/17216081 finish (~17:30 BST):
   if both fit in 8h, submit full `--array=0-259` for each. If k=24 OOMs or times
   out, rule it out (add to Things ruled out).

3. **Extend `process_pilot.R`** to loop over all 5 matrices × 3 models once long
   runs accumulate ≥50 trees/run (earliest: ~24h from now, so ~09:00 BST 2026-05-19).
   Output: CID table with columns `matrix | model | n | CID_mean | CID_sd`.

4. **Write analysis script** to aggregate CID across all 6 matrices (including
   07203 pilot) and compare t_9v vs t_kv and by_nt_9v vs `by_nt_kv` (the
   existing baseline — pull from neotrans `inst/results/`).

5. **Pull neotrans `by_nt_kv` baseline CID** for the 6 matrices, so we can do the
   (i) and (ii) comparisons from the original plan.

6. **Pre-pub `_9i` (informative coding) variants**: deferred, after pilots confirm
   M9 direction.

## Technical pointers

- **SSH**: PowerShell ssh to `hamilton8.dur.ac.uk` works; Bash-tool ssh does NOT
  (different config). Always route Hamilton commands through PowerShell.
- **PowerShell quoting**: use single-quoted strings for ssh commands containing `$`
  to avoid PowerShell variable interpolation. PowerShell 5.1 lacks `&&` / `||`;
  chain with `;` only when failure is acceptable.
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
- **RB .trees format**: NOT NEXUS. Tab-separated: `Iteration Posterior Likelihood Prior phylogeny`. Read with `read.delim(f)$phylogeny` then `ape::read.tree(text = s)` per row.
- **WCT pruning**: `KeepTip(wcTree, intersect(sampleTips, wcTree$tip.label))` before CID. Matrix→WCT mapping: `paste0("0720", 0:6)` → `wellCorroboratedTrees.nwk` lines 1–7. syab07203→asher (outgroup: Didelphis, Macropus).
- **OOM threshold**: fnJC(9) trans-only model needs 64G for matrices with ≥214 trans.nex lines (syab07204/07205). Matrices with ≤125 lines (07200, 07202, 07206) fit in 32G.
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
  6-way CID pilot. More flexibility always helps in this sim, up to k=9. Whether
  the ramp continues past k=9 is what k=15/k=24 arms test.

## Suggested first action

```powershell
# Check whether the 07203 pilots have finished and, if so, pull their trees:
ssh pjjg18@hamilton8.dur.ac.uk 'sacct -j 17216901,17216902 --format=JobID,State,ExitCode,Elapsed'
# If COMPLETED, pull and run process_pilot.R:
scp 'pjjg18@hamilton8.dur.ac.uk:/nobackup/pjjg18/m9-pilot/syab07203/by_nt_9v_run_*.trees' dev/m9-pilot/syab07203/
scp 'pjjg18@hamilton8.dur.ac.uk:/nobackup/pjjg18/m9-pilot/syab07203/t_9v_run_*.trees'     dev/m9-pilot/syab07203/
Rscript dev/m9-pilot/process_pilot.R
```
