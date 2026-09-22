# Prior validation campaign — log

Operational record of the empirical campaign comparing priors on k' for the mkp validation work. Started 2026-05-12 with commit 9332a54 (empirical_geometric default). Companion to the code-level audit in `dev/red-team/`.

Entries are dated, newest at top. Each entry: what was run, what came back, what was decided. Job IDs, file paths, and finding IDs cross-reference the red-team trail.

---

## 2026-05-16 — Summarization arrays + EG-001 empirical test: **u_post anchor is prior shape, not truncation bug**

**Summarization 17185790 (mk, 260) + 17185791 (geo, 26):**
First submission failed all 286 tasks instantly due to a stray `}` in `summarize_array.slurm` — the `${1:?usage: ... {mk|mkp_geo}}` parameter expansion's first `}` closed the `${...}`, leaving the trailing `}` as a literal suffix on `$ARM` so the case dispatch read `mk}` / `mkp_geo}`. Fix: removed inner braces from error message. Resubmitted; all 286 COMPLETED, 117 MB of summary RDS files pulled to `dev/pilots/2026-05-12-prior-validation/summary/`.

**EG-001 u_post comparison (geo arm vs EG arm, same 26 trees × rep_01):**

| prior | n_chars | mean(u_post) | median | p(u_post<0.5) | bias=mean(u_post−u_true) | Spearman(u_post,u_true) |
|-------|---------|--------------|--------|---------------|---------------------------|-------------------------|
| EG    | 1295    | **1.31**     | 1.16   | 22%           | **+0.07**                 | **−0.08** (per-task −0.11±0.15) |
| geo   | 1295    | **0.14**     | 0.11   | 81%           | **−1.11**                 | **−0.10** (per-task −0.12±0.15) |

**Key finding — EG-001's hypothesis is refuted as a cause:** EG-001 claimed the truncation bug under-penalises high p and over-shrinks u, producing the u_post≈1 anchor. The geo prior has no such truncation issue. Empirically: **geo u_post is MORE anchored toward 0** (mean 0.14, 81% < 0.5) than EG (mean 1.31, 22% < 0.5). Both arms show Spearman(u_post, u_true) ≈ −0.1, i.e. u is **unidentified by 50-character likelihoods** under either prior. The posterior shape is dominated by the prior in both cases:
- geo with p_post ≈ 0.86 → steep geometric (1−p)^k decay → tiny u_post.
- EG with body decay 0.67, 0.21, 0.067 → moderate decay → u_post ≈ 1.

**Status changes:**
- **EG-001** marked REFUTED-AS-CAUSE. The bug exists as a coding defect; it is not the cause of the live u-anchor observation.
- **EG-003** upgraded MED → HIGH and marked CONFIRMED. "Prior shape dominates u_post" is now the headline scientific issue, not a side concern.

Outputs: `dev/pilots/2026-05-12-prior-validation/analysis/eg001_upost_compare.{R,rds,png}` + `eg001_upost_vs_utrue.png`.

**CID three-prior comparison (tree recovery):**

| comparison | n=26 common (rep_01) | n=260 all tasks |
|------------|----------------------|------------------|
| EG  | mean=0.249, sd=0.050 | mean=0.267, sd=0.057 |
| geo | mean=0.256, sd=0.051 | (26 tasks only) |
| mk  | mean=0.260, sd=0.051 | mean=0.275, sd=0.057 |

Paired Wilcoxon on common 26 tasks:
- geo − EG = +0.0067 (p=6.6e−4) — EG slightly better than geo
- mk  − EG = +0.0101 (p=5.7e−7) — EG slightly better than mk floor
- geo − mk = −0.0034 (p=0.063) — geo marginally better than mk

**Interpretation:** Tree recovery is robust across priors — all three give CID ≈ 0.25–0.27 — but the ordering is statistically clean: **EG < geo < mk floor**. Having a k′ parameter at all (EG or geo) recovers slightly better trees than fixing k′ = kObs, even though u itself is unidentified (see EG-001/EG-003 above). Differences are small (≈ 1% CID) but reproducible across tasks. The prior on k′ shapes the **u marginal** strongly (above) but shapes **tree topology** only weakly. Outputs: `dev/pilots/2026-05-12-prior-validation/analysis/cid_three_prior.{R,rds,png}`.

Next: file HARNESS-001 investigation as a separate task (R maxTime=7.5h did not fire on Hamilton despite fix being deployed; all 286 production tasks TIMEOUT'd at 8h SLURM wall — data fully recovered via streaming, but clean resume semantics broken).

## 2026-05-15 (night) — STREAM-001 empirically confirmed FIXED; mk resubmission unblocked

**Geo pilot 17167921 streaming confirmed at ~40 min elapsed:**
- t01_r01: 62,501 lines in `mkp_geo_trees.nwk`
- t05_r01: 57,501 lines
- t10_r01: 77,001 lines

All tasks RUNNING. Trees streaming normally across the board. **STREAM-001 status updated to FIXED** in `dev/red-team/findings.md`.

**Side observation:** `[DIAG] fopen failed iter=XXXXXX` appears every iteration in `.err` logs — this is `src/mcmc.cpp:4902` failing to open `pcl_diag.txt` in the SLURM cwd. Not affecting R-side tree streaming; tree files are written correctly. The fopen failure is a C++ debug artefact (Hamilton scratch may be read-only from the job cwd). No action required for the campaign.

**Results directory structure clarified:** per-task output is at `/nobackup/pjjg18/mkp-study/results/t{NN}_r{MM}/` (e.g., `t01_r01/`), not `t{NN}/rep_{MM}/` as assumed in earlier log entries.

**mk resubmission:** submitted as job **17170442** (260 tasks, `mk_array.slurm`, `--mem=3G`, `--time=08:00:00`). All SLURM mem values recalibrated to ≤50% headroom over observed MaxRSS: mk→3G (peak 2.05 GB × 1.5), geo/mkp→2G (EG proxy peak 914 MB × 1.5 = 1.37 GB; will tighten after geo pilot MaxRSS available).

**Mid-run check at ~1h40m elapsed (15:10 UTC, 2026-05-15):**
- mk array 17170442: 260 tasks RUNNING (some still in warmup at iter ~85k, e.g. task 124). Lower-index tasks already streaming: t01_r01=26k lines, t05_r01=6.8k, t10_r01=13.5k. **STREAM-001 fixed for mk arm too** — empirical confirmation; previous-run 1-line files (May 14 timestamps) are stale leftovers, not new-run output.
- geo pilot 17167921: 26 tasks RUNNING, all streaming at 135–156k lines.

**Mid-run check #2 at ~2h40m elapsed (16:13 UTC):** mk 0 COMPLETED / 260 RUNNING — all 7 sampled tasks now streaming (74k–163k lines, mtimes 16:13 today). geo 0 COMPLETED / 26 RUNNING (184k–210k lines). STREAM-001 comprehensively confirmed for both arms.

**Mid-run check #3 at ~3h45m elapsed (17:15 UTC):** mk 0 COMPLETED / 260 RUNNING (235k/246k/218k for t01/t13/t26 — adding ~80k lines/hr). geo 0 COMPLETED / 26 RUNNING (226k/220k/263k — ~50k/hr). Trees growing healthily; graceful-exit expected ~21:00 UTC when maxTime=7.5h fires.

## 2026-05-16 (morning) — Both jobs hit 8h SLURM wall; STREAM-001 fully validated, HARNESS-001 reopened

**Final state:** all 260 mk + 26 geo tasks ended with state TIMEOUT (SLURM kill at 8h00m00s wall). 0 partial `.rds` files written — maxTime=7.5h graceful exit did **not** fire. mk task 0's last Sample-phase line showed `ETA: ~113h` — the mk floor benchmark cannot converge in 8h for this dataset, as expected. geo arm tasks were closer (early-run ETA ~7.6h, may have approached convergence).

**Recoverable output:**
- mk arm: 260/260 `mk_trees.nwk` files, all >100 lines (vs previous run: 232/260 empty). **STREAM-001 fully confirmed FIXED.**
- geo arm: 26/26 `mkp_geo_trees.nwk` files, all >100 lines.
- All scalar samples streamed to `*_run_{1,2}.log` files (post-Round-4 checkpoint fix).
- Periodic checkpoints written to `*_checkpoint.rds` (post-Round-4 STREAM-002 fix).
- **Missing:** 0/260 mk + 0/26 geo partial summary RDS files (HARNESS-001 regression).

**HARNESS-001 reopened** (`dev/red-team/findings.md`): the R-level maxTime check in `RunMkPrime`'s Sample-phase orchestrator does not appear to fire between SLURM provision and SLURM kill. Investigation deferred — the streamed data is sufficient for the campaign analysis.

**Next:** proceed with three-prior comparison (EG / geometric / mk floor) from the streamed scalar `.log` files and `.nwk` tree files. Pull to mkprime sister repo.

---

## 2026-05-15 (night) — Package reinstalled on Hamilton; geo pilot submitted

**Install job 17167855** — COMPLETED (exit 0:0, 1m51s). MkPrime_0.0.0.9000.tar.gz installed under `/nobackup/pjjg18/mkp-study/lib` with post-Round-4 fixes (STREAM-002, STREAM-003, STREAM-001 partial, HARNESS-001).

**Geo pilot submitted: job 17167921** (`mkp_geo_array.slurm`, 26 tasks, `--time=08:00:00`). This is the empirical test of:
1. STREAM-001 (mk-arm tree streaming) — if `mkp_geo_trees.nwk` have >1 line for non-converged tasks, the brColStart + resume-treeFile fixes are sufficient.
2. EG-003 prior shape (geometric-prior arm as control for u_post anchor under EG).

**Next check:** poll job state + first `mkp_geo_trees.nwk` line counts after ~1h. If trees are streaming (≥100 lines per non-converged task by hour 2), STREAM-001 closes empirically and the mk resubmission (260 tasks) can proceed.

---

## 2026-05-15 (evening) — Streaming bugs investigated; fixes integrated; HARNESS-001 closed

**Opus subagent (worktree `agent-ab5d3f56687d5fa3e`) report integrated into main.** Round 4 entry in the red-team log archive (discussion #124) has the full diagnostic narrative; finding-status updates:

- STREAM-002 → **FIXED**. Root cause: `.RunSerialRuns` was passing `checkpointFile = NULL` to `.RunMkPrimeSingleRun` at `R/RunMkPrime.R:1486` and `:1572`, suppressing every per-batch save under the standard `nRuns >= 2 && maxRhat` orchestrator. Only the iter=0 init snapshot and post-Phase-1/per-epoch saves ever fired. Fix passes `mcmc$checkpointFile` through and rewrites the in-batch save call sites to checkpoint `shared$runs` (all runs' coherent state) rather than `list(r)` (current run only). Regression: `tests/testthat/test-checkpoint.R::"Cross-run R-hat orchestrator saves checkpoint mid-run"`.
- STREAM-003 **NEW + FIXED**. `brColStart` calc in both `RunMkPrime` (:277) and `ResumeMkPrime` (:2369) omitted the two diagnostic columns C++ inserts at every saved row (`swap_cold` and `topo_hash`, written at `src/mcmc.cpp:4864-4867`). Reconstructed trees had 1–2 garbage edge lengths (one of them `topo_hash`, ~1e16) and dropped the last 1–2 true `br_*` values. Affects **both** arms, not just mk. Trees on disk were parseable Newick but unusable for any branch-length-dependent analysis. Regression: `tests/testthat/test-tree-thin.R::"mk arm ... writes trees with finite edge lengths"`.
- STREAM-001 → **PARTIAL** (subagent's honest call, kept). Two causal-chain defects fixed (STREAM-003 above + `ResumeMkPrime` passing `treeFile = NULL` to inner calls, dropping all post-resume tree writes). But local reproducer at `dev/red-team/repros/stream-001-002.R` does NOT reproduce empty-tree-file symptom even pre-fix, so the Hamilton-specific factor producing 232 empty `mk_trees.nwk` files remains unexplained. Most likely remaining hypothesis: NFS write semantics under SIGKILL, or a Hamilton R 4.5.1 vs local R-devel difference. Will be re-tested empirically by the mk arm rerun — if trees now stream as expected, STREAM-001 transitions to FIXED; if not, Hamilton-side `message()` instrumentation per Round 3's "next reviewer" note is required.
- STREAM-004 **NEW + OPEN** (LOW). brColStart off-by-1 for `kPrimePrior = "beta_geometric"` — that prior contributes 2 hyperparam cols, not 1. Out of scope for this campaign (no beta_geometric arm planned); filed for the next reviewer.
- HARNESS-001 → **FIXED**. `data-raw/hamilton/run_one.R:171` `maxTime` lowered from `10 * 3600` to `7.5 * 3600` so the R graceful-stop fires inside the 8h SLURM wall.

**Implication for prior-validation campaign:**

The Hamilton checkpoints already on disk (260 mk + 260 EG, all iter=0) were written by the buggy code path and are **unrecoverable** for resume. Both arms need to rerun from scratch with the fixed code installed on Hamilton. Path forward:

1. Build the package locally and push the source tarball (or git-pull on Hamilton).
2. Reinstall under `/nobackup/pjjg18/mkp-study/lib`.
3. Submit `mkp_geo_array.slurm` (26 tasks) as the first test — if trees stream correctly and checkpoints advance past iter=0, STREAM-001 closes empirically.
4. Resubmit mk array (260 tasks) once geo confirms streaming works.



**New artefacts in repo, not yet on Hamilton:**

- `data-raw/hamilton/run_one.R` — added `mkp_geo` arm (lines 305–348). Pattern matches `mkp_eg` block but passes `kPrimePrior = "geometric"` to `MkPrimeModel(...)`. Output goes to `mkp_geo_<tag>.rds`, checkpoint prefix `mkp_geo`. `match.arg` updated to accept the new arm.
- `data-raw/hamilton/mkp_geo_array.slurm` — 26-task pilot (1 rep per tree, `TREE=SLURM_ARRAY_TASK_ID+1`, `REP=1`). Job name `mkp-geo`, 8h wall (same trap as the failed mk arm — flagged inline as HARNESS-001).

**`kPrimePrior` value confirmed** against `R/MkPrimeModel.R:132` (`match.arg(c("empirical_geometric", "geometric", "beta_geometric", "logseries"))`).

**Structural finding worth recording:** `mkp_eg_array.slurm` does **not** exist in the repo — only on Hamilton (at `/nobackup/pjjg18/mkp-study/mkp_eg_array.slurm`). The repo has `mk_`, `mkp_`, `combine_`, `mkp_check`, `mkp_test`, and now `mkp_geo_array.slurm`. The eg SLURM file was presumably authored on Hamilton and never copied back. Reproducibility gap: if the Hamilton scratch is lost, the eg arm submission is not bit-identical reproducible from this repo alone. Worth pulling a copy back at some point.

**Smoke test not run.** Driver requires `tree_NN/rep_MM/chr*.nex` + `tree_NN/tree.nwk` data layout under `data_root`; the only local `.nex` is the flat `project3832.nex`, not in harness format. R-level dispatch (`match.arg`, branch routing) verified by reading the diff; runtime not verified.

**Not yet submitted.** Pilot is gated on STREAM-001 fix landing — if the streaming bug also affects this arm, the pilot will produce empty `mkp_geo_trees.nwk` for any task killed mid-chain. Once the fix lands and is verified, the path forward is: copy `run_one.R` + `mkp_geo_array.slurm` to Hamilton, `sbatch mkp_geo_array.slurm`.

## 2026-05-15 — mkp-mk arm post-mortem; bugs filed; pivot to bug-fix + geometric pilot

**Jobs reviewed:** 17154895 (mkp-mk array, 260 tasks).

**Outcome:** Full SLURM TIMEOUT for all 260 tasks at 2026-05-14T21:07:14 (8h wall, exit 0:0). Resume context's "~4h remaining at 01:23 UK on 2026-05-15" was incorrect — the array had been dead for ~28h.

**Recoverable output, per-task inventory of `/nobackup/pjjg18/mkp-study/results/t*/`:**

| arm | `*_trees.nwk` ≤1 line | ≥100 lines | per-task tree count (when written) |
|---|---|---|---|
| `mkp_eg_trees.nwk` | 0 / 260 | 260 / 260 | 138k–262k |
| `mk_trees.nwk` | **232 / 260** | 28 / 260 | 720k–846k (those that wrote) |

Per-task `mk_run_1.log` is fully populated for all tasks (Sample-phase confirmed — streaming log writes fired). Per-task `mk_<tag>.rds` (the `saveRDS` partial after `RunMkPrime` returns) does not exist for any of the 260 tasks — R never reached the post-call save block under either timeout or convergence stop.

**Three bugs filed in `dev/red-team/findings.md`:**
- STREAM-001 HIGH — mk-arm tree streaming fails for unfinished runs ([R/RunMkPrime.R:925–942](../../../R/RunMkPrime.R))
- STREAM-002 HIGH — periodic checkpoint flush never reaches disk ([R/RunMkPrime.R:1229–1251](../../../R/RunMkPrime.R)); resume-from-checkpoint is functionally a restart for the 232 unfinished tasks
- HARNESS-001 LOW — `maxTime=10h` in `run_one.R` exceeds `--time=08:00:00` in `mk_array.slurm`; no graceful exit window

**Decisions taken (user-confirmed):**
1. Fix STREAM-001 before any mk resubmission. Investigation thread: instrument the Sample-phase tree write locally, no Hamilton round-trip needed.
2. Fix STREAM-002 in parallel (separable bug; same buffer-flush mechanism).
3. Launch a small geometric-prior pilot to test EG-001 directly. EG-001 is logically independent of mk-vs-EG performance comparison: it concerns whether the EG prior is correctly normalised, not whether EG is better than fixed-k'=kObs. The u_post≈1 anchor diagnostic is EG-specific and would remain unexplained even with mk in hand. Pilot scope TBD but on the order of 26 trees × 1 rep.
4. Defer full mk re-run until STREAM-001 fix lands.

**What survives for paper-level reporting from this round:**
- EG arm complete (260/260): `report-data/mkp-eg-260/post_means.rds`, `report-data/mkp-eg-260/cid_compare_260.rds` (the latter currently with EG=260 / mk=28 paired).
- mk arm only 28 paired tasks for CID comparison; biased toward early-finishing (likely smaller-data) tasks.
- Preliminary tree-recovery direction (n=28): EG closer to truth in 24/28, mean ΔCID −0.006. Not paper-strength on its own.

---

## 2026-05-14 — Red-team rounds 1 + 2; EG-001 surfaced

Two Opus-Xhigh red-team passes on EG prior math (round 1, area #1) and the case-30 logit-MH p sampler (round 2, area #2). 10 findings + 1 info note filed; details in `dev/red-team/{findings,log}.md`. Headline: **EG-001 HIGH** — missing per-character truncation normaliser in `LogPrior` (`R/MkPrimeModel.R:489–497` + `src/mcmc.cpp:269–315`). Sampler math (case 30) is correct and **cannot** be the cause of the u-anchor; the leading remaining hypothesis is EG-001.

## 2026-05-13 to -14 — mkp-eg arm completed; u_post anchor surfaces

**Jobs:** 17152099 + 17152211–17152219 (mkp-eg array, 260 tasks). All converged.

Posterior summaries pulled to `report-data/mkp-eg-260/post_means.rds`. Diagnostic: pooled Σu_post/Σu_true = 0.952 (encouraging at the aggregate level) but per-character u_post is anchored near 1 regardless of u_true (Spearman ≈ −0.06, MAE 1.53). This anchor is the empirical anomaly that motivated the red-team rounds.

## 2026-05-12 — Campaign opened

Commit 9332a54: empirical_geometric prior on k' made default. Campaign design: three priors compared on the same 26-tree × 10-rep validation set — mkp-eg (above), mkp (native logseries), mk (no u, fixed k'=kObs). The mkp native arm was later cancelled to reduce scope; mk was intended as the floor benchmark and is the one now blocked on STREAM-001.

---

## Cross-references

- Red-team findings table: [dev/red-team/findings.md](../../red-team/findings.md)
- Red-team session log: [discussion #124](https://github.com/agent-issues/MkPrime/discussions/124) (archive of the retired `dev/red-team/log.md`)
- Red-team focus rotation: [dev/red-team/focus-areas.md](../../red-team/focus-areas.md)
- Hamilton HPC conventions: `~/.claude/skills/hamilton-hpc/SKILL.md`
- Posterior summaries (EG, complete): `report-data/mkp-eg-260/` in the sister `mkprime` repo
