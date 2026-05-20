# mkp — hand-off 2026-05-20 (prior-vs-fixed-k report patched to n=260 mk_tlshrink)

## What this repo is

MkPrime — R package for Bayesian phylogenetic inference under the Mk′ model
for discrete morphological characters. The current focus is a 2026-05-12
prior-validation pilot comparing arms on 26 simulated trees × 10 reps to find
which k-specification best recovers the true tree (CID-to-truth). A 13th arm,
*mk_tlshrink*, runs as an interventional test of the regularisation-via-saturation
hypothesis. A parallel real-data pilot (in neotrans) tests whether the winning
arm from simulation transfers to six empirical matrices compared against
well-corroborated reference trees.

## Where we left off

This session (continuation):

- **Completed mk_tlshrink array** (Hamilton job 17234909, all 260 tasks).
  Submitted summariser job 17245451 for the 205 outstanding indices; rebuilt
  the digest CSV at n=260; refreshed `mkprime/report-data/cid_thirteen_arm_digest.csv`;
  re-rendered `mk-prime-prior-vs-fixed-k.qmd` (mkprime commit dc65348, pushed).
  Patched the qmd to remove the "partial" TODO and update fig caption /
  panel title to "250 of 260" (was "every"). Headline at n=260:

  | Arm | Mean CID (n = 260 paired) | Mean TL |
  |---|---:|---:|
  | mk_k40       | 0.2385 | 1.19 |
  | mk           | 0.2754 | 1.39 |
  | mk_tlshrink  | 0.2723 | 0.78 |

  mk_tlshrink loses to mk_k40 on 250/260 paired tasks (binomial sign test,
  p = 3.67 × 10⁻⁶¹), and beats kObs-Mk in 165/260 (p = 1.68 × 10⁻⁵) with a
  margin of 0.0031 CID — well below the 0.0369 mk vs mk_k40 gap. Mean TL =
  0.775 (stable vs the 0.78 interim value). The n=55 first batch was
  favourable to mk_tlshrink (then ranked fifth on the leaderboard); at full
  n=260 it sits second-worst, just above mk. **Conclusion unchanged**: TL
  shrinkage alone is insufficient to replicate the mk_k40 advantage; the
  saturation effect on the topology likelihood is the missing piece.

Previous session:

- **Pulled mk_tlshrink interim results** (n=55/260; first batch favourable).

- **Drafted the co-author write-up** at
  [mkprime/mk-prime-prior-vs-fixed-k.qmd](../mkprime/mk-prime-prior-vs-fixed-k.qmd)
  (mkprime commit 11f7f5d, pushed to ms609/mkprime). Three findings,
  themed-not-chronological, in `/method-voice`:
    1. Mk′'s data-aware prior is a passenger — per-character k′ posterior is
       prior-dominated; aggregate CID is flat across Mk′ family; even the
       oracle (mk_ktrue) loses to mk_k40.
    2. Mk with fixed large k monotonically wins (the ramp).
    3. Mechanism is saturation regularisation, not TL inflation
       (TL_HYPOTHESIS_TEST.md) and not explicit TL shrinkage (mk_tlshrink).
  Supporting CSV digest at `mkprime/report-data/cid_thirteen_arm_digest.csv`
  (2941 rows across 13 arms; pulled from
  `/nobackup/pjjg18/mkp-study/summary/*.rds`).
  Bib expanded at `mkprime/inst/REFERENCES.bib` from 1 to 7 entries
  (added Wright2014, Harrison2015, GelmanRubin1992, Felsenstein1981,
  Geyer1992, SmithCID).

- **db4c348 (pre-session)** — TreeESS exported (Option A: drop dot prefix).
  This was committed before the conversation compact and is in the package now.

## 13-arm CID results (full n=260 mk_tlshrink)

| Arm | n | Mean CID | Family |
|---|---:|---:|---|
| **mk_k40** | 260 | **0.2385** | fixed-k (winner) |
| mk_k24 | 260 | 0.2398 | fixed-k |
| mk_k15 | 260 | 0.2418 | fixed-k |
| mk_k9 | 260 | 0.2451 | fixed-k |
| mk_kp2 | 260 | 0.2535 | fixed-k |
| mkp_geo | 26 | 0.2573 | Mk′ geometric |
| mk_kp1 | 260 | 0.2595 | fixed-k |
| mkp_logs | 260 | 0.2595 | Mk′ logseries(c=0.95) |
| mk_ktrue | 260 | 0.2624 | oracle (k=k_true per char) |
| mkp_highk | 260 | 0.2639 | Mk′ geom Beta(1,20) |
| mkp_eg | 260 | 0.2677 | Mk′ empirical_geom |
| **mk_tlshrink** | 260 | **0.2723** | Gamma(20, 20/0.7) on TL, k = kObs |
| mk | 260 | 0.2754 | kObs |

At n=260, mk_tlshrink sits second-worst — modestly above the kObs baseline
mk (Δ = 0.0031, p = 1.7e-5) but well below every fixed-k arm and all Mk′
variants except mkp_eg. Its mean posterior TL (0.775) is the shortest of any
arm — about half of truth (1.40) and well below mk_k40's 1.19 — yet that
shrinkage buys almost nothing on CID. Saturation does the work.

## Pending jobs

| Type | ID | Status | ETA | On completion |
|------|----|--------|-----|---------------|
| HPC | 17217093-17217107 + 17217660-17217662 | 14 M9 long-form jobs, RUNNING (2d 5h of 3d walltime) | ~20 h (~04:30 BST 2026-05-21) | scp `.trees` from `/nobackup/pjjg18/m9-long/<matrix>/`; run `dev/m9-pilot/process_pilot.R <pid>` per matrix; build 6-matrix × 3-model CID table. Compare against 07203 (asher). |

## Open items / next steps

1. ~~**Patch mk-prime-prior-vs-fixed-k.qmd when mk_tlshrink finishes.**~~
   **Done 2026-05-20** (mkprime commit dc65348). Digest rebuilt at n=260,
   scp'd, qmd re-rendered with updated narrative and figure captions.

2. ~~**Resolve Felsenstein1981 vs 1978 citation choice**~~ **Done 2026-05-20**
   (mkprime commit 6036a49). Added Felsenstein1978 (Syst Zool 27:401–410,
   "Cases in which parsimony or compatibility methods will be positively
   misleading") to `inst/REFERENCES.bib`; swapped `[@Felsenstein1981]` →
   `[@Felsenstein1978]` in §3.3; re-rendered HTML.

3. **Flip the AI-callout to reviewed** in the qmd when the user has read it
   (template note already in the callout comment block).

4. ~~**Diagnose the spurious R parser error** in mk_ktrue / mk_tlshrink~~
   **Workaround applied 2026-05-20** (mkp commit 75dd825; live on Hamilton
   at `/nobackup/pjjg18/mkp-study/run_one.R`). Added `options(warn = 1L)`
   at the top of `run_one.R` so warnings print immediately rather than
   queuing for R's exit-time `Warning messages:` formatter, which was
   triggering a srcref-into-comment-block parser error on the only two
   arms that emit warnings during MCMC. **Caveat**: symptom-targeted, not
   root-caused — the agent could not reproduce the parser error in
   isolation, and the em-dash theory does not unify mk_ktrue and
   mk_tlshrink. Verify on the next array submission that tasks show
   COMPLETED rather than FAILED; if not, the diagnosis is wrong.

5. **Long-form M9 collect** (~2026-05-21 morning local): scp `.trees`, run
   `process_pilot.R <pid>` for 07200, 07202, 07204, 07205, 07206; build the
   full 6-matrix × 3-model CID table. Compare against 07203 (asher).

6. **Decide on 17217659 retry** at 128 G vs skip syab07204 by_nt_9v.

7. **Re-frame the paper around regularisation-via-saturation** — the new
   findings invalidate the original "high-k beats kObs because the data have
   hidden multistate" framing. The new write-up in mkprime is the first
   draft of that re-frame.

8. ~~**Commit canonical run_one.R + auxiliary scripts**~~ **Done 2026-05-20**
   (mkp commits f72fddc, 75dd825). All Hamilton dispatch scripts now under
   VC in `data-raw/hamilton/`: `run_one.R`, `summarize_streamed.R` (already
   matched Hamilton byte-for-byte), `summarize_array.slurm`, the per-arm
   `*_array.slurm` files including newly-fetched `mkp_eg_array.slurm`, plus
   previously-untracked `install_dt.R`, `summarize_inspect.sh`,
   `verify_summary.R`. `run_one.R` includes the item-#4 `options(warn = 1L)`
   fix.

9. **Direct test of saturation hypothesis**: re-run a single task with mk_k40
   but with a *very tight* prior on tree length forcing TL = 1.4 (the truth).
   If CID degrades back toward Mk's level, the mechanism is confirmed as
   branch-length-mediated.

10. **The mk_no_relabel arm**: take Mk′ (with kPrimePrior) but bypass the
    relabel correction (treat characters as type="known" while still
    inferring k′). Tests whether the relabel correction is the gap between
    mkp_highk and mk_k40.

11. **Pull neotrans `by_nt_kv` baseline CID** for the 6 matrices.

12. **Clean 35 contaminated rows** in ground_truth.csv files (orig_idx ≤ 0,
    rbind leakage). Per-character analyses only.

13. **Install new MkPrime package on Hamilton** — currently-running jobs (mk_tlshrink
    17234909, M9 long-form) use the OLD package loaded at job start. New per-run
    tree files + parallel ResumeMkPrime are committed locally but not yet
    installed on Hamilton. Defer until current runs finish to avoid mid-run
    behaviour drift.

## Technical pointers

- **SSH**: PowerShell ssh to `hamilton8.dur.ac.uk` works; Bash-tool ssh does NOT
  (different config). Always route Hamilton commands through PowerShell.
- **PowerShell tool history**: The Claude Code PowerShell tool failed silently
  (exit 1, no output) in all sessions prior to 2026-05-20 because pwsh was
  installed only via the Windows Store (app execution alias — can't be spawned
  from Electron). Fixed by installing the MSI via `choco install powershell-core`.
  The MSI puts a real exe at `C:\Program Files\PowerShell\7\pwsh.exe`. Hamilton
  SSH was working via the `hamilton-hpc` skill, not the PowerShell tool directly.
- **PowerShell quoting**: use single-quoted strings for ssh commands containing
  `$` to avoid PowerShell variable interpolation. PowerShell 5.1 lacks `&&` /
  `||`; chain with `;` only when failure is acceptable.
- **PowerShell heredocs for ssh**: PowerShell `@'...'@` then piped to `ssh ... 'cat > /tmp/X.R && Rscript /tmp/X.R'` works cleanly for multi-line R scripts on Hamilton. Avoid embedded heredoc-in-heredoc (broke 2026-05-20).
- **Background SSH on Windows can hang** if the ssh session is initiated as a
  Bash background task; route through PowerShell instead. Even PowerShell ssh
  can be slow to first-handshake (~30s) when the cluster is under load.
- **Hamilton paths (mkp sim study)**:
  - Run scripts: `/nobackup/pjjg18/mkp-study/{run_one.R, summarize_streamed.R, *_array.slurm}`
  - Raw streamed logs: `/nobackup/pjjg18/mkp-study/results/t{NN}_r{MM}/{arm}_run_1.log`
  - Tree files: `/nobackup/pjjg18/mkp-study/results/t{NN}_r{MM}/{arm}_trees.nwk` (legacy) or `{arm}_trees_run_*.nwk` (new code)
  - Summaries: `/nobackup/pjjg18/mkp-study/summary/{arm}_t{NN}_r{MM}.rds`
  - Data root: `/nobackup/pjjg18/mkprime-files/tree-inference/tree_{NN}/rep_{MM}/chr*.nex`
- **Hamilton paths (M9 real-data pilot)**:
  - Pilot (4h): `/nobackup/pjjg18/m9-pilot/syab07203/{model}_run_{1,2}.trees`
  - Long runs (3d): `/nobackup/pjjg18/m9-long/<matrix>/{model}_run_{1,2}.trees`
  - Logs: `/nobackup/pjjg18/m9-long/logs/m9L-<matrix>-<model>_<jobid>.{out,err}`
  - Source nex files: `/nobackup/pjjg18/07200_by_t_kv/`, `/nobackup/pjjg18/07202_ns_n_kv/`, etc.
- **Sibling repo for write-ups**: `../mkprime` is served via GitHub Pages.
  New reports go there, not in mkp itself. Pattern: report YAML matches
  `mk-prime-simulation-report.qmd`; data digests in `report-data/`; bib in
  `inst/REFERENCES.bib`.
- **Local tree files are gitignored**: `dev/m9-pilot/syab*/` is in `.gitignore`
  (intentional — big outputs); CID summaries are too. Record numerical results
  in commits / RESUME.md, not files.
- **RB .trees format**: NOT NEXUS. Tab-separated: `Iteration Posterior Likelihood Prior phylogeny`. Read with `read.delim(f)$phylogeny` then `ape::read.tree(text = s)` per row.
- **WCT pruning**: `KeepTip(wcTree, intersect(sampleTips, wcTree$tip.label))` before CID. Matrix→WCT mapping: `paste0("0720", 0:6)` → `wellCorroboratedTrees.nwk` lines 1–7. syab07203→asher (outgroup: Didelphis, Macropus).
- **OOM threshold (Mk′ M9 real-data)**: fnJC(9) trans-only model needs ≥64G for
  matrices with ≥214 trans.nex lines (syab07204/07205). 64G itself is **not
  always enough** — see 17217659 (syab07204 by_nt_9v at 64G OOM'd at 1m43s).
  by_nt models may need 128G+ for these larger matrices.
- **mk_k24 convergence behaviour**: uses built-in convergence criterion; most tasks
  finish in 1–4h not 8h. Tree counts per task: 750–66k (wide range by difficulty).
  The earlier "43× slower / ~2,350 samples" estimate was a pilot artefact — disregard.
- **Per-arm MCMC speed (mkp sim)**: rough iter/sec on t01_r01 (from tree counts):
  - mk_k9:  ~120k trees in 8h (hits walltime, thin=100)
  - mk_k15: ~105k trees, converges; ~74 tasks hit 8h walltime
  - mk_k24: converges early (~750–66k trees); 172/260 done within 4h of submission
- **kObs alignment**: `chr*.nex` filenames are unpadded. `run_one.R` uses plain `sort()` → lex order (chr1, chr10, chr11, …). Analysis must match.
- **EG arm result (corrected, 2026-05-17)**: median u_post = 1.01 (binary chars, n=9917); EG posterior on k′ is calibrated near kObs+1. mk_kp1 collapses this to a point estimate, beating EG on CID.
- **No oversampling**: thin streamed MCMC to ~30k samples/run (thin_iters ≥ 16 for 8h job). See `feedback_no_oversample.md`.
- **Summariser production vs canonical**: `tmp_summarize_streamed.R` (repo root) is the production working copy on Hamilton; `data-raw/hamilton/summarize_streamed.R` is stale (missing mk_k15/mk_k24 in match.arg). Reconcile before next major summariser change.
- **mkprime/inst/REFERENCES.bib**: as of 2026-05-20 contains Lewis2001, Wright2014, Harrison2015, GelmanRubin1992, Felsenstein1981, Geyer1992, SmithCID. Felsenstein1978 (LBA) is NOT in there — add before next revision.
- **mk_tlshrink TL prior**: Gamma(20, 20/0.7) = shape 20, rate 20/0.7 ≈ 28.57 → mean = shape/rate = 0.7, sd = √(shape)/rate = √20/28.57 ≈ 0.157. Coded in run_one.R; not yet committed to canonical.

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
  6-way CID pilot. Extended in 7-way: mk_k15 > mk_k9 (p=4.3e-22). Extended in
  8-way: mk_k24 > mk_k15 (p=2.0e-14). Ramp is monotonic up to k=24 with strongly
  diminishing returns (Δk15→k24 = −0.0020 vs Δkp2→k9 = −0.0084).
- **mk_k24 infeasibility**: earlier "43× slower" pilot estimate was wrong. mk_k24
  converges via built-in criterion, finishing within the 8h walltime for all tasks.
- **High-k forces longer branches (the natural intuition)**: rejected
  2026-05-19. TL goes DOWN monotonically with k (mk_k40 TL=1.19, mk TL=1.39),
  not up. Spearman ρ(ΔTL, ΔCID) = −0.37, p=1.9e-9.
- **Pure TL-shrinkage explanation of mk_k40 advantage**: rejected 2026-05-20.
  mk_tlshrink with Gamma(20, 20/0.7) prior on T achieves TL=0.775 (shorter
  than mk_k40's 1.19) but CID=0.2723 vs mk_k40's 0.2385 (n=260 paired,
  10/260 wins for tlshrink, p=3.67e-61). TL shrinkage is necessary but not
  sufficient; the saturation effect on the topology likelihood is the
  missing piece.

## Suggested first action

Items #2, #4, #8 cleared 2026-05-20 by parallel subagents. Remaining
priorities: #3 (flip AI-callout to reviewed once the user has read
the report); #7 (re-frame paper around regularisation-via-saturation);
#5 (M9 long-form collect when the 14 RUNNING jobs land, ~04:30 BST
2026-05-21); #9 (direct saturation test — mk_k40 with very tight TL
prior pinning TL = 1.4). #7 is the heaviest piece; #3 is a one-line
edit; #9 needs a new Hamilton array but is the cleanest mechanism
confirmation.
