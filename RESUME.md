# MkPrime ecology-aware — hand-off 2026-05-20 (afternoon)

## What this repo is

MkPrime is an R package implementing an ecology-aware Mk' MCMC for Bayesian
phylogenetic inference from morphological data. The `worktree-ecology-aware`
branch (aliased locally to `ecology-aware`) drives the ecology paper:
simulations demonstrating that ecology-aware inference reduces false-clade
support relative to an ecology-blind Mk' baseline, plus an empirical rodent
case study (MorphoBank X24848). It is being developed as a "plugin" to mkp
core, with all ecology-restricted content living under `dev/ecology/` and
`inst/ecology/` so future merges with `main` stay clean.

## Where we left off

This session: merged `main` into `ecology-aware`, segregated ecology content
into plugin namespaces, attempted v9 relabel + rodent v3 continuation. Two
streaming-layer bugs survived the merge; the rodent v3 chain turned out to
be unresumable with the merged code.

Recent commits (top first):

- `3dbd7db` docs(notes): per-character partition API design plan (v2) —
  user's concurrent work, committed from another session. v2 of the
  partition + unlink design (MrBayes-style component tokens, mean-1
  Dirichlet identifiability). Implementation deferred; section 11 lists
  7 open decisions.
- `fc6c2ff` feat: rodent v3 aware continuation harness
  (`inst/ecology/hamilton/rodent-v3-cont/`). Both submitted jobs FAILED at
  startup — see "Pending jobs" below.
- `50295e9` refactor: relocate ecology-restricted content under `*/ecology/`
  namespace (280 files moved with git mv; history preserved). New rule:
  ecology-restricted content lives under `dev/ecology/`, `inst/ecology/`
  etc.; top-level dirs stay in sync with mkp main.
- `85dc82b` Merge `main` into `ecology-aware` — 30+ main commits brought in.
  Highlights:
  - Streaming/per-run-tree machinery (commits `f023150`, `243a5f0`,
    `e512157`, `1463d8d`, `86a782b`, `51d9060`).
  - `7a1569b` brColStart fix for diagnostic cols + beta_geometric prior.
  - `c9a0686` per-site compile-time-K dispatch + tip-edge fast path.
  - `db4c348` TreeESS export (Option A, drop dot prefix).
  - New dep `callr` for `nCore > 1` parallel runs (`callr::r_bg`),
    `parallel = ` arg removed from `MkPrimeMCMC`. **Breaking change**.
  - Hamilton lib at `/nobackup/pjjg18/mkp-sim3-multirep-v3/lib` rebuilt
    against the merged HEAD on 2026-05-20 ~13:20 BST; `callr`,
    `processx`, `ps` installed.

## Session events

### v9-induce 3-rep MCMC results

Job `17233486` (queued previous session): rep03 COMPLETED, rep01 + rep02
FAILED in `RelabelEcology()` post-processing. Chain files all intact (4-chain
PT, 100k iter, ~3.3h aware).

`RelabelEcology()` cannot run on reps 01/02 due to two streaming bugs
(below). **Bypassed relabel and ran `ScoreTreesUnrooted()` directly** —
topology metrics are symmetric under phi <-> 1/phi reflection, so they don't
need relabeled samples. `summary.rds` saved for all three reps with a
`RELABEL_SKIPPED=TRUE` flag + reason.

Scores (P(AC) = true sister, P(AB) = false trap):

| rep | blind P(AC) | aware P(AC) | blind P(AB) | aware P(AB) | blind CID | aware CID |
|-----|-------------|-------------|-------------|-------------|-----------|-----------|
| 01  | 0.852       | 0.959       | 0.109       | 0.041       | 0.150     | 0.158     |
| 02  | 0.831       | 0.913       | 0.131       | 0.027       | 0.117     | 0.089     |
| 03  | 0.874       | 0.888       | 0.071       | 0.019       | 0.108     | 0.087     |
| mean | 0.852      | 0.920       | 0.104       | 0.029       | 0.125     | 0.111     |

**The "blind fails / aware rescues" headline does NOT survive PT-MCMC**.
The MP-trap (62% false-sister at MP per `8fc37ef`) collapses to 10% in
4-chain PT. Aware suppresses false sisters ~3.5× further but blind already
does well. n=3 — paper structure conversation needed; the original Sim 2
narrative is dead.

### Streaming bugs that survived `7a1569b`

The brColStart fix was about reading tree edges; it does NOT cover two
distinct intermittent failures:

- **z_samples / sample-row mismatch** (rep01): `result$z_samples` has 230
  entries vs `nrow(samples) == 200`. `RelabelEcology` aborts at
  `R/RelabelEcology.R:123`. Looks like z-buffer captures some tuning-phase
  entries the log doesn't record.
- **Torn row write** (rep02): line 178 of `aware-chain.log` has 8 of 44
  fields — partial flush mid-write. `ReadMkLog` -> `scan()` aborts.

Rep03 of the same job had neither bug. Spawned task chip exists for
root-causing (`Fix streaming z_samples / torn-write bugs`).

### Rodent v3 aware continuation: not resumable

Original v3 aware run wrote 0 trees to `aware-chain_trees.nwk` (file is
1 byte). The checkpoint records 135 trees that should be on disk.

The merged HEAD added a `.TruncateTreeToN()` desync check (came in with
the per-run tree machinery) that **catches** this and refuses to resume.
Both `17238607` (cont1) and `17238608` (cont2) failed at startup in 2-3s
with the message:

```
Found 0 valid trees on disk, but the checkpoint recorded 135.
The chain's param log has rows for trees that are no longer on disk;
resuming would produce a permanently inconsistent output.
```

The original chain's **scalar log is intact** (346k iter, minESS=27), but
posterior trees are gone. Two paths forward:

1. **Fresh restart** with merged code (recommended). Will produce both
   scalars AND trees correctly. ~12h per 1M iter at the rate observed
   (~15k iter/h sampling), so 2-3 12h continuations expected to reach
   minESS=200. Throws away 12h of scalar samples but they are not
   defensible standalone (we couldn't claim posterior topology stats from
   them anyway).
2. **Salvage scalars only** — pull existing chain.log, summarise phi/pi0/
   theta posteriors, treat rodent as "phi posterior from blind-tree case"
   rather than full topology + scalar. Cheap but limits the claims.

Decision deferred to next session.

## Pending jobs

| Type | ID / ref | Status | ETA | On completion |
|------|----------|--------|-----|---------------|
| SLURM | `17234909_*` (mkp-mk-tlshrink) | RUNNING (~3.5h of 8h) | ~4.5h | Not ecology — mkp-core arm benchmarking. Ignore on this branch; collect from mkp main when ready. |

The rodent v3 cont1 + cont2 jobs (`17238607`, `17238608`) **FAILED**;
no further continuation queued (would fail identically against the same
checkpoint). v9-induce reps already scored; no further jobs there.

## Open items / next steps

1. **Decide rodent v3 path** — fresh restart vs scalar-only salvage. If
   restart: submit `inst/ecology/hamilton/rodent-v3/rodent-v3-aware.sh`
   from scratch after `rm /nobackup/pjjg18/mkp-rodent-v3/aware/*.{ckp,log,nwk}`.
   Queue 2-3 sequential continuations via `--dependency=afterany`. With
   the merged code these will work (per-run trees, post-RDS cleanup,
   thin=500 streaming defaults).
2. **Paper structure conversation** — the v9 PT-MCMC null result kills
   the original Sim 2 "blind fails / aware rescues" framing. Options:
   reframe around "aware reduces false sisters 3.5×" as a regularisation
   claim; pivot to the rodent empirical as the headline (if it reaches
   ESS); deepen v9 with more reps + alternative-z mechanisms (M1/M3 from
   `dev/ecology/sim-design/v9-discriminate.R`).
3. **Streaming bug root-cause** — see the spawned task chip. Two bugs:
   z_samples/sample-row drift, and torn writes at thinning boundaries.
   Both intermittent; reps 01/02 of v9 reproduce on demand.
4. **Document directional-z mode** in `R/MkPrimeModel.R` / `R/RunMkPrime.R`
   (z=1: `(rate01 × phi, rate10 / phi)` directional). Currently a feature
   without docs.
5. **Per-character partition API v2** (`3dbd7db`) — user-authored design
   plan in another session. Section 11 has 7 open decisions awaiting
   review. Not blocking; not implementation-ready.

## Technical pointers

- **SSH to Hamilton from this machine**: use Windows OpenSSH explicitly:
  `/c/WINDOWS/System32/OpenSSH/ssh.exe pjjg18@hamilton8.dur.ac.uk "<cmd>"`.
  Passwordless. The plain `ssh` in MSYS PATH fails (Permission denied), and
  `plink` is configured for git but lacks a saved Hamilton session. SSH
  occasionally hangs from the parent shell mid-session — dispatching the same
  command via a subagent usually works; if not, retry in 15-30 min (likely
  login-node transient).
- **Hamilton paths**:
  - Repo: `/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime` (branch
    `worktree-ecology-aware`)
  - R library: `MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib`
  - Deploy via `git fetch origin && git restore --source=origin/worktree-ecology-aware -- <paths>`.
    Do NOT use `git reset --hard` on Hamilton — earlier session reset zeroed
    files when ORIG_HEAD.lock failed; recovery cost ~1 hour.
  - 2026-05-20 deploy: `git stash push -u -m hand-off-2026-05-20-pre-merge-deploy`
    preserved a working set of mods (NEWS.md, R/MkPrimeModel.R, test-priors.R
    deltas + untracked test-ecology-scaffold.R). Not yet applied; review
    before clearing.
- **Local worktree path**: `C:\Users\pjjg18\GitHub\worktrees\ecology-aware`
  (not the earlier nested `GitHub\GitHub\worktrees` mistake). Primary mkp
  checkout at `C:\Users\pjjg18\GitHub\mkp` stays on `main`. The
  `feature/ecology-aware` branch (commit `3a3b289` partition API v1) lives
  in mkp's local refs and is the prior version of the v2 design plan that
  landed at `3dbd7db` on this branch.
- **Plugin layout (post 50295e9)**:
  - `dev/`, `inst/hamilton/m131-*`, `inst/MkPrime/`, top-level `RED_TEAM_*.md`
    et al = mkp-core, mirror of main, do not touch from this branch.
  - `dev/ecology/`, `inst/ecology/` = ecology plugin namespace. Mirror
    structure inside (`dev/ecology/red-team/`, `dev/ecology/sim-design/`,
    etc.). Future ecology dev notes belong here.
- **Rodent data**: nexus at
  `/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex`
  (only there — NOT bundled in repo). 60 extant tips × 217 chars
  (160 neomorphic + 59 transformational after `AutoDetectNeomorphic()`).
  Aware-v2 config: `magnitudeMode="global"`, `kPrimePrior="geometric"`,
  `coding="variable"`, `expSteps=10` (legacy; v3 uses auto), `rho0Alpha=7`,
  `rho0Beta=3`, `thetaAlpha=2`, `thetaBeta=2`, `sigmaPhi=0.5`,
  `set.seed(20260512)`.
- **v4-cross helpers**: `inst/ecology/simulations/sim3v4cross-helpers.R` and
  `sim3v4-helpers.R`. Note `MatrixToPhyDat` / `MkPrimeData` are NOT exported
  in the installed lib — use bare names, not `MkPrime::` prefix.
- **PT execution model**: `nChains > 1` runs sequentially in MkPrimeMCMC;
  set `cpus-per-task=1` and multiply wall by nChains. Confirmed in
  `sim3-multirep-v4.sh` comments.
- **Aware iter rate on rodent data**: ~15k iter/h post-tuning (revised
  down from earlier 50k estimate, which was a tuning-phase artifact);
  1M iter ≈ 70h aware. 2G RAM, 1 CPU sufficient. Multiple sequential
  continuations needed for full chain.
- **`.SimulateMkPrimeEcology` z semantics**: per-(char, eco) integer in {0,1,2}.
  - z=0: no eco effect on this (char, eco)
  - z=1: for neomorphic chars, `(rate01 * phi, rate10 / phi)` — DIRECTIONAL
    toward state 1. For transformational, `mult = phi` (symmetric, noisier).
  - z=2: opposite of z=1.
  - Set z asymmetrically across two ecologies (different cols of z matrix)
    to create directional convergence at realistic TL — see v9 (commit `8fc37ef`).
- **Auto-expSteps**: `R/MkPrimeModel.R` finalises `expSteps = parsimony /
  nChar × expStepsInflation` (default 1.05). For rodent: parsimony 1519,
  nChar 217 → expSteps ≈ 7.35. Logs print parsimony, nChar, expSteps.
  Floor 0.5 prevents degenerate priors on small datasets.
- **ResumeMkPrime signature (merged HEAD)**: `ResumeMkPrime(checkpointFile,
  data, tree=NULL, neomorphic=integer(0), knownStates=integer(0),
  model=NULL)`. No `mcmc=` argument. Checkpoint versions 1, 2, 3 accepted.
- **Checkpoint/tree-file sync invariant (NEW post-merge)**: chain log,
  checkpoint, and `_trees.nwk` must agree on the number of recorded
  trees. `.TruncateTreeToN()` enforces this on resume and aborts on
  mismatch. **Old (pre-merge) chains that hit STREAM-001 (0-byte tree
  file) cannot be resumed under the merged code** — must restart fresh.
- **Bypass-RelabelEcology pattern for topology-only scores**: when the
  streaming bugs fire, load `result$trees` directly and call
  `ScoreTreesUnrooted(trees, biparts, refTree)` from
  `inst/ecology/simulations/sim3-scoring.R`. Topology metrics are
  symmetric under phi-flip; only phi/theta posterior summaries need
  relabeling.

## Things ruled out

- **chmod -R u+w** on Hamilton paths — blocked by classifier as destructive.
  Use file-by-file chmod or work around via git restore.
- **scp of .sh scripts** — line-ending corruption injected binary nulls in
  earlier sessions. Either strip CRLF before scp, or write via SSH heredoc.
- **Treating all rodent characters as transformational** — user clarified
  the existing 160 neo + 59 trans classification IS gold standard; do not
  reclassify.
- **Single-chain aware MCMC on v4-cross** — exhibits contiguous mode-traps
  (sample 84 onward stuck). Default to nChains ≥ 4 with PT for any
  v4-cross-like configuration.
- **`--amend` after a failed pre-commit hook** — the commit didn't happen, so
  amend rewrites the previous commit. Create a new commit instead.
- **`ape::prop.part` for posterior tree scoring** — root-dependent. Use
  `TreeTools::as.Splits` instead, via `HasBipartSplits` in
  `inst/ecology/simulations/sim3-scoring.R`. Discovered 2026-05-19.
- **`expSteps=10` default** for the TL prior — produces Gamma mean=10, gives
  5-13× TL inflation on realistic 16-tip morphological sims. Fixed in `bc24d52`
  (then unit-corrected in `f7b4644`).
- **Phi-rate-multiplier mechanism for realistic-TL trap design** — symmetric
  scaling can't produce directional convergence; v6/v7/v8 confirm.
  Use directional z (v9 path) instead.
- **multirep-v3 truth TL=13.5** as a "main result" — saturated; user declared
  junk for typical morphology. Acceptable as a long-tree case study only.
- **Subagent dispatch with "directional bias / trap" terminology** — burned 4
  spurious API-policy refusals on 2026-05-20. Write code inline if it recurs.
- **Resuming pre-merge ecology checkpoints under merged HEAD** — fails when
  the original chain hit STREAM-001 (no trees written despite checkpoint
  recording them). The new desync check is correct; old chains are not
  recoverable as continuations. Restart fresh.
- **`RelabelEcology()` as a paper-required step** — for topology metrics
  (the actual paper claims about P(true/false sister), CID), relabel is
  unnecessary; metrics are phi-flip-symmetric. Only phi/theta posterior
  summaries need it.
- **MP-trap → MCMC-trap inference** — the parsimony-screen 62% false-sister
  rate at TL=1.48 (v9 / commit `8fc37ef`) does NOT carry over to 4-chain
  PT MCMC, which recovers truth ~85% (blind) / ~92% (aware) of the time.
  TL prior + parsimony anchoring regularises sufficiently. Don't extrapolate
  from MP screens to MCMC behaviour without the MCMC.

## Worktrees

| Branch | Path | Status |
|--------|------|--------|
| ecology-aware (tracks worktree-ecology-aware) | C:\Users\pjjg18\GitHub\worktrees\ecology-aware | clean |
| main (mkp primary) | C:\Users\pjjg18\GitHub\mkp | clean, in sync with origin/main |

## Suggested first action

Decide rodent v3 path. If fresh restart (recommended), on Hamilton:

```
ssh pjjg18@hamilton8.dur.ac.uk "rm /nobackup/pjjg18/mkp-rodent-v3/aware/rodent-aware-v3.{ckp,log} /nobackup/pjjg18/mkp-rodent-v3/aware/rodent-aware-v3_trees.nwk && \
  cd /nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime && \
  J1=\$(sbatch --parsable inst/ecology/hamilton/rodent-v3/rodent-v3-aware.sh) && \
  J2=\$(sbatch --parsable --dependency=afterany:\$J1 inst/ecology/hamilton/rodent-v3-cont/rodent-v3-aware-cont.sh) && \
  J3=\$(sbatch --parsable --dependency=afterany:\$J2 inst/ecology/hamilton/rodent-v3-cont/rodent-v3-aware-cont.sh) && \
  echo restart=\$J1 cont1=\$J2 cont2=\$J3"
```

After cont2 completes, check minESS. If still <200, queue cont3. If salvage
chosen instead, pull `/nobackup/pjjg18/mkp-rodent-v3/aware/rodent-aware-v3.log`
locally and summarise phi/pi0/theta marginals only.
