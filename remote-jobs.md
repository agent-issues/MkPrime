# Pending Remote Jobs

Track asynchronous jobs (Hamilton SLURM, long-running GHA, etc.) that
produce results an agent needs to retrieve later.

## How this works

- **Add a row** when you submit a remote job whose results won't be
  consumed in the same conversation turn.
- **Delete the row** once results have been retrieved and acted on
  (committed to repo, written up in coordination.md, etc.).
- Agents check this file at `/assign` time, after triaging `u.*`
  files but before claiming from `to-do.md`. If a retrievable
  job is listed, retrieving and processing it takes priority.

## Jobs

| Submitted | Type | Job ID | Branch | Description | Retrieve how | Owner |
|-----------|------|--------|--------|-------------|-------------|-------|
| 2026-03-31 | SLURM | 16632292 | main | M-131: Warmup stabilisation validation (8 datasets × 4 seeds, 200k iter each). Results in `/nobackup/pjjg18/m131/results/`. | `ssh::scp_download(session, "/nobackup/pjjg18/m131/results/", "inst/hamilton/")` | B |
| 2026-05-12 | SLURM | 17140607 | main | mkp-study eg-arm (260 tasks = 26 trees × 10 reps, empirical_geometric prior, heat=0.1, maxTime=10h). Symptom: median minESS=27, max 73, 0/26 reach 100; logP swings 30–50. Results in `/nobackup/pjjg18/mkp-study/results/` as `mkp_eg_t##_r##.rds` + `t##_r##/mkp_eg_run.log`. | `source("inst/hamilton/mkp-study-retrieve.R")` (interactive — prompts for ssh passphrase) | — |

<!-- Example row:
| 2026-03-29 | SLURM | 16622483 | main | Benchmark: Mk' vs RevBayes on Sun2018 (5 seeds) | `scp hamilton:scratch/mkp_bench/*.csv benchmark/` | E |
-->
