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

<!-- Example row:
| 2026-03-29 | SLURM | 16622483 | main | Benchmark: Mk' vs RevBayes on Sun2018 (5 seeds) | `scp hamilton:scratch/mkp_bench/*.csv benchmark/` | E |
-->
