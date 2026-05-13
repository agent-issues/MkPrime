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
| ~~2026-03-31~~ | ~~SLURM~~ | ~~16632292~~ | ~~main~~ | ~~M-131: Warmup stabilisation validation (8 datasets × 4 seeds, 200k iter each).~~ **CLOSED 2026-05-13:** all 32 array tasks COMPLETED with ExitCode 0:0 but ran 4–11 s each; `/nobackup/pjjg18/m131/results/` is empty. The job exited cleanly without producing output — likely an early-exit guard in the R script. Nothing to retrieve. If M-131 validation is still wanted, resubmit (see [inst/hamilton/m131-warmup-validation.slurm](inst/hamilton/m131-warmup-validation.slurm)). | — | B |

<!-- Example row:
| 2026-03-29 | SLURM | 16622483 | main | Benchmark: Mk' vs RevBayes on Sun2018 (5 seeds) | `scp hamilton:scratch/mkp_bench/*.csv benchmark/` | E |
-->
