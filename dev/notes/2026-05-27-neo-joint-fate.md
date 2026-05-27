# neo_joint fate (A2 verdict)

**Decision: drop `neo_joint` (case 18).**

## Measurements

Source: `mkprime_635_by_nt_9v.rds` pooled posterior, n = 850 post-burnin
samples from a 29-min smoke run on Hamilton (jobs 17302689/17302690,
2026-05-27). pid 635 is a mixed dataset (neo + trans), `by_nt_9v` model.

| Correlation pair | ρ |
|---|---|
| cor(log T,         log rate_neo)  | **+0.072** |
| cor(log T,         log rate_loss) | **−0.140** |
| cor(log rate_neo,  log rate_loss) | **−0.124** |

All three are weak (|ρ| < 0.15). This is outcome 3 from the A2 decision
tree in `dev/notes/2026-05-27-rate-neo-ridge-and-joint-moves.md`.

## Why the ridge is weak here

The theoretical post-fix ridge `T ∝ (1+r)` is strongest when the posterior
is unconstrained on the trans side. pid 635 (`by_nt_9v`) has 36 characters
(~22 trans, ~14 neo per the matrix metadata) and a 9-state neomorphic model
— enough data to anchor T independently of r, giving a weak effective ridge.

The (T, r) ridge from Issue 1 is real in principle but is beaten into
flatness by the data likelihood for these mid-sized matrices. For very small
matrices or extreme neo/trans ratios it may be stronger; SBC (A3) will check.

## Decision

**Drop `neo_joint` (move type 18, case 18 in `src/mcmc.cpp`).**

Rationale:
- The (rate_loss, rate_neo) correlation it was designed to ride is −0.12:
  no better than statistical noise at this sample size.
- `joint_tl_rn` (A1) covers the (T, r) direction — also weak here but
  still the right direction to have a proposal for.
- Keeping neo_joint adds weight scheduling overhead for a move that is
  neither helping mixing nor addressing a real posterior ridge.

### What to remove

| Location | What |
|---|---|
| `src/mcmc.cpp` | case 18 (`neo_joint`) |
| `R/RunMkPrime.R` | `neo_joint` entry in `.BuildMoves` |
| `R/RunMkPrime.R` | `neo_joint` in `.moveTypes` |
| `R/MkPrimeMCMC.R` | `neo_joint` in valid move names list |
| `R/MkPrimeMCMC.R` | `scale_neo_joint` default tuning, if present |
| test files | update any `multi_dim_moves` or allowlist tests that reference `neo_joint` |

### Before removing

Check `squeue` / any long runs on Hamilton that have a checkpoint using
`neo_joint` — the resume path will error if the move name is gone. The
smoke run cells (950, 635) finished, so no open checkpoints use it.

## Open: A3 (SBC on mixed data)

The rb-equivalence scalar pass for pid 635 (all 4 params, rhat < 1.025,
ESS > 128) is encouraging but is not a substitute for SBC. A3 remains open.
