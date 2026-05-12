# MkPrime — performance notes

## Overall bottleneck (S-PROF round 2, 2026-03-28)

C++ Felsenstein pruning is ~90% of wall time. OPP-1–6 optimizations achieved
1.80× cumulative speedup. Diminishing returns on further pruning optimization
without algorithmic restructuring.

## Gibbs kPrime sweep (M-155, 2026-03-31; commit `224bda8`)

**Pre-optimization** (per-character individual traversals):

| Move | Cost | Coverage |
|------|------|----------|
| Gibbs kPrime sweep (moveType 25) | 340 ms | All 99 trans chars |
| 99 × int_walk kPrime (moveType 7) | 33 ms  | All 99 trans chars |
| Block kPrime shift (moveType 26)  | 0.66 ms | All 99 (collective mode) |

Root cause of 10× overhead: `single_char_loglik_jc()` did per-character,
per-candidate-k' individual tree traversals with per-call heap allocation.

**Post-optimization:** Batched partition-level pruning with progressive
early termination.

| Config | Chars | nCat | Old (est.) | New | Speedup |
|--------|-------|------|------------|-----|---------|
| all-trans | 225 | 1 | ~70 ms  | 35.7 ms | ~2×   |
| all-trans | 225 | 6 | ~770 ms | 214 ms  | ~3.6× |
| with-neo  | 31  | 1 | ~70 ms  | 13.0 ms | ~5.4× |
| with-neo  | 31  | 6 | ~107 ms | 71 ms   | ~1.5× |

nCat scaling remains roughly linear. Speedup is larger with more characters
(better amortization of batched traversals).

## VTune hotspot profile (S-PROF round 4, 2026-03-31)

Sun2018 (54 taxa, 225 chars, all transformational, nCat=6 ACRV). 15k
iterations, ~90s CPU time. VTune 2025.10, user-mode sampling.

| Function | Source | CPU Time | % |
|----------|--------|----------|---|
| `pruning_jc_acrv_persite` | mcmc_likelihood.cpp | 53.0s  | 59.1% |
| `_expl_internal` (exp)    | compiler runtime    | 10.5s  | 11.7% |
| `constant_site_prob_jc`   | ascertainment.cpp   |  5.9s  |  6.6% |
| `pruning_jc_acrv_flat`    | mcmc_likelihood.cpp |  3.5s  |  3.9% |
| `jc_transition`           | gibbs_partial_cl.h  |  1.3s  |  1.4% |
| `std::vector` copies      | stl_vector.h        | ~1.3s  | ~1.5% |
| Rcpp bounds checks        | traits.h            | ~1.0s  | ~1.1% |
| Everything else           |                     | ~13s   | ~14.5% |

### Key findings

1. **Gibbs kPrime sweep** (`persite` + ascertainment) dominates at ~65–70% of
   MkPrime CPU. Regular MH proposals (`flat`) are only 3.9%.
2. **`exp()` calls at 11.7%** — one per edge × rate category per traversal
   (`mcmc_likelihood.cpp` line ~370). Already amortized across characters;
   cannot factor across candidate k' values (eigenvalue depends on k). Fast
   approximate `exp()` could save ~5–8% of total.
3. **Ascertainment at 6.6%** — `constant_site_prob_jc` does separate tree
   traversals. Could batch like the main pruning.
4. **Vector copies at 1.5%** — heap allocation in the hot path.
5. **Rcpp bounds checks at 1.1%** — `check_index` on `operator[]`.

## VTune workflow specifics

- Override `DLLFLAGS` via `MAKEFLAGS` env var (not `src/Makevars.win`).
- Add `-g -fno-omit-frame-pointer` to `PKG_CXXFLAGS` in `src/Makevars.win`.
- **Remove profiling flags after collection.** `src/Makevars.win` must never
  be committed — it is in `.gitignore`.
- VTune installed at `C:/Program Files (x86)/Intel/oneAPI/vtune/latest/bin64/`.
- This PC (Intel i7-10700, 10th gen) supports hardware sampling.
- Use the `r-package-profiling` skill (`Skill r-package-profiling`) for the
  full driver-script workflow.
