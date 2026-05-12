# M-125: Block Dirichlet Branch-Length Proposal

## Problem

The current branch-length proposal (`beta_simplex`, case 4) updates only **2 edges at a time** from 44 relative branch lengths. On Vinther2008 (23 taxa, fixed topology, 20k iter), branch-length ESS is 2–27 (median ~10) while `tree_length` and `rate_log_sd` reach ESS 115–177. The log-likelihood ESS (~5) is completely bottlenecked by branch lengths.

## Design: `mvDirichletSimplex`-style Proposal

RevBayes uses two complementary simplex moves together:
- **`mvBetaSimplex`** — update 1 element, rescale others (target accept: 0.44)
- **`mvDirichletSimplex`** — update K elements simultaneously from a Dirichlet distribution, scale the rest (target accept: 0.234)

Their `mvDirichletSimplex` has a `numCats` parameter controlling how many elements to update per proposal. This is the missing piece in MkPrime.

### Algorithm

Given relative branch lengths `r = (r₁, ..., rₙ)` with `sum(r) = nEdge`:

1. **Choose K indices** uniformly at random (K = `nCats`, configurable)
2. **Form reduced simplex**: `x = (r[i₁]/S, ..., r[iₖ]/S, (S_rest)/S)` where `S = sum(r)`, `S_rest = S - sum(selected)`
3. **Propose** `z ~ Dirichlet(α·x₁ + κ, ..., α·xₖ + κ, α·x_{K+1} + κ)` where α is the concentration (tuning parameter) and κ is a small offset (default 0)
4. **Update selected**: `r[iⱼ] = z[j] · S` for j = 1..K
5. **Scale unselected** by factor `f = z[K+1] / x[K+1]`: `r[i] *= f` for all unselected i
6. **Hastings ratio**: `log H = log Dir(x | α·z + κ) - log Dir(z | α·x + κ) + (n - K - 1)·log(f)`

The `(n - K - 1)·log(f)` term is the Jacobian for the scaling transformation on unselected elements.

When K = n (all elements), steps 5–6 simplify: no scaling, no Jacobian.

**Cost: 1 full likelihood evaluation per proposal** (same as current beta_simplex).

### Default K

For branch lengths, a sensible default is `nCats = min(nEdge, 10)`. This updates ~10 edges simultaneously — enough to make coordinated moves without the curse of dimensionality that a full 44-dimensional Dirichlet would face.

The adaptive scheduler can discover the right weight balance between the existing 2-edge beta_simplex (cheap local moves) and the new K-edge Dirichlet (more expensive but more global moves).

## Implementation Plan

### Step 1: C++ proposal function (`src/proposals.cpp`)

New function `dirichlet_simplex_impl()`:

```cpp
bool dirichlet_simplex_impl(
    NumericVector& x,           // relative branch lengths (modified in-place)
    int nCats,                  // K: number of elements to update
    double alpha,               // concentration parameter (tuning)
    double kappa,               // offset (default 0)
    double& logHastings,        // output: log Hastings ratio
    IntegerVector& rollbackIdx, // output: indices modified (for rollback)
    NumericVector& rollbackVal  // output: old values (for rollback)
);
```

Implementation details:
- Use `rgamma` to sample from Dirichlet (draw K+1 Gamma variates, normalize)
- Clamp minimum element values to 1e-8 to avoid numerical issues
- Store rollback info: all n old values (since unselected are also scaled)
- Actually, simpler: save the full `x` vector snapshot before modification for O(n) rollback

### Step 2: New move type in C++ engine (`src/mcmc.cpp`)

- New case `23` (dirichlet_branch) in the `do_move_impl` switch
- Pre-move: snapshot `relBrLengths` for rollback
- Call `dirichlet_simplex_impl()`
- On rejection: restore from snapshot
- Likelihood evaluation uses the updated `relBrLengths` × `treeLength`

### Step 3: R-side move construction (`R/RunMkPrime.R`)

In `.BuildMoves()`:
```r
list(name = "dirichlet_branch", type = "dirichlet_simplex",
     target = "rel_br_lengths",
     weight = max(1, nEdge / 4),
     dim = as.integer(min(nEdge, 10L)),  # K edges per proposal
     nCats = as.integer(min(nEdge, 10L)))
```

In `.kMoveTypes`:
```r
dirichlet_branch = 23L
```

### Step 4: Tuning integration

- New tuning parameter `dirichlet_alpha` in the tuning object
- Initial value: `10.0` (conservative — high concentration = small perturbation)
- Target acceptance rate: `0.234` (multi-dimensional proposal)
- Adaptation: same multiplicative scheme as other moves, but divide by adjustment (like beta_simplex — higher alpha = more conservative)
- Cap: `max(1.0, alpha)` floor, `alpha <= 1000` ceiling

### Step 5: Rollback in C++

Two options for rollback efficiency:

**Option A — Full vector snapshot (simple, O(n)):**
Save entire `relBrLengths` before the move. On rejection, copy back. Cost: ~44 doubles = 352 bytes. Trivial.

**Option B — Selected-index rollback (complex, O(K)):**
Save only modified indices + old values + the scaling factor for unselected. More complex logic for marginal memory savings.

**Decision: Option A.** The cost is negligible compared to a likelihood evaluation, and the code is much simpler. Store the snapshot in a pre-allocated vector in `McmcState`.

### Step 6: Tests

1. **Proposal correctness** (`test-proposals.R`):
   - Simplex constraint maintained (sum preserved to machine epsilon)
   - All elements positive after proposal
   - Hastings ratio is finite
   - For K=n, no unselected elements remain (pure Dirichlet)
   - For K=2, behavior similar to beta_simplex

2. **Detailed balance** (`test-proposals.R`):
   - Round-trip test: propose forward, swap old/new, propose reverse
   - Check `log H_fwd + log H_rev ≈ 0` (not exactly, but statistically)

3. **Mixing improvement** (slow test):
   - Compare branch-length ESS with and without dirichlet_branch on Vinther2008
   - Expect ≥ 3× ESS improvement for branch lengths

### Step 7: Integration with existing moves

The new move **complements** (does not replace) the existing beta_simplex:
- `beta_simplex`: fast, local, 2-edge pairwise redistribution
- `dirichlet_branch`: moderate, K-edge coordinated updates

The adaptive scheduler will discover the optimal weight split. Initial weights:
- `branch_lengths` (beta_simplex): `nEdge / 3` (unchanged)
- `dirichlet_branch`: `nEdge / 4` (new)

## Files to Modify

| File | Change |
|------|--------|
| `src/proposals.cpp` | Add `dirichlet_simplex_impl()` |
| `src/mcmc.cpp` | Add case 23, rollback logic, snapshot buffer in McmcState |
| `R/RunMkPrime.R` | `.BuildMoves()` entry, `.kMoveTypes` entry, tuning in `.AdaptTuning()` |
| `tests/testthat/test-proposals.R` | Proposal correctness tests |
| `tests/testthat/test-m125-dirichlet-branch.R` | Integration tests |

## Risk Assessment

- **Low risk**: This is an additive change — new move alongside existing ones. If the Dirichlet proposal has low acceptance, the scheduler will downweight it automatically.
- **Numerical concern**: Dirichlet proposals can produce very small values for some elements. Clamping to 1e-8 and checking for `sum > 0` prevents degenerate proposals.
- **Tuning sensitivity**: The initial alpha value matters. Starting too low (bold proposals) → near-zero acceptance. Starting too high (timid proposals) → no mixing benefit. 10.0 is a safe conservative start; warmup adaptation will find the right level.
