# Per-character partition API for `RunMkPrime` — design plan (v4)

**Status:** draft v4, awaiting implementation sign-off.
**Driver:** AutoPart Casali production pipeline (T0–T4 treatments).
**Branch:** `feature/partition-api`, off `main` at `db4c348`.
**Cross-references:** MrBayes `unlink` and `Ratepr`
(`NBISweden/MrBayes` `src/command.c`, `model.c`, `proposal.c`),
RevBayes partition tutorial, Casali et al. (2023)
*Syst. Biol.* `10.1093/sysbio/syac046`.

This document proposes an extension to `RunMkPrime` that lets a caller
supply (a) a per-character partition assignment, (b) a per-class model
type (neomorphic / transformational), and (c) a list of model components
to *unlink* across classes. It is deliberately a design conversation:
every section ends with the decisions a senior reviewer would want to
weigh in on.

**v4 changes from v3** (Martin's feedback, 2026-05-20 session 2):

- Mixed-type classes: **SUPPORTED** (§3.2 reversed — a class may contain both neomorphic
  and transformational characters; Q1 resolved).
- Token name corrected: `"ratemultiplier"` throughout (§3.3, §4.2, §1 table; Q4 resolved).
- §5.2 reparameterised: within a class the geometric mean of neomorphic forward (gain) and
  reverse (loss) rates equals `class_rate[c]`, matching the transformational rate. This
  decouples the mean-1 constraint from `eta_neo` and simplifies the formula (Q5 resolved).
- `knownStates` kept unchanged; Q2 resolved.
- Token vocabulary `c("shape", "ratemultiplier", "brlens")` for v1 confirmed; Q3 resolved.

**v3 changes from v2** (Martin's feedback, 2026-05-20):

- Plan migrated off `ecology-aware` onto its own feature branch.
  Ecology-aware is treated as a hypothetical future merge, not a design
  driver (§5.3 reframed; see end of §1).
- Partition flagged neomorphic or transformational at the **class
  level**, not per-character (§3.3). Adds `partitionType` arg.
- `rate_neo` stays as a **single global scalar** shared across all
  neomorphic characters, multiplied by `class_rate[c]` for class `c`
  (§5.2). Parameterisation is the one-free-parameter `eta_neo` form
  that satisfies char-weighted mean = 1 by construction.
- For Casali, all classes are transformational, so `rate_neo` is never
  exercised in production.
- Decision Q6 (numeric-tolerance equivalence in the §7(b) contract):
  accepted as proposed.

---

## 1. What the caller needs

The AutoPart pipeline ([SPEC §Partition file](../../../auto-part/dev/benchmarks/casali/SPEC.md))
needs `RunMkPrime` to accept five treatments. All Casali classes are
flagged transformational (`partitionType` defaults are sufficient; no
neomorphic chars in the Casali corpus).

| Treatment | `nClasses` | What is unlinked across classes |
|-----------|-----------:|----------------------------------|
| T0 (unpartitioned)        | 1     | — (nothing to unlink)                       |
| T1 (anatomical)           | K≈2   | `shape`, `ratemultiplier`                   |
| T2a (AutoPart, isolated)  | K+1   | `shape`, `ratemultiplier`                   |
| T2b (AutoPart, merged)    | K     | `shape`, `ratemultiplier`                   |
| T3 (AutoPart unlinked)    | K+1   | `shape`, `ratemultiplier`, `brlens`         |
| T4 (random control)       | K+1   | `shape`, `ratemultiplier`                   |

**Ecology-aware composability note** (informational only, **not a design
driver**): the `ecology-aware` branch adds a per-character `z` and
per-edge `phi` mixture model. If that work lands later, the multiplicative
composition is straightforward (§5.3) provided each rate-modifier factor
independently satisfies char-weighted mean = 1. We adopt that invariant
here, but make no other accommodation for `ecology-aware` and write no
code that depends on it. The branch may never merge; if it does, this
plan must remain valid as-is.

---

## 2. Current state of the engine

Survey findings (`R/RunMkPrime.R`, `R/partition.R`, `R/likelihood.R`,
`src/mcmc_state.h`, `src/mcmc_likelihood.cpp`):

- **`MkPrimeData$partitions`** is built by `.BuildPartitions()` and
  groups characters by `(type, kObs)`. The (type, kObs) grouping is
  load-bearing for the C++ pruning loop — JC dispatch needs a fixed `k`
  per group, and MkN vs JC are different code paths. The new user-class
  layer wraps around this; it does not replace it.
- **Per-character type today**: `mkd$type[i]` records each character's
  model. In the new design (v4), types remain per-character — set by
  the existing `neomorphic` and `knownStates` args, unchanged. Classes
  may mix types freely (§3.2).
- **MCMC state is one-of-each globally**: `tree_length`,
  `rel_br_lengths`, `rate_log_sd`, `rate_loss`, `rate_neo`. Each
  becomes a vector of length 1 or nClasses depending on `unlink`,
  except `rate_neo` and `rate_loss` which stay as global scalars
  (§5.2, §5.3).

---

## 3. Proposed `RunMkPrime` signature

```r
RunMkPrime <- function(data, tree = NULL,
                       neomorphic = integer(0),
                       knownStates = integer(0),
                       model = NULL,
                       mcmc = NULL,
                       fixTopology = FALSE,
                       overwrite = FALSE,
                       partition = NULL,                # NEW
                       unlink = character(0),           # NEW
                       ...)
```

`partitionType` is **removed** from the v4 design. Per-character model types are already
carried by `mkd$type[i]`, set via the existing `neomorphic` and `knownStates` args.
Supporting mixed-type classes makes a per-class type label redundant.

### 3.1 `partition`

`NULL` (default) → behaves exactly as today (legacy code path; see §7a).
Or `integer(nChar)` with values in `1:nClasses`, no NAs. Validated
against `mkd$nChar` after invariant-character drop.

### 3.2 Per-character type assignment (v4: `partitionType` removed)

Per-character model types (`neomorphic` / `transformational` / `known`)
remain per-character, not per-class. They are set via the existing
`neomorphic = integer(0)` and `knownStates` args, exactly as today.
Mixed-type classes are **supported** [Q1 resolved]: `eta_neo` (§5.2)
applies to neomorphic chars in a class; transformational chars in the
same class use `class_rate[c]` directly.

Backward-compat: the `neomorphic` and `knownStates` args are unchanged
and work whether or not `partition` is supplied.

**Known characters:** the legacy `knownStates` mechanism survives
unchanged [Q2 resolved]. Known-state characters can be members of any
user class; their per-character `k` is set from `knownStates`
independently of class membership.

### 3.3 `unlink`

Character vector of component tokens. Default `character(0)`:
everything linked. Token vocabulary borrowed from MrBayes
(`Help_Unlink`, `command.c` L12798):

| Token                | Effect                                                            | MrBayes analogue   |
|----------------------|-------------------------------------------------------------------|--------------------|
| `"shape"`            | per-class `rate_log_sd[c]` (ACRV Γ shape)                         | `shape`            |
| `"ratemultiplier"`   | per-class `class_rate[c]` multiplier (mean-1 Dirichlet)           | `ratemultiplier`   |
| `"brlens"`           | per-class branch lengths under shared topology (subParam idiom)   | `brlens`           |

Out of scope for v1: `pinvar`, `statefreq`, `tratio`, etc. — MkPrime
doesn't expose these as user-tunable.

**Matching rule:** case-insensitive **with partial-prefix matching**
and a `cli::cli_warn` when a prefix match resolves (`"shap" → "shape"`).
Ambiguous prefixes hard-error. Unknown tokens hard-error with an
`agrep`-driven "did you mean …?" suggestion.

`unlink` is silently coerced to `character(0)` with a `cli_alert_info`
when `partition = NULL` or `nClasses == 1` — lets the AutoPart
dispatcher pass `unlink` uniformly across treatments.

---

## 4. Internal data model

### 4.1 User-class layer atop existing grouping

`.BuildPartitions()` keeps its job, operating **within each user
class**. Each `PartInfo` (R-side list and C++ struct) gains a
`classIdx` field. With `partition = NULL`, all chars get `classIdx =
1L` and the layout is identical to today.

### 4.2 Per-class state fields

| State field             | Length when "linked" | Length when "unlinked" | Token                |
|-------------------------|----------------------|------------------------|----------------------|
| `class_rate_log_sd`     | 1                    | nClasses               | `"shape"`            |
| `class_rate`            | implicit ≡ 1         | nClasses (on simplex)  | `"ratemultiplier"`   |
| `class_rel_br_lengths`  | 1 (length-nEdge)     | nClasses (each length-nEdge) | `"brlens"` |
| `class_tree_length`     | 1                    | nClasses               | `"brlens"`           |

`rate_neo` and `rate_loss` stay as **global scalars** regardless of
`unlink` — they're shared across all neomorphic classes (§5.2). Not
listed in the unlink vocabulary in v1; if a user wants per-class
neomorphic asymmetry, that's a future extension.

The state-field layout is **uniform regardless of `unlink`**: only
vector lengths change. Move-list construction is uniform: every move
has weight = (per-class weight) × (length of that class-vector). A
linked-`shape` chain has one `rate_log_sd` move (length-1 vector); an
unlinked-`shape` chain has nClasses `rate_log_sd` moves.

### 4.3 Brlens unlink: MrBayes' subParam-of-topology idiom

When `"brlens"` is unlinked, the **single shared topology** has
nClasses brlens vectors hanging off it (MrBayes' BRLENS subParams of
the TOPOLOGY param; `mcmc.c` L11349). Each per-class brlens vector is
a free Dirichlet simplex × tree-length scalar.

Topology proposals (NNI / SPR / etc.) operate on the shared topology,
and the edge-set re-mapping is applied **to every class's brlens
vector consistently** in the same C++ call (MrBayes `proposal.c`
L836–921, L1438–1456). Our existing C++ topology proposals return
parent/child arrays for the new edge configuration; the R-side
wrapper iterates over `c = 1..nClasses` and remaps each
`class_rel_br_lengths[[c]]` against the new edge IDs.

---

## 5. Priors and identifiability

### 5.1 `class_rate` — char-weighted mean-1 Dirichlet

**Constraint:** `sum(nChar_c * class_rate[c]) / nChar = 1` (char-weighted
arithmetic mean), matching MrBayes (`Help_Prset`, `command.c` L11385–
L11422) and RevBayes (`part_rate := part_rate_mult * n_sites /
num_sites_subset[i]`). Keeps `tree_length` interpretable as the expected
number of changes per character.

**Parameterisation:** sample `(w_1, …, w_nClasses) ~ Dirichlet(α)` where
`w_c = class_rate[c] * nChar_c / nChar` is the share of expected total
changes contributed by class `c`. Recover `class_rate[c] = w_c * nChar
/ nChar_c`.

**Prior:** Dirichlet(α₁, …, α_nClasses) with default α = 1 (flat).
User-tunable via a new `MkPrimeModel(classRateConcentration = 1)` arg.

**Move:** Dirichlet-simplex proposal on `w`, reusing the existing
`beta_simplex` move infrastructure (currently used for
`rel_br_lengths`).

### 5.2 `rate_neo` → `eta_neo` reparameterisation [Q5 resolved]

`eta_neo` is a **single global scalar** controlling gain→loss vs
loss→gain asymmetry in neomorphic characters. The constraint is:

> Within any class `c`, the **geometric mean** of the neomorphic
> forward (gain) rate and reverse (loss) rate equals the
> transformational rate, i.e. `class_rate[c]`.

- Free parameter: `eta_neo > 0` ~ LogNormal(0, `rateNeoSdlog`),
  default `rateNeoSdlog = 1`.
- Within class `c`, for each neomorphic char:
  ```
  r_gain_c = class_rate[c] * sqrt(eta_neo)
  r_loss_c = class_rate[c] / sqrt(eta_neo)
  ```
  Geometric mean: `sqrt(r_gain_c * r_loss_c) = class_rate[c]` ✓.
  Ratio: `r_gain_c / r_loss_c = eta_neo` (same in every class).
- Transformational chars in class `c`: rate = `class_rate[c]` (unchanged).
- At `eta_neo = 1`: `r_gain_c = r_loss_c = class_rate[c]` (no
  asymmetry; neomorphic and transformational chars in the same class
  have the same rate).
- `eta_neo > 1`: gain dominates (acquisition more common than loss).

**Mean-1 constraint** is carried entirely by `class_rate` (§5.1),
independently of `eta_neo`:
```
sum_c (nChar_c * class_rate[c]) / nChar = 1
```
All characters (neomorphic and transformational) are weighted equally;
the geometric-mean formulation ensures this without coupling `eta_neo`
into the normalisation.

Edge cases: when `nNeo_c = 0` (pure-trans class), `r_gain_c` and
`r_loss_c` are undefined and unused. When `nTrans_c = 0` (pure-neo
class), the formulas apply directly.

**For Casali**: `nNeo = 0` globally (all transformational chars), so
`eta_neo` is never sampled — the production chains are unaffected by
this change.

**Behavioural change**: existing chains with non-default `rate_neo`
will produce a different posterior under the new parameterisation
because the prior shape on `eta_neo` is the same LogNormal but the
quantity scaling tip evolution is now derived (`r_neo`,
`r_trans`) rather than (`rate_neo`, 1). At `eta_neo = 1` the chains
agree exactly. Document in NEWS.md.

### 5.3 Composability with hypothetical `ecology-aware`

(Informational; not load-bearing for this work.)

The ecology-aware likelihood applies per-(character, ecology) `z` and
per-edge `phi` rate multipliers. If that work lands, the effective
per-(char, edge) rate becomes a product:

```
effective_rate(char, edge) =
  class_rate[ class(char) ] * (eta_neo-derived r_neo_c or r_trans_c)
  * phi(z(char, eco), edge_has_eco)
```

Provided each factor independently satisfies char-weighted mean = 1
across the chars it acts on, the product does too, and `tree_length`
stays interpretable. The ecology-aware model's `phi` already satisfies
this (it's centred at 1 under the `pi_0` no-effect mass). Nothing in
the partition API forces accommodation of `ecology-aware`; this is a
forward-compatibility note only.

### 5.4 Other priors

| Parameter                | Prior                                                                | Notes |
|--------------------------|----------------------------------------------------------------------|-------|
| `class_rate_log_sd[c]`   | `Gamma(rateLogSdShape, rateLogSdRate)`, i.i.d. per class             | hierarchical pooling deferred |
| `class_tree_length[c]`   | `Gamma(treeLengthShape, treeLengthRate)`, i.i.d. per class           | unlinked-brlens only |
| `class_rel_br_lengths[[c]]` | flat Dirichlet on the edge-length simplex                         | unlinked-brlens only |

---

## 6. MCMC moves

Move construction iterates over each component:

| Component            | Linked                                          | Unlinked                                           |
|----------------------|-------------------------------------------------|----------------------------------------------------|
| `shape`              | one `rate_log_sd` move (current)                | nClasses `rate_log_sd` moves                       |
| `ratemultiplier`     | none (vector is implicit ≡ 1)                   | one Dirichlet-simplex move on `w`                  |
| `brlens`             | one `tree_length` + one `rel_br_lengths` move (current) | nClasses copies of each, plus topology moves that update all per-class simplexes |

Topology moves change the edge set across all classes simultaneously
(MrBayes idiom).

### C++ surface change

```cpp
double cpp_log_likelihood(
    const McmcData& data,
    IntegerVector parent, IntegerVector child,
    NumericMatrix edgeLen,           // nEdge × nClasses (1 column when brlens linked)
    const IntegerVector& kPrime,
    double rateLoss,                 // scalar (still global in v1)
    NumericVector rateLogSd,         // length nClasses (1 when shape linked)
    NumericVector classRate,         // length nClasses (1 when rate linked)
    double etaNeo,                   // scalar (still global)
    double betaScale,
    ClWorkspace* ws);
```

For each `PartInfo`, look up `classIdx`, pick the right column of
`edgeLen`, the right `rateLogSd`, and the right `classRate`; derive
`r_neo_c`/`r_trans_c` from `etaNeo` and per-class counts; compute
per-partition likelihood as today. When all `unlink`-vectors are
length 1 (default behaviour), this reduces to the current single-value
path.

---

## 7. Backward compatibility contract

**(a) Bit-identity for `partition = NULL`.** When `partition` is
`NULL`, `RunMkPrime` routes through the **unchanged** legacy code
path. New per-class code is a sibling, not a replacement. Hard
guarantee; regression test checks `identical()` against a stored
reference sample matrix on `Lobo.phy` × 1000 iter × fixed seed.

**(b) Numeric-tolerance equivalence for `partition = rep(1L, nChar),
unlink = character(0)`.** Routes through the new code path with all
class-vectors length-1. Likelihood and log-prior agree with HEAD to
~1e-10 absolute on the initial state. Sample matrix is not required
to match bit-for-bit (RNG ordering may differ). **Confirmed
acceptable by Martin.**

**(c) `rate_neo` → `eta_neo` reparameterisation (§5.2)** is a
deliberate behavioural change. Bit-identity guarantee (a) holds only
when the user keeps default `rate_neo = 1` (i.e. `eta_neo = 1`).
Document in NEWS.md.

---

## 8. Risks and open design questions

1. **`rate_neo` → `eta_neo` change of variable.** Necessary for
   consistency with the mean-1 rule across the entire model. Martin
   confirmed nervousness about the current non-mean-1 behaviour; v3
   §5.2 fixes it with a single free scalar. Open question is just the
   default `rateNeoSdlog`.
2. **Topology moves under `"brlens"` unlinked.** MrBayes proves this
   is tractable, but it requires editing every C++ topology move
   (~8 functions in `src/proposals.cpp` and `src/tree_moves.cpp`) to
   accept and update a `NumericMatrix` of per-class branch lengths.
   Estimate ~1 week + tests. **Clean wedge** — defer to Layer 2.
3. **Mixed-type classes [resolved].** v4 supports classes mixing
   neomorphic and transformational chars. `eta_neo` applies only to
   neomorphic chars within each class; the geometric-mean formulation
   (§5.2) keeps the per-class mean-1 constraint clean.
4. **Token vocabulary scope.** v4 uses `"shape"`, `"ratemultiplier"`,
   `"brlens"`. Defer `rateloss`, `betascale`, `pinvar`, etc. to v2.
5. **Composability with `ecology-aware`.** Informational only (§5.3);
   no code dependency. If it lands later, retest.
6. **Log-file column count.** Each `unlink` component adds 1 or
   nEdge columns per class. T3 on a 30-tip tree (~57 edges) × 5 classes
   = 285 brlens columns + 5 shape + 5 `w_c`. AutoPart's
   `03_diagnostics.R` and `04_metrics.R` read by name; the column
   scheme `class<c>_rate_log_sd`, `class<c>_br_<j>`, `w_<c>` is stable.

---

## 9. Test plan

Under `tests/testthat/`:

- `test-partition-bitcompat-null.R`: `partition = NULL` chain matches
  stored reference matrix on a small dataset. Hard guarantee.
- `test-partition-numequiv-singleclass.R`: `partition = rep(1L, nChar),
  unlink = character(0)` reproduces HEAD log-likelihood and log-prior
  to 1e-10 on initial state.
- `test-partition-unlink-shape.R`: 2-class chain with `unlink =
  "shape"`, asserts per-class `rate_log_sd` columns exist, chain
  runs, no errors.
- `test-partition-unlink-rate.R`: same with `unlink = "rate"`,
  asserts `sum(nChar_c * class_rate[c]) == nChar` exactly every
  sample (by construction).
- `test-partition-unlink-brlens.R`: same with `unlink = "brlens"`,
  asserts per-class brlens columns exist, topology is shared
  (matching topology hash across classes), chain runs.
- `test-partition-unlink-all.R`: T3 equivalent — `unlink =
  c("shape", "rate", "brlens")`.
- `test-partition-type-neomorphic.R`: a dataset with neomorphic chars
  in a partitioned run triggers `eta_neo` sampling; mixed-type class
  (neo + trans chars in same class) runs without error.
- `test-partition-validation.R`: bad inputs produce clean errors
  (NAs in partition, out-of-range class IDs, unknown unlink token).
- `test-partition-unlink-matching.R`: partial-prefix match warns;
  ambiguous prefix errors; unknown token suggests via `agrep`.
- `test-partition-rate-neo-reparam.R`: `eta_neo = 1` reproduces the
  default `rate_neo = 1` likelihood; `eta_neo > 1` shifts neo gain and
  loss rates symmetrically around `class_rate[c]` (geometric mean
  preserved).

Smoke test (manual): dispatch all five Casali treatments against a
small matrix via the AutoPart 02 script.

---

## 10. Implementation phasing

1. **Layer 1 — `partition`, `partitionType`, `unlink = "shape"` and
   `unlink = "rate"`.** Plus `rate_neo` → `eta_neo` reparameterisation.
   Lands T0/T1/T2a/T2b/T4 (5 of 6 treatments). No C++ topology
   refactor. ~1 week.
2. **Layer 2 — `unlink = "brlens"`.** C++ topology refactor for
   per-class brlens subParams. Lands T3. ~1 week.
3. **Layer 3 — AutoPart cross-update**: edit
   `auto-part/dev/benchmarks/casali/02_run_treatments.R`'s
   `.CallMkPrime` to pass `partition`, `partitionType`, `unlink`
   instead of `branchModel`; re-run `--dry-run`. ~half a day.

The split between Layers 1 and 2 is the clean wedge: 5 of 6 treatments
ship after Layer 1; T3 follows.

---

## 11. Decisions

All design decisions are now resolved. Layer 1 can start.

1. **Mixed-type classes.** [**resolved 2026-05-20: SUPPORTED.**]
   Classes may freely mix neomorphic and transformational characters;
   `partitionType` arg removed from the signature (§3.2).
2. **`knownStates` handling.** [**resolved 2026-05-20: keep unchanged.**]
   `knownStates` mechanism survives unmodified; `partitionType` removal
   makes this a non-question — per-char types come from existing args.
3. **Token vocabulary scope.** [**resolved 2026-05-20: confirmed
   `c("shape", "ratemultiplier", "brlens")` for v1; defer `rateloss`,
   `betascale`, etc.**]
4. **Token name for rate multiplier.** [**resolved 2026-05-20:
   `"ratemultiplier"` (exact MrBayes match).**]
5. **`rate_neo` → `eta_neo` reparameterisation.** [**resolved
   2026-05-20: geometric-mean formulation (§5.2).** Within each class,
   `sqrt(r_gain_c * r_loss_c) = class_rate[c]`; `eta_neo = r_gain /
   r_loss` is the asymmetry ratio; mean-1 constraint is independent of
   `eta_neo`.]
6. **Numeric-tolerance equivalence contract (§7).** [**resolved
   2026-05-20: accepted as proposed.**]
7. **Scope of v1.** [**resolved 2026-05-20: Layer 1 ships 5/6
   treatments; Layer 2 adds T3. Acceptable.**]
8. **Partial-prefix matching with warning (§3.3).** [**resolved
   2026-05-20: implement as proposed.**]
9. **Ecology-aware composability test.** [**resolved 2026-05-20:
   out of scope; §5.3 is informational only.**]

---

## Appendix A. MrBayes / RevBayes survey crib sheet

- **MrBayes `unlink` parameter list** (`command.c` L12798–12822):
  `Tratio, Revmat, Omega, Statefreq, Shape, Pinvar, Correlation,
  Ratemultiplier, Switchrates, Topology, Brlens, Popsize, Growthrate,
  Aamodel, Cpprate, Cppmultdev, Cppevents, TK02var, WNvar, IGRvar,
  ILNvar`.
- **MrBayes `Ratepr = variable`** (`Help_Prset`, `command.c` L11385–
  L11422): "the rate is allowed to vary across partitions subject to
  the constraint that the **average rate of substitution across the
  partitions is 1** … Dirichlet(1,…,1) prior on the **weighted
  rates**".
- **MrBayes brlens unlink** (`mcmc.c` L11349, `proposal.c` L504,
  L836–921): brlens vectors are subParams of the topology; topology
  moves walk subParams and rewrite each partition's brlens
  consistently.
- **RevBayes partition tutorial**: `part_rate := part_rate_mult *
  n_sites / num_sites_subset[i]` — same char-weighted mean-1
  convention.
- **Casali 2023**: uses MrBayes' `prset ratepr = variable` directly.

Convention is unambiguous: **char-weighted arithmetic mean = 1**,
Dirichlet(1) prior on the weighted simplex. v3 adopts this.
