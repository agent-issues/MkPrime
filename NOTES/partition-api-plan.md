# Per-character partition API for `RunMkPrime` — design plan (v3)

**Status:** draft v3, awaiting user discussion.
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
| T0 (unpartitioned)        | 1     | — (nothing to unlink)             |
| T1 (anatomical)           | K≈2   | `shape`, `rate`                   |
| T2a (AutoPart, isolated)  | K+1   | `shape`, `rate`                   |
| T2b (AutoPart, merged)    | K     | `shape`, `rate`                   |
| T3 (AutoPart unlinked)    | K+1   | `shape`, `rate`, `brlens`         |
| T4 (random control)       | K+1   | `shape`, `rate`                   |

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
  model. In the new design, type is assigned at the **class** level
  (§3.3), and the per-character `mkd$type` is derived as `type[i] =
  partitionType[ partition[i] ]`.
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
                       partitionType = NULL,            # NEW
                       unlink = character(0),           # NEW
                       ...)
```

### 3.1 `partition`

`NULL` (default) → behaves exactly as today (legacy code path; see §7a).
Or `integer(nChar)` with values in `1:nClasses`, no NAs. Validated
against `mkd$nChar` after invariant-character drop.

### 3.2 `partitionType`

`NULL` (default) → all classes treated as `"transformational"`. Or
`character(nClasses)` with values in `{"transformational",
"neomorphic"}`, one per class. Per-character types are derived as
`mkd$type[i] = partitionType[ partition[i] ]`.

Backward-compat with the legacy `neomorphic = integer(0)` arg:

- If `partition = NULL`, the legacy `neomorphic` arg works as today.
- If `partition` is supplied **and** `partitionType = NULL`, default
  to all-transformational. Pass `partitionType` explicitly to opt in
  to neomorphic classes.
- If both `partition` and the legacy `neomorphic` arg are supplied,
  hard-error: the type spec must come from one source or the other,
  not both.

**Open question (Q1 below):** can a single class contain both
neomorphic and transformational characters? Martin's clarification was
"a partition is flagged as either neomorphic or transformational" —
reads as *no*. I'm going with that constraint.

**Known characters:** the legacy `knownStates` mechanism survives
unchanged. Known-state characters can be members of any user class;
their per-character `k` is still set from `knownStates`. The class's
`partitionType` is then a constraint on the *other* (non-known) chars
in the same class. Acceptable, or do we want a third `partitionType`
value `"known"`? I'd vote no (it makes the class type orthogonal to
the existence of known chars), but flagging.

### 3.3 `unlink`

Character vector of component tokens. Default `character(0)`:
everything linked. Token vocabulary borrowed from MrBayes
(`Help_Unlink`, `command.c` L12798):

| Token         | Effect                                                            | MrBayes analogue   |
|---------------|-------------------------------------------------------------------|--------------------|
| `"shape"`     | per-class `rate_log_sd[c]` (ACRV Γ shape)                         | `shape`            |
| `"rate"`      | per-class `class_rate[c]` multiplier (mean-1 Dirichlet)           | `ratemultiplier`   |
| `"brlens"`    | per-class branch lengths under shared topology (subParam idiom)   | `brlens`           |

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

| State field             | Length when "linked" | Length when "unlinked" | Token        |
|-------------------------|----------------------|------------------------|--------------|
| `class_rate_log_sd`     | 1                    | nClasses               | `"shape"`    |
| `class_rate`            | implicit ≡ 1         | nClasses (on simplex)  | `"rate"`     |
| `class_rel_br_lengths`  | 1 (length-nEdge)     | nClasses (each length-nEdge) | `"brlens"` |
| `class_tree_length`     | 1                    | nClasses               | `"brlens"`   |

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

### 5.2 `rate_neo` reparameterisation under mean-1

`rate_neo` stays as a **single global scalar** (per Martin's
clarification). For mean-1 consistency we drop the current LogNormal(0, 2)
*free* parameterisation and instead derive `(r_neo, r_trans)` from one
free positive scalar `eta_neo > 0`:

- Free parameter: `eta_neo` ~ LogNormal(0, `rateNeoSdlog`), default
  `rateNeoSdlog = 1` (current default).
- Derived: `r_neo = nChar * eta_neo / (nNeo * eta_neo + nTrans)`,
  `r_trans = nChar / (nNeo * eta_neo + nTrans)`.
- Identity: `(nNeo * r_neo + nTrans * r_trans) / nChar = 1` by
  construction. At `eta_neo = 1`, `r_neo = r_trans = 1` (no
  asymmetry).

Where `nNeo`, `nTrans` are dataset-global counts.

**Composition with `class_rate`:** within class `c`, derive per-class
effective rates:

```
r_trans_c = class_rate[c] * nChar_c / (nNeo_c * eta_neo + nTrans_c)
r_neo_c   = eta_neo * r_trans_c
```

so the per-class char-weighted mean of effective rates equals
`class_rate[c]`. The ratio `r_neo_c / r_trans_c = eta_neo` is the same
in every class — that is what Martin's "single rate_neo shared by all
neomorphic characters" means in this formulation. The global
char-weighted mean across all chars in the dataset is exactly 1 by
class_rate's mean-1 constraint.

Edge cases: when `nNeo_c = 0` (pure trans class), `r_neo_c` is
undefined but unused; when `nTrans_c = 0` (pure neo class), `r_trans_c
= class_rate[c] * nChar_c / (nNeo_c * eta_neo)` reduces cleanly.

**For Casali**: all classes have `partitionType = "transformational"`,
so `nNeo = 0` globally and `eta_neo` is never sampled — the production
chains are unaffected by this change.

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

| Component | Linked                                          | Unlinked                                           |
|-----------|-------------------------------------------------|----------------------------------------------------|
| `shape`   | one `rate_log_sd` move (current)                | nClasses `rate_log_sd` moves                       |
| `rate`    | none (vector is implicit ≡ 1)                   | one Dirichlet-simplex move on `w`                  |
| `brlens`  | one `tree_length` + one `rel_br_lengths` move (current) | nClasses copies of each, plus topology moves that update all per-class simplexes |

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
3. **Pure-typed classes constraint.** v3 §3.2 forbids classes mixing
   neomorphic and transformational chars. Acceptable for Casali (all
   trans) and for foreseeable use; if a user needs mixed classes,
   they can simply put neomorphic chars in their own class.
4. **Token vocabulary scope.** v3 proposes `shape`, `rate`, `brlens`.
   Defer `rateloss`, `betascale`, `pinvar`, etc. to v2.
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
- `test-partition-type-neomorphic.R`: a class flagged
  `partitionType = "neomorphic"` runs without error, `eta_neo` is
  sampled.
- `test-partition-validation.R`: bad inputs produce clean errors
  (NAs in partition, out-of-range class IDs, mismatched
  `partitionType` length, unknown unlink token, mixed-type class).
- `test-partition-unlink-matching.R`: partial-prefix match warns;
  ambiguous prefix errors; unknown token suggests via `agrep`.
- `test-partition-rate-neo-reparam.R`: `eta_neo = 1` reproduces the
  default `rate_neo = 1` likelihood; `eta_neo > 1` shifts neo and
  trans rates symmetrically.

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

## 11. Decisions for Martin

Numbered for easy reply. Items Martin has already resolved are marked
[**resolved**].

1. **Mixed-type classes (§3.2).** Forbid (current v3 default), or
   allow? Forbid is simpler and Martin's wording suggests forbid.
2. **`partitionType` value for known-state-space classes (§3.2).**
   Keep `knownStates` mechanism unchanged and treat
   `partitionType = "transformational"` as a no-op for known chars
   in the class (current v3), or add a third
   `partitionType = "known"` value? I lean keep it as-is.
3. **Token vocabulary scope (§3.3).** Confirm `c("shape", "rate",
   "brlens")` for v1; defer `rateloss`, `betascale`, etc.
4. **Token name for "rate"** — terse `"rate"` (my proposal) or
   `"ratemultiplier"` (exact MrBayes match)? Partial-prefix matching
   resolves both.
5. **`rate_neo` → `eta_neo` reparameterisation (§5.2).** Martin
   confirmed the mean-1 nervousness and the "single global scalar"
   scope; my v3 §5.2 is the minimal mean-1-respecting form. Confirm
   this is what you meant.
6. **Numeric-tolerance equivalence contract (§7).** [**resolved
   2026-05-20: accepted as proposed.**]
7. **Scope of v1.** Layer 1 ships 5/6 treatments standalone; Layer 2
   adds T3. Acceptable? (§10)
8. **Partial-prefix matching with warning (§3.3).** [**resolved
   2026-05-20: implement as proposed.**]
9. **Ecology-aware composability test.** Out of scope while
   `ecology-aware` is hypothetical (§5.3). Re-evaluate if it merges.
   [**effectively resolved 2026-05-20.**]

Once 1–5 and 7 are settled, Layer 1 can start.

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
