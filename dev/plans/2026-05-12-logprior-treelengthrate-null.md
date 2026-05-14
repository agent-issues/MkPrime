# LogPrior returns numeric(0) when treeLengthRate is unset

## Problem

`MkPrimeModel()` leaves `treeLengthRate = NULL` unless either
`treeLengthRate` or `expSteps` is supplied. The rate is normally derived
from `expSteps` inside `.FinalizeModel()` (rate = 2 / expSteps), so any
caller that goes through `RunMkPrime()` is fine.

But `LogPrior(state, model, mkd)` is also called directly from tests
and from `.InitState()`, and it does:

```r
lp <- lp + dgamma(state$tree_length,
                  shape = model$treeLengthShape,
                  rate = model$treeLengthRate,
                  log = TRUE)
```

`dgamma(x, rate = NULL)` returns `numeric(0)`, which silently poisons
`lp`. The function then returns `numeric(0)` (not `-Inf`, not a scalar),
which propagates into `state$log_prior` and breaks the MH ratio.

## Reproducer

```r
mat <- matrix(c(0, 1, 0, 1, 2, 1, 0, 1, 2, 0),
              nrow = 5, ncol = 2,
              dimnames = list(paste0("t", 1:5), NULL))
pd  <- TreeTools::MatrixToPhyDat(mat)
mkd <- MkPrimeData(pd)
state <- list(
  tree_length = 1, rate_log_sd = 0.5,
  rel_br_lengths = rep(1/7, 7),
  rate_loss = 1, kPrime = mkd$kObs, p = 0.5
)
LogPrior(state, MkPrimeModel(), mkd)
#> numeric(0)
```

Caught while writing tests for the ecology-aware NT model: any test that
constructs an `MkPrimeModel` without `expSteps` and then evaluates
`LogPrior` directly produces a length-0 result. Workaround so far has
been to always pass `expSteps = 10` in test fixtures.

## Fix options

1. **Validate in `LogPrior`**: error / return `-Inf` if
   `model$treeLengthRate` is NULL. Surfaces the misuse immediately and
   makes the contract explicit.
2. **Default in `MkPrimeModel`**: if neither `treeLengthRate` nor
   `expSteps` is supplied, fill a sensible placeholder (e.g.
   `treeLengthRate = 2 / 10`). Hides the contract — any model that
   isn't finalised would silently use this default.
3. **Mark `MkPrimeModel` output as incomplete until finalised**: add a
   class flag, error in `LogPrior` if the flag is unset. Cleanest but
   touches everywhere LogPrior is called.

Option 1 is the smallest fix and matches the existing `-Inf` boundary
behaviour of the rest of `LogPrior`. Probably the right call.

## Out of scope

- Ecology-aware extensions (separate worktree).
- Other unfinalised hyperparameters that may have similar issues
  (`rateLossMeanlog`/`Sdlog`, etc.) — audit alongside whichever fix
  lands.
