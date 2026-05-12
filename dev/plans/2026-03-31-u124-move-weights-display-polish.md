# u.124: Move weights display polish

## Summary

Restyle the "Move weights frozen" CLI output to use a tidyverse-inspired
colour scheme. Scope is limited to `.FormatMoveWeights()` in
`R/RunMkPrime.R`.

## Changes

In `.FormatMoveWeights()` (~line 3440):

1. **Category headings** (Topology:, Branches:, etc.): change from
   `cli::col_silver()` to **white** (`cli::col_white()`).
2. **Move names** (nni, spr, etc.): change from `cli::col_silver()` to
   **blue** (`cli::col_blue()`).
3. **Separator**: replace `=` with `:` (e.g. `nni:1.0%` instead of
   `nni=1.0%`).
4. **Value colours**: keep green (≥10%) and yellow (5–10%); change <5%
   from `cli::col_white()` to **silver** (`cli::col_silver()`).

Also update `.FormatMoveWeightsPlain()` to use `:` separator for
consistency (plain text for log files).

No changes to `.PrintMoveWeights()` or `.LogMoveWeights()`.

## Testing

Visual — rebuild, load, and print a sample set of weights to confirm the
styling looks right.
