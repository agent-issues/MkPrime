# The suite should report test outcomes and nothing else, so run MkPrime in
# silent mode.  A test that asserts on console output must raise the level for
# itself -- `withr::local_options(MkPrime.verbosity = 1)`, or a `verbosity`
# argument on the call under test -- and set an expectation.
options(MkPrime.verbosity = 0L)

# `ape` and `TreeTools` are attached once here rather than inside individual
# test files, so the search path is the same for every file.  A file that must
# attach them itself uses `library(quietly = TRUE)`, which suppresses `ape`'s
# startup banner (re-emitted via TreeTools' Depends); re-attaching MkPrime
# itself errors outright.
# `ape` goes last so that it masks `TreeTools`, as it did when each file
# attached them for itself.
suppressPackageStartupMessages({
  library("TreeTools")
  library("ape")
})

# Suggested packages announce themselves the first time they load: `phangorn`
# reports re-registering `[.phyDat` over TreeTools' method, and
# `shiny::testServer()` attaches `shiny` via `require()`.  Do both here, so the
# announcement does not surface mid-test and the search path is the same for
# every file rather than changing once test-bayesian-module.R has run.
suppressMessages(suppressPackageStartupMessages({
  requireNamespace("phangorn", quietly = TRUE)
  if (requireNamespace("shiny", quietly = TRUE)) library("shiny")
}))
