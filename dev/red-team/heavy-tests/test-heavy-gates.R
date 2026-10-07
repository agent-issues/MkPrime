#!/usr/bin/env Rscript
# Checks that the gibbs-spr-db, funnel-stress and subtree-swap-db harnesses and
# the Hamilton install script cannot read green when they should not (#390).
#
# Run from the repository root:
#   Rscript dev/red-team/heavy-tests/test-heavy-gates.R
#
# Nothing here loads MkPrime or runs a harness: the verdict logic is lifted out
# of gibbs-spr-db.R by name, and the other checks read the scripts' source.

ok <- function(label, cond) {
  cat(sprintf("%-66s %s\n", label, if (isTRUE(cond)) "PASS" else "*** FAIL"))
  if (!isTRUE(cond)) stop(label)
}

hdir <- "dev/red-team/heavy-tests"

## ---- gibbs-spr-db verdict and exit status -----------------------------------
src <- parse(file.path(hdir, "gibbs-spr-db.R"))
env <- new.env(parent = globalenv())
for (e in src) {
  if (is.call(e) && identical(as.character(e[[1L]]), "<-") &&
      is.name(e[[2L]]) &&
      as.character(e[[2L]]) %in% c(".Verdict", ".VerdictStatus")) {
    eval(e, envir = env)
  }
}
ok("verdict logic found in the script",
   is.function(env$.Verdict) && is.function(env$.VerdictStatus))

V <- function(genOk = TRUE, ctrlOk = TRUE, fixedOk = TRUE, q1Fail = FALSE,
              q2Fail = FALSE, powerFloorC = 0.5)
  env$.Verdict(genOk, ctrlOk, fixedOk, q1Fail, q2Fail, powerFloorC)

ok("clean run with a power floor is PASS", V() == "PASS")
ok("a failing control is WARN", V(ctrlOk = FALSE) == "WARN")
ok("a failing generator self-check is WARN", V(genOk = FALSE) == "WARN")
ok("no power floor is INCONCLUSIVE", V(powerFloorC = NA_real_) == "INCONCLUSIVE")
ok("gibbs_spr failing is FAIL", V(q1Fail = TRUE) == "FAIL")
ok("a composite br_* arm failing is FAIL, not PASS", V(q2Fail = TRUE) == "FAIL")
ok("PASS exits 0", env$.VerdictStatus("PASS") == 0L)
ok("FAIL exits non-zero", env$.VerdictStatus("FAIL") != 0L)
ok("WARN exits non-zero", env$.VerdictStatus("WARN") != 0L)
ok("INCONCLUSIVE exits non-zero", env$.VerdictStatus("INCONCLUSIVE") != 0L)

gsrc <- readLines(file.path(hdir, "gibbs-spr-db.R"))
ok("gibbs-spr-db loads through load-mkprime.R",
   any(grepl("LoadMkPrime(", gsrc, fixed = TRUE)) &&
     !any(grepl("^\\s*(pkgload|devtools)::load_all", gsrc)))

## ---- funnel-stress pins what its contrast does not vary ---------------------
fsrc <- parse(file.path(hdir, "funnel-stress-hyperprior-sigma.R"))
Calls <- function(x, name) {
  if (is.call(x)) {
    c(if (identical(x[[1L]], as.name(name))) list(x),
      unlist(lapply(as.list(x), Calls, name), recursive = FALSE))
  }
}
modelCalls <- unlist(lapply(fsrc, Calls, "MkPrimeModel"), recursive = FALSE)
ok("funnel-stress builds exactly one model", length(modelCalls) == 1L)
pinned <- names(as.list(modelCalls[[1L]]))
ok("it pins nCat, expSteps and priorVariant",
   all(c("nCat", "expSteps", "priorVariant") %in% pinned))
ok("it no longer prints a headline claimed for NEWS.md",
   !any(grepl("NEWS", readLines(file.path(hdir,
                                          "funnel-stress-hyperprior-sigma.R")))))

## ---- subtree-swap-db is marked retired --------------------------------------
swapR  <- readLines(file.path(hdir, "subtree-swap-db.R"), n = 12L)
swapSh <- readLines(file.path(hdir, "submit-swap-db-n6.sh"))
ok("subtree-swap-db.R header says RETIRED", any(grepl("RETIRED", swapR)))
ok("submit-swap-db-n6.sh says RETIRED", any(grepl("RETIRED", swapSh)))
sh <- file.path(hdir, "submit-swap-db-n6.sh")
out <- suppressWarnings(system2("bash", shQuote(sh), stdout = TRUE,
                                stderr = TRUE,
                                env = "MKP_RUN_RETIRED="))
ok("the submit script refuses to run unless MKP_RUN_RETIRED=1",
   !is.null(attr(out, "status")) && attr(out, "status") != 0L &&
     any(grepl("retired", out)))

## ---- data.table floor in install_dt.R ---------------------------------------
dt <- paste(readLines("data-raw/hamilton/install_dt.R"), collapse = "\n")
ok("install_dt.R reinstalls a data.table older than 1.16.0",
   grepl("1.16.0", dt, fixed = TRUE) && grepl("packageVersion", dt, fixed = TRUE))

cat("\nALL CHECKS PASSED\n")
