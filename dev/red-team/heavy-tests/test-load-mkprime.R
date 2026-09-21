#!/usr/bin/env Rscript
# Checks for the race-safe loader in load-mkprime.R (#15).
#
# Run from the repository root:
#   Rscript dev/red-team/heavy-tests/test-load-mkprime.R
#
# Nothing here compiles or loads MkPrime: the point is exactly that the guard
# decides BEFORE any compilation would start, so it can be exercised against
# fabricated src/ trees.

source("dev/red-team/heavy-tests/load-mkprime.R")

ok <- function(label, cond) {
  cat(sprintf("%-58s %s\n", label, if (isTRUE(cond)) "PASS" else "*** FAIL"))
  if (!isTRUE(cond)) stop(label)
}

# A fake package root with a src/ directory in a chosen state.
MakeTree <- function(shlib, sourceNewer) {
  root <- tempfile("pkg"); dir.create(file.path(root, "src"), recursive = TRUE)
  writeLines("int main() { return 0; }", file.path(root, "src", "mcmc.cpp"))
  writeLines("PKG_CXXFLAGS = -O2",       file.path(root, "src", "Makevars"))
  if (shlib) {
    writeLines("binary",                 file.path(root, "src", "MkPrime.dll"))
    if (sourceNewer) {
      # Touch a source so it is strictly newer than the shared object.
      Sys.sleep(1.1)
      writeLines("int main() { return 1; }", file.path(root, "src", "mcmc.cpp"))
    }
  }
  root
}

built    <- MakeTree(shlib = TRUE,  sourceNewer = FALSE)
stale    <- MakeTree(shlib = TRUE,  sourceNewer = TRUE)
unbuilt  <- MakeTree(shlib = FALSE, sourceNewer = FALSE)

ok("a compiled, current tree is current",  .MkpCompiledAndCurrent(built))
ok("a tree with a newer source is stale", !.MkpCompiledAndCurrent(stale))
ok("a tree with no shared object is not built",
   !.MkpCompiledAndCurrent(unbuilt))
ok("a nonexistent src/ is not built",
   !.MkpCompiledAndCurrent(tempfile("nothing")))

# The guard's decision, without ever reaching pkgload.
WouldRefuse <- function(root, arrayTask = "3", allowCompile = "") {
  old <- Sys.getenv(c("SLURM_ARRAY_TASK_ID", "MKP_ALLOW_COMPILE"),
                    names = TRUE, unset = NA)
  on.exit({
    for (nm in names(old)) {
      if (is.na(old[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv,
                                                          setNames(list(old[[nm]]), nm))
    }
  }, add = TRUE)
  Sys.setenv(SLURM_ARRAY_TASK_ID = arrayTask, MKP_ALLOW_COMPILE = allowCompile)
  res <- tryCatch(LoadMkPrime(root), error = function(e) conditionMessage(e))
  is.character(res) && grepl("Refusing to compile", res)
}

ok("refuses an unbuilt tree inside an array task", WouldRefuse(unbuilt))
ok("refuses a stale tree inside an array task",    WouldRefuse(stale))
ok("names the pre-build command in the message",
   grepl("pkgload::load_all",
         tryCatch({
           Sys.setenv(SLURM_ARRAY_TASK_ID = "3", MKP_ALLOW_COMPILE = "")
           LoadMkPrime(unbuilt)
         }, error = conditionMessage)))
ok("points at the dependency submission too",
   grepl("--dependency=afterok",
         tryCatch({
           Sys.setenv(SLURM_ARRAY_TASK_ID = "3", MKP_ALLOW_COMPILE = "")
           LoadMkPrime(unbuilt)
         }, error = conditionMessage)))
ok("MKP_ALLOW_COMPILE=1 overrides the refusal",
   !WouldRefuse(unbuilt, allowCompile = "1"))

Sys.unsetenv("SLURM_ARRAY_TASK_ID")
Sys.unsetenv("MKP_ALLOW_COMPILE")

# Outside an array — a login node, or a single task — nothing is refused,
# because there is no sibling to race with.
ok("outside an array an unbuilt tree is not refused",
   !isTRUE(tryCatch({
     res <- tryCatch(LoadMkPrime(unbuilt), error = function(e) conditionMessage(e))
     is.character(res) && grepl("Refusing to compile", res)
   }, error = function(e) FALSE)))

cat("\nALL CHECKS PASSED\n")
