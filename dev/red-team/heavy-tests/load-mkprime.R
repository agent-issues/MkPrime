# Load MkPrime for a heavy test without racing on the shared src/ tree.
#
# `pkgload::load_all()` compiles into `src/` **in place**. Under a SLURM array
# every task does that against the same tree, so they race on the same `.o` and
# `.dll`: SBC v9 lost 3 of 6 tasks to `MkPrime.so: file too short` (#15).
#
# The mitigation has been a convention — pre-build once on the login node before
# `sbatch` — and the failure mode when someone forgets is a partial array with
# no obvious cause, easily misread as a code failure. This makes it structural:
# under an array task, a tree that is not already compiled and current is a hard
# stop, so nothing ever starts a compile that a sibling task could corrupt.
#
# `load_all` is kept in preference to an installed package on purpose. A heavy
# test exists to exercise the branch's sources, and silently running an
# installed copy instead is SBC-HARNESS-006 all over again — the failure where
# six weeks of runs used a stale out-of-tree script without anyone noticing.
#
# Escape hatch: `MKP_ALLOW_COMPILE=1` permits compiling even inside an array
# task, for the case where you know only one task is running.

.MkpCompiledAndCurrent <- function(pkgRoot = ".") {
  src <- file.path(pkgRoot, "src")
  if (!dir.exists(src)) return(FALSE)

  shlib <- list.files(src, "\\.(dll|so|dylib)$", full.names = TRUE)
  if (!length(shlib)) return(FALSE)

  sources <- list.files(src, "\\.(cpp|h|hpp|c)$", full.names = TRUE)
  makevars <- list.files(src, "^Makevars", full.names = TRUE)
  inputs <- c(sources, makevars)
  if (!length(inputs)) return(FALSE)

  max(file.mtime(shlib)) >= max(file.mtime(inputs))
}

LoadMkPrime <- function(pkgRoot = ".") {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    # No pkgload: an installed package is the only option, and it cannot race.
    library(MkPrime)
    return(invisible("installed"))
  }

  inArray <- nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  allowed <- nzchar(Sys.getenv("MKP_ALLOW_COMPILE"))

  if (inArray && !allowed && !.MkpCompiledAndCurrent(pkgRoot)) {
    stop(
      "Refusing to compile inside SLURM array task ",
      Sys.getenv("SLURM_ARRAY_TASK_ID"), ".\n",
      "  src/ under ", normalizePath(pkgRoot, mustWork = FALSE),
      " is not built, or a source file is newer than the shared object.\n",
      "  Every task in this array would compile into the same directory and\n",
      "  corrupt each other's objects (#15; SBC v9 lost 3 of 6 tasks this way).\n",
      "\n",
      "  Pre-build once on the login node, then resubmit:\n",
      "    cd <src> && Rscript -e 'pkgload::load_all(getwd())'\n",
      "\n",
      "  Or submit the build as a dependency, which SLURM enforces:\n",
      "    BUILD=$(sbatch --parsable dev/red-team/heavy-tests/submit-build.sh)\n",
      "    sbatch --dependency=afterok:$BUILD dev/red-team/heavy-tests/submit-sbc.sh\n",
      "\n",
      "  Set MKP_ALLOW_COMPILE=1 to override (only when one task is running).",
      call. = FALSE
    )
  }

  pkgload::load_all(pkgRoot, quiet = TRUE)
  invisible("load_all")
}
