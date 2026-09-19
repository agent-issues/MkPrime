#!/bin/bash
#SBATCH --job-name=rt-build
#SBATCH --partition=shared
#SBATCH --time=00:45:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/nobackup/%u/mkp-study/red-team/logs/build_%j.out
#SBATCH --error=/nobackup/%u/mkp-study/red-team/logs/build_%j.err
#
# Compile the shared source tree ONCE, as a job the heavy-test arrays can
# depend on.
#
# `pkgload::load_all()` compiles into `src/` in place, so N concurrent array
# tasks race on the same `.o` and `.dll` files; SBC v9 lost 3 of 6 tasks to
# `MkPrime.so: file too short` (#15). The standing mitigation was a convention
# — pre-build on the login node before `sbatch` — and a convention has no way
# to fail loudly when it is forgotten.
#
# Submit this first and make the array depend on it. SLURM then enforces what
# the convention only asked for:
#
#   BUILD=$(sbatch --parsable dev/red-team/heavy-tests/submit-build.sh)
#   sbatch --dependency=afterok:$BUILD dev/red-team/heavy-tests/submit-sbc.sh
#
# `afterok` means the array does not start at all if the build fails, rather
# than starting and failing task by task.

set -euo pipefail

module load r/4.5.1
module load gcc/14.2 || true

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

PROJECT=/nobackup/${USER}/mkp-study
RT=${PROJECT}/red-team
SRC=${RT}/mkp-source
export R_LIBS_USER="${RT}/lib:${PROJECT}/lib"
export R_LIBS="${RT}/lib:${PROJECT}/lib"

cd "${SRC}"

echo "[$(date)] building ${SRC}"
Rscript -e 'pkgload::load_all(getwd(), quiet = FALSE)'

# The arrays' loader checks exactly this: a shared object at least as new as
# every source file under src/. Assert it here so a build that silently did
# nothing fails the dependency rather than the array.
Rscript -e '
  source("dev/red-team/heavy-tests/load-mkprime.R")
  if (!.MkpCompiledAndCurrent(getwd())) {
    stop("build finished but src/ has no current shared object")
  }
  cat("src/ is compiled and current\n")
'

echo "[$(date)] build complete"
