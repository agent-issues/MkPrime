#!/bin/bash
#SBATCH --job-name=mkp-check
#SBATCH --partition=shared
#SBATCH --time=0:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/nobackup/pjjg18/mkp-study/logs/mkp_check_%j.out
#SBATCH --error=/nobackup/pjjg18/mkp-study/logs/mkp_check_%j.err

# CI-style R CMD check + testthat (including slow detailed-balance tests).
# Usage:
#   scp MkPrime_0.0.0.9000.tar.gz hamilton:/nobackup/pjjg18/mkp-study/
#   sbatch mkp_check.slurm

module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKPRIME_SLOW_TESTS=true

STUDY=/nobackup/pjjg18/mkp-study
LIB=$STUDY/lib
TARBALL=$STUDY/MkPrime_0.0.0.9000.tar.gz

mkdir -p "$LIB" "$STUDY/logs"

echo "=== Installing MkPrime from tarball ==="
Rscript -e "
.libPaths(c('$LIB', .libPaths()))
install.packages('$TARBALL', lib = '$LIB', repos = NULL, type = 'source')
cat('Installed:', as.character(packageVersion('MkPrime', lib.loc = '$LIB')), '\n')
"

echo ""
echo "=== Running testthat (MKPRIME_SLOW_TESTS=true) ==="
Rscript -e "
.libPaths(c('$LIB', .libPaths()))
library(testthat)
library(MkPrime, lib.loc = '$LIB')
cat('MKPRIME_SLOW_TESTS =', Sys.getenv('MKPRIME_SLOW_TESTS'), '\n')
results <- test_package('MkPrime', reporter = 'summary')
q(status = any(as.data.frame(results)\$failed > 0))
"

echo ""
echo "=== Done (exit $?) ==="
