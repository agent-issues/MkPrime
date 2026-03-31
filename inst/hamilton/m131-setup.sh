#!/bin/bash
# M-131: Hamilton setup script
# Run on Hamilton after uploading files to /nobackup/$USER/m131/

set -e

BASEDIR=/nobackup/$USER/m131
LIB=$BASEDIR/lib
cd $BASEDIR

mkdir -p lib logs results data-raw

module load r/4.5.1
module load gcc/14.2

# Install MkPrime from tarball (and dependencies)
echo "Installing dependencies..."
Rscript -e "
  lib <- '$LIB'
  .libPaths(c(lib, .libPaths()))
  install.packages(c('Rcpp', 'ape', 'TreeTools', 'cli'),
                   lib = lib, repos = 'https://cloud.r-project.org',
                   Ncpus = 4L)
"

echo "Installing MkPrime..."
R CMD INSTALL --library=$LIB MkPrime_*.tar.gz

echo "Setup complete. Submit with: sbatch m131-warmup-validation.slurm"
