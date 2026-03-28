#!/bin/bash
#SBATCH --job-name=mkp-install-tt
#SBATCH --partition=shared
#SBATCH --time=0:20:00
#SBATCH --ntasks=2
#SBATCH --mem=4G
#SBATCH --output=/nobackup/pjjg18/mkp-study/logs/install_tt_%j.out
#SBATCH --error=/nobackup/pjjg18/mkp-study/logs/install_tt_%j.err

module load r/4.5.1 gcc/14.2
LIB=/nobackup/pjjg18/mkp-study/lib
Rscript -e "
  .libPaths(c('$LIB', .libPaths()))
  pkgs <- c('testthat', 'xml2', 'brio', 'desc', 'praise', 'processx', 'ps', 'R6', 'waldo')
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  cat('Installing:', paste(to_install, collapse=', '), '\n')
  if (length(to_install))
    install.packages(to_install, lib = '$LIB',
      repos = 'https://cloud.r-project.org', Ncpus = 2L)
  cat('Done. testthat:', as.character(packageVersion('testthat', lib.loc='$LIB')), '\n')
"
