#!/bin/bash
#SBATCH --job-name=mkp-install
#SBATCH --partition=shared
#SBATCH --time=0:45:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# `Ncpus = 4L` below parallelises one R process, so it needs four CORES.
# `--ntasks=4` allocated four separate tasks instead, which on the shared
# partition can leave that one process oversubscribing a single core (#103).
#SBATCH --output=/nobackup/pjjg18/mkp-study/logs/install_%j.out
#SBATCH --error=/nobackup/pjjg18/mkp-study/logs/install_%j.err

module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

LIB=/nobackup/pjjg18/mkp-study/lib
mkdir -p "$LIB"

Rscript - <<'EOF'
lib <- "/nobackup/pjjg18/mkp-study/lib"
.libPaths(c(lib, .libPaths()))

# data.table is what summarize_streamed.R reads the logs with; it was absent
# here and installed only by the separate install_dt.R, so a fresh library
# satisfied every arm and then failed at summarisation (#103).
pkgs <- c("Rcpp", "rlang", "cli", "Rdpack", "ape", "TreeTools", "TreeDist",
          "data.table")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) {
  install.packages(to_install,
    lib = lib,
    repos = "https://cloud.r-project.org",
    Ncpus = 4L,
    quiet = FALSE)
}

# Install mkp from tarball
install.packages(
  "/nobackup/pjjg18/mkp-study/MkPrime_0.0.0.9000.tar.gz",
  lib = lib,
  repos = NULL,
  type = "source")

# Verify
library(MkPrime, lib.loc = lib)
cat("MkPrime loaded OK\n")
cat("Version:", as.character(packageVersion("MkPrime", lib.loc = lib)), "\n")
EOF
