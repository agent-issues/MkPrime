#!/bin/bash
set -e
module load r/4.5.1 >/dev/null 2>&1
LIB=/nobackup/pjjg18/mkp-study/red-team/lib
export R_LIBS_USER="$LIB:/nobackup/pjjg18/mkp-study/lib"
Rscript -e '
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (!requireNamespace("goftest", quietly = TRUE)) {
  options(Ncpus = 4)
  install.packages("goftest",
                   repos = "https://cloud.r-project.org",
                   lib = "/nobackup/pjjg18/mkp-study/red-team/lib")
}
cat("goftest:", as.character(packageVersion("goftest")), "\n")
'
