#!/bin/bash
set -e
module load r/4.5.1
module load gcc/14.2 || true
export PROJECT=/nobackup/$USER/mkp-study
export RT=$PROJECT/red-team
export SRC=$RT/mkp-source
export R_LIBS_USER="$RT/lib:$PROJECT/lib"
export R_LIBS="$RT/lib:$PROJECT/lib"
cd "$SRC"
echo "=== precompile via pkgload (post-patch f9ea652) ==="
Rscript -e 'pkgload::load_all(".", quiet = TRUE); cat("OK\n")'
echo "=== compiled .so present ==="
ls -la src/*.so 2>/dev/null | head -3
