#!/bin/bash
# Submit T9 Model A beta_geometric SBC to Hamilton.
#
# Usage (from login node, after git pull + pre-build):
#   bash ${SRC}/inst/hamilton/submit-t9-bg-sbc.sh
#
# Pre-build requirement (race-avoidance; SBC-BUILD-RACE-001):
#   module load r/4.5.1 && (module load gcc/14.2 || true)
#   export R_LIBS_USER=/nobackup/pjjg18/mkp-study/red-team/lib:/nobackup/pjjg18/mkp-study/lib
#   export R_LIBS=$R_LIBS_USER
#   Rscript -e 'pkgload::load_all(".")'

set -euo pipefail

SRC=/nobackup/pjjg18/mkp-study/red-team/mkp-source
RT=/nobackup/pjjg18/mkp-study/red-team
LIB=${RT}/lib:/nobackup/pjjg18/mkp-study/lib

cat > /tmp/t9-bg-sbc.slurm << 'SLURM'
#!/bin/bash
#SBATCH --job-name=mkp-t9-bg-sbc
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=/nobackup/pjjg18/mkp-study/red-team/t9-bg-sbc-%j.log
#SBATCH --error=/nobackup/pjjg18/mkp-study/red-team/t9-bg-sbc-%j.err

set -euo pipefail

SRC=/nobackup/pjjg18/mkp-study/red-team/mkp-source
RT=/nobackup/pjjg18/mkp-study/red-team
LIB=${RT}/lib:/nobackup/pjjg18/mkp-study/lib

module load r/4.5.1
module load gcc/14.2 || true

export R_LIBS_USER=${LIB}
export R_LIBS=${LIB}

cd ${SRC}
Rscript ${SRC}/dev/red-team/heavy-tests/kprime-viability/T9-bg-model-a-sbc.R
SLURM

JOB=$(sbatch /tmp/t9-bg-sbc.slurm | awk '{print $NF}')
echo "Submitted T9 beta_geo SBC as job ${JOB}"
echo "Log: /nobackup/pjjg18/mkp-study/red-team/t9-bg-sbc-${JOB}.log"
