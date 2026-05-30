#!/bin/bash
#SBATCH --job-name=rt-marg-k-agg
#SBATCH --partition=shared
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --gres=tmp:2G
#SBATCH --output=/nobackup/%u/mkp-study/red-team/logs/marg-k-agg_%j.out
#SBATCH --error=/nobackup/%u/mkp-study/red-team/logs/marg-k-agg_%j.err
#
# Aggregation step for the marginal-k SBC array (submit-marginal-k-sbc.sh).
# Submit with an afterany dependency on the array job:
#
#     ARR=$(sbatch --parsable submit-marginal-k-sbc.sh)
#     sbatch --dependency=afterany:${ARR} submit-marginal-k-agg.sh
#
# afterany (not afterok): a single dead shard then costs its 5 sims, not the
# whole verdict — the harness' min_good filter tolerates missing shards.
#
# Runs the SAME harness in AGGREGATE mode: it merges shards/sims-shard-*.rds
# into the full sims list and runs the identical AD / verdict code, so the
# verdict is identical-by-construction to what a monolith would have produced.

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

mkdir -p "${RT}/logs"
mkdir -p "${SRC}/dev/red-team/heavy-tests/marginal-k/sbc-results-hamilton"

echo "[$(date)] marginal-k SBC AGGREGATE (job ${SLURM_JOB_ID})"
cd "${SRC}"

export MARGINAL_K_SBC_NSHARD=40
export MARGINAL_K_SBC_AGGREGATE=1
Rscript "${SRC}/dev/red-team/heavy-tests/marginal-k/T-SBC-marginal-geometric.R"

# Mirror artefacts into the Hamilton-results dir for easier collection.
RESULTS_SRC="${SRC}/dev/red-team/heavy-tests/marginal-k/sbc-results"
RESULTS_DST="${SRC}/dev/red-team/heavy-tests/marginal-k/sbc-results-hamilton"
if [ -d "${RESULTS_SRC}" ]; then
  cp -r "${RESULTS_SRC}"/* "${RESULTS_DST}/" || true
fi

echo "[$(date)] marginal-k SBC aggregate complete"
echo "--- verdict.txt ---"
cat "${RESULTS_SRC}/verdict.txt" || true
