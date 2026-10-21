#!/usr/bin/env bash
#SBATCH --job-name=mwm-auto
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=mwm-auto-%j.out
set -euo pipefail
project="${WORK:?WORK must be set}/jacob/awn_pipeline"
export MPLBACKEND=Agg
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
mkdir -p "$WORK/jacob/mwm_results"
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/mwm-mpl-${SLURM_JOB_ID}"
mkdir -p "$MPLCONFIGDIR"
# Optional first argument selects one experiment; no argument processes all.
selection=()
if [[ $# -gt 0 ]]; then selection=(--experiment "$1"); fi
"$project/.venv/bin/python" "$project/tools/auto_mwm.py" \
  --root "$WORK/jacob/mwm_exports" \
  --out "$WORK/jacob/mwm_results/batch_${SLURM_JOB_ID}" \
  "${selection[@]}"
