#!/usr/bin/env bash
#SBATCH --job-name=mwm-all
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=mwm-all-%j.out

set -euo pipefail

# Process every immediate experiment folder under $WORK/jacob/mwm_exports.
# Each experiment should contain metadata/*_total.csv and trajectories/*.csv,
# or the corresponding flat export files. Source exports are read only.

WORK_ROOT="${WORK:?WORK must be set}"
PROJECT="$WORK_ROOT/jacob/awn_pipeline"
EXPORT_ROOT="$WORK_ROOT/jacob/mwm_exports"
RESULT_ROOT="$WORK_ROOT/jacob/mwm_results"
PYTHON="$PROJECT/.venv/bin/python"

export MPLBACKEND=Agg
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/mwm-mpl-${SLURM_JOB_ID:-manual}"

mkdir -p "$RESULT_ROOT" "$MPLCONFIGDIR"

if [[ ! -d "$PROJECT" ]]; then
  echo "ERROR: Pipeline directory not found: $PROJECT" >&2
  exit 2
fi
if [[ ! -x "$PYTHON" ]]; then
  echo "ERROR: Pipeline Python environment not found or not executable: $PYTHON" >&2
  echo "Install/repair the environment before submitting this job." >&2
  exit 2
fi
if [[ ! -d "$EXPORT_ROOT" ]]; then
  echo "ERROR: MWM export directory not found: $EXPORT_ROOT" >&2
  exit 2
fi

echo "Started: $(date)"
echo "Node: ${SLURMD_NODENAME:-unknown}"
echo "Export root: $EXPORT_ROOT"
echo "Results root: $RESULT_ROOT"
echo "Experiment folders detected: $(find "$EXPORT_ROOT" -mindepth 1 -maxdepth 1 -type d ! -name '.*' | wc -l)"
echo "Each experiment's config.json will be used when present."
echo ""

# No --experiment argument means auto_mwm.py processes every experiment folder.
"$PYTHON" "$PROJECT/tools/auto_mwm.py" \
  --root "$EXPORT_ROOT" \
  --out "$RESULT_ROOT/batch_${SLURM_JOB_ID:-manual}"

echo "Finished: $(date)"
