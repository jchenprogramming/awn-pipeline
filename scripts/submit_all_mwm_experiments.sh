#!/usr/bin/env bash
# Submit one independent Slurm job per exported MWM experiment folder.
set -euo pipefail

WORK_ROOT="${WORK:?WORK must be set}"
PROJECT="$WORK_ROOT/jacob/awn_pipeline"
EXPORT_ROOT="$WORK_ROOT/jacob/mwm_exports"
SBATCH_SCRIPT="$PROJECT/scripts/auto_mwm.slurm"

if [[ ! -f "$SBATCH_SCRIPT" ]]; then
  echo "ERROR: Per-experiment Slurm script not found: $SBATCH_SCRIPT" >&2
  exit 2
fi
if [[ ! -d "$EXPORT_ROOT" ]]; then
  echo "ERROR: Export directory not found: $EXPORT_ROOT" >&2
  exit 2
fi

experiments=()
for dir in "$EXPORT_ROOT"/*/; do
  [[ -d "$dir" ]] || continue
  name="$(basename "$dir")"
  [[ "$name" == .* ]] && continue
  # Only submit folders containing an ANY-maze total export, either in
  # metadata/ or directly in the experiment folder.
  if compgen -G "$dir"'metadata/*_total.csv' >/dev/null || \
     compgen -G "$dir"'*_total.csv' >/dev/null; then
    experiments+=("$name")
  else
    echo "SKIP $name: no *_total.csv found at expected level"
  fi
done

if [[ ${#experiments[@]} -eq 0 ]]; then
  echo "ERROR: No experiment folders with *_total.csv found under $EXPORT_ROOT" >&2
  exit 2
fi

echo "Submitting ${#experiments[@]} experiment jobs from $EXPORT_ROOT"
for experiment in "${experiments[@]}"; do
  echo -n "$experiment: "
  sbatch --parsable "$SBATCH_SCRIPT" "$experiment"
done
