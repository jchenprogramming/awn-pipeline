#!/usr/bin/env bash
#SBATCH --job-name=mwm-plots
#SBATCH --time=02:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --output=mwm-plots-%j.out
set -euo pipefail
project="${WORK:?WORK must be set}/jacob/awn_pipeline"
export MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/mwm-plots-${SLURM_JOB_ID:-manual-$$}"
mkdir -p "$MPLCONFIGDIR"
export INPUT_CSV="${INPUT_CSV:-$WORK/jacob/mwm_results/analysis_336083/all_trials_with_metadata.csv}"
export PLOT_OUT="${PLOT_OUT:-$WORK/jacob/mwm_results/other_experiment_plots_${SLURM_JOB_ID:-manual-$(date +%Y%m%d_%H%M%S)}}"
"$project/.venv/bin/python" - <<'PY'
import os
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

source = Path(os.environ['INPUT_CSV'])
out = Path(os.environ['PLOT_OUT'])
if not source.is_file():
    raise SystemExit(f'Analysis CSV not found: {source}')
out.mkdir(parents=True, exist_ok=False)
trials = pd.read_csv(source, keep_default_na=False)
required = {'session_id','subject_id','stage','day','trial_id','awn','qc_status','qc_flags'}
missing = required - set(trials.columns)
if missing:
    raise SystemExit(f'Missing columns in {source}: {sorted(missing)}')
trials = trials[trials.session_id.ne('MWM_220207')].copy()
trials['awn'] = pd.to_numeric(trials.awn, errors='coerce')
trials['day'] = pd.to_numeric(trials.day, errors='coerce')
trials = trials[trials.qc_status.isin(['pass','review']) & trials.awn.notna()]

factors = {
    'diet': 'metadata_Diet',
    'genotype': 'metadata_Genotype',
    'treatment': 'metadata_Treatment',
    'sex': 'metadata_Sex',
}
status, summaries = [], []
for experiment, exp_data in trials.groupby('session_id', sort=True):
    exp_dir = out / str(experiment)
    exp_dir.mkdir()
    for factor, col in factors.items():
        if col not in exp_data:
            status.append(dict(experiment=experiment,factor=factor,status='skipped',detail='Metadata column absent'))
            continue
        data = exp_data.copy()
        data['factor_group'] = data[col].astype(str).str.strip()
        data = data[data.factor_group.ne('')].copy()
        # Group labels must be consistent across a mouse's repeated trials.
        if data.groupby('subject_id').factor_group.nunique().gt(1).any():
            status.append(dict(experiment=experiment,factor=factor,status='skipped',detail='Inconsistent label within mouse'))
            continue
        if data.factor_group.nunique() < 2:
            status.append(dict(experiment=experiment,factor=factor,status='skipped',detail='Fewer than two labeled groups'))
            continue
        # For diet/genotype, show sexes separately to avoid pooling them.
        # Treatment labels may already encode sex; retain their recorded labels.
        strata = ['metadata_Sex'] if factor in ('diet','genotype') and 'metadata_Sex' in data else []
        if strata:
            data['plot_group'] = data.factor_group + ' | ' + data.metadata_Sex.astype(str).str.strip()
        else:
            data['plot_group'] = data.factor_group
        factor_dir = exp_dir / factor
        factor_dir.mkdir()
        mouse_day = (data.groupby(['subject_id','stage','day','factor_group','plot_group'], dropna=False)
                     .agg(awn=('awn','mean'), n_trials=('trial_id','size')).reset_index())
        mouse_day.insert(0, 'experiment', experiment)
        mouse_day.to_csv(factor_dir / 'mouse_day_awns.csv', index=False)
        for stage, stage_data in mouse_day.groupby('stage', sort=True):
            fig, ax = plt.subplots(figsize=(9, 5.5))
            for i, (group, values) in enumerate(stage_data.groupby('plot_group', sort=True)):
                color = f'C{i % 10}'
                for _, mouse in values.groupby('subject_id'):
                    mouse = mouse.sort_values('day')
                    ax.plot(mouse.day, mouse.awn, color=color, alpha=.18, linewidth=.85)
                mean = values.groupby('day').awn.agg(['mean','sem']).reset_index()
                ax.errorbar(mean.day, mean['mean'], yerr=mean['sem'], marker='o', capsize=3,
                            color=color, linewidth=2, label=str(group))
                count = values.groupby('day').subject_id.nunique().rename('n_mice').reset_index()
                count['experiment'], count['factor'], count['stage'], count['plot_group'] = experiment, factor, stage, group
                summaries.append(count)
            ax.set(title=f'{experiment} | {factor} | {stage}', xlabel='Study day', ylabel='Mean trial AWN per mouse')
            ax.legend(title='Recorded group (faint lines = individual mice)', fontsize=8)
            ax.grid(alpha=.2)
            fig.tight_layout()
            fig.savefig(factor_dir / f'{stage}_awn.png', dpi=160)
            plt.close(fig)
        status.append(dict(experiment=experiment,factor=factor,status='plotted',
                           detail=f'{data.factor_group.nunique()} labels; pass/review trials; failed/error excluded'))

pd.DataFrame(status).to_csv(out / 'plot_status.csv', index=False)
if summaries:
    pd.concat(summaries, ignore_index=True).to_csv(out / 'group_day_mouse_counts.csv', index=False)
qc = trials.groupby(['session_id','qc_status']).size().rename('n_trials').reset_index()
qc.to_csv(out / 'qc_counts.csv', index=False)
print(pd.DataFrame(status).groupby('status').size().to_string())
print(f'Plots and tables saved to {out}')
print('Review QC flags and small group sizes before interpreting; review trials are included for exploration.')
PY
