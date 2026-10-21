#!/usr/bin/env bash
#SBATCH --job-name=mwm-analysis
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=mwm-analysis-%j.out
set -euo pipefail
project="${WORK:?WORK must be set}/jacob/awn_pipeline"
cd "$project"
export MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/awn-analysis-${SLURM_JOB_ID:-manual-$$}"
mkdir -p "$MPLCONFIGDIR"
# Pin the successful submission series; avoid mixing in older or newer reruns.
# Override via: sbatch --export=ALL,FIRST_JOB=123,LAST_JOB=154 ...
export FIRST_JOB="${FIRST_JOB:-336051}" LAST_JOB="${LAST_JOB:-336082}"
export ANALYSIS_OUT="${ANALYSIS_OUT:-$WORK/jacob/mwm_results/analysis_${SLURM_JOB_ID:-manual-$(date +%Y%m%d_%H%M%S)}}"
"$project/.venv/bin/python" - <<'PY'
"""Exploratory AWN associations; never infer genotype or intervention from names.
One experiment per model, daily mean trial AWN, random intercept per mouse.
Pass + review trials are used for exploratory summaries; fail/error excluded.
Models additionally require declared units and gap threshold. No causal claims.
"""
import html
import json
import os
import traceback
from pathlib import Path
import numpy as np
import pandas as pd
from winding_pipeline.pipeline import summaries
from winding_pipeline.analysis import model
from winding_pipeline.io import sha256, versions

root = Path(os.environ['WORK']) / 'jacob/mwm_results'
out = Path(os.environ['ANALYSIS_OUT'])
out.mkdir(parents=True, exist_ok=False)
statuses, exports, sources, figures, coefficients = [], [], [], [], []
seen = set()
factors = {'Diet': 'diet', 'Genotype': 'genotype', 'Treatment': 'treatment', 'Sex': 'sex'}

def status(experiment, factor, stage, state, detail):
    statuses.append(dict(experiment=experiment, factor=factor, stage=stage,
                         status=state, detail=str(detail)))
    pd.DataFrame(statuses).to_csv(out / 'analysis_status.csv', index=False)
    print(f'{experiment} / {factor} / {stage}: {state}: {detail}', flush=True)

def read(path):
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    d.columns = d.columns.str.strip()
    if d.columns.duplicated().any():
        raise ValueError('Duplicate metadata headers')
    return d.apply(lambda col: col.str.strip())

for job in range(int(os.environ['FIRST_JOB']), int(os.environ['LAST_JOB']) + 1):
    batch = root / f'batch_{job}'
    paths = sorted(batch.glob('*/run/trial_metrics.csv'))
    if not paths:
        status(f'batch_{job}', '', '', 'no_trial_output', 'Check batch_status.csv and Slurm log')
    for path in paths:
        experiment = path.parent.parent.name
        try:
            if experiment in seen:
                raise ValueError('Experiment appears in multiple selected jobs; resolve batch selection')
            seen.add(experiment)
            destination = out / experiment
            destination.mkdir()
            meta_path = path.parent.parent / 'source_total.csv'
            meta = read(meta_path)
            trials = read(path)
            if trials.trial_id.duplicated().any():
                raise ValueError('Duplicate trial IDs')
            meta['test_key'] = pd.to_numeric(meta['Test'], errors='raise')
            if meta.test_key.duplicated().any() or not np.isfinite(meta.test_key).all():
                raise ValueError('Missing/duplicate Test metadata')
            trials['test_key'] = pd.to_numeric(trials.trial_id.str.extract(r'_Test(\d+)$')[0], errors='raise')
            # Retain every source field with a prefix, including raw treatment labels.
            meta = meta.rename(columns={c: 'metadata_' + c for c in meta if c != 'test_key'})
            trials = trials.merge(meta, on='test_key', how='left', validate='one_to_one', indicator=True)
            expected = experiment + '::' + trials['metadata_Animal'].fillna('')
            if not trials['_merge'].eq('both').all() or not trials.subject_id.eq(expected).all():
                raise ValueError('Test/Animal identity mismatch')
            trials = trials.drop(columns='_merge')
            trials['source_job'] = job
            for col in ['day', 'trial', 'awn']:
                trials[col] = pd.to_numeric(trials[col], errors='raise')
            if not trials.session_id.eq(experiment).all():
                raise ValueError('Unexpected session ID')
            # Use explicit metadata only. Normalize recorded sex spelling.
            for source, target in [('Sex','sex'), ('Genotype','genotype'), ('Age (months)','age_months')]:
                name = 'metadata_' + source
                trials[target] = trials[name].fillna('') if name in trials else ''
            trials['sex'] = trials.sex.replace({'Female':'F','Male':'M','female':'F','male':'M'})
            trials.to_csv(destination / 'trials_with_metadata.csv', index=False)
            exports.append(trials)
            config_path = path.parent / 'config.json'
            cfg = json.loads(config_path.read_text())
            sources.append(dict(experiment=experiment, job=job, trials=str(path),
                                trials_sha256=sha256(path), metadata_sha256=sha256(meta_path),
                                config=cfg))
            pd.DataFrame([dict(qc_status=k, n_trials=v) for k,v in trials.qc_status.value_counts().items()]).to_csv(destination / 'qc_counts.csv', index=False)
            for source, factor in factors.items():
                column = 'metadata_' + source
                if column not in trials:
                    status(experiment, factor, '', 'skipped', 'Metadata column absent')
                    continue
                d = trials.copy()
                d['group'] = d[column].fillna('')
                if source == 'Sex':
                    d['group'] = d['sex']
                # Group assignment must be stable within each mouse.
                if d.groupby('subject_id').group.nunique().gt(1).any():
                    status(experiment, factor, '', 'skipped', 'Missing/inconsistent group labels within mouse')
                    continue
                missing = d.group.eq('')
                d.loc[missing].to_csv(destination / f'{factor}_missing_labels.csv', index=False)
                d = d.loc[~missing].copy()
                if d.empty:
                    status(experiment, factor, '', 'skipped', 'All group labels blank')
                    continue
                # Record metadata overlap; composite Treatment labels are never decomposed.
                overlap_cols = ['subject_id','group','sex','genotype','age_months']
                for extra in ['metadata_Diet','metadata_Treatment']:
                    if extra in d:
                        overlap_cols.append(extra)
                d[overlap_cols].drop_duplicates().to_csv(destination / f'{factor}_mouse_metadata.csv', index=False)
                summary_dir = destination / factor
                summary_dir.mkdir()
                day = summaries(d, summary_dir, include_review=True)
                for figure in sorted((summary_dir / 'summary_figures').glob('*.png')):
                    figures.append((f'{experiment}: {source}', figure.relative_to(out).as_posix()))
                status(experiment, factor, '', 'plots_and_summary', 'Exploratory pass + review; failed/error excluded')
                for stage in sorted(day.stage.unique()):
                    usable = day[(day.stage == stage) & (day.n_trials_included > 0)].copy()
                    if usable.group.nunique() < 2:
                        status(experiment, factor, stage, 'model_skipped', 'Fewer than two groups with usable AWN')
                        continue
                    if cfg.get('coordinate_unit', 'unknown') in ('unknown', '', None) or cfg.get('max_time_gap_s') is None:
                        status(experiment, factor, stage, 'model_skipped', 'Units/gap threshold unresolved; plots remain preliminary')
                        continue
                    # Adjust for recorded sex/genotype/age when they vary. Rank checks
                    # in the pipeline block fully confounded or overparameterized designs.
                    covars = []
                    for cov in ['sex','genotype','age_months']:
                        if cov == factor:
                            continue
                        if usable[cov].replace('', np.nan).dropna().nunique() > 1:
                            covars.append(cov)
                    if any(d.groupby('subject_id')[c].nunique().gt(1).any() for c in ['sex','genotype']):
                        status(experiment, factor, stage, 'model_skipped', 'Sex/genotype inconsistent within mouse')
                        continue
                    try:
                        result = model(summary_dir / 'subject_day.csv', summary_dir / f'model_{stage}',
                                       experiment, stage, outcome='awn', covariates=covars)
                        status(experiment, factor, stage, result['status'],
                               'Exploratory association; covariates=' + ','.join(covars) +
                               '; inspect mouse metadata for age/diet/genotype overlap')
                        if result['status'] == 'fit_ok_exploratory':
                            coef = pd.read_csv(summary_dir / f'model_{stage}/fixed_effects.csv')
                            coef['experiment'], coef['factor'], coef['stage'] = experiment, factor, stage
                            coefficients.append(coef)
                    except Exception as error:
                        status(experiment, factor, stage, 'model_skipped_or_failed', error)
        except Exception as error:
            status(experiment, '', '', 'experiment_error', error)
            (out / f'ERROR_{experiment}.txt').write_text(traceback.format_exc())

if exports:
    pd.concat(exports, ignore_index=True).to_csv(out / 'all_trials_with_metadata.csv', index=False)
if coefficients:
    coefs = pd.concat(coefficients, ignore_index=True)
    # BH across all finite, non-intercept coefficient tests from successful fits.
    # These are individual contrasts, NOT omnibus group or interaction tests.
    mask = coefs.term.ne('Intercept') & np.isfinite(coefs.p_wald_unadjusted)
    p = coefs.loc[mask, 'p_wald_unadjusted'].to_numpy()
    if len(p):
        order = np.argsort(p)
        q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
        adjusted = np.empty(len(p))
        adjusted[order] = np.minimum(q, 1)
        coefs.loc[mask, 'q_bh_all_nonintercept_contrasts'] = adjusted
    coefs.to_csv(out / 'exploratory_model_coefficients.csv', index=False)
(out / 'provenance.json').write_text(json.dumps(dict(
    sources=sources, versions=versions(), first_job=os.environ['FIRST_JOB'], last_job=os.environ['LAST_JOB'],
    outcome='awn', aggregation='Mean trial AWN within mouse/stage/day',
    inclusion='pass + review; fail/error excluded; preliminary exploratory analysis',
    factors=factors, model='group * categorical day + varying recorded covariates; mouse random intercept',
    limitation='Not design approval or causal inference. Missing genotype/age cannot be reconstructed. Treatment is a composite label. No pooling of experiments. No omnibus tests.'
), indent=2))
intro = ('<h1>Exploratory AWN progress report</h1><p>Daily mean trial AWN per mouse. '
         'Pass and review trials included; fail/error excluded. Faint lines are mice; '
         'group curves show mean ± SEM. Acquisition, probe, and reversal are separate. '
         'Recorded labels are not proof of comparable groups or causal effects. '
         'Treatment labels may combine genotype, age, sex, and intervention. '
         'Models require resolved extraction settings and numerical checks. '
         'Read analysis_status.csv and mouse metadata before interpreting contrasts.</p>')
body = pd.DataFrame(statuses).to_html(index=False, escape=True)
for label, relative in figures:
    body += f'<h2>{html.escape(label)}</h2><img src="{html.escape(relative, quote=True)}" style="max-width:100%;width:850px">'
(out / 'report.html').write_text('<!doctype html><meta charset="utf-8"><title>AWN progress</title>'
    '<style>body{font:16px system-ui;margin:2rem;max-width:1200px}td,th{padding:.4rem;text-align:left}table{border-collapse:collapse}tr{border-bottom:1px solid #ddd}</style>' + intro + body)
print(f'Finished. Download the entire folder: {out}', flush=True)
PY
