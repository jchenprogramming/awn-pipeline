#!/usr/bin/env python3
"""Match ANY-maze Test/Animal exports and batch the existing AWN calculator."""
import argparse
import csv
import json
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from winding_pipeline.io import REQUIRED, OPTIONAL, read_csv_strict, load_settings, sha256, write_json
from winding_pipeline.pipeline import extract


def integer(value):
    value = str(value).strip()
    if not re.fullmatch(r'\d+(?:\.0+)?', value):
        raise ValueError(f'Expected nonnegative integer, got {value!r}')
    return str(int(float(value)))


def stage_day(stage, overrides):
    if stage in overrides:
        return overrides[stage]['stage'], int(overrides[stage]['day'])
    m = re.fullmatch(r'Day[ _-]*(\d+)', stage, re.I)
    if m:
        return 'acquisition', int(m[1])
    m = re.fullmatch(r'(Probe|Acquisition|Reversal|Habituation)[ _-]*(?:D|Day)[ _-]*(\d+)', stage, re.I)
    if m:
        return m[1].lower(), int(m[2])
    raise ValueError(f'Unrecognized Stage {stage!r}; supply --stage-map JSON')


def write_csv(path, rows, fields):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, help='Archive or copied exports; one experiment per subfolder')
    p.add_argument('--out', required=True, help='New output directory; existing targets are rejected')
    p.add_argument('--experiment', help='Process only this experiment folder')
    p.add_argument('--config', help='Global config; experiment/config.json takes precedence')
    p.add_argument('--stage-map', help='JSON mapping unfamiliar stage labels to {stage, day}')
    p.add_argument('--group-column', help='Explicit metadata column for group; otherwise left blank')
    p.add_argument('--no-plots', action='store_true')
    args = p.parse_args()
    root = Path(args.root).resolve()
    out = Path(args.out).resolve()
    if not root.is_dir():
        p.error(f'No export root: {root}')
    if out == root or root in out.parents:
        p.error('Output must be outside the export root')
    folders = [root / args.experiment] if args.experiment else sorted(x for x in root.iterdir() if x.is_dir() and not x.name.startswith('.'))
    if args.experiment and not folders[0].is_dir():
        p.error(f'Experiment not found: {folders[0]}')
    overrides = json.loads(Path(args.stage_map).read_text()) if args.stage_map else {}
    out.mkdir(parents=True, exist_ok=False)
    summaries = []
    combined = []
    for folder in folders:
        target = out / folder.name
        target.mkdir()
        result = dict(experiment=folder.name, status='error', matched=0, mapping_issues=0, pass_trials=0, review_trials=0, fail_trials=0, error_trials=0, detail='')
        try:
            # Copied exports use metadata/ and trajectories/; raw exports are flat.
            meta_dir = folder / 'metadata' if (folder / 'metadata').is_dir() else folder
            traj_dir = folder / 'trajectories' if (folder / 'trajectories').is_dir() else folder
            totals = sorted(x for x in meta_dir.iterdir() if x.is_file() and x.name.lower().endswith('_total.csv') and not x.name.startswith('.'))
            if len(totals) != 1:
                raise ValueError(f'Expected one *_total.csv; found {len(totals)}')
            total = totals[0]
            df = read_csv_strict(total)
            required = {'Test', 'Animal', 'Stage', 'Trial'}
            if not required.issubset(df.columns):
                raise ValueError(f'Total missing columns: {sorted(required - set(df.columns))}')
            records = df.to_dict('records')
            tests = [integer(r['Test']) for r in records]
            if len(tests) != len(set(tests)):
                raise ValueError('Duplicate Test numbers in total CSV')
            candidates = defaultdict(list)
            issues = []
            for f in sorted(traj_dir.iterdir()):
                if not f.is_file() or f.suffix.lower() != '.csv' or f.name.startswith('.') or f.name.lower().endswith('_total.csv'):
                    continue
                # Supported formats: <Animal>_<Test>.csv or <experiment> - Test <Test>.csv.
                match = re.fullmatch(r'.+ - Test\s+(\d+)', f.stem, re.I)
                animal = None
                if not match:
                    match = re.fullmatch(r'(.+)_(\d+)', f.stem)
                    if match:
                        animal = match[1]
                if not match:
                    issues.append(dict(test='', animal='', file=str(f), reason='unrecognized_filename'))
                    continue
                test = integer(match[1] if animal is None else match[2])
                candidates[test].append((f, animal))
            rows = []
            for record, test in zip(records, tests):
                animal = record['Animal'].strip()
                options = candidates.pop(test, [])
                reason = None
                if not animal:
                    reason = 'blank_Animal'
                elif len(options) != 1:
                    reason = f'expected_one_trajectory_found_{len(options)}'
                elif options[0][1] is not None and options[0][1] != animal:
                    reason = f'filename_animal_mismatch:{options[0][1]}'
                try:
                    stage, day = stage_day(record['Stage'].strip(), overrides)
                    trial = integer(record['Trial'])
                except ValueError as e:
                    reason = str(e)
                if reason:
                    issues.append(dict(test=test, animal=animal, file=';'.join(str(x[0]) for x in options), reason=reason))
                    continue
                f = options[0][0]
                row = dict.fromkeys(REQUIRED + OPTIONAL, '')
                row.update(trial_id=f'{folder.name}_Test{test}', subject_id=f'{folder.name}::{animal}',
                           session_id=folder.name, stage=stage, day=day, trial=trial,
                           trajectory_file=str(f.resolve()), sex=record.get('Sex', '').strip(),
                           age_months=record.get('Age (months)', '').strip(),
                           notes=f'Animal={animal}; Test={test}; source Stage={record["Stage"]}')
                if args.group_column:
                    row['group'] = record.get(args.group_column, '').strip()
                rows.append(row)
            for test, options in candidates.items():
                for f, animal in options:
                    issues.append(dict(test=test, animal=animal or '', file=str(f), reason='test_not_in_total'))
            write_csv(target / 'mapping_issues.csv', issues, ['test', 'animal', 'file', 'reason'])
            write_csv(target / 'manifest.csv', rows, REQUIRED + OPTIONAL)
            shutil.copyfile(total, target / 'source_total.csv')
            cfg_path = folder / 'config.json'
            if not cfg_path.is_file():
                cfg_path = Path(args.config) if args.config else Path(__file__).resolve().parents[1] / 'examples/config.json'
            cfg = load_settings(cfg_path)
            write_json(target / 'mapping_provenance.json', dict(total_path=str(total), total_sha256=sha256(total),
                       total_rows=len(records), matched_rows=len(rows), issues=len(issues),
                       group_column=args.group_column, subject_id_policy='experiment::Animal; cross-experiment identity not inferred',
                       stage_overrides=overrides, config_source=str(cfg_path)))
            result.update(matched=len(rows), mapping_issues=len(issues))
            if not rows:
                raise ValueError('No unambiguous mapped trials; see mapping_issues.csv')
            trials = extract(target / 'manifest.csv', cfg, target / 'run', plots=not args.no_plots)
            counts = trials.qc_status.value_counts().to_dict()
            for status in ('pass', 'review', 'fail', 'error'):
                result[f'{status}_trials'] = counts.get(status, 0)
            result['status'] = 'complete_with_issues' if issues or counts.get('error', 0) or counts.get('fail', 0) else 'complete'
            for r in trials.to_dict('records'):
                r['experiment'] = folder.name
                combined.append(r)
        except Exception as e:
            result['detail'] = f'{type(e).__name__}: {e}'
            (target / 'ERROR.txt').write_text(result['detail'] + '\n')
        summaries.append(result)
        print(f'{folder.name}: {result["status"]}; matched={result["matched"]}; mapping issues={result["mapping_issues"]}; {result["detail"]}', flush=True)
        write_csv(out / 'batch_status.csv', summaries, list(result))
    if combined:
        fields = list(dict.fromkeys(k for row in combined for k in row))
        write_csv(out / 'all_trial_metrics.csv', combined, fields)
    print(f'Results: {out}', flush=True)
    print('Review mapping_issues.csv and run/qc_report.html. Unknown units/gap settings remain QC review flags.')
    return 2 if not summaries or any(r['status'] != 'complete' for r in summaries) else 0

if __name__ == '__main__':
    raise SystemExit(main())
