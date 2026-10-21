# Automatic ANY-maze batch processing

Unzip this update inside your existing awn_pipeline directory. It adds
`tools/auto_mwm.py` and `scripts/auto_mwm.slurm`.

## Run on the cluster

The existing `.venv` must contain the pipeline's numpy, pandas, matplotlib dependencies.

```bash
cd "$WORK/jacob/awn_pipeline"
sbatch scripts/auto_mwm.slurm MWM_080326  # first experiment
sbatch scripts/auto_mwm.slurm             # all copied experiments
```

Source: `$WORK/jacob/mwm_exports/<experiment>/metadata/*_total.csv` and
`trajectories/*.csv`. A flat original archive can also be supplied to the Python
runner using `--root "$WORK/jacob/MWM092226"`. Files are read in place.

Outputs: `$WORK/jacob/mwm_results/batch_<jobid>/` with `batch_status.csv`,
`all_trial_metrics.csv`, and a folder for each experiment containing
`manifest.csv`, `mapping_issues.csv`, `source_total.csv`, provenance, and
`run/trial_metrics.csv`, diagnostics, plots and `run/qc_report.html`.

## Matching rules

Requires one total CSV and columns Test, Animal, Stage, Trial per experiment.
Supported individual names: `<Animal>_<Test>.csv` (animal must match exactly)
or `<experiment> - Test <Test>.csv` (matches unique Test number).
Duplicate, absent, unfamiliar, or conflicting filenames are reported.
Duplicate Test numbers reject that experiment. Files never match by row order.

Day1, Day2, etc. map to acquisition with the stated day. Probe_D5, Probe_D8,
Reversal_D1, Habituation_D1, and analogous acquisition/day labels are parsed.
Unknown stages are reported, not guessed. Add explicit mappings with
`--stage-map path.json`, e.g. `{"Training 1":{"stage":"acquisition","day":1}}`.

Subject IDs are experiment::Animal so records from different studies are not
silently treated as the same mouse. Linking repeated mice across studies
requires a verified crosswalk. Original animal IDs remain in manifest notes
and the copied master metadata.

## Config and interpretation

Uses each experiment's `config.json` when present, otherwise the pipeline's
example config. Optional `--config` supplies another global fallback.
Unknown coordinate units and unset gap thresholds allow initial calculations
but remain review flags. Default pass-only summaries can therefore be empty;
all calculated AWN values still appear in trial_metrics.csv. Confirm calibration,
gap settings, plot quality and mapping before final biological comparisons.
Rows with missing X/Y before the first measured coordinate (for example, an
empty ANY-maze startup row at 0 s) are ignored for trajectory calculations and
counted as `n_ignored_leading_xy_rows`. Missing X/Y after tracking starts remain
in the series as gaps and trigger `invalid_xy_rows`; they are not bridged.
Repeated timestamps on rows without valid X/Y are treated as missing-coordinate
event rows: they break the trajectory and remain QC flagged, while timestamp
monotonicity is checked on tracked XY samples. Repeated/reversed timestamps on
tracked XY samples still fail QC.
No group is guessed from missing Treatment. Choose a known metadata column
explicitly with `--group-column Diet` (or Treatment) if appropriate to the study.
Different historical files may need their own exact XY/time columns, delimiter,
encoding, skiprows and stage mappings. These are recorded as errors/issues.

Each invocation requires a new output directory. Exit 2 indicates at least one
mapping issue, failed/error trial, or failed experiment; partial outputs remain
available and the runner continues to other experiments. Review-only QC flags
do not force exit 2. This script performs extraction/QC, not final statistical
inference or automatic approval of biological results.
