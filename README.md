# QIAL AWN pipeline — Jacob's working version

Version 0.1.0 • Built from an audit of the supplied `awn-main(1).zip` • 2026-10-06

## What this does

Turns explicitly identified ANY-maze trajectory CSVs into trial-level **absolute
tangent turning (AWN)**, tracking QC, reviewable PNGs, subject-day summaries and
optional exploratory longitudinal models. It preserves original exports and
never infers mouse identity from file order. No GUI, XQuartz, GPU, or MATLAB is
needed for Python extraction.

**Status:** research software tested on mathematical fixtures and synthetic
trajectories. NOT yet validated against real exports or execution in MATLAB.
The original ZIP has summary tables, not raw trajectories or paired numerical
MATLAB validation outputs. The optional statsmodels fitting test could not run
in the build environment. See `docs/AUDIT.md` and `docs/VALIDATION.md`.

## Quick start on your cluster

Place this new directory inside your own folder, for example:

```text
/mnt/newStor/paros/paros_WORK/jacob/winding_project/awn_pipeline/
```

Keep Jimin's and Aurora's folders unchanged. Follow your lab's compute-node/job
policy for batch work. Use the Python environment your lab recommends. In this
directory, with Python 3.10 or newer:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[stats]'
python -m unittest discover -s tests -v
```

Installation needs package access; ask for a lab-provided environment if package
downloads are blocked. Do not change system Python or shared installations.
Core extraction needs numpy, pandas and matplotlib only. Without installation,
`python -m winding_pipeline` works from this folder if dependencies exist.

### 1. Run a synthetic practice batch first

```bash
python -m winding_pipeline demo --out practice_inputs
python -m winding_pipeline extract --manifest practice_inputs/manifest.csv --config practice_inputs/config.json --out practice_run
```

The practice data are synthetic and must never be reported as mouse results.
Open `practice_run/qc_report.html` locally after downloading the whole run folder
(including `figures/`). The HTML links use relative paths. PNGs work without X11.
Use a new output directory for every run; existing directories are rejected.

### 2. Inventory a small, known export folder

```bash
python -m winding_pipeline inventory --directory /absolute/path/to/your/trajectory_exports --out export_inventory
```

Only that folder is listed, not the broader shared storage. This creates
`export_inventory/manifest_to_complete.csv`; all animal/trial fields are blank
by design. Match them using the ANY-maze test schedule, not a presumed 27-animal
rotation. Filter out non-trajectory CSVs before running extraction.

### 3. Prepare a verified manifest and config

Each row is exactly one trajectory/test. Required fields:

| Field | Meaning |
|---|---|
| trial_id | Unique verified test identifier; can combine cohort and ANY-maze test ID |
| subject_id | Stable globally unique mouse ID, kept as text |
| session_id | Explicit study/testing session, not inferred from a filename |
| stage | `acquisition`, `probe`, `reversal`, or `habituation` |
| day | Nonnegative study-day integer; acquisition day 5 and probe day 5 remain separate |
| trial | Nonnegative trial number within the stage/day |
| trajectory_file | Exact path; relative paths resolve against the manifest's directory |

Optional fields: `group,cohort,sex,genotype,age_months,notes`. Missing treatment
does NOT become control. Copy optional metadata explicitly from your verified
master sheet. Global subject IDs must not be reused for distinct mice across
cohorts. The supported design assumes group/sex/genotype/cohort stable within
mouse; crossover or changing-treatment designs need a different schema/model.

Edit a copy of `examples/config.json`:

- Set exact time/X/Y column names, delimiter, decimal character and encoding.
- `skiprows` is the number of physical preamble lines before the CSV header.
- Times must be elapsed seconds or elapsed `MM:SS` / `HH:MM:SS`; calendar dates,
  midnight rollover, milliseconds-as-numeric and Excel day fractions are NOT
  inferred. Convert those explicitly if needed.
- Set `coordinate_unit` from the export/calibration, not from the data's scale.
  `coordinate_scale` applies equally to both axes; unit names describe values
  AFTER scaling. Default 1 means no conversion. No automatic factor-of-ten fix.
- Set `max_time_gap_s` after inspecting the export cadence and protocol. Null
  disables that filter and marks the trial for review. There is no universal
  safe value. The synthetic value 0.5 s is NOT a biological recommendation.
- `max_speed` is optional, in coordinate units/second. Choose from protocol and
  tracking validation, not to force a treatment effect. Null means unchecked.
- `min_step` defaults to zero: exactly stationary steps cannot define headings.
  A positive value breaks short-step turns; it is not a validated jitter filter.
- `resample_step` optionally adds a separately named spatial-resolution
  sensitivity metric. It does not replace `awn` or silently smooth the data.
- `y_axis` is the plotting convention only. Reflections do not change absolute AWN.

Export only needed time/centre-coordinate fields where feasible. The software
cannot distinguish an event-only empty row from an untracked frame; it treats
both as gaps and flags them. Do not drop those rows until their meaning is known.

### 4. Extract a small verified subset

```bash
python -m winding_pipeline extract --manifest export_inventory/manifest_to_complete.csv --config your_config.json --out first_real_run
```

Start with 3–5 diverse trials: clean movement, a pause, a missing section and a
probe. Compare their tracks and AWN with MATLAB on the SAME rows and columns.
Do not launch the full study until mapping and numerical agreement are checked.

### 5. Review, then summarize

`pass` = numerical checks passed, NOT visual/biological approval.
`review` = calculable but requires review of listed flags.
`fail` = non-increasing timestamps or no usable turns.
`error` = file/parser/processing problem. All records remain in trial_metrics.csv.

Default summaries use only pass trials and keep rows for empty subject-days.
This is a conservative computational default, NOT a validated exclusion rule.
Complete `review_template.csv` with `include`/`exclude`, reviewer and reason for
EVERY trial after checking mapping and plots, then run:

```bash
python -m winding_pipeline summarize --trials first_real_run/trial_metrics.csv --review completed_review.csv --out reviewed_summary
```

Failed/error trials cannot be included; fix the input/config and rerun. A
separate `--include-review` option is provided for explicitly labeled sensitivity
analysis, never as a substitute for the final review log.

### 6. Optional exploratory biological model

After agreeing on endpoints, trial completeness and covariates with Dr. Badea:

```bash
python -m winding_pipeline model --subject-day reviewed_summary/subject_day.csv --session-id YOUR_SESSION --stage acquisition --outcome awn --covariates age_months sex genotype cohort --out learning_model
```

Include only prespecified, measured, varying covariates. Constant or missing
covariates trigger errors. This fits `AWN ~ group * categorical_day + covariates`
with a mouse random intercept, after trial averaging within mouse/day/stage.
It is one working model, not an automatic final analysis. Random slopes, litter,
matched sets, session nesting and other design terms are not implemented.

Review residuals, sample sizes, confounding and `model_status.json`. Rank-deficient
designs are blocked. Warnings/nonconvergence/near-zero random-effect variance
suppress p-values and confidence intervals. No fallback to ordinary independent
OLS is made. Asymptotic Wald inference can be weak for small cohorts. P-values
are explicitly UNADJUSTED; multiple endpoints/contrasts require a prespecified
testing family and correction. `--log1p` is an explicitly requested sensitivity
model, not an outcome-selected transformation. The command expects repeated days;
one probe endpoint needs a separate subject-level model agreed with your mentor.

### 7. Prepare eventual diffusion MRI integration

Use `examples/crosswalk_template.csv` to map behavioral IDs to MRI IDs such as
`B25122902`. Never infer this link from similar dates or row order. The MRI CSV
must have one row per `imaging_id`; select scan/session explicitly first.

```bash
python -m winding_pipeline join-imaging --subject-day reviewed_summary/subject_day.csv --crosswalk verified_crosswalk.csv --imaging diffusion_metrics.csv --session-id YOUR_SESSION --stage probe --day 5 --out behavior_mri_join
```

This selects ONE declared behavioral endpoint per mouse, checks unique IDs,
preserves unmatched records, and prefixes imaging columns with `mri_`. It does
not automatically compute correlations or call them causal effects. Review
imaging masks, metric provenance and scan protocol differences before analysis.

## Output guide

- `trial_metrics.csv`: every trial, QC, raw AWN, optional sensitivity AWN, path,
  duration, timing and sample-retention information, input hash and error details.
- `diagnostics/`: original row positions, scaled XY/time, accepted edges and turns.
- `figures/` and `qc_report.html`: headless visual review; discarded edges are not
  drawn as trajectory lines. Figures show pre-duplicate-check per-trial status;
  the final table/report is authoritative for duplicate-content flags.
- `subject_day.csv`: trial-average metrics per mouse/session/stage/day, total and
  included trial counts. These are NOT pooled across acquisition/probe stages.
- `group_day_descriptive.csv`: group means and SEM across mice, not across trials.
- `summary_figures/`: descriptive mouse trajectories and group means ± SEM;
  see the matching counts in the descriptive table. Single-mouse SEM stays undefined.
- `analysis_inclusion.csv`: every trial and its current inclusion decision.
- `provenance.json`: config, versions, manifest/code hashes and QC counts.

`awn_per_moving_s` and `awn_per_moving_unit` use accepted moving time/path
denominators; they are descriptive alternatives, not automatic confound removal.
Spatial resampling has its own heading support and can omit sharp details.
The tool provides no platform-centered winding, entropy, thigmotaxis or validated
memory score. An AWN difference alone is not evidence of cognitive improvement.

## Preserve and reproduce

Keep original exports, manifest, config, review log, output directory, software
version and environment together. On the cluster, record the resolved environment:

```bash
python -m pip freeze > installed_environment.txt
```

Exit code 2 means a validation error or at least one failed/error trial (outputs
may still exist). Review flags alone do not force nonzero exit. Source files are
opened read-only; output files are confined to the newly created run directory.

See `docs/AUDIT.md`, `docs/VALIDATION.md`, and `docs/BIOLOGICAL_ANALYSIS_PLAN.md`
for the detailed reasoning and remaining limitations.

`synthetic_example/` contains the fully generated practice inputs and output
from this build. Its resolved file paths/provenance identify the build workspace;
generate a fresh demo on your cluster for a directly runnable local example.
