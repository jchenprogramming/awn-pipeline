# From verified trajectories to biological analysis

## Milestone 1 — export and identity

Ask BJ for representative per-test raw time/centre-X/centre-Y CSVs, their ANY-maze
test schedule, export settings, and a few original MATLAB AWN values. Obtain the
master animal/group/cohort table. Retain raw inputs unchanged. Confirm how blank
rows, event records, repeated timestamps and untracked frames are encoded.

Deliverable: verified trial manifest and an input dictionary. Do not reuse
the legacy round-robin mapping without checking Test IDs against source records.

## Milestone 2 — computation validation

Start with a clean trial, stationary stretches, a gappy trial, a probe and an
unusually high-AWN trajectory. Compare BOTH the formula and the preprocessing
support (same rows, units and axes). Run `examples/matlab_validation.m` using the
original function and retain MATLAB version and outputs.

Use `awn_legacy_python` only to diagnose differences from the historical Python
policy. `awn_matlab_clean_reference` is a mathematical reference on clean data,
not proof that MATLAB was run. Reproducing a historical handling error is not the
same as selecting a sound primary estimator.

Deliverable: per-trial validation table, discrepancy explanations, reviewed plots.

## Milestone 3 — prespecified QC and sensitivity

Agree with Dr. Badea on whether conservative breaks at pauses/gaps match the
intended endpoint; decide permitted missingness and minimum temporal/spatial
support BEFORE comparing groups. Inspect export cadence, tracking video, jumps,
retained-turn fraction, unit calibration and missing-data patterns by cohort.

Run separately named settings for plausible gap/short-step/resampling choices.
Changing `min_step` breaks turns rather than reconstructing a noise-free path.
Spatial resampling creates a different scale-specific estimator. Never choose
the setting with the smallest treatment p-value. Keep settings and exclusions
in the record even when the numerical result appears unfavorable.

Deliverable: trial-by-trial inclusion log and sensitivity comparison, plus sample
counts/missingness by group/session/day/stage. For unequal included trials within
a mouse/day, specify the minimum acceptable completeness; means from one versus
four trials need not be equally reliable. The pipeline records counts but does
not choose a universal biological completeness cutoff.

## Milestone 4 — mouse-level biological endpoints

Separate acquisition, reversal and probes. Decide a primary endpoint and a
small, justified set of secondary outcomes before testing. Possibilities for
discussion, not final choices: acquisition AWN trajectory, a fixed-day probe AWN,
and prespecified within-mouse probe change. Duration, path length and normalized
AWN can be descriptive companion measures. Whether to adjust for speed/duration
depends on the scientific question; they may be treatment-mediated outcomes.

Use mice as independent biological units, not frames or trials. A mouse-level
random-intercept model is a starting point for repeated days, not an automatic
final model. Inspect whether random slopes, litter, matched set, session, cage or
other dependence needs representation. These extensions are not in v0.1.

Check group overlap in age/sex/genotype/cohort. If all treated mice are in one
cohort and all controls another, the treatment contrast is confounded; no
regression can separate a perfectly confounded design. Historical comparisons
should be clearly labeled exploratory, with acquisition/method differences stated.

Deliverable: prespecified analysis plan, reviewed model fit, effect estimates and
intervals, group sample sizes, residual diagnostics and a declared multiplicity
strategy. Do not treat many nominal p<0.05 values as independent discoveries.

## Milestone 5 — behavior and imaging

Create a verified behavioral-ID to MRI-ID crosswalk. Choose one behavioral endpoint
per mouse and one specified scan/session before joining. Do not multiply MRI rows
by joining every trial as an independent record. Keep unmatched IDs visible.

Verify MRI metric provenance and masks. Matching image dimensions/header codes
alone is not proof of anatomical registration quality. Whole-image nonzero means
can have different voxel support between FA/MD/AD/RD, particularly after clipping;
they need not be clean common-mask biological endpoints. Use reviewed common
support or clearly defined regional metrics as appropriate to the study.

Preselect biologically justified imaging regions/metrics, or define correction
over the complete planned set. Inspect cohort, genotype, sex, age and scan protocol
confounding. Brain–behavior correlations are associations, not causal effects.

Deliverable: analysis-ready, one-row-per-mouse endpoint/MRI table with a crosswalk
audit. The pipeline supplies the safe join; confirmatory correlations and multiple-
comparison procedures should follow an agreed endpoint plan rather than an
automatic sweep of hundreds of regions.
