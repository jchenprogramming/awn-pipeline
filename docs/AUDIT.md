# Audit of the uploaded MATLAB / Python winding-number code

Date: 2026-10-06. Scope: the actual `awn-main(1).zip` supplied in this conversation.
The core calculator, MATLAB wrapper, sanity script, file renamers, plotting code,
representative statistical functions across Python variants, and the R learning
models were inspected. This was a source review, not execution of all original
scripts: those scripts contain hard-coded paths and top-level file writes.

## Executive conclusion

The central MATLAB geometry is a reasonable implementation of **total absolute
tangent turning** for clean, sampled trajectories. The primary Python calculator
largely reproduces that idea, but input mapping and missing-data policies can
change results without sufficient warning. A second Python script calculates a
different metric. The existing statistical scripts provide useful exploratory
plots but are not a single, sufficiently validated end-to-end biological pipeline.

Build a new, isolated package rather than overwriting the legacy scripts. Keep
historical outputs and methods identifiable. Do not promise exact real-data
MATLAB parity until representative paired files have been checked.

## 1. What the MATLAB AWN actually measures

For points p_i=(x_i,y_i), let v_i=p_(i+1)-p_i. The function computes

    AWN = sum_i acos(clip(dot(v_i,v_(i+1))/(|v_i| |v_(i+1)|), -1, 1)) / (2*pi)

There are N-2 angle opportunities for N samples. Time is ignored by this function.
Each angle is between 0 and pi. The result is nonnegative, generally not an
integer, and includes left and right turns without cancellation. It measures
how much the travel direction changes; it is NOT the number of times the animal
encircles a platform or arena center, nor the animal's body rotation.

The new implementation uses atan2(abs(cross), dot) on normalized adjacent
displacements. This is mathematically equivalent for clean moving data and
avoids some acos roundoff sensitivity near straight motion. Very small
floating-point differences are expected, so comparisons use tolerances.

| Clean polyline | Original-style AWN | Reason |
|---|---:|---|
| Three collinear points going forward | 0 | One zero-degree turn |
| Three points forming a right angle | 0.25 | 90 degrees / 360 |
| Out and back | 0.5 | 180-degree reversal |
| Square with start repeated once at the end | 0.75 | Three interior turns; no cyclic heading closure |
| 100-segment sampled circle, including repeated endpoint | 0.99 | 99 measured interior heading changes |

Translation, rotation, reflection and positive uniform scaling do not change
this geometric AWN. Changing sample cadence, tracking noise, gaps, smoothing,
axis calibration nonuniformly, or trial duration CAN change the result. Reflection
invariance means a downward image Y axis is fine for absolute turning.

## 2. Findings and their consequences

### Critical — two different meanings of winding

`winding_sanity_check.py:86–105` uses atan2 of positions relative to a center,
unwraps polar angles, and sums polar-angle changes. `AbsWinding.m:16–23` uses
changes between displacement directions. The sanity script's default center is
the trajectory mean (`CENTER_METHOD`, and the processing loop), so its reference
point also changes with each path. It is neither the measured pool center nor a
verified platform coordinate. Its `Net_turning_degrees` is polar, not tangent,
turning. For a closed square around its center, polar absolute winding is 1,
but MATLAB AWN on five corner samples is 0.75.

**Action:** do not mix or compare these columns as equivalent endpoints. The
replacement implements tangent AWN only. Platform-centered winding would require
its own name, fixed reference coordinates, center-crossing handling and validation.

### Critical — animal/trial mapping can silently misassign data

`winding_calc_051226.py` hard-codes 27 animals and calculates file suffix as
`position + 27 * within_animal_counter`. The counter is based on metadata row
order, not a verified Test ID/day/trial join. If records are shuffled, missing or
duplicated, the link can fail or target the wrong existing file. That can attach
a valid-looking number to the wrong behavioral record.

Both `filerenamer*.sh` scripts rename originals in place using a CSV row counter;
they do not use a real CSV parser for quoted commas, do not guard destination
collisions, and assume the table/file sequence correspondence. The sanity script
also extracts `44` from a cohort-like directory name and calls it Animal; this
cannot establish individual mouse identity and can pool different mice.

**Action:** explicit unique manifest keys, exact filenames, no destructive
renaming, no round-robin inference, duplicate path/content flags and input hashes.
Human verification against the export schedule is still necessary.

### High — missing samples become invented connections

`winding_calc_051226.py:38–44` drops nonfinite XY rows before differencing. Both
polar winding and path length in the sanity script do the same. This draws a
straight connection between points on opposite sides of a tracking gap and
creates headings that were never observed. Long timestamp gaps are also unchecked.

**Action:** preserve row order and break turns at invalid rows/gaps. Report
retention and observed duration. Do not impute unknown turns. A segmented sum
has reduced observational support; it is NOT a repaired whole-trial value.
Event-only missing rows need export-format review rather than automatic deletion.

### High — stationary samples and jitter need explicit policies

`AbsWinding.m` divides by the product of two step lengths without checking zero.
The NaN behavior of the subsequent MATLAB min/max operations should be checked
in the lab's MATLAB release; the code does not explicitly define its intended
handling. The primary Python calculator skips zero-length adjacent pairs, which
is a policy choice rather than demonstrated MATLAB parity for degenerate inputs.
An entirely stationary valid trajectory becomes zero in that Python function,
which could be confused with valid straight swimming.

**Action:** zero steps break heading pairs; no usable turns produces missing AWN
and failed QC. A valid straight moving trajectory still produces zero. Optional
step thresholds and spatial resampling are sensitivity settings, not silent
changes to the primary endpoint. Raw jitter can still inflate AWN, even when
numerical QC passes. Inspect videos/plots and calibration before group analysis.

### High — ambiguous coordinate columns and time parsing

The preferred XY aliases in `winding_calc_051226.py:80–100` do not include the
exact `Centre position X`/`Centre position Y` headers used by the sanity script.
Fallback chooses columns by numeric dtype and position. Adding a speed/event
column or switching time to numeric can change the selected fields. It may
work accidentally for a three-column export, but is not a safe general parser.

The MATLAB wrapper skips the first data-table row unconditionally. It labels
datevec's six outputs as if they were individual time digits, although datevec
returns year/month/day/hour/minute/second components. This timestamp conversion
is not correct as written, but time is not used in the current AWN function.
The wrapper also assumes two special directory entries, preallocates N-2 after
already subtracting two, processes arbitrary directory entries, and its minimum
nonzero selection indexes a filtered vector without mapping back to original
file indices. Suppressing all warnings makes these problems harder to notice.

**Action:** explicit column names, delimiter, elapsed-time parser, duplicate/ragged
header checks and visible errors. No automatic first-row deletion or sorting.

### High — ordinary OLS on repeated days ignores dependence

`GLP1_Controls_Stats.py:469–490`, `GLP1_Controls_Stats2.py:539–558`,
`winding_matlab_stats.py:436–459`, and `winding_matlab_stats_covars.py:750–768`
fit group-by-day ordinary least squares models after mouse/day averaging.
Averaging trials first is good, but days from the SAME mouse remain dependent.
Adding age, sex and genotype does not account for that correlation.

In contrast, `MWM_Example.R:131,150,183` uses `(1|Animal)` for learning data. The
R approach is closer to the repeated-measures structure. Earlier Python
`GLP1vsCtrlAnalysis*.py` scripts collapse to mouse-level means, which avoids
trial-level pseudoreplication for those endpoints but can mix stages/days and
lose the learning trajectory. These are different analyses, not interchangeable
versions of one validated result.

**Action:** separate acquisition/probe/session; retain trial counts; use a
prespecified repeated-measures model for learning. New optional model uses a
mouse random intercept and categorical day. It still needs design review, fit
diagnostics, and adequate independent mice; it is not universally sufficient.

### High — unknown metadata and cohort confounding

Several scripts map missing treatment to `None`, or classify anything other
than GLP1 as untreated when constructing Group. Unknown is not evidence of no
treatment. Multiple routines pick the first metadata row without checking
conflicting values. `winding_plots.py` converts IDs/sex/genotype to strings before
dropna; missing values can become the literal string `nan`. It then drops
duplicate IDs keeping the first row and joins AnimalID to BadeaID without a
separate verified crosswalk.

Historical-control versus newly treated comparisons can be confounded by
recording cohort, apparatus, export cadence, preprocessing and age. Finding an
age-near control does not remove those differences. Some later matching
functions choose nearest age without a caliper; earlier versions expose a
caliper/sex/genotype options. Matched-set dependence is not then modeled.

**Action:** no inferred controls; preserve text IDs; check metadata consistency;
document group-by-cohort overlap; reject rank-deficient models. Perfect group–
cohort confounding cannot be repaired by simply adding cohort as a covariate.

### Medium/high — unit conversion is hard-coded, not verified

Several scripts multiply old distances by 10. The R example divides distances
by 10, but may refer to differently prepared data, so this is not proof either
conversion is wrong. The ZIP does not include the actual hard-coded GLP1/control
input tables or calibration records needed to establish correctness. Uniform
coordinate scale does not change tangent AWN; it does change path/speed and any
thresholds based on physical distance.

**Action:** verify native units, store explicit conversion and provenance, never
infer conversion from group distribution. Raw exported path length need not
equal ANY-maze's reported distance, which applies movement processing.

### Medium/high — inference needs a defined family and successful fit

Some scripts use Welch tests and BH-FDR, which is useful. However correction is
generally within each function/output (often one metric), not necessarily across
all metrics, days, subsets and contrasts examined. Many scripts suppress all
warnings. Residualizing outcomes on covariates and plotting them is not a
substitute for estimating adjusted group effects in the correct joint model.

**Action:** predeclare primary endpoint and comparisons; keep sample accounting;
define multiplicity across the planned family; display warnings and uncertainty.
New model is explicitly exploratory and does NOT present its unadjusted Wald
p-values as multiplicity-controlled results. No p-values are generated by default
extraction, summaries or MRI linking.

### Engineering — scripts are exploratory rather than reusable packages

Most Python modules execute on import, write to hard-coded `/Volumes/...` paths,
depend on absent local datasets, and mix parsing, analysis and plotting. The
uploaded `glp1_tests.py` is a descriptive-analysis script, not a unit-test suite.
No original automated numerical regression tests are supplied.

**Action:** import-safe modules, CLI, new output directories, explicit configs,
checksums, preserved error rows, test fixtures and headless plots. Original
scripts remain untouched. The new package is a working v0.1, not a certified tool.

## 3. What was and was not validated

Executed: clean mathematical fixtures, transformation invariance, randomized
agreement with an independent dot/acos reference; parity with the original Python
function extracted by AST (without its top-level writes); timing/gap/pause/parser/
identity/output tests; headless plotting and a 72-trial synthetic extraction.

Not executed: MATLAB itself, the original R analyses, real ANY-maze trials,
actual GLP1/control comparisons, or mixed-model fitting in this environment.
See the exact test record in `docs/VALIDATION.md`.

## 4. Why these checks matter biologically

Higher AWN means more total absolute direction change at the analyzed sampling
scale. It does not alone mean worse memory or a pathological search strategy.
A mouse with a longer trial, longer path, more tracking noise, different export
cadence or more slow movement can differ in AWN for non-cognitive reasons.
Normalized AWN/time or AWN/path are alternative estimands, not automatic cures.
In acquisition, escape latency can change observation duration as learning
improves. In a fixed-duration probe, the support is more comparable but speed,
tracking quality, strategy and coverage still matter.

The first scientific output should be a validated trial dataset and QC report,
not a treatment p-value. Only then define learning/probe endpoints and link to
MRI through verified animal IDs and agreed imaging QC.

## Sources checked

1. Supplied source archive, especially `AbsWinding.m`, `extract_data.m`,
   `winding_calc_051226.py`, `winding_sanity_check.py`, `winding_matlab_stats*.py`,
   `GLP1*Stats*.py`, `GLP1vsCtrlAnalysis*.py`, `winding_plots.py`, `MWM_Example.R`.
   Original source is included under `legacy_reference/` for traceability; do
   not run those scripts casually because of hard-coded paths and writes.
2. ANY-maze, “Exporting XY tracking coordinates”:
   https://www.any-maze.com/support/guides/exporting-xy-tracking-coordinates/
   Checked 2026-10-06. Export format/units are configurable; reported movement
   metrics need not equal a simple sum/average over raw exported samples.
3. statsmodels, “Linear Mixed Effects Models”:
   https://www.statsmodels.org/stable/mixed_linear.html
   Checked 2026-10-06. Supports dependent longitudinal observations through
   random effects; does not remove the need for design and fit diagnostics.

The geometric comparisons and code-specific defects above are derived directly
from the supplied implementation and tests, not inferred from those web pages.
