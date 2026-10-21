# Validation record — v0.1.0

Run date: 2026-10-06. Full logs and machine-readable results are in `validation/`.

## Executed successfully

- 45 unittest cases ran: **44 passed, 1 skipped, 0 errors, 0 failures**.
- The skipped case is the actual optional statsmodels mixed-model fit. The build
  environment lacks statsmodels and the package was not available through its
  configured installer. Input/model-data checks were tested; fitted coefficients,
  inference and convergence behavior are NOT runtime-validated here.
- 100 additional deterministic random trajectories, 500 points each, matched an
  independent clean dot/acos implementation with maximum absolute AWN difference
  **6.821210263296962e-13**. This is formula agreement, not execution in MATLAB.
- The original uploaded Python `abs_winding_from_xy` function was isolated from
  its AST and run on finite, gappy and stationary fixtures. The new
  `awn_legacy_python` diagnostic matched it. No legacy top-level scripts were run.
- Mathematical fixtures: straight line, right angle, reversal, open square,
  sampled circle, translation/rotation/reflection/scale/time reversal invariance.
- Robustness: missing rows, all-stationary data, too-short data, time gaps,
  non-increasing timestamps, jumps, short steps, and segmentation-preserving
  spatial resampling.
- CSV/schema: elapsed-time formats, decimal commas, semicolon delimiter, preamble,
  blank records, ragged records, ambiguous/missing coordinate headers, duplicate
  headers, text/leading-zero IDs, duplicate trials/files and metadata conflicts.
- Batch/QC: input hash unchanged, headless QC image, new-directory-only output,
  missing/error records retained, all-error summary, duplicated content flags,
  review gating and acquisition/probe separation.
- Imaging join: keyed linkage, unmatched records and duplicate-scan rejection.
- End-to-end synthetic practice: 12 fictional mice × 3 days × 2 trials = 72
  trajectories; extraction and subject-day summary completed. These are not
  biological observations or evidence of a treatment effect.

## Environment used for core tests

Python 3.12.14; NumPy 2.3.5; pandas 2.2.3; matplotlib 3.10.8; SciPy 1.17.0.
statsmodels: not installed. See `validation/build_environment.txt` for the full
environment record. Dependency ranges in pyproject.toml are not a claim that all
permitted version combinations were tested.

## Not yet validated

1. Actual MATLAB execution: run `examples/matlab_validation.m` with the original
   `AbsWinding.m`, then compare representative real tracks on identical support.
2. Real ANY-maze export variants, frame/event semantics, empirical thresholds,
   calibration and ID mapping. No raw trajectories were in the supplied archive.
3. Real GLP1/control cohort results, treatment labels or associations with MRI.
4. Whether the chosen pause/gap policy best answers the lab's research question.
5. Optional statistical fitting in a fully provisioned environment and review of
   its assumptions by the lab. Its integration test is included for rerunning.
6. Large-cohort performance, all encodings/export variants and Windows-specific
   installation. Individual CSVs are currently read into memory.

## How to repeat tests

From the package root after installing dependencies:

```bash
AWN_ORIGINAL_SOURCE="$PWD/legacy_reference" python -m unittest discover -s tests -v
```

The release validation builder writes a NEW `validation/` folder and deliberately
refuses to overwrite an existing one. To rerun the builder, use a fresh checkout
or preserve the previous validation directory under a different name first.

## Interpretation boundary

Passing software tests demonstrates tested behavior, not validity as a cognitive
biomarker. Do not use these outputs as final biological evidence until real-data
validation, protocol review, inclusion criteria and study-design checks are done.
