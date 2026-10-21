Leading empty X/Y row fix

This patch ignores only consecutive leading rows with missing X/Y before the
first measured coordinate. Those rows are counted in
n_ignored_leading_xy_rows. Missing X/Y rows after tracking begins remain in the
trajectory as gaps and trigger invalid_xy_rows.

Apply from the root of your existing awn_pipeline folder:
  unzip -o leading_xy_startup_fix.zip

Then rerun the analysis into a NEW output folder. Existing results are not
modified. This change does not set coordinate units or the time-gap threshold;
those remain review flags until confirmed from lab/ANY-maze documentation.

Validation in this environment: 49 unit tests passed, 2 optional tests skipped.
On uploaded sample 260622_01_651.csv: 1 leading row ignored, 0 invalid XY rows,
AWN unchanged at 15.6131865658; coordinate_unit_unknown and
gap_threshold_unset remain review flags.
