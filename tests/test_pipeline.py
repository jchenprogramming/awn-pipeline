import ast
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
from winding_pipeline.core import Settings, calculate, matlab_clean_reference, legacy_python_reference, turn_angles, spatial_resample
from winding_pipeline.io import parse_time, read_trajectory, read_manifest, sha256
from winding_pipeline.pipeline import extract, select_trials, summaries
from winding_pipeline.analysis import join_imaging, prepare_model_data

CFG = Settings(coordinate_unit="pixels", max_time_gap_s=2)


def calc(points, times=None, cfg=CFG):
    p = np.asarray(points, float).reshape(-1, 2)
    t = np.arange(len(p), dtype=float) if times is None else times
    return calculate(t, p, cfg)[0]


class GeometryTests(unittest.TestCase):
    def test_straight(self):
        self.assertEqual(calc([[0, 0], [1, 0], [2, 0]])["awn"], 0)

    def test_right_angle(self):
        self.assertAlmostEqual(calc([[0, 0], [1, 0], [1, 1]])["awn"], .25)

    def test_reversal(self):
        self.assertAlmostEqual(calc([[0, 0], [1, 0], [0, 0]])["awn"], .5)

    def test_open_square_not_closed(self):
        self.assertAlmostEqual(calc([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]])["awn"], .75)

    def test_circle_finite_sample(self):
        a = np.linspace(0, 2 * np.pi, 101)
        self.assertAlmostEqual(calc(np.column_stack([np.cos(a), np.sin(a)]))["awn"], .99)

    def test_clean_matlab_agreement(self):
        rng = np.random.default_rng(541)
        for n in (3, 10, 100, 500):
            for _ in range(10):
                p = np.cumsum(rng.normal(size=(n, 2)), axis=0)
                self.assertAlmostEqual(calc(p)["awn"], matlab_clean_reference(p), places=10)

    def test_rigid_scale_reversal_invariance(self):
        p = np.cumsum(np.random.default_rng(3).normal(size=(40, 2)), axis=0)
        w = calc(p)["awn"]
        for q in (p + [100, -400], p * 15, p @ np.array([[0, -1], [1, 0]]), p * [1, -1], p[::-1]):
            self.assertAlmostEqual(w, calc(q)["awn"], places=10)

    def test_stationary_is_undefined_not_zero(self):
        m = calc(np.zeros((10, 2)))
        self.assertTrue(np.isnan(m["awn"]))
        self.assertEqual(m["qc_status"], "fail")
        self.assertEqual(m["awn_legacy_python"], 0)

    def test_empty_short(self):
        for p in ([], [[0, 0]], [[0, 0], [1, 0]]):
            self.assertTrue(np.isnan(calc(p)["awn"]))

    def test_missing_breaks_not_bridges(self):
        p = [[0, 0], [1, 0], [np.nan, np.nan], [1, 1], [2, 1], [2, 2]]
        m = calc(p)
        self.assertAlmostEqual(m["awn"], .25)
        self.assertGreater(m["awn_legacy_python"], m["awn"])
        self.assertEqual(m["n_turns"], 1)
        self.assertIn("invalid_xy_rows", m["qc_flags"])

    def test_leading_missing_xy_rows_are_ignored_but_counted(self):
        points = [[np.nan, np.nan], [np.nan, np.nan],
                  [0, 0], [1, 0], [1, 1], [2, 1]]
        metrics, diag = calculate([0, .1, .348, .448, .548, .648], points, CFG)
        self.assertNotIn("invalid_xy_rows", metrics["qc_flags"])
        self.assertEqual(metrics["n_ignored_leading_xy_rows"], 2)
        self.assertEqual(metrics["n_rows"], 6)
        self.assertEqual(metrics["n_analyzed_rows"], 4)
        self.assertEqual(metrics["n_invalid_xy"], 0)
        self.assertAlmostEqual(metrics["awn"], .5)
        self.assertEqual(diag["row"].tolist(), [3, 4, 5, 6])

    def test_all_missing_xy_rows_are_not_silently_ignored(self):
        metrics = calc([[np.nan, np.nan], [np.nan, np.nan]])
        self.assertEqual(metrics["n_ignored_leading_xy_rows"], 0)
        self.assertEqual(metrics["n_invalid_xy"], 2)
        self.assertIn("invalid_xy_rows", metrics["qc_flags"])
        self.assertIn("no_usable_turns", metrics["qc_flags"])

    def test_valid_zero_time_coordinate_is_retained(self):
        metrics = calc([[0, 0], [1, 0], [1, 1]], [0, .348, .448])
        self.assertEqual(metrics["n_ignored_leading_xy_rows"], 0)
        self.assertEqual(metrics["n_analyzed_rows"], 3)
        self.assertAlmostEqual(metrics["awn"], .25)

    def test_pause_breaks_turn(self):
        m = calc([[0, 0], [1, 0], [1, 0], [1, 1], [2, 1]])
        self.assertEqual(m["n_turns"], 1)
        self.assertAlmostEqual(m["awn"], .25)
        self.assertEqual(m["qc_status"], "review")

    def test_time_gap(self):
        m = calc([[0, 0], [1, 0], [1, 1], [2, 1], [2, 2]], [0, 1, 9, 10, 11])
        self.assertEqual(m["n_time_gaps"], 1)
        self.assertAlmostEqual(m["awn"], .25)
        self.assertEqual(m["observed_s"], 3)

    def test_nonmonotone_fails(self):
        for t in ([0, 1, 1, 2], [0, 2, 1, 3], [2, np.nan, 1, 3]):
            self.assertEqual(calc([[0, 0], [1, 0], [1, 1], [2, 1]], t)["qc_status"], "fail")

    def test_repeated_time_on_missing_xy_event_rows_does_not_fail(self):
        points = [[np.nan, np.nan], [0, 0], [np.nan, np.nan],
                  [np.nan, np.nan], [1, 0], [1, 1], [2, 1]]
        times = [0, 1, 1, 1, 2, 3, 4]
        metrics = calc(points, times)
        self.assertEqual(metrics["qc_status"], "review")
        self.assertNotIn("non_increasing_time", metrics["qc_flags"])
        self.assertIn("invalid_xy_rows", metrics["qc_flags"])
        self.assertEqual(metrics["n_turns"], 1)
        self.assertAlmostEqual(metrics["awn"], .25)

    def test_speed_jump(self):
        m = calc([[0, 0], [1, 0], [100, 100], [2, 0], [3, 0], [4, 0]], cfg=Settings(coordinate_unit="pixels", max_time_gap_s=2, max_speed=5))
        self.assertEqual(m["n_speed_jumps"], 2)
        self.assertEqual(m["awn"], 0)

    def test_threshold_no_bridging(self):
        m = calc([[0, 0], [1, 0], [1.01, 0], [1.01, 1]], cfg=Settings(coordinate_unit="pixels", max_time_gap_s=2, min_step=.1))
        self.assertTrue(np.isnan(m["awn"]))

    def test_resampling_separate(self):
        cfg = Settings(coordinate_unit="pixels", max_time_gap_s=2, resample_step=.5)
        m = calc([[0, 0], [1, 0], [1, 1]], cfg=cfg)
        self.assertAlmostEqual(m["awn_spatial_resampled"], .25)
        self.assertEqual(m["awn"], .25)

    def test_resample_gap_not_connected(self):
        cfg = Settings(coordinate_unit="pixels", max_time_gap_s=2, resample_step=.2)
        m = calc([[0, 0], [1, 0], [np.nan, np.nan], [1, 2], [1, 3]], cfg=cfg)
        self.assertEqual(m["awn_spatial_resampled"], 0)
        self.assertTrue(np.isnan(m["awn"]))

    def test_polar_is_different(self):
        p = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]])
        angles = np.unwrap(np.arctan2(p[:, 1] - .5, p[:, 0] - .5))
        polar = np.abs(np.diff(angles)).sum() / (2 * np.pi)
        self.assertAlmostEqual(polar, 1)
        self.assertNotAlmostEqual(polar, calc(p)["awn"])

    def test_unknown_settings_flagged(self):
        self.assertEqual(calc([[0, 0], [1, 0], [1, 1]], cfg=Settings())["qc_status"], "review")

    def test_invalid_config(self):
        for kwargs in ({"coordinate_scale": 0}, {"max_time_gap_s": -1}, {"min_step": -1}, {"max_speed": float("nan")}, {"x_column": "Time"}):
            with self.assertRaises(ValueError):
                Settings(**kwargs)

    def test_original_python_parity_if_source_available(self):
        root = os.environ.get("AWN_ORIGINAL_SOURCE")
        if not root:
            self.skipTest("Set AWN_ORIGINAL_SOURCE to the supplied legacy source folder.")
        source = Path(root, "winding_calc_051226.py").read_text()
        tree = ast.parse(source)
        func = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "abs_winding_from_xy")
        # Isolate this pure function only; NEVER execute legacy top-level writes.
        scope = {"np": np, "pd": pd}
        exec(compile(ast.Module(body=[func], type_ignores=[]), "original_function", "exec"), scope)
        for p in (np.zeros((4, 2)), np.array([[0, 0], [1, 0], [np.nan, 0], [1, 1], [2, 1]]), np.cumsum(np.random.default_rng(2).normal(size=(200, 2)), axis=0)):
            expected = scope["abs_winding_from_xy"](pd.Series(p[:, 0]), pd.Series(p[:, 1]))
            self.assertAlmostEqual(legacy_python_reference(p), expected, places=12)


class IOAndBatchTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.raw = self.root / "track.csv"
        self.raw.write_text("Time,Centre position X,Centre position Y\n0,0,0\n1,1,0\n2,1,1\n")
        self.manifest = self.root / "manifest.csv"
        self.rows = [dict(trial_id="t1", subject_id="001", session_id="s1", stage="acquisition", day=1, trial=1, trajectory_file="track.csv", group="control")]
        self.save_manifest()

    def save_manifest(self):
        pd.DataFrame(self.rows).to_csv(self.manifest, index=False)

    def tearDown(self):
        self.tmp.cleanup()

    def test_time_formats(self):
        for text, expected in (("1.5", 1.5), ("01:02.5", 62.5), ("01:02:03.5", 3723.5), ("00:00:00", 0)):
            self.assertEqual(parse_time(text), expected)
        for text in ("", "NA", "12:65", "-1", "inf", "2026-01-01", "1.2:03"):
            self.assertTrue(np.isnan(parse_time(text)))

    def test_csv_decimal_and_preamble(self):
        self.raw.write_text('experiment\nTime;Centre position X;Centre position Y\n00:00:00;0;0\n00:00:01;1,5;0\n00:00:02;1,5;2\n')
        t, p = read_trajectory(self.raw, Settings(skiprows=1, delimiter=";", decimal=","))
        self.assertEqual(p[1, 0], 1.5)
        self.assertEqual(t[2], 2)

    def test_no_column_guessing(self):
        self.raw.write_text("Time,Speed,Event\n0,2,1\n")
        with self.assertRaises(ValueError):
            read_trajectory(self.raw, CFG)

    def test_duplicate_headers_rejected(self):
        self.raw.write_text("Time,Time,Y\n0,2,1\n")
        with self.assertRaises(ValueError):
            read_trajectory(self.raw, CFG)

    def test_ragged_rows_rejected(self):
        self.raw.write_text("Time,Centre position X,Centre position Y\n0,2,1,5\n")
        with self.assertRaises(ValueError):
            read_trajectory(self.raw, CFG)

    def test_malformed_csv_preserved_as_error(self):
        self.raw.write_text('Time,Centre position X,Centre position Y\n0,0,"unterminated\n')
        result = extract(self.manifest, CFG, self.root / "run", plots=False)
        self.assertEqual(result.qc_status.iloc[0], "error")

    def test_blank_row_preserved(self):
        self.raw.write_text("Time,Centre position X,Centre position Y\n0,0,0\n\n2,1,1\n")
        t, p = read_trajectory(self.raw, CFG)
        self.assertEqual(len(t), 3)
        self.assertTrue(np.isnan(p[1]).all())

    def test_manifest_leading_zero_id(self):
        self.assertEqual(read_manifest(self.manifest).subject_id.iloc[0], "001")

    def test_manifest_duplicate_trial(self):
        self.rows.append(dict(self.rows[0]))
        self.save_manifest()
        with self.assertRaises(ValueError):
            read_manifest(self.manifest)

    def test_manifest_duplicate_file(self):
        self.rows.append(dict(self.rows[0], trial_id="t2", trial=2))
        self.save_manifest()
        with self.assertRaises(ValueError):
            read_manifest(self.manifest)

    def test_metadata_conflict(self):
        self.rows.append(dict(self.rows[0], trial_id="t2", trial=2, trajectory_file="other.csv", group="treated"))
        self.save_manifest()
        with self.assertRaises(ValueError):
            read_manifest(self.manifest)

    def test_batch_and_input_unchanged(self):
        before = sha256(self.raw)
        out = self.root / "run"
        trials = extract(self.manifest, CFG, out, plots=True)
        self.assertEqual(before, sha256(self.raw))
        self.assertEqual(trials.awn.iloc[0], .25)
        self.assertTrue((out / "figures/trial_00001.png").exists())
        self.assertTrue((out / "qc_report.html").exists())
        self.assertEqual(pd.read_csv(out / "subject_day.csv").n_trials_included.iloc[0], 1)
        with self.assertRaises(FileExistsError):
            extract(self.manifest, CFG, out, plots=False)

    def test_missing_file_preserved(self):
        self.rows.append(dict(self.rows[0], trial_id="t2", subject_id="002", trajectory_file="missing.csv"))
        self.save_manifest()
        result = extract(self.manifest, CFG, self.root / "run", plots=False)
        self.assertEqual(len(result), 2)
        self.assertEqual(result.qc_status.tolist(), ["pass", "error"])
        days = pd.read_csv(self.root / "run/subject_day.csv")
        self.assertEqual(days.n_trials_included.tolist(), [1, 0])

    def test_all_error_does_not_crash_summary(self):
        self.raw.unlink()
        result = extract(self.manifest, CFG, self.root / "run", plots=False)
        self.assertEqual(result.qc_status.iloc[0], "error")

    def test_duplicate_contents_review(self):
        other = self.root / "other.csv"
        other.write_bytes(self.raw.read_bytes())
        self.rows.append(dict(self.rows[0], trial_id="t2", subject_id="002", trajectory_file="other.csv"))
        self.save_manifest()
        result = extract(self.manifest, CFG, self.root / "run", plots=False)
        self.assertTrue(result.qc_status.eq("review").all())

    def test_review_decisions(self):
        result = extract(self.manifest, Settings(), self.root / "run", plots=False)
        self.assertTrue(select_trials(result).empty)
        self.assertEqual(len(select_trials(result, include_review=True)), 1)
        review = pd.DataFrame([dict(trial_id="t1", decision="include", reviewer="J", reason="Verified export and track")])
        self.assertEqual(len(select_trials(result, review=review)), 1)
        result["qc_status"] = "fail"
        with self.assertRaises(ValueError):
            select_trials(result, review=review)

    def test_stage_separation(self):
        self.raw2 = self.root / "probe.csv"
        self.raw2.write_text("Time,Centre position X,Centre position Y\n0,0,0\n1,1,0\n2,2,0\n")
        self.rows.append(dict(self.rows[0], trial_id="t2", stage="probe", trajectory_file="probe.csv"))
        self.save_manifest()
        extract(self.manifest, CFG, self.root / "run", plots=False)
        self.assertEqual(len(pd.read_csv(self.root / "run/subject_day.csv")), 2)

    def test_join_and_missing_crosswalk(self):
        pd.DataFrame([dict(subject_id=x, session_id="s1", stage="probe", day=5, awn=.25, n_trials_included=1) for x in ["001", "002", "003"]]).to_csv(self.root / "day.csv", index=False)
        (self.root / "cross.csv").write_text("subject_id,imaging_id\n001,B001\n")
        (self.root / "mri.csv").write_text("imaging_id,FA\nB001,.4\nB004,.5\n")
        result = join_imaging(self.root / "day.csv", self.root / "cross.csv", self.root / "mri.csv", self.root / "joined", "s1", "probe", 5)
        self.assertEqual(len(result), 3)
        self.assertEqual(result.imaging_match.eq("both").sum(), 1)

    def test_join_duplicate_imaging_rejected(self):
        pd.DataFrame([dict(subject_id="001", session_id="s1", stage="probe", day=5, awn=.25, n_trials_included=1)]).to_csv(self.root / "day.csv", index=False)
        (self.root / "cross.csv").write_text("subject_id,imaging_id\n001,B001\n")
        (self.root / "mri.csv").write_text("imaging_id,FA\nB001,.4\nB001,.5\n")
        with self.assertRaises(ValueError):
            join_imaging(self.root / "day.csv", self.root / "cross.csv", self.root / "mri.csv", self.root / "joined", "s1", "probe", 5)


class AnalysisTests(unittest.TestCase):
    def data(self):
        return pd.DataFrame([dict(subject_id=f"s{i}", day=d, group="A" if i < 6 else "B", awn=i + d + .1, age_months=10 + i) for i in range(12) for d in [1, 2, 3]])

    def test_prepare_model(self):
        d, excluded = prepare_model_data(self.data(), "awn", ["age_months"], True)
        self.assertEqual(len(d), 36)
        self.assertEqual(len(excluded), 0)
        np.testing.assert_allclose(d.response, np.log1p(d.awn))

    def test_model_missing_covariates_explicit(self):
        d = self.data()
        d.loc[0, "age_months"] = np.nan
        clean, excluded = prepare_model_data(d, "awn", ["age_months"])
        self.assertEqual(len(excluded), 1)
        self.assertEqual(len(clean), 35)

    def test_model_no_group_guess(self):
        with self.assertRaises(ValueError):
            prepare_model_data(self.data().drop(columns="group"), "awn", [])

    def test_model_too_few_subjects(self):
        with self.assertRaises(ValueError):
            prepare_model_data(self.data().iloc[:6], "awn", [])

    def test_optional_mixed_model_integration(self):
        if importlib.util.find_spec("statsmodels") is None:
            self.skipTest("statsmodels not installed; mixed-model fitting unverified in this environment")
        from winding_pipeline.analysis import model
        rng = np.random.default_rng(341)
        d = self.data()
        intercept = dict(zip(d.subject_id.unique(), rng.normal(0, 1, 12)))
        d["awn"] = [10 + intercept[s] - .3 * day + rng.normal(0, .2) for s, day in zip(d.subject_id, d.day)]
        d["session_id"] = "s1"
        d["stage"] = "acquisition"
        d["n_trials_included"] = 4
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "day.csv"
            d.to_csv(path, index=False)
            result = model(path, Path(tmp) / "model", "s1", "acquisition")
            self.assertIn(result["status"], ["fit_ok_exploratory", "review_fit_no_pvalues"])


if __name__ == "__main__":
    unittest.main()
