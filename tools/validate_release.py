"""Run release tests and document reproducible geometry evidence.

Run from the package root: python tools/validate_release.py
This writes only validation/ (new directory required), never raw study inputs.
"""
from pathlib import Path
import io
import json
import os
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from winding_pipeline.core import Settings, calculate, matlab_clean_reference
from winding_pipeline.io import versions, write_json


def main():
    output = ROOT / "validation"
    output.mkdir(exist_ok=False)
    os.environ["AWN_ORIGINAL_SOURCE"] = str(ROOT / "legacy_reference")
    log = io.StringIO()
    suite = unittest.defaultTestLoader.discover(str(ROOT / "tests"))
    result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
    (output / "test_results.txt").write_text(log.getvalue())
    print(log.getvalue())
    rng = np.random.default_rng(712)
    errors = []
    cfg = Settings(coordinate_unit="pixels", max_time_gap_s=2)
    for _ in range(100):
        p = np.cumsum(rng.normal(size=(500, 2)), axis=0)
        w, _ = calculate(np.arange(len(p)), p, cfg)
        errors.append(abs(w["awn"] - matlab_clean_reference(p)))
    evidence = {"versions": versions(), "tests_run": result.testsRun,
                "tests_passed": result.testsRun - len(result.failures) - len(result.errors) - len(result.skipped),
                "failures": len(result.failures), "errors": len(result.errors),
                "skipped": [dict(test=str(t), reason=r) for t, r in result.skipped],
                "random_reference_trials": 100, "samples_per_random_trial": 500,
                "maximum_absolute_reference_difference": max(errors),
                "MATLAB_executed": False, "real_trajectories_tested": False}
    write_json(output / "validation_results.json", evidence)
    fixtures = {
        "straight": [[0, 0], [1, 0], [2, 0]],
        "right_angle": [[0, 0], [1, 0], [1, 1]],
        "square_open": [[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]],
        "missing_row": [[0, 0], [1, 0], [float('nan'), float('nan')], [1, 1], [2, 1], [2, 2]],
        "stationary": [[0, 0], [0, 0], [0, 0]],
    }
    import pandas as pd
    records = []
    for name, p in fixtures.items():
        metrics, _ = calculate(np.arange(len(p)), p, cfg)
        records.append(dict(fixture=name, **metrics))
    pd.DataFrame(records).to_csv(output / "fixture_results.csv", index=False)
    from winding_pipeline.pipeline import plot_trial
    p = np.array(fixtures["missing_row"])
    m, diag = calculate(np.arange(len(p)), p, cfg)
    plot_trial(diag, output / "synthetic_gap_qc.png", "SYNTHETIC missing-row fixture | observed AWN=0.25", cfg)
    # Capture a matching runtime package list without requiring pip.
    from importlib.metadata import distributions
    packages = sorted({f"{d.metadata['Name']}=={d.version}" for d in distributions() if d.metadata.get('Name')})
    (output / "build_environment.txt").write_text("\n".join(packages) + "\n")
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
