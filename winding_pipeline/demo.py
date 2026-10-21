"""Deterministic synthetic data: never report these as mouse results."""
import numpy as np
import pandas as pd
from .io import new_output, write_json
from .core import Settings


def make_demo(output):
    out = new_output(output)
    (out / "trajectories").mkdir()
    rng = np.random.default_rng(847)
    rows = []
    for mouse in range(12):
        offset = rng.normal(0, .08)
        for day in range(1, 4):
            for trial in range(1, 3):
                key = f"SYNTH_{mouse:02d}_D{day}_T{trial}"
                t = np.arange(151) / 10
                heading = np.cumsum(rng.normal(0, .10 + .02 * (3 - day), len(t))) + offset
                xy = np.cumsum(np.column_stack([np.cos(heading), np.sin(heading)]), axis=0)
                pd.DataFrame({"Time": t, "Centre position X": xy[:, 0], "Centre position Y": xy[:, 1]}).to_csv(out / "trajectories" / f"{key}.csv", index=False)
                rows.append(dict(trial_id=key, subject_id=f"SYNTH_{mouse:02d}", session_id="SYNTHETIC",
                    stage="acquisition", day=day, trial=trial, trajectory_file=f"trajectories/{key}.csv",
                    group="A" if mouse < 6 else "B", cohort="SYNTHETIC", sex="F" if mouse % 2 else "M", age_months=10 + mouse / 10))
    pd.DataFrame(rows).to_csv(out / "manifest.csv", index=False)
    write_json(out / "config.json", Settings(coordinate_unit="pixels", max_time_gap_s=.5, max_speed=50, resample_step=2).to_dict())
    (out / "SYNTHETIC_ONLY.txt").write_text("These are generated trajectories, not experimental mouse data.\n")

