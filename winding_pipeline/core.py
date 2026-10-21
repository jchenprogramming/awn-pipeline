"""Geometry only. Never infer animal identities or close an open trajectory."""
from dataclasses import dataclass, asdict
import numpy as np


@dataclass(frozen=True)
class Settings:
    x_column: str = "Centre position X"
    y_column: str = "Centre position Y"
    time_column: str = "Time"
    delimiter: str = ","
    decimal: str = "."
    encoding: str = "utf-8-sig"
    skiprows: int = 0
    coordinate_unit: str = "unknown"
    coordinate_scale: float = 1.0
    y_axis: str = "down"
    max_time_gap_s: float | None = None
    max_speed: float | None = None
    min_step: float = 0.0
    resample_step: float | None = None

    def __post_init__(self):
        if len(self.delimiter) != 1 or self.decimal not in (".", ","):
            raise ValueError("Use a one-character delimiter and decimal '.' or ','.")
        if self.y_axis not in ("up", "down"):
            raise ValueError("y_axis must be up or down.")
        if self.coordinate_unit not in ("pixels", "mm", "cm", "m", "unknown"):
            raise ValueError("coordinate_unit must be pixels, mm, cm, m, or unknown.")
        if not isinstance(self.skiprows, int) or self.skiprows < 0:
            raise ValueError("skiprows must be a nonnegative integer.")
        for key in ("coordinate_scale", "max_time_gap_s", "max_speed", "resample_step"):
            value = getattr(self, key)
            if value is not None and (not np.isfinite(value) or value <= 0):
                raise ValueError(f"{key} must be finite and positive.")
        if not np.isfinite(self.min_step) or self.min_step < 0:
            raise ValueError("min_step must be finite and nonnegative.")
        if len({self.x_column, self.y_column, self.time_column}) != 3:
            raise ValueError("Time, X and Y must be distinct columns.")

    def to_dict(self):
        return asdict(self)


def turn_angles(points):
    """Nonnegative heading changes on finite, nonzero adjacent displacements."""
    v = np.diff(np.asarray(points, dtype=float), axis=0)
    lengths = np.linalg.norm(v, axis=1)
    if len(v) < 2:
        return np.empty(0)
    if not np.all(np.isfinite(v)) or np.any(lengths <= 0):
        raise ValueError("turn_angles requires finite points and nonzero steps.")
    u = v / lengths[:, None]
    dot = np.einsum("ij,ij->i", u[:-1], u[1:])
    cross = u[:-1, 0] * u[1:, 1] - u[:-1, 1] * u[1:, 0]
    return np.arctan2(np.abs(cross), dot)


def matlab_clean_reference(points):
    """Direct dot/acos translation; only defined here for >=3 clean moving points.

    Does not claim MATLAB runtime equivalence for NaN/zero-step inputs, whose
    min/max NaN behavior can depend on MATLAB version and calling semantics.
    """
    p = np.asarray(points, dtype=float)
    v = np.diff(p, axis=0)
    lengths = np.linalg.norm(v, axis=1)
    if len(p) < 3 or not np.all(np.isfinite(p)) or np.any(lengths == 0):
        return np.nan
    u = v / lengths[:, None]
    cos = np.einsum("ij,ij->i", u[:-1], u[1:])
    return float(np.arccos(np.clip(cos, -1, 1)).sum() / (2 * np.pi))


def legacy_python_reference(points):
    """Replicate winding_calc_051226.py's drop-invalid / skip-zero policy.

    Diagnostic ONLY: this can connect across missing rows and ignores time.
    """
    p = np.asarray(points, dtype=float)
    p = p[np.isfinite(p).all(axis=1)]
    if len(p) < 3:
        return np.nan
    total = 0.0
    for a, b in zip(np.diff(p, axis=0)[:-1], np.diff(p, axis=0)[1:]):
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        if na == 0 or nb == 0:
            continue
        total += np.arccos(np.clip(np.dot(a, b) / (na * nb), -1, 1))
    return float(total / (2 * np.pi))


def edge_runs(edges):
    """Return point-index arrays for contiguous true edges, including 2-point runs."""
    indices = np.flatnonzero(edges)
    if not len(indices):
        return []
    groups = np.split(indices, np.flatnonzero(np.diff(indices) != 1) + 1)
    return [np.arange(g[0], g[-1] + 2) for g in groups]


def spatial_resample(points, step):
    """Piecewise-linear, equal arc-length samples within ONE valid moving run.

    This changes the estimator; never used for the primary raw-vertex AWN.
    """
    p = np.asarray(points, dtype=float)
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
    if s[-1] <= 0:
        return p[:1]
    if s[-1] / step > 1_000_000:
        raise ValueError("resample_step would create >1 million samples per segment.")
    grid = np.arange(0, s[-1], step)
    if len(grid) and np.isclose(grid[-1], s[-1], rtol=1e-12, atol=1e-12):
        grid[-1] = s[-1]
    else:
        grid = np.r_[grid, s[-1]]
    return np.column_stack([np.interp(grid, s, p[:, j]) for j in (0, 1)])


def calculate(time, xy, cfg=Settings()):
    """Return trial metrics plus per-sample diagnostics, without dropping rows.

    Missing rows, excessive time gaps, speed jumps and short steps break turns.
    Non-increasing timestamps on rows with valid XY fail QC. Rows with missing
    coordinates are treated as tracking gaps and never connected; timestamps on
    those rows do not define tracking-sample order (ANY-maze may repeat them on
    event-only rows).
    """
    t = np.asarray(time, dtype=float)
    raw = np.asarray(xy, dtype=float)
    if raw.shape != (len(t), 2):
        raise ValueError("Expected Nx2 coordinates and N timestamps.")
    n_input_rows = len(t)
    p = raw * cfg.coordinate_scale
    input_finite_xy = np.isfinite(p).all(axis=1)
    # ANY-maze may emit an empty pre-tracking/event row (often at t=0).
    # Ignore only the contiguous prefix before the first measured XY point.
    # Missing coordinates after tracking starts remain in the series as gaps.
    first_xy = int(np.flatnonzero(input_finite_xy)[0]) if input_finite_xy.any() else 0
    n_ignored_leading_xy = first_xy if input_finite_xy.any() else 0
    t = t[first_xy:]
    p = p[first_xy:]
    original_rows = np.arange(first_xy + 1, n_input_rows + 1)
    n = len(t)
    finite_xy = np.isfinite(p).all(axis=1)
    finite_t = np.isfinite(t)
    valid = finite_xy & finite_t
    v = np.diff(p, axis=0)
    lengths = np.linalg.norm(v, axis=1)
    dt = np.diff(t)
    base = valid[:-1] & valid[1:]
    nonpositive = base & (dt <= 0)
    tracking_times = t[valid]
    time_order_bad = bool(np.any(np.diff(tracking_times) <= 0))
    positive = base & (dt > 0)
    gap = positive & (dt > cfg.max_time_gap_s) if cfg.max_time_gap_s else np.zeros(len(dt), bool)
    speed = np.divide(lengths, dt, out=np.full(len(dt), np.nan), where=positive)
    jump = positive & (speed > cfg.max_speed) if cfg.max_speed else np.zeros(len(dt), bool)
    observed = positive & ~gap & ~jump
    moving = observed & (lengths > cfg.min_step)
    usable_turn = moving[:-1] & moving[1:]
    radians = np.full(n, np.nan)
    indices = np.flatnonzero(usable_turn) + 1
    if len(indices):
        a = v[indices - 1] / lengths[indices - 1, None]
        b = v[indices] / lengths[indices, None]
        cross = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
        dot = np.einsum("ij,ij->i", a, b)
        radians[indices] = np.arctan2(np.abs(cross), dot)
    nturn = len(indices)
    awn = float(np.nansum(radians) / (2 * np.pi)) if nturn else np.nan
    observed_s = float(dt[observed].sum())
    moving_s = float(dt[moving].sum())
    moving_path = float(lengths[moving].sum())
    observed_path = float(lengths[observed].sum())
    elapsed = float(t[-1] - t[0]) if n > 1 and finite_t[[0, -1]].all() and not time_order_bad else np.nan
    flags = []
    for condition, name in (
        (not finite_xy.all(), "invalid_xy_rows"),
        (not finite_t.all(), "invalid_time_rows"),
        (time_order_bad, "non_increasing_time"),
        (np.any(gap), "time_gaps"),
        (np.any(jump), "speed_jumps"),
        (np.any(observed & (lengths == 0)), "zero_length_steps"),
        (np.any(observed & (lengths > 0) & (lengths <= cfg.min_step)), "short_steps_removed"),
        (cfg.coordinate_unit == "unknown", "coordinate_unit_unknown"),
        (cfg.max_time_gap_s is None, "gap_threshold_unset"),
        (nturn == 0, "no_usable_turns"),
    ):
        if condition:
            flags.append(name)
    positive_dt = dt[positive]
    metrics = {
        "n_rows": n_input_rows, "n_analyzed_rows": n,
        "n_ignored_leading_xy_rows": n_ignored_leading_xy,
        "n_valid_rows": int(valid.sum()),
        "n_invalid_xy": int((~finite_xy).sum()), "n_invalid_time": int((~finite_t).sum()),
        "n_zero_steps": int((observed & (lengths == 0)).sum()),
        "n_short_steps": int((observed & (lengths <= cfg.min_step)).sum()),
        "n_time_gaps": int(gap.sum()), "n_speed_jumps": int(jump.sum()),
        "n_nonpositive_adjacent_dt": int(nonpositive.sum()),
        "n_observed_edges": int(observed.sum()), "n_moving_edges": int(moving.sum()),
        "n_turns": nturn, "n_moving_segments": len(edge_runs(moving)),
        "usable_turn_fraction": nturn / (n - 2) if n > 2 else np.nan,
        "valid_row_fraction": float(valid.mean()) if n else 0.0,
        "elapsed_s": elapsed, "observed_s": observed_s, "moving_s": moving_s,
        "observed_time_fraction": observed_s / elapsed if elapsed > 0 else np.nan,
        "median_dt_s": float(np.median(positive_dt)) if len(positive_dt) else np.nan,
        "max_dt_s": float(positive_dt.max()) if len(positive_dt) else np.nan,
        "dt_cv": float(np.std(positive_dt) / np.mean(positive_dt)) if len(positive_dt) else np.nan,
        "observed_path_length": observed_path, "moving_path_length": moving_path,
        "coordinate_unit": cfg.coordinate_unit,
        "mean_observed_speed": observed_path / observed_s if observed_s > 0 else np.nan,
        "awn": awn,
        "awn_per_moving_s": awn / moving_s if moving_s > 0 else np.nan,
        "awn_per_moving_unit": awn / moving_path if moving_path > 0 else np.nan,
        "awn_legacy_python": legacy_python_reference(p),
        "awn_matlab_clean_reference": matlab_clean_reference(p),
        "qc_status": "fail" if time_order_bad or not nturn else ("review" if flags else "pass"),
        "qc_flags": ";".join(flags),
    }
    if cfg.resample_step is not None:
        resampled_angles = [turn_angles(spatial_resample(p[g], cfg.resample_step)) for g in edge_runs(moving)]
        counts = sum(len(a) for a in resampled_angles)
        metrics["awn_spatial_resampled"] = float(sum(a.sum() for a in resampled_angles) / (2 * np.pi)) if counts else np.nan
        metrics["n_resampled_turns"] = counts
    diagnostic = {
        "row": original_rows, "time_s": t, "x": p[:, 0], "y": p[:, 1],
        "valid_row": valid,
        "accepted_edge_to_next": np.r_[observed, False] if n else np.array([], bool),
        "moving_edge_to_next": np.r_[moving, False] if n else np.array([], bool),
        "turn_radians": radians,
        "cumulative_awn": np.cumsum(np.nan_to_num(radians)) / (2 * np.pi),
    }
    return metrics, diagnostic

