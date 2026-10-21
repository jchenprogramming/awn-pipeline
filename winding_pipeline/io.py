"""Explicit CSV schema handling and immutable-input provenance."""
import csv
import hashlib
import json
import platform
from pathlib import Path
from importlib.metadata import version, PackageNotFoundError
import numpy as np
import pandas as pd
from . import __version__
from .core import Settings

REQUIRED = ["trial_id", "subject_id", "session_id", "stage", "day", "trial", "trajectory_file"]
OPTIONAL = ["group", "cohort", "sex", "genotype", "age_months", "notes"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def versions():
    found = {"qial-awn": __version__, "python": platform.python_version()}
    for name in ("numpy", "pandas", "matplotlib", "scipy", "statsmodels"):
        try:
            found[name] = version(name)
        except PackageNotFoundError:
            found[name] = "not installed"
    return found


def load_settings(path):
    with open(path, encoding="utf-8") as f:
        return Settings(**json.load(f))


def read_csv_strict(path, delimiter=",", encoding="utf-8-sig", skiprows=0):
    """Reject duplicate headers/ragged rows rather than silently reindex a CSV."""
    with open(path, encoding=encoding, newline="") as f:
        for _ in range(skiprows):
            next(f, None)
        reader = csv.reader(f, delimiter=delimiter, strict=True)
        header = next(reader, None)
        if not header:
            raise ValueError(f"Empty CSV/header: {path}")
        header = [c.strip() for c in header]
        if any(not c for c in header) or len(set(header)) != len(header):
            raise ValueError(f"Blank or duplicate headers: {header}")
        rows = []
        for number, row in enumerate(reader, start=skiprows + 2):
            if not row:
                row = [""] * len(header)  # preserve a blank record as a gap
            if len(row) != len(header):
                raise ValueError(f"Ragged CSV at record {number}: expected {len(header)} fields, got {len(row)}")
            rows.append(row)
    return pd.DataFrame(rows, columns=header)


def parse_time(value, decimal="."):
    s = str(value).strip().replace(decimal, ".")
    try:
        parts = s.split(":")
        if len(parts) == 1:
            result = float(s)
        elif len(parts) in (2, 3):
            numbers = [float(x) for x in parts]
            if not all(np.isfinite(numbers)) or any(x < 0 for x in numbers):
                return np.nan
            if numbers[-1] >= 60 or (len(parts) == 3 and numbers[-2] >= 60):
                return np.nan
            if any(x != int(x) for x in numbers[:-1]):
                return np.nan
            result = sum(x * 60 ** i for i, x in enumerate(reversed(numbers)))
        else:
            return np.nan
        return result if np.isfinite(result) and result >= 0 else np.nan
    except (ValueError, OverflowError):
        return np.nan


def read_trajectory(path, cfg):
    df = read_csv_strict(path, cfg.delimiter, cfg.encoding, cfg.skiprows)
    required = [cfg.time_column, cfg.x_column, cfg.y_column]
    if not set(required).issubset(df.columns):
        raise ValueError(f"Required columns {required}; found {df.columns.tolist()}. Configure exact names; no numeric-column guessing.")
    time = np.array([parse_time(v, cfg.decimal) for v in df[cfg.time_column]])
    coords = []
    for col in (cfg.x_column, cfg.y_column):
        text = df[col].str.strip().str.replace(cfg.decimal, ".", regex=False)
        coords.append(pd.to_numeric(text, errors="coerce").to_numpy(float))
    return time, np.column_stack(coords)


def read_manifest(path):
    path = Path(path).resolve()
    df = read_csv_strict(path)
    missing = set(REQUIRED) - set(df.columns)
    unknown = set(df.columns) - set(REQUIRED + OPTIONAL)
    if missing or unknown:
        raise ValueError(f"Manifest missing columns {sorted(missing)}; unsupported columns {sorted(unknown)}.")
    for c in df:
        df[c] = df[c].str.strip()
    if df.empty or (df[REQUIRED] == "").any().any():
        raise ValueError("Manifest must contain trials and every required field must be populated.")
    if not df.stage.isin(["acquisition", "probe", "reversal", "habituation"]).all():
        raise ValueError("stage must explicitly be acquisition, probe, reversal, or habituation.")
    for c in ("day", "trial"):
        nums = pd.to_numeric(df[c], errors="coerce")
        if not (np.isfinite(nums) & (nums >= 0) & (nums == np.floor(nums))).all():
            raise ValueError(f"{c} must be a nonnegative integer.")
        df[c] = nums.astype(int)
    if "age_months" in df:
        nums = pd.to_numeric(df.age_months.replace("", np.nan), errors="coerce")
        if ((df.age_months != "") & (~np.isfinite(nums) | (nums < 0))).any():
            raise ValueError("age_months must be nonnegative numeric or blank.")
        df["age_months"] = nums
    if df.trial_id.duplicated().any() or df.duplicated(["subject_id", "session_id", "stage", "day", "trial"]).any():
        raise ValueError("Duplicate trial_id or subject/session/stage/day/trial key.")
    for col in ("group", "cohort", "sex", "genotype"):
        if col in df and df.groupby("subject_id")[col].nunique().gt(1).any():
            raise ValueError(f"Inconsistent {col} for a subject. Do not silently choose first metadata row.")
    resolved = [(path.parent / name).resolve() for name in df.trajectory_file]
    if len(set(resolved)) != len(resolved):
        raise ValueError("A trajectory file is assigned to more than one trial.")
    df["trajectory_file"] = [str(p) for p in resolved]
    return df


def new_output(path):
    out = Path(path).resolve()
    # Reject any existing target: no accidental overwrite of earlier analyses.
    out.mkdir(parents=True, exist_ok=False)
    return out


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")

