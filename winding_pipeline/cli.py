import argparse
from pathlib import Path
import sys
import shutil
import pandas as pd
from .io import load_settings, new_output, versions, sha256, write_json
from .pipeline import extract, summaries


def main(argv=None):
    parser = argparse.ArgumentParser(description="Auditable ANY-maze absolute tangent turning (AWN).")
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("extract", help="Extract trajectories using explicit manifest and config.")
    p.add_argument("--manifest", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--no-plots", action="store_true")
    p = commands.add_parser("summarize", help="Create reviewed subject-day analysis tables.")
    p.add_argument("--trials", required=True)
    p.add_argument("--out", required=True)
    g = p.add_mutually_exclusive_group()
    g.add_argument("--review", help="Completed review_template.csv, all trials explicitly decided.")
    g.add_argument("--include-review", action="store_true", help="Exploratory sensitivity only; includes review but never failed/error trials.")
    p = commands.add_parser("model", help="Optional exploratory random-intercept mixed model.")
    p.add_argument("--subject-day", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--session-id", required=True)
    p.add_argument("--stage", required=True, choices=["acquisition", "probe", "reversal"])
    p.add_argument("--outcome", default="awn")
    p.add_argument("--covariates", nargs="*", default=[])
    p.add_argument("--log1p", action="store_true")
    p = commands.add_parser("join-imaging", help="Link one behavioral endpoint per mouse to an MRI table.")
    for flag in ("subject-day", "crosswalk", "imaging", "out", "session-id", "stage"):
        p.add_argument(f"--{flag}", required=True)
    p.add_argument("--day", type=int, required=True)
    p = commands.add_parser("demo", help="Write synthetic test trajectories; no real animals.")
    p.add_argument("--out", required=True)
    p = commands.add_parser("inventory", help="List CSVs directly in one folder; create a manifest with IDs left blank.")
    p.add_argument("--directory", required=True)
    p.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "extract":
            trials = extract(args.manifest, load_settings(args.config), args.out, not args.no_plots)
            print(trials.qc_status.value_counts().to_string())
            print(f"Results: {Path(args.out).resolve()}")
            return 2 if trials.qc_status.isin(["fail", "error"]).any() else 0
        if args.command == "summarize":
            d = pd.read_csv(args.trials, dtype={"trial_id": str, "subject_id": str, "session_id": str}, keep_default_na=False)
            review = pd.read_csv(args.review, dtype=str, keep_default_na=False) if args.review else None
            # Validate decisions before creating output folder.
            from .pipeline import select_trials
            select_trials(d, args.include_review, review)
            out = new_output(args.out)
            summaries(d, out, args.include_review, review)
            if args.review:
                shutil.copyfile(args.review, out / "applied_review.csv")
            write_json(out / "summary_provenance.json", {"versions": versions(), "trial_metrics_sha256": sha256(args.trials),
                "include_review": args.include_review, "review_sha256": sha256(args.review) if args.review else None})
        elif args.command == "model":
            from .analysis import model
            result = model(args.subject_day, args.out, args.session_id, args.stage, args.outcome, args.covariates, args.log1p)
            print(result["status"])
            return 0 if result["status"] == "fit_ok_exploratory" else 2
        elif args.command == "join-imaging":
            from .analysis import join_imaging
            join_imaging(args.subject_day, args.crosswalk, args.imaging, args.out, args.session_id, args.stage, args.day)
        elif args.command == "demo":
            from .demo import make_demo
            make_demo(args.out)
        elif args.command == "inventory":
            from .io import REQUIRED, OPTIONAL
            folder = Path(args.directory).resolve()
            if not folder.is_dir():
                raise ValueError("Trajectory directory does not exist or is inaccessible.")
            files = sorted(p for p in folder.iterdir() if p.is_file() and p.suffix.lower() == ".csv" and not p.name.startswith("."))
            if not files:
                raise ValueError("No CSV files directly in this folder. No recursive search was attempted.")
            out = new_output(args.out)
            rows = []
            for i, path in enumerate(files, start=1):
                r = dict.fromkeys(REQUIRED + OPTIONAL, "")
                r.update(trial_id=f"file_{i:05d}", trajectory_file=str(path), notes="VERIFY identity against ANY-maze Test ID/export schedule")
                rows.append(r)
            pd.DataFrame(rows).to_csv(out / "manifest_to_complete.csv", index=False)
            print(f"Inventoried {len(files)} files. Subject/session/stage/day/trial MUST be filled before extraction.")
        print(f"Saved to {Path(args.out).resolve()}")
        return 0
    except (OSError, ValueError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
