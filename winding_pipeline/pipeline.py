"""Batch extraction, headless QC, review gating and provenance."""
from pathlib import Path
from datetime import datetime, timezone
import shutil
import html
import csv
import numpy as np
import pandas as pd
from .core import calculate, edge_runs
from .io import read_manifest, read_trajectory, new_output, sha256, versions, write_json


def plot_trial(diag, path, title, cfg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    p = np.column_stack([diag["x"], diag["y"]])
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for run in edge_runs(np.asarray(diag["accepted_edge_to_next"][:-1], bool)):
        axes[0].plot(p[run, 0], p[run, 1], color="#177e89", linewidth=0.9)
    finite = np.isfinite(p).all(axis=1)
    axes[0].scatter(p[finite, 0], p[finite, 1], s=2, color="#177e89", alpha=0.35)
    axes[0].set_aspect("equal", adjustable="datalim")
    if cfg.y_axis == "down":
        axes[0].invert_yaxis()
    axes[0].set(xlabel=f"X ({cfg.coordinate_unit})", ylabel=f"Y ({cfg.coordinate_unit})", title="Observed track; gaps not connected")
    axes[1].plot(diag["row"], diag["cumulative_awn"], color="#955196")
    axes[1].set(xlabel="Input row (not elapsed time)", ylabel="Cumulative absolute turning / 2π", title="Only usable turns contribute")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def extract(manifest, cfg, output, plots=True):
    df = read_manifest(manifest)
    out = new_output(output)
    (out / "diagnostics").mkdir()
    if plots:
        (out / "figures").mkdir()
    shutil.copyfile(manifest, out / "input_manifest.csv")
    write_json(out / "config.json", cfg.to_dict())
    df.to_csv(out / "resolved_manifest.csv", index=False)
    rows = []
    for i, meta in enumerate(df.to_dict("records"), start=1):
        row = dict(meta, trajectory_sha256="", error="", diagnostic_file="", figure_file="")
        try:
            path = Path(meta["trajectory_file"])
            row["trajectory_sha256"] = sha256(path)
            time, xy = read_trajectory(path, cfg)
            metrics, diag = calculate(time, xy, cfg)
            row.update(metrics)
            stem = f"trial_{i:05d}"
            row["diagnostic_file"] = f"diagnostics/{stem}.csv"
            pd.DataFrame(diag).to_csv(out / row["diagnostic_file"], index=False)
            if plots:
                row["figure_file"] = f"figures/{stem}.png"
                plot_trial(diag, out / row["figure_file"], f"{meta['trial_id']} | {metrics['qc_status']} | AWN={metrics['awn']:.4g}", cfg)
        except (OSError, ValueError, UnicodeError, csv.Error) as error:
            row.update(qc_status="error", qc_flags="processing_error", error=str(error), awn=np.nan)
        rows.append(row)
    trials = pd.DataFrame(rows)
    hashes = trials.trajectory_sha256
    duplicates = hashes.ne("") & hashes.duplicated(keep=False)
    for idx in trials.index[duplicates]:
        trials.at[idx, "qc_flags"] = ";".join(filter(None, [trials.at[idx, "qc_flags"], "duplicate_trajectory_content"]))
        if trials.at[idx, "qc_status"] == "pass":
            trials.at[idx, "qc_status"] = "review"
    trials.to_csv(out / "trial_metrics.csv", index=False)
    review_cols = ["trial_id", "subject_id", "session_id", "stage", "day", "trial", "qc_status", "qc_flags"]
    review = trials[review_cols].copy()
    review["decision"] = ""
    review["reviewer"] = ""
    review["reason"] = ""
    review.to_csv(out / "review_template.csv", index=False)
    summaries(trials, out, include_review=False)
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(), "versions": versions(),
        "manifest_sha256": sha256(manifest), "config": cfg.to_dict(),
        "n_trials": len(trials), "status_counts": trials.qc_status.value_counts().to_dict(),
        "interpretation": "QC pass means numerical checks only, not manual validation or biological validity.",
        "code_sha256": {p.name: sha256(p) for p in sorted(Path(__file__).parent.glob('*.py'))},
    }
    write_json(out / "provenance.json", provenance)
    cells = []
    for r in trials.to_dict("records"):
        label = html.escape(str(r["trial_id"]))
        if r.get("figure_file"):
            label = f'<a href="{r["figure_file"]}">{label}</a>'
        cells.append(f"<tr><td>{label}</td><td>{html.escape(str(r['subject_id']))}</td><td>{r['qc_status']}</td><td>{html.escape(str(r['qc_flags']))}</td><td>{r.get('awn', np.nan):.5g}</td><td>{html.escape(str(r['error']))}</td></tr>")
    report = '<!doctype html><meta charset="utf-8"><title>AWN trial QC</title><style>body{font:16px system-ui;margin:2em}td,th{padding:.6em;border-bottom:1px solid #ddd;text-align:left}table{border-collapse:collapse}</style><h1>AWN trial QC</h1><p>AWN is absolute tangent turning, not platform-centered winding. All trials remain in trial_metrics.csv; default summaries use numerical pass only. Review figures and input identity before interpretation.</p><table><tr><th>Trial / figure</th><th>Subject</th><th>Status</th><th>Flags</th><th>AWN</th><th>Error</th></tr>' + "".join(cells) + '</table>'
    (out / "qc_report.html").write_text(report, encoding="utf-8")
    return trials


def select_trials(trials, include_review=False, review=None):
    allowed = trials.qc_status.eq("pass") | (include_review & trials.qc_status.eq("review"))
    if review is not None:
        required = {"trial_id", "decision", "reviewer", "reason"}
        if not required.issubset(review.columns) or review.trial_id.duplicated().any():
            raise ValueError("Review requires unique trial_id, decision, reviewer, reason.")
        if set(review.trial_id) != set(trials.trial_id):
            raise ValueError("Review must cover exactly the extraction's trial IDs.")
        r = review.set_index("trial_id").reindex(trials.trial_id).fillna("")
        if not r.decision.isin(["include", "exclude"]).all():
            raise ValueError("Complete every review decision as include or exclude.")
        if r.reviewer.str.strip().eq("").any() or r.reason.str.strip().eq("").any():
            raise ValueError("Document reviewer and reason for every decision.")
        allowed = r.decision.eq("include").to_numpy()
        if np.any(allowed & trials.qc_status.isin(["fail", "error"]).to_numpy()):
            raise ValueError("Failed/error trials cannot be included; correct inputs and rerun.")
    awn = pd.to_numeric(trials.awn, errors="coerce")
    return trials.loc[allowed & np.isfinite(awn)].copy()


def summaries(trials, out, include_review=False, review=None):
    selected = select_trials(trials, include_review, review)
    metric_cols = [c for c in ("awn", "awn_per_moving_s", "awn_per_moving_unit", "observed_path_length", "elapsed_s", "observed_s", "moving_s", "awn_spatial_resampled") if c in trials]
    for c in metric_cols:
        selected[c] = pd.to_numeric(selected[c], errors="coerce")
    keys = ["subject_id", "session_id", "stage", "day"]
    aggregations = {c: (c, "mean") for c in metric_cols}
    aggregations["n_trials_included"] = ("trial_id", "size")
    for c in ("group", "cohort", "sex", "genotype"):
        if c in trials:
            aggregations[c] = (c, "first")
    if "age_months" in trials:
        selected["age_months"] = pd.to_numeric(selected.age_months, errors="coerce")
        aggregations["age_months"] = ("age_months", "mean")
    total = trials.groupby(keys, dropna=False).size().rename("n_trials_total").reset_index()
    day = selected.groupby(keys, dropna=False).agg(**aggregations).reset_index()
    day = total.merge(day, on=keys, how="left", validate="one_to_one")
    day["n_trials_included"] = day.n_trials_included.fillna(0).astype(int)
    day.to_csv(out / "subject_day.csv", index=False)
    # Retain every trial and its inclusion decision, including empty subject-days.
    decisions = trials[["trial_id", "subject_id", "qc_status", "qc_flags"]].copy()
    decisions["included"] = decisions.trial_id.isin(selected.trial_id)
    decisions.to_csv(out / "analysis_inclusion.csv", index=False)
    if "group" in day:
        group = day.groupby(["session_id", "stage", "day", "group"], dropna=False).awn.agg(["mean", "std", "count"]).reset_index()
        group = group.rename(columns={"count": "n_subjects"})
        group["sem"] = group["std"] / np.sqrt(group.n_subjects)
        group.to_csv(out / "group_day_descriptive.csv", index=False)
        plot_subject_days(day, out)
    return day


def plot_subject_days(day, out):
    """Descriptive curves with mice, not trials, as the group-summary unit."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    usable = day[np.isfinite(day.awn) & day.group.notna() & day.group.ne("")].copy()
    if usable.empty:
        return
    target = out / "summary_figures"
    target.mkdir(exist_ok=True)
    index = []
    for number, ((session, stage), data) in enumerate(usable.groupby(["session_id", "stage"]), start=1):
        fig, ax = plt.subplots(figsize=(8, 5))
        for j, (group, values) in enumerate(data.groupby("group")):
            color = f"C{j % 10}"
            for _, mouse in values.groupby("subject_id"):
                mouse = mouse.sort_values("day")
                ax.plot(mouse.day, mouse.awn, color=color, alpha=.15, linewidth=.8)
            aggregate = values.groupby("day").awn.agg(["mean", "sem", "count"]).reset_index()
            ax.plot(aggregate.day, aggregate["mean"], "o-", color=color, label=str(group))
            has_sem = aggregate["sem"].notna()
            ax.errorbar(aggregate.loc[has_sem, "day"], aggregate.loc[has_sem, "mean"],
                        yerr=aggregate.loc[has_sem, "sem"], fmt="none", color=color, capsize=3)
        ax.set(xlabel="Study day", ylabel="Mean trial AWN per mouse", title=f"{session} / {stage}\nDescriptive group mean ± SEM; faint lines = mice")
        ax.set_xticks(sorted(data.day.unique()))
        ax.legend(title="Recorded group")
        fig.tight_layout()
        filename = f"summary_figures/curve_{number:03d}.png"
        fig.savefig(out / filename, dpi=140)
        plt.close(fig)
        index.append(dict(session_id=session, stage=stage, figure_file=filename))
    pd.DataFrame(index).to_csv(out / "summary_figure_index.csv", index=False)
