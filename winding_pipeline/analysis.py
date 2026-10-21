"""Optional exploratory repeated-measures analysis and safe MRI linking.

These functions do not decide the lab's endpoint, exclusions, or causal model.
"""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from .io import new_output, read_csv_strict, sha256, versions, write_json


def subject_day_data(path, session_id, stage):
    d = pd.read_csv(path, dtype={"subject_id": str, "session_id": str}, keep_default_na=False)
    required = {"subject_id", "session_id", "stage", "day", "awn", "n_trials_included"}
    if not required.issubset(d):
        raise ValueError(f"Missing subject-day columns: {required - set(d.columns)}")
    d = d[(d.session_id == session_id) & (d.stage == stage)].copy()
    if d.empty:
        raise ValueError("No records match the selected session and stage.")
    if d.duplicated(["subject_id", "day"]).any():
        raise ValueError("Expected at most one observation per subject/day.")
    return d


def prepare_model_data(d, outcome, covariates, log1p=False):
    allowed_covars = {"age_months", "sex", "genotype", "cohort"}
    allowed_outcomes = {"awn", "awn_per_moving_s", "awn_per_moving_unit", "observed_path_length", "awn_spatial_resampled"}
    if outcome not in allowed_outcomes or set(covariates) - allowed_covars:
        raise ValueError("Unsupported outcome/covariate; use documented column names.")
    if len(set(covariates)) != len(covariates):
        raise ValueError("Duplicate covariate.")
    needed = [outcome, "subject_id", "day", "group", *covariates]
    if not set(needed).issubset(d):
        raise ValueError(f"Missing required model columns: {set(needed) - set(d.columns)}")
    d = d.copy()
    for c in (outcome, "day", "age_months"):
        if c in d:
            d[c] = pd.to_numeric(d[c], errors="coerce").replace([np.inf, -np.inf], np.nan)
    d[needed] = d[needed].replace("", np.nan)
    excluded = d[d[needed].isna().any(axis=1)].copy()
    d = d.dropna(subset=needed).copy()
    if d.empty:
        raise ValueError("No complete model records remain.")
    if d.subject_id.nunique() < 6 or d.group.nunique() < 2:
        raise ValueError("Exploratory model guard: need >=6 mice and >=2 groups; not a power guarantee.")
    if d.groupby("group").subject_id.nunique().min() < 3:
        raise ValueError("Need >=3 mice in each group; not a power guarantee.")
    if d.groupby("subject_id").size().min() < 2 or d.day.nunique() < 2:
        raise ValueError("Each retained mouse needs >=2 days for this longitudinal model.")
    if d.groupby("subject_id").group.nunique().gt(1).any():
        raise ValueError("Group changes within mouse: crossover design needs a different model.")
    for c in covariates:
        if d[c].nunique() < 2:
            raise ValueError(f"Covariate {c} is constant; explicitly remove it from the model.")
    if log1p and (d[outcome] < 0).any():
        raise ValueError("log1p analysis requires nonnegative outcomes.")
    d["response"] = np.log1p(d[outcome]) if log1p else d[outcome]
    if d.response.nunique() < 2:
        raise ValueError("Outcome has no variation.")
    return d, excluded


def model(subject_day, output, session_id, stage, outcome="awn", covariates=(), log1p=False):
    try:
        import statsmodels.formula.api as smf
        import statsmodels.api as sm
        import patsy
    except ImportError as e:
        raise RuntimeError("Model dependency missing. Install with: python -m pip install '.[stats]' ") from e
    original = subject_day_data(subject_day, session_id, stage)
    d, excluded = prepare_model_data(original, outcome, covariates, log1p)
    out = new_output(output)
    d.to_csv(out / "model_rows.csv", index=False)
    excluded.to_csv(out / "excluded_missing_model_fields.csv", index=False)
    # Day categorical: no unverified assumption of linear learning across days.
    formula = "response ~ C(group) * C(day)"
    for c in covariates:
        formula += " + " + (c if c == "age_months" else f"C({c})")
    provenance = {"versions": versions(), "subject_day_sha256": sha256(subject_day),
                  "formula": formula, "outcome": outcome, "log1p": log1p,
                  "stage": stage, "session_id": session_id,
                  "n_subjects": int(d.subject_id.nunique()), "n_rows": len(d),
                  "n_excluded_rows": len(excluded), "interpretation": "Exploratory, noncausal. Wald p-values unadjusted; no automatic endpoint search."}
    if "cohort" in d:
        pd.crosstab(d.drop_duplicates("subject_id").group,
                    d.drop_duplicates("subject_id").cohort).to_csv(out / "group_by_cohort.csv")
        if d.cohort.nunique() > 1:
            provenance["cohort_warning"] = "Multiple cohorts: inspect overlap. Age matching or covariates cannot identify a treatment effect completely confounded with cohort."
    y, x = patsy.dmatrices(formula, d, return_type="dataframe", NA_action="raise")
    rank = int(np.linalg.matrix_rank(x.to_numpy()))
    provenance.update(design_rank=rank, design_columns=x.shape[1])
    if rank < x.shape[1] or len(d) <= x.shape[1] + 2:
        provenance["status"] = "blocked_rank_or_residual_df"
        write_json(out / "model_status.json", provenance)
        raise ValueError("Rank-deficient/confounded or overparameterized design; see model_status.json.")
    caught = []
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            fit = smf.mixedlm(formula, d, groups=d.subject_id, re_formula="1", missing="raise").fit(reml=True, method="lbfgs")
        provenance["warnings"] = [str(w.message) for w in caught]
        variance = float(fit.cov_re.iloc[0, 0])
        stable = (bool(fit.converged) and not caught and np.isfinite(variance)
                  and variance > 1e-8 * max(float(fit.scale), 1e-12)
                  and np.all(np.isfinite(fit.bse_fe)))
        provenance.update(status="fit_ok_exploratory" if stable else "review_fit_no_pvalues",
                          converged=bool(fit.converged),
                          random_intercept_variance=variance if np.isfinite(variance) else None)
        k = len(fit.fe_params)
        ci = fit.conf_int().iloc[:k]
        coef = pd.DataFrame({"term": fit.fe_params.index, "estimate": fit.fe_params.values,
                             "standard_error": fit.bse_fe.values,
                             "ci95_low": ci.iloc[:, 0].to_numpy() if stable else np.nan,
                             "ci95_high": ci.iloc[:, 1].to_numpy() if stable else np.nan,
                             "p_wald_unadjusted": fit.pvalues.iloc[:k].to_numpy() if stable else np.nan})
        coef.to_csv(out / "fixed_effects.csv", index=False)
        # Fixed-effects residuals still include subject heterogeneity; label them.
        predicted = np.asarray(x @ fit.fe_params)
        diagnostic = d[["subject_id", "day", "group"]].copy()
        diagnostic["response"] = d.response
        diagnostic["fixed_effects_prediction"] = predicted
        diagnostic["marginal_residual"] = d.response.to_numpy() - predicted
        diagnostic.to_csv(out / "model_diagnostics.csv", index=False)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(9, 4))
        ax[0].scatter(predicted, diagnostic.marginal_residual, s=14)
        ax[0].axhline(0, color="black", linewidth=.7)
        ax[0].set(xlabel="Fixed-effects prediction", ylabel="Marginal residual")
        sm.qqplot(diagnostic.marginal_residual.to_numpy(), line="s", ax=ax[1])
        fig.suptitle("Exploratory model: marginal residual checks")
        fig.tight_layout()
        fig.savefig(out / "model_diagnostics.png", dpi=140)
        plt.close(fig)
    except (ValueError, np.linalg.LinAlgError) as e:
        provenance.update(status="fit_failed", error=str(e), warnings=[str(w.message) for w in caught])
    write_json(out / "model_status.json", provenance)
    return provenance


def join_imaging(subject_day, crosswalk, imaging, output, session_id, stage, day):
    """One prespecified behavioral endpoint per mouse; one MRI row per scan ID."""
    b = subject_day_data(subject_day, session_id, stage)
    b = b[pd.to_numeric(b.day, errors="coerce") == day].copy()
    if b.empty:
        raise ValueError("No behavioral rows for the selected day.")
    c = read_csv_strict(crosswalk)
    m = read_csv_strict(imaging)
    if set(c.columns) != {"subject_id", "imaging_id"} or "imaging_id" not in m:
        raise ValueError("Crosswalk needs exactly subject_id,imaging_id; MRI CSV needs imaging_id.")
    for table, ids in ((b, ["subject_id"]), (c, ["subject_id", "imaging_id"]), (m, ["imaging_id"])):
        for key in ids:
            table[key] = table[key].astype(str).str.strip()
            if table[key].eq("").any() or table[key].duplicated().any():
                raise ValueError(f"Empty or duplicate {key}; resolve repeated scans/IDs explicitly.")
    # Prefix all MRI columns except the key, preventing unnoticed metadata clashes.
    m = m.rename(columns={col: f"mri_{col}" for col in m if col != "imaging_id"})
    if (set(b) | set(c)) & (set(m) - {"imaging_id"}):
        raise ValueError("MRI-prefixed columns collide with existing behavioral fields.")
    joined = b.merge(c, on="subject_id", how="left", validate="one_to_one", indicator="crosswalk_match")
    # Several unmapped mice may have null imaging IDs. Non-null IDs were already
    # checked one-to-one in the crosswalk; permit these unmatched rows to remain.
    joined = joined.merge(m, on="imaging_id", how="left", validate="many_to_one", indicator="imaging_match")
    out = new_output(output)
    joined.to_csv(out / "behavior_imaging_join.csv", index=False)
    c[~c.subject_id.isin(b.subject_id)].to_csv(out / "crosswalk_without_behavior.csv", index=False)
    m[~m.imaging_id.isin(joined.imaging_id.dropna())].to_csv(out / "imaging_without_behavior.csv", index=False)
    write_json(out / "join_provenance.json", {"versions": versions(), "stage": stage, "day": day,
        "session_id": session_id, "input_sha256": {"behavior": sha256(subject_day), "crosswalk": sha256(crosswalk), "imaging": sha256(imaging)},
        "note": "No ID inference, correlations, group comparisons or causal claims are made by this join."})
    return joined
