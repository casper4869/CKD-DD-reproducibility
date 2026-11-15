"""Single prespecified country-average age-profile sensitivity (offline).

Usage: python run_age_reference.py --ckd PATH --dd PATH --output PATH
Optional --replication-comparator PATH audits the original stored 408 results.
No network calls or data downloads occur in this script.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd

AGES = [f"{a}-{a+4} years" for a in range(5, 95, 5)] + ["95+ years"]
SEED = 20261115
N_BOOTSTRAP = 2000


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_rates(path, disease):
    df = pd.read_csv(path)
    assert set(df["year"]) == {2021}, "Input must contain only 2021"
    assert set(df["metric_name"]) == {"Rate"}
    assert set(df["measure_name"]) == {"Incidence"}
    expected_cause = {"CKD": "Chronic kidney disease", "DD": "Depressive disorders"}[disease]
    assert set(df["cause_name"]) == {expected_cause}, "Disease label does not match input role"
    d = df.loc[df.age_name.isin(AGES)].copy()
    assert set(d.sex_name) == {"Female", "Male"}
    assert d.location_id.nunique() == 204
    assert len(d) == 204 * 2 * 19
    assert not d.duplicated(["location_id", "sex_name", "age_name"]).any()
    assert np.isfinite(d.val).all() and (d.val >= 0).all()
    assert not d.location_id.eq(1).any(), "Contract assumes no official Global aggregate"
    return d[["location_id", "location_name", "sex_name", "age_name", "val"]].rename(columns={"val": disease})


def pearson(x, y):
    a, b = np.asarray(x) - np.mean(x), np.asarray(y) - np.mean(y)
    return float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))


def residualize(y, reference):
    # Centered OLS with intercept, 19 observations, 2 parameters, 17 residual df.
    x0 = reference - reference.mean()
    y0 = y - y.mean()
    slope = np.dot(x0, y0) / np.dot(x0, x0)
    intercept = y.mean() - slope * reference.mean()
    e = y0 - slope * x0
    r2 = 1 - np.dot(e, e) / np.dot(y0, y0)
    # Relative orthogonality checks avoid dependence on absolute disease-rate scale.
    tol = 1e-10
    assert abs(e.sum()) / max(1.0, np.linalg.norm(y)) < tol
    assert abs(np.dot(e, x0)) / max(1.0, np.linalg.norm(y) * np.linalg.norm(x0)) < tol
    return e, float(intercept), float(slope), float(r2)


def interval(values):
    return [float(v) for v in np.quantile(values, [0.025, 0.975])]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ckd", required=True, type=Path)
    ap.add_argument("--dd", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--replication-comparator", type=Path)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    contract_path = Path(__file__).with_name("analysis_contract.md")
    assert contract_path.is_file(), "Prespecified analysis contract must accompany script"

    ckd, dd = read_rates(args.ckd, "CKD"), read_rates(args.dd, "DD")
    # Assert complete identity of analysis rows and country names before one-to-one join.
    keys = ["location_id", "location_name", "sex_name", "age_name"]
    data = ckd.merge(dd, on=keys, how="outer", validate="one_to_one", indicator=True)
    assert data._merge.eq("both").all()
    data = data.drop(columns="_merge")
    records, residual_records = [], []
    for sex in ["Female", "Male"]:
        sub = data.loc[data.sex_name.eq(sex)]
        ids = sorted(sub.location_id.unique())
        names = sub.groupby("location_id").location_name.first().to_dict()
        mat = {k: sub.pivot(index="location_id", columns="age_name", values=k).loc[ids, AGES].to_numpy() for k in ["CKD", "DD"]}
        assert all(m.shape == (204, 19) for m in mat.values())
        for i, location_id in enumerate(ids):
            ref = {k: (m.sum(axis=0) - m[i, :]) / 203 for k, m in mat.items()}
            fit = {k: residualize(m[i, :], ref[k]) for k, m in mat.items()}
            raw = pearson(mat["CKD"][i], mat["DD"][i])
            adjusted = pearson(fit["CKD"][0], fit["DD"][0])
            item = {"location_id": int(location_id), "location_name": names[location_id], "sex_name": sex,
                    "n_age": 19, "n_reference_countries": 203, "residual_df_each_model": 17,
                    "profile_r_original": raw, "profile_r_reference_adjusted": adjusted,
                    "paired_difference": adjusted - raw}
            for k in ["CKD", "DD"]:
                item.update({f"{k}_reference_intercept": fit[k][1], f"{k}_reference_slope": fit[k][2], f"{k}_reference_R2": fit[k][3]})
            records.append(item)
            for j, age in enumerate(AGES):
                residual_records.append({"location_id": int(location_id), "location_name": names[location_id], "sex_name": sex, "age_name": age,
                                         "CKD": mat["CKD"][i,j], "DD": mat["DD"][i,j],
                                         "CKD_reference": ref["CKD"][j], "DD_reference": ref["DD"][j],
                                         "CKD_residual": fit["CKD"][0][j], "DD_residual": fit["DD"][0][j]})
    profiles = pd.DataFrame(records)
    residuals = pd.DataFrame(residual_records)
    assert profiles.shape[0] == 408 and residuals.shape[0] == 7752
    max_diff = None
    if args.replication_comparator:
        prior = pd.read_csv(args.replication_comparator)
        check = profiles.merge(prior, on=["location_id", "location_name", "sex_name", "n_age"], validate="one_to_one")
        assert len(check) == 408
        max_diff = float(np.max(np.abs(check.profile_r_original - check.profile_r)))
        assert max_diff < 1e-12

    summaries = []
    # Paired bootstrap uses country indices identically for original/adjusted profiles.
    rng = np.random.default_rng(SEED)
    for sex, d in profiles.groupby("sex_name", sort=True):
        raw = d.profile_r_original.to_numpy()
        adj = d.profile_r_reference_adjusted.to_numpy()
        ii = rng.integers(0, len(d), size=(N_BOOTSTRAP, len(d)))
        sr = {"sex_name": sex, "n_country": len(d), "n_age": 19,
              "n_positive_original": int((raw > 0).sum()), "n_positive_adjusted": int((adj > 0).sum()),
              "n_positive_to_nonpositive": int(((raw > 0) & (adj <= 0)).sum()),
              "n_nonpositive_to_positive": int(((raw <= 0) & (adj > 0)).sum()),
              "median_original": float(np.median(raw)), "q25_original": float(np.quantile(raw, .25)), "q75_original": float(np.quantile(raw, .75)),
              "median_adjusted": float(np.median(adj)), "q25_adjusted": float(np.quantile(adj, .25)), "q75_adjusted": float(np.quantile(adj, .75)),
              "min_adjusted": float(adj.min()), "max_adjusted": float(adj.max()),
              "median_paired_difference": float(np.median(adj-raw)),
              "median_CKD_reference_R2": float(d.CKD_reference_R2.median()),
              "median_DD_reference_R2": float(d.DD_reference_R2.median())}
        for name, arr in [("median_original", np.median(raw[ii], axis=1)), ("median_adjusted", np.median(adj[ii], axis=1)), ("median_paired_difference", np.median((adj-raw)[ii], axis=1))]:
            low, high = interval(arr)
            sr[f"{name}_bootstrap_lower"] = low
            sr[f"{name}_bootstrap_upper"] = high
        summaries.append(sr)
    summary = pd.DataFrame(summaries)
    profiles.to_csv(args.output / "S3_age_reference_country.csv", index=False, float_format="%.17g")
    residuals.to_csv(args.output / "S3_age_reference_residuals.csv", index=False, float_format="%.17g")
    summary.to_csv(args.output / "S3_age_reference_summary.csv", index=False, float_format="%.17g")
    report = {"status": "calculated_pending_independent_review", "contract_sha256": sha256(contract_path),
              "script_sha256": sha256(__file__), "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
              "input_sha256": {args.ckd.name: sha256(args.ckd), args.dd.name: sha256(args.dd)},
              "replication_max_abs_error": max_diff,
              "reference": "sex/disease-specific equal-country-weighted leave-one-country-out mean (203 reference countries); no official GBD Global aggregate used",
              "method": "Pearson correlation of OLS residuals after separate intercept-plus-reference regressions for CKD and DD across 19 age groups",
              "bootstrap_seed": SEED, "bootstrap_replicates": N_BOOTSTRAP, "bootstrap_unit": "country",
              "summary": summaries,
              "outputs_sha256": {p.name: sha256(p) for p in sorted(args.output.glob("S3_age_reference*.csv"))}}
    (args.output / "age_reference_verification.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
