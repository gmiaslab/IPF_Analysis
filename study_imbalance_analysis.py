#!/usr/bin/env python3

from __future__ import annotations

import math
import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "study_imbalance_results"
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"
TARGET_STUDY = "GSE213001"

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS_DIR / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def ensure_dirs() -> None:
    for path in (RESULTS_DIR, TABLES_DIR, FIGURES_DIR, RESULTS_DIR / ".mplconfig"):
        path.mkdir(parents=True, exist_ok=True)


def load_retained_metadata() -> pd.DataFrame:
    metadata = pd.read_csv(BASE_DIR / "combined_metadata.csv")
    counts_header = pd.read_csv(BASE_DIR / "combined_study_matrix.csv", nrows=1).columns[1:]

    required = ["Run", "Condition", "Sex", "Age", "GEO", "LibraryLayout", "Subject"]
    metadata = metadata.loc[metadata["Run"].isin(counts_header)].copy()
    metadata = metadata.dropna(subset=required)
    metadata = metadata.loc[
        metadata["Condition"].isin(["Healthy", "IPF"]) & metadata["Sex"].isin(["female", "male"])
    ].copy()

    metadata["Age"] = pd.to_numeric(metadata["Age"])
    metadata["Group"] = metadata["Condition"] + "_" + metadata["Sex"]
    metadata["Set"] = np.where(metadata["GEO"] == TARGET_STUDY, TARGET_STUDY, "Other studies")

    subject_counts = metadata["Subject"].value_counts()
    metadata["RepeatedSubjectSample"] = metadata["Subject"].map(subject_counts).gt(1)
    return metadata.sort_values("Run").reset_index(drop=True)


def pooled_sd(x1: pd.Series, x0: pd.Series) -> float:
    n1 = x1.notna().sum()
    n0 = x0.notna().sum()
    if n1 < 2 or n0 < 2:
        return float("nan")
    var1 = x1.var(ddof=1)
    var0 = x0.var(ddof=1)
    pooled_var = ((n1 - 1) * var1 + (n0 - 1) * var0) / (n1 + n0 - 2)
    if pooled_var <= 0:
        return float("nan")
    return math.sqrt(pooled_var)


def smd_continuous(x1: pd.Series, x0: pd.Series) -> float:
    sd = pooled_sd(x1, x0)
    if np.isnan(sd) or sd == 0:
        return float("nan")
    return (x1.mean() - x0.mean()) / sd


def smd_binary(x1: pd.Series, x0: pd.Series) -> float:
    p1 = x1.mean()
    p0 = x0.mean()
    denom = math.sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
    if denom == 0:
        return float("nan")
    return (p1 - p0) / denom


def build_study_summary(metadata: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for geo, df in metadata.groupby("GEO", sort=True):
        rows.append(
            {
                "GEO": geo,
                "n_samples": len(df),
                "n_subjects": df["Subject"].nunique(),
                "n_repeated_subjects": int(df["Subject"].value_counts().gt(1).sum()),
                "n_repeated_subject_samples": int(df["RepeatedSubjectSample"].sum()),
                "n_healthy": int((df["Condition"] == "Healthy").sum()),
                "n_ipf": int((df["Condition"] == "IPF").sum()),
                "n_female": int((df["Sex"] == "female").sum()),
                "n_male": int((df["Sex"] == "male").sum()),
                "n_healthy_female": int((df["Group"] == "Healthy_female").sum()),
                "n_healthy_male": int((df["Group"] == "Healthy_male").sum()),
                "n_ipf_female": int((df["Group"] == "IPF_female").sum()),
                "n_ipf_male": int((df["Group"] == "IPF_male").sum()),
                "age_mean": round(df["Age"].mean(), 2),
                "age_sd": round(df["Age"].std(ddof=1), 2),
                "age_median": round(df["Age"].median(), 2),
                "age_min": round(df["Age"].min(), 2),
                "age_max": round(df["Age"].max(), 2),
            }
        )
    return pd.DataFrame(rows)


def build_set_summary(metadata: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for set_name, df in metadata.groupby("Set", sort=False):
        rows.append(
            {
                "Set": set_name,
                "n_samples": len(df),
                "n_subjects": df["Subject"].nunique(),
                "n_repeated_subjects": int(df["Subject"].value_counts().gt(1).sum()),
                "n_repeated_subject_samples": int(df["RepeatedSubjectSample"].sum()),
                "n_healthy": int((df["Condition"] == "Healthy").sum()),
                "n_ipf": int((df["Condition"] == "IPF").sum()),
                "n_female": int((df["Sex"] == "female").sum()),
                "n_male": int((df["Sex"] == "male").sum()),
                "n_healthy_female": int((df["Group"] == "Healthy_female").sum()),
                "n_healthy_male": int((df["Group"] == "Healthy_male").sum()),
                "n_ipf_female": int((df["Group"] == "IPF_female").sum()),
                "n_ipf_male": int((df["Group"] == "IPF_male").sum()),
                "age_mean": round(df["Age"].mean(), 2),
                "age_sd": round(df["Age"].std(ddof=1), 2),
                "age_median": round(df["Age"].median(), 2),
                "age_min": round(df["Age"].min(), 2),
                "age_max": round(df["Age"].max(), 2),
            }
        )
    return pd.DataFrame(rows)


def build_binary_imbalance_rows(metadata: pd.DataFrame) -> pd.DataFrame:
    target = metadata.loc[metadata["Set"] == TARGET_STUDY].copy()
    other = metadata.loc[metadata["Set"] != TARGET_STUDY].copy()

    specs = [
        ("Condition", "IPF", metadata["Condition"].eq("IPF")),
        ("Sex", "female", metadata["Sex"].eq("female")),
        ("Group", "Healthy_female", metadata["Group"].eq("Healthy_female")),
        ("Group", "Healthy_male", metadata["Group"].eq("Healthy_male")),
        ("Group", "IPF_female", metadata["Group"].eq("IPF_female")),
        ("Group", "IPF_male", metadata["Group"].eq("IPF_male")),
        ("LibraryLayout", "PAIRED", metadata["LibraryLayout"].eq("PAIRED")),
        ("LibraryLayout", "SINGLE", metadata["LibraryLayout"].eq("SINGLE")),
        ("Repeated sampling", "RepeatedSubjectSample", metadata["RepeatedSubjectSample"]),
    ]

    rows = []
    for domain, level, indicator in specs:
        target_indicator = indicator.loc[target.index].astype(int)
        other_indicator = indicator.loc[other.index].astype(int)
        rows.append(
            {
                "variable_type": "binary",
                "domain": domain,
                "level": level,
                f"{TARGET_STUDY}_n": int(target_indicator.sum()),
                f"{TARGET_STUDY}_denom": int(len(target_indicator)),
                f"{TARGET_STUDY}_prop": round(float(target_indicator.mean()), 4),
                "other_n": int(other_indicator.sum()),
                "other_denom": int(len(other_indicator)),
                "other_prop": round(float(other_indicator.mean()), 4),
                "difference_in_proportion": round(float(target_indicator.mean() - other_indicator.mean()), 4),
                "smd": round(float(smd_binary(target_indicator, other_indicator)), 4),
                "abs_smd": round(abs(float(smd_binary(target_indicator, other_indicator))), 4),
            }
        )
    return pd.DataFrame(rows).sort_values(["abs_smd", "domain", "level"], ascending=[False, True, True])


def build_age_imbalance_rows(metadata: pd.DataFrame) -> pd.DataFrame:
    target = metadata.loc[metadata["Set"] == TARGET_STUDY].copy()
    other = metadata.loc[metadata["Set"] != TARGET_STUDY].copy()

    rows = [
        {
            "subset": "Overall",
            "subset_value": "All samples",
            f"{TARGET_STUDY}_n": len(target),
            f"{TARGET_STUDY}_mean_age": round(target["Age"].mean(), 2),
            f"{TARGET_STUDY}_sd_age": round(target["Age"].std(ddof=1), 2),
            "other_n": len(other),
            "other_mean_age": round(other["Age"].mean(), 2),
            "other_sd_age": round(other["Age"].std(ddof=1), 2),
            "mean_difference": round(float(target["Age"].mean() - other["Age"].mean()), 2),
            "smd": round(float(smd_continuous(target["Age"], other["Age"])), 4),
            "abs_smd": round(abs(float(smd_continuous(target["Age"], other["Age"]))), 4),
        }
    ]

    for subset_col in ["Condition", "Sex", "Group"]:
        for subset_value in sorted(metadata[subset_col].unique()):
            target_subset = target.loc[target[subset_col] == subset_value, "Age"]
            other_subset = other.loc[other[subset_col] == subset_value, "Age"]
            if len(target_subset) < 2 or len(other_subset) < 2:
                continue
            smd = smd_continuous(target_subset, other_subset)
            rows.append(
                {
                    "subset": subset_col,
                    "subset_value": subset_value,
                    f"{TARGET_STUDY}_n": len(target_subset),
                    f"{TARGET_STUDY}_mean_age": round(target_subset.mean(), 2),
                    f"{TARGET_STUDY}_sd_age": round(target_subset.std(ddof=1), 2),
                    "other_n": len(other_subset),
                    "other_mean_age": round(other_subset.mean(), 2),
                    "other_sd_age": round(other_subset.std(ddof=1), 2),
                    "mean_difference": round(float(target_subset.mean() - other_subset.mean()), 2),
                    "smd": round(float(smd), 4),
                    "abs_smd": round(abs(float(smd)), 4),
                }
            )

    return pd.DataFrame(rows).sort_values(["abs_smd", "subset", "subset_value"], ascending=[False, True, True])


def save_group_tables(metadata: pd.DataFrame) -> None:
    for column, stem in [
        ("Condition", "condition_counts"),
        ("Sex", "sex_counts"),
        ("Group", "group_counts"),
        ("LibraryLayout", "library_layout_counts"),
    ]:
        table = pd.crosstab(metadata["Set"], metadata[column])
        table.to_csv(TABLES_DIR / f"{stem}_gse213001_vs_others.csv")


def write_summary_text(
    metadata: pd.DataFrame,
    binary_imbalance: pd.DataFrame,
    age_imbalance: pd.DataFrame,
) -> None:
    top_binary = binary_imbalance.head(5)
    top_age = age_imbalance.head(5)
    target_n = int((metadata["Set"] == TARGET_STUDY).sum())
    other_n = int((metadata["Set"] != TARGET_STUDY).sum())

    lines = [
        "Study imbalance summary",
        f"Target study: {TARGET_STUDY}",
        f"Retained-sample comparison: {TARGET_STUDY} (n={target_n}) vs other studies combined (n={other_n})",
        "",
        "Top binary/categorical imbalances by absolute standardized mean difference:",
    ]
    for _, row in top_binary.iterrows():
        lines.append(
            f"- {row['domain']} = {row['level']}: {TARGET_STUDY} proportion={row[f'{TARGET_STUDY}_prop']:.3f}, "
            f"other proportion={row['other_prop']:.3f}, diff={row['difference_in_proportion']:.3f}, abs_smd={row['abs_smd']:.3f}"
        )
    lines.append("")
    lines.append("Top age imbalances by absolute standardized mean difference:")
    for _, row in top_age.iterrows():
        lines.append(
            f"- {row['subset']} | {row['subset_value']}: {TARGET_STUDY} mean age={row[f'{TARGET_STUDY}_mean_age']:.2f}, "
            f"other mean age={row['other_mean_age']:.2f}, diff={row['mean_difference']:.2f}, abs_smd={row['abs_smd']:.3f}"
        )

    (RESULTS_DIR / "imbalance_summary.txt").write_text("\n".join(lines) + "\n")


def make_love_plot(binary_imbalance: pd.DataFrame, age_imbalance: pd.DataFrame) -> None:
    plot_rows = []
    for row in binary_imbalance.itertuples(index=False):
        plot_rows.append(
            {
                "label": f"{row.domain}: {row.level}",
                "abs_smd": row.abs_smd,
                "family": "Categorical",
            }
        )
    for row in age_imbalance.itertuples(index=False):
        plot_rows.append(
            {
                "label": f"Age | {row.subset}: {row.subset_value}",
                "abs_smd": row.abs_smd,
                "family": "Age",
            }
        )

    plot_df = pd.DataFrame(plot_rows).sort_values("abs_smd", ascending=True).tail(12)

    plt.figure(figsize=(10, 6))
    colors = plot_df["family"].map({"Categorical": "#b53a2f", "Age": "#2e6f9e"})
    plt.barh(plot_df["label"], plot_df["abs_smd"], color=colors)
    plt.axvline(0.1, color="gray", linestyle="--", linewidth=1)
    plt.axvline(0.2, color="gray", linestyle=":", linewidth=1)
    plt.xlabel("Absolute standardized mean difference")
    plt.ylabel("")
    plt.title(f"{TARGET_STUDY} vs other studies")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "gse213001_vs_others_love_plot.pdf")
    plt.close()


def main() -> None:
    ensure_dirs()
    metadata = load_retained_metadata()
    metadata.to_csv(TABLES_DIR / "retained_metadata_with_sets.csv", index=False)

    study_summary = build_study_summary(metadata)
    study_summary.to_csv(TABLES_DIR / "study_level_demographics.csv", index=False)

    set_summary = build_set_summary(metadata)
    set_summary.to_csv(TABLES_DIR / "gse213001_vs_others_summary.csv", index=False)

    save_group_tables(metadata)

    binary_imbalance = build_binary_imbalance_rows(metadata)
    binary_imbalance.to_csv(TABLES_DIR / "categorical_imbalance_metrics.csv", index=False)

    age_imbalance = build_age_imbalance_rows(metadata)
    age_imbalance.to_csv(TABLES_DIR / "age_imbalance_metrics.csv", index=False)

    write_summary_text(metadata, binary_imbalance, age_imbalance)
    make_love_plot(binary_imbalance, age_imbalance)

    print(f"Wrote study imbalance analysis to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
