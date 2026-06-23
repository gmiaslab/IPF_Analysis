#!/usr/bin/env python3

from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

PANEL_SIZES = [10, 25, 50]
C_GRID = [0.01, 0.1, 1.0, 10.0, 100.0]
ELASTIC_NET_L1_GRID = [0.25, 0.5, 0.75]
LOG_CPM_PSEUDOCOUNT = 0.5
POSITIVE_CLASS = "IPF"
MODEL_LABEL = "Elastic Net Logistic Regression"


def get_param_grid() -> list[dict[str, float]]:
    return [
        {"C": c_value, "l1_ratio": l1_ratio}
        for l1_ratio in ELASTIC_NET_L1_GRID
        for c_value in C_GRID
    ]


def format_param_value(value: float) -> str:
    if isinstance(value, float):
        return f"{value:g}"
    return str(value)


def params_to_label(params: dict[str, float]) -> str:
    return ", ".join(f"{key}={format_param_value(value)}" for key, value in params.items())


def expand_param_record(prefix: str, params: dict[str, float]) -> dict[str, float | str]:
    return {
        f"{prefix}_c": params["C"],
        f"{prefix}_l1_ratio": params["l1_ratio"],
        f"{prefix}_param_label": params_to_label(params),
    }


def safe_divide(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return math.nan
    return numerator / denominator


def compute_binary_metrics(y_true: np.ndarray, y_prob: np.ndarray, threshold: float = 0.5) -> dict[str, float]:
    y_true = np.asarray(y_true, dtype=int)
    y_prob = np.asarray(y_prob, dtype=float)
    y_pred = (y_prob >= threshold).astype(int)

    tp = int(np.sum((y_true == 1) & (y_pred == 1)))
    tn = int(np.sum((y_true == 0) & (y_pred == 0)))
    fp = int(np.sum((y_true == 0) & (y_pred == 1)))
    fn = int(np.sum((y_true == 1) & (y_pred == 0)))

    sensitivity = safe_divide(tp, tp + fn)
    specificity = safe_divide(tn, tn + fp)
    precision = safe_divide(tp, tp + fp)
    fpr = safe_divide(fp, fp + tn)
    fnr = safe_divide(fn, fn + tp)
    accuracy = safe_divide(tp + tn, len(y_true))
    balanced_accuracy = np.nanmean([sensitivity, specificity])

    metrics = {
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "n": int(len(y_true)),
        "n_positive": int(np.sum(y_true == 1)),
        "n_negative": int(np.sum(y_true == 0)),
        "accuracy": accuracy,
        "balanced_accuracy": balanced_accuracy,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "precision": precision,
        "fpr": fpr,
        "fnr": fnr,
    }

    if len(np.unique(y_true)) == 2:
        metrics["roc_auc"] = roc_auc_score(y_true, y_prob)
        metrics["average_precision"] = average_precision_score(y_true, y_prob)
    else:
        metrics["roc_auc"] = math.nan
        metrics["average_precision"] = math.nan

    return metrics


def fit_elastic_net_model(
    x_train: pd.DataFrame,
    y_train: pd.Series,
    x_test: pd.DataFrame,
    params: dict[str, float],
) -> tuple[StandardScaler, LogisticRegression, np.ndarray]:
    scaler = StandardScaler()
    x_train_scaled = scaler.fit_transform(x_train)
    x_test_scaled = scaler.transform(x_test)

    model = LogisticRegression(
        penalty="elasticnet",
        C=float(params["C"]),
        l1_ratio=float(params["l1_ratio"]),
        class_weight="balanced",
        solver="saga",
        max_iter=20000,
        random_state=0,
    )

    model.fit(x_train_scaled, y_train)
    y_prob = model.predict_proba(x_test_scaled)[:, 1]
    return scaler, model, y_prob


def choose_elastic_net_params(
    x_train: pd.DataFrame,
    y_train: pd.Series,
    train_studies: pd.Series,
) -> tuple[dict[str, float], pd.DataFrame]:
    inner_records: list[dict[str, float | str]] = []
    unique_studies = list(dict.fromkeys(train_studies.tolist()))

    for params in get_param_grid():
        fold_scores: list[float] = []
        for held_out_study in unique_studies:
            inner_train_mask = train_studies != held_out_study
            inner_val_mask = train_studies == held_out_study

            _, model, y_prob = fit_elastic_net_model(
                x_train.loc[inner_train_mask],
                y_train.loc[inner_train_mask],
                x_train.loc[inner_val_mask],
                params,
            )
            metrics = compute_binary_metrics(y_train.loc[inner_val_mask].to_numpy(), y_prob)
            fold_scores.append(metrics["balanced_accuracy"])
            record = {
                "held_out_inner_study": held_out_study,
                "balanced_accuracy": metrics["balanced_accuracy"],
                "sensitivity": metrics["sensitivity"],
                "specificity": metrics["specificity"],
                "precision": metrics["precision"],
                "roc_auc": metrics["roc_auc"],
                "n_validation": metrics["n"],
            }
            record.update(expand_param_record("candidate", params))
            inner_records.append(record)

        mean_score = float(np.nanmean(fold_scores))
        mean_record = {
            "held_out_inner_study": "mean",
            "balanced_accuracy": mean_score,
            "sensitivity": math.nan,
            "specificity": math.nan,
            "precision": math.nan,
            "roc_auc": math.nan,
            "n_validation": int(np.sum(train_studies.isin(unique_studies))),
        }
        mean_record.update(expand_param_record("candidate", params))
        inner_records.append(mean_record)

    inner_df = pd.DataFrame(inner_records)
    mean_df = inner_df[inner_df["held_out_inner_study"] == "mean"].copy()
    best_row = mean_df.sort_values(
        ["balanced_accuracy", "candidate_l1_ratio", "candidate_c"],
        ascending=[False, True, True],
    ).iloc[0]
    best_params = {
        "C": float(best_row["candidate_c"]),
        "l1_ratio": float(best_row["candidate_l1_ratio"]),
    }
    return best_params, inner_df


def aggregate_subject_counts(metadata: pd.DataFrame, counts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = metadata.copy()
    metadata["unit_id"] = (
        metadata["GEO"].astype(str)
        + "::"
        + metadata["Subject"].astype(str)
        + "::"
        + metadata["Condition"].astype(str)
        + "::"
        + metadata["Sex"].astype(str)
    )

    grouped = metadata.groupby("unit_id", sort=False)
    unit_meta = grouped.agg(
        GEO=("GEO", "first"),
        Subject=("Subject", "first"),
        Condition=("Condition", "first"),
        Sex=("Sex", "first"),
        n_runs=("Run", "size"),
        runs=("Run", lambda x: ";".join(x)),
    )

    counts_by_unit = counts.T.assign(unit_id=metadata["unit_id"].to_numpy())
    counts_by_unit = counts_by_unit.groupby("unit_id", sort=False).sum().T
    counts_by_unit = counts_by_unit[unit_meta.index]

    return unit_meta, counts_by_unit


def compute_log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    library_sizes = counts.sum(axis=0)
    cpm = counts.div(library_sizes, axis=1) * 1e6
    return np.log2(cpm + LOG_CPM_PSEUDOCOUNT)


def build_fixed_panels(de_table: pd.DataFrame) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    ranked = de_table.loc[
        (de_table["adj.P.Val"] < 0.05) & (de_table["logFC"].abs() >= 1)
    ].copy()
    ranked["panel_rank"] = np.arange(1, len(ranked) + 1)

    panels: dict[str, pd.DataFrame] = {}
    records: list[pd.DataFrame] = []

    for panel_size in PANEL_SIZES:
        panel_name = f"top_{panel_size}"
        panel_df = ranked.head(panel_size).copy()
        panel_df["panel_name"] = panel_name
        panels[panel_name] = panel_df
        records.append(panel_df)

    return panels, pd.concat(records, ignore_index=True)


def render_metric_heatmaps(metrics_df: pd.DataFrame, output_path: Path) -> None:
    plot_metrics = ["balanced_accuracy", "sensitivity", "specificity", "precision"]
    panel_order = [f"top_{n}" for n in PANEL_SIZES]
    study_order = metrics_df["held_out_study"].drop_duplicates().tolist()

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    cmap = plt.get_cmap("YlOrRd")

    for ax, metric in zip(axes.flat, plot_metrics):
        pivot = (
            metrics_df.pivot(index="held_out_study", columns="panel_name", values=metric)
            .reindex(index=study_order, columns=panel_order)
        )
        data = pivot.to_numpy(dtype=float)
        im = ax.imshow(data, vmin=0, vmax=1, cmap=cmap, aspect="auto")
        ax.set_title(metric.replace("_", " ").title())
        ax.set_xticks(np.arange(len(panel_order)))
        ax.set_xticklabels(panel_order, rotation=35, ha="right")
        ax.set_yticks(np.arange(len(study_order)))
        ax.set_yticklabels(study_order)

        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                value = data[i, j]
                label = "NA" if math.isnan(value) else f"{value:.2f}"
                ax.text(j, i, label, ha="center", va="center", fontsize=9)

        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle(f"Leave-One-Study-Out IPF Classifier Performance by Fixed Panel\n{MODEL_LABEL}", fontsize=14)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def render_macro_metric_plot(summary_df: pd.DataFrame, output_path: Path) -> None:
    panel_order = [f"top_{n}" for n in PANEL_SIZES]
    metric_order = ["balanced_accuracy", "sensitivity", "specificity", "precision"]
    colors = {
        "balanced_accuracy": "#0072B2",
        "sensitivity": "#D55E00",
        "specificity": "#009E73",
        "precision": "#E69F00",
    }

    fig, ax = plt.subplots(figsize=(11, 6))
    x = np.arange(len(panel_order))
    width = 0.18

    for idx, metric in enumerate(metric_order):
        values = (
            summary_df.set_index("panel_name")
            .reindex(panel_order)[metric]
            .to_numpy(dtype=float)
        )
        ax.bar(x + (idx - 1.5) * width, values, width=width, label=metric.replace("_", " ").title(), color=colors[metric])

    ax.set_xticks(x)
    ax.set_xticklabels(panel_order, rotation=35, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Metric value")
    ax.set_title(f"Macro-Averaged LOSO Metrics by Fixed Panel\n{MODEL_LABEL}")
    ax.legend(frameon=False, ncol=2)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def summarize_confusion_cells(metrics_df: pd.DataFrame) -> pd.DataFrame:
    cell_order = ["tn", "fp", "fn", "tp"]
    records: list[dict[str, float | int | str]] = []
    denominator_map = {
        "tn": "n_negative",
        "fp": "n_negative",
        "fn": "n_positive",
        "tp": "n_positive",
    }

    for panel_name, panel_metrics in metrics_df.groupby("panel_name", sort=False):
        for cell_name in cell_order:
            denominators = panel_metrics[denominator_map[cell_name]]
            study_percentages = 100 * panel_metrics[cell_name] / denominators
            records.append(
                {
                    "panel_name": panel_name,
                    "cell": cell_name.upper(),
                    "mean_study_pct": float(study_percentages.mean()),
                    "sd_study_pct": float(study_percentages.std(ddof=1)),
                    "n_studies": int(panel_metrics.shape[0]),
                }
            )

    return pd.DataFrame(records)


def render_confusion_summary(confusion_df: pd.DataFrame, output_path: Path) -> None:
    panel_order = [f"top_{n}" for n in PANEL_SIZES]
    row_labels = ["True Healthy", "True IPF"]
    col_labels = ["Predicted Healthy", "Predicted IPF"]
    cell_positions = {
        "TN": (0, 0),
        "FP": (0, 1),
        "FN": (1, 0),
        "TP": (1, 1),
    }

    fig, axes = plt.subplots(1, len(panel_order), figsize=(5 * len(panel_order), 4.8), constrained_layout=True)
    if len(panel_order) == 1:
        axes = [axes]
    cmap = plt.get_cmap("Blues")
    vmax = float(confusion_df["mean_study_pct"].max()) if not confusion_df.empty else 100.0

    for ax, panel_name in zip(axes, panel_order):
        panel_df = confusion_df.loc[confusion_df["panel_name"] == panel_name].copy()
        matrix = np.full((2, 2), np.nan)
        annotation = np.full((2, 2), "", dtype=object)

        for _, row in panel_df.iterrows():
            i, j = cell_positions[row["cell"]]
            matrix[i, j] = row["mean_study_pct"]
            annotation[i, j] = (
                f"{row['cell']}\n"
                f"{row['mean_study_pct']:.1f}% +/- {row['sd_study_pct']:.1f}%"
            )

        im = ax.imshow(matrix, vmin=0, vmax=vmax, cmap=cmap)
        ax.set_xticks(np.arange(2))
        ax.set_xticklabels(col_labels)
        ax.set_yticks(np.arange(2))
        ax.set_yticklabels(row_labels)
        ax.set_title(panel_name.replace("_", " ").title())

        for i in range(2):
            for j in range(2):
                value = matrix[i, j]
                text_color = "white" if value >= (vmax * 0.55) else "#0b2545"
                ax.text(
                    j,
                    i,
                    annotation[i, j],
                    ha="center",
                    va="center",
                    fontsize=10.5,
                    fontweight="bold",
                    color=text_color,
                )

    fig.colorbar(im, ax=axes, fraction=0.035, pad=0.03, label="Mean row-normalized held-out study percentage")
    fig.suptitle(f"Confusion Summary Across Held-Out Studies\nMean row percentage +/- study-level SD\n{MODEL_LABEL}", fontsize=14)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def write_summary_text(
    output_path: Path,
    unit_meta: pd.DataFrame,
    panel_summary: pd.DataFrame,
    pooled_summary: pd.DataFrame,
) -> None:
    best_panel = panel_summary.sort_values("balanced_accuracy", ascending=False).iloc[0]
    lines = [
        f"Fixed-signature leave-one-study-out IPF classifier summary ({MODEL_LABEL})",
        "",
        f"Subject-level units analyzed: {len(unit_meta)}",
        "Held-out studies: " + ", ".join(panel_summary['held_out_studies'].iloc[0].split(";")),
        f"Best macro balanced-accuracy panel: {best_panel['panel_name']} ({best_panel['balanced_accuracy']:.3f})",
        f"Corresponding macro sensitivity/specificity/precision: {best_panel['sensitivity']:.3f} / {best_panel['specificity']:.3f} / {best_panel['precision']:.3f}",
        "",
        "Macro-average metrics by panel:",
    ]

    for _, row in panel_summary.iterrows():
        lines.append(
            f"- {row['panel_name']}: balanced_accuracy={row['balanced_accuracy']:.3f}, "
            f"sensitivity={row['sensitivity']:.3f}, specificity={row['specificity']:.3f}, "
            f"precision={row['precision']:.3f}, roc_auc={row['roc_auc']:.3f}"
        )

    lines.extend(["", "Pooled held-out predictions by panel:"])
    for _, row in pooled_summary.iterrows():
        lines.append(
            f"- {row['panel_name']}: tp={int(row['tp'])}, fp={int(row['fp'])}, tn={int(row['tn'])}, fn={int(row['fn'])}, "
            f"balanced_accuracy={row['balanced_accuracy']:.3f}, roc_auc={row['roc_auc']:.3f}"
        )

    output_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    results_dir = script_dir / "ipf_classifier_results"
    tables_dir = results_dir / "tables"
    figures_dir = results_dir / "figures"
    for path in (results_dir, tables_dir, figures_dir):
        path.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(script_dir / "ipf_analysis_results" / "qc" / "analysis_metadata_used.csv")
    counts = pd.read_csv(script_dir / "combined_study_matrix.csv", index_col=0)
    de_table = pd.read_csv(script_dir / "ipf_analysis_results" / "tables" / "de_disease_average.csv")

    metadata = metadata.loc[metadata["Run"].isin(counts.columns)].copy()
    metadata = metadata.sort_values("Run").reset_index(drop=True)
    counts = counts.loc[:, metadata["Run"]]

    unit_meta, unit_counts = aggregate_subject_counts(metadata, counts)
    expr = compute_log_cpm(unit_counts)

    panels, panel_definitions = build_fixed_panels(de_table)
    panel_definitions.to_csv(tables_dir / "fixed_panel_definitions.csv", index=False)
    unit_meta.to_csv(tables_dir / "subject_level_metadata.csv")

    study_order = unit_meta["GEO"].drop_duplicates().tolist()
    y_all = (unit_meta["Condition"] == POSITIVE_CLASS).astype(int)

    outer_metrics_records: list[dict[str, float | int | str]] = []
    prediction_records: list[pd.DataFrame] = []
    coefficient_records: list[pd.DataFrame] = []
    tuning_records: list[pd.DataFrame] = []

    for panel_name, panel_df in panels.items():
        gene_ids = [gene_id for gene_id in panel_df["gene_id"] if gene_id in expr.index]
        x_panel = expr.loc[gene_ids].T.copy()
        gene_labels = panel_df.set_index("gene_id").loc[gene_ids, "gene_label"]

        for held_out_study in study_order:
            train_mask = unit_meta["GEO"] != held_out_study
            test_mask = unit_meta["GEO"] == held_out_study

            x_train = x_panel.loc[train_mask]
            x_test = x_panel.loc[test_mask]
            y_train = y_all.loc[train_mask]
            y_test = y_all.loc[test_mask]
            train_studies = unit_meta.loc[train_mask, "GEO"]

            best_params, tuning_df = choose_elastic_net_params(x_train, y_train, train_studies)
            tuning_df["panel_name"] = panel_name
            tuning_df["outer_held_out_study"] = held_out_study
            tuning_records.append(tuning_df)

            scaler, model, y_prob = fit_elastic_net_model(x_train, y_train, x_test, best_params)
            metrics = compute_binary_metrics(y_test.to_numpy(), y_prob)
            metrics.update(
                {
                    "panel_name": panel_name,
                    "panel_size": len(gene_ids),
                    "held_out_study": held_out_study,
                    "n_train_units": int(train_mask.sum()),
                    "n_test_units": int(test_mask.sum()),
                }
            )
            metrics.update(expand_param_record("chosen", best_params))
            outer_metrics_records.append(metrics)

            pred_df = unit_meta.loc[test_mask, ["GEO", "Subject", "Condition", "Sex", "n_runs", "runs"]].copy()
            pred_df["unit_id"] = pred_df.index
            pred_df["panel_name"] = panel_name
            pred_df["y_true"] = y_test.to_numpy()
            pred_df["predicted_probability_ipf"] = y_prob
            pred_df["predicted_class"] = np.where(y_prob >= 0.5, "IPF", "Healthy")
            for key, value in expand_param_record("chosen", best_params).items():
                pred_df[key] = value
            prediction_records.append(pred_df.reset_index(drop=True))

            coef_df = pd.DataFrame(
                {
                    "panel_name": panel_name,
                    "held_out_study": held_out_study,
                    "gene_id": gene_ids,
                    "gene_label": gene_labels.to_numpy(),
                    "coefficient": model.coef_.ravel(),
                }
            )
            for key, value in expand_param_record("chosen", best_params).items():
                coef_df[key] = value
            coefficient_records.append(coef_df)

    outer_metrics = pd.DataFrame(outer_metrics_records)
    predictions = pd.concat(prediction_records, ignore_index=True)
    coefficients = pd.concat(coefficient_records, ignore_index=True)
    tuning = pd.concat(tuning_records, ignore_index=True)

    macro_summary = (
        outer_metrics.groupby("panel_name", as_index=False)
        .agg(
            panel_size=("panel_size", "first"),
            balanced_accuracy=("balanced_accuracy", "mean"),
            sensitivity=("sensitivity", "mean"),
            specificity=("specificity", "mean"),
            precision=("precision", "mean"),
            accuracy=("accuracy", "mean"),
            roc_auc=("roc_auc", "mean"),
            average_precision=("average_precision", "mean"),
        )
    )
    macro_summary["held_out_studies"] = ";".join(study_order)

    pooled_records: list[dict[str, float | int | str]] = []
    for panel_name in predictions["panel_name"].drop_duplicates():
        pred_subset = predictions.loc[predictions["panel_name"] == panel_name]
        pooled = compute_binary_metrics(
            pred_subset["y_true"].to_numpy(),
            pred_subset["predicted_probability_ipf"].to_numpy(),
        )
        pooled["panel_name"] = panel_name
        pooled["panel_size"] = int((outer_metrics.loc[outer_metrics["panel_name"] == panel_name, "panel_size"]).iloc[0])
        pooled_records.append(pooled)

    pooled_summary = pd.DataFrame(pooled_records)
    confusion_summary = summarize_confusion_cells(outer_metrics)

    outer_metrics.to_csv(tables_dir / "loso_metrics_by_study.csv", index=False)
    predictions.to_csv(tables_dir / "loso_predictions_by_subject.csv", index=False)
    coefficients.to_csv(tables_dir / "loso_coefficients.csv", index=False)
    tuning.to_csv(tables_dir / "loso_inner_tuning.csv", index=False)
    macro_summary.to_csv(tables_dir / "loso_macro_summary.csv", index=False)
    pooled_summary.to_csv(tables_dir / "loso_pooled_summary.csv", index=False)
    confusion_summary.to_csv(tables_dir / "loso_confusion_summary.csv", index=False)

    render_metric_heatmaps(outer_metrics, figures_dir / "loso_metric_heatmaps.pdf")
    render_macro_metric_plot(macro_summary, figures_dir / "loso_macro_metrics.pdf")
    render_confusion_summary(confusion_summary, figures_dir / "loso_confusion_summary.pdf")
    write_summary_text(results_dir / "classifier_summary.txt", unit_meta, macro_summary, pooled_summary)

    print(f"Wrote classifier analysis to {results_dir}")


if __name__ == "__main__":
    main()
