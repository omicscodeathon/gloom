"""
step2b_qc_report.py — Sample-level QC and cohort-artefact diagnostics
=====================================================================
Runs after step 2 (preprocessing) and produces the diagnostics a reviewer needs to judge
whether a tumor-vs-normal contrast reflects biology or a cohort / pipeline artefact.

Computed from the step-2 matrices (log2(x + 1)), restricted to the genes present in both groups:

  1. Per-sample summary (median, IQR, mean, fraction of zeros, detected genes) and a
     per-group summary of those statistics.
  2. Sample PCA (2 PCs) on the most variable genes — PC scores and explained variance saved as
     CSV, plus a figure; the separation of the groups along PC1 is reported.
  3. Global differential-expression diagnostics from a quick Welch test (a diagnostic only —
     step4 produces the official DE table): median log2FC over all genes, fraction of genes
     that are DE-significant, number up / down and the up:down ratio.
  4. WARNING when any of the following holds (thresholds in config.QC_*):
        |median log2FC| > 1, more than 80 % of genes DE-significant, up:down ratio > 10 or < 0.1.
     The warning is logged and written to results/qc_cohort_warning.txt (section SAMPLE_QC).

Outputs (results/):
  qc_sample_summary.csv, qc_group_summary.csv, qc_sample_pca.csv,
  qc_pca_explained_variance.csv, qc_global_de_summary.csv, qc_cohort_warning.txt (if warned)
Figures: qc_sample_pca.png, qc_sample_medians.png

This step never stops the pipeline (optional).
"""
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
from gloom_utils import (auroc_score, benjamini_hochberg, evaluate_qc_warnings,
                         update_report_section)

config.create_output_dirs()
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)


def _load_matrix(path) -> pd.DataFrame:
    """Read a genes x samples CSV, upper-case the symbols and collapse duplicate genes."""
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).str.strip().str.upper()
    if df.index.duplicated().any():
        df["_mean"] = df.mean(axis=1)
        df = df.sort_values("_mean", ascending=False)
        df = df[~df.index.duplicated(keep="first")].drop(columns=["_mean"])
    return df


def sample_summary(tumor: pd.DataFrame, normal: pd.DataFrame) -> pd.DataFrame:
    """Per-sample median / IQR / mean / fraction of zeros / number of detected genes."""
    log_thr = np.log2(config.MIN_EXPRESSION_VALUE + 1)
    rows = []
    for group, df in (("tumor", tumor), ("normal", normal)):
        q25 = df.quantile(0.25, axis=0)
        q75 = df.quantile(0.75, axis=0)
        frac_zero = (df == 0).sum(axis=0) / df.shape[0]
        detected = (df > log_thr).sum(axis=0)
        for s in df.columns:
            rows.append({
                "sample_id": s, "group": group,
                "median": float(df[s].median()), "iqr": float(q75[s] - q25[s]),
                "mean": float(df[s].mean()), "fraction_zero": float(frac_zero[s]),
                "n_detected_genes": int(detected[s]),
            })
    return pd.DataFrame(rows)


def group_summary(per_sample: pd.DataFrame) -> pd.DataFrame:
    """Median and IQR of the per-sample medians / IQRs, by group."""
    rows = []
    for group, sub in per_sample.groupby("group"):
        rows.append({
            "group": group,
            "n_samples": len(sub),
            "median_of_sample_medians": float(sub["median"].median()),
            "iqr_of_sample_medians": float(sub["median"].quantile(0.75) - sub["median"].quantile(0.25)),
            "median_of_sample_iqrs": float(sub["iqr"].median()),
            "iqr_of_sample_iqrs": float(sub["iqr"].quantile(0.75) - sub["iqr"].quantile(0.25)),
            "median_fraction_zero": float(sub["fraction_zero"].median()),
        })
    return pd.DataFrame(rows)


def sample_pca(tumor: pd.DataFrame, normal: pd.DataFrame, top_genes: int, n_components: int = 2):
    """
    PCA of the samples on the ``top_genes`` most variable genes (genes centred; missing
    values filled with the gene mean).  Returns (scores, explained_variance) DataFrames.
    """
    combined = pd.concat([tumor, normal], axis=1)
    arr = combined.to_numpy(dtype=float)
    row_mean = np.nanmean(np.where(np.isfinite(arr), arr, np.nan), axis=1)
    arr = np.where(np.isnan(arr), row_mean[:, None], arr)
    combined = pd.DataFrame(arr, index=combined.index, columns=combined.columns).dropna(how="any")
    var = combined.var(axis=1)
    keep = var.sort_values(ascending=False).index[:min(top_genes, len(var))]
    X = combined.loc[keep].to_numpy(dtype=float).T            # samples x genes
    X = X - X.mean(axis=0, keepdims=True)
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    explained = (S ** 2) / np.sum(S ** 2)
    k = min(n_components, len(S))
    groups = ["tumor"] * tumor.shape[1] + ["normal"] * normal.shape[1]
    scores = pd.DataFrame(U[:, :k] * S[:k], columns=[f"PC{i + 1}" for i in range(k)])
    scores.insert(0, "group", groups)
    scores.insert(0, "sample_id", list(tumor.columns) + list(normal.columns))
    ev = pd.DataFrame({
        "component": [f"PC{i + 1}" for i in range(k)],
        "explained_variance_ratio": explained[:k],
    })
    return scores, ev


def global_de_diagnostics(tumor: pd.DataFrame, normal: pd.DataFrame) -> dict:
    """Quick vectorised Welch test + BH, summarised globally (diagnostic only)."""
    genes = tumor.index.intersection(normal.index)
    t, n = tumor.loc[genes], normal.loc[genes]
    log2fc = t.mean(axis=1) - n.mean(axis=1)
    _, p = stats.ttest_ind(t.to_numpy(), n.to_numpy(), axis=1, equal_var=False, nan_policy="omit")
    p = np.where(np.isnan(p), 1.0, p)
    padj = benjamini_hochberg(p)
    sig = (padj <= config.DE_PVALUE_THRESHOLD) & (log2fc.abs().to_numpy() >= config.DE_LOG2FC_THRESHOLD)
    up = int((sig & (log2fc.to_numpy() > 0)).sum())
    down = int((sig & (log2fc.to_numpy() < 0)).sum())
    return {
        "n_genes": int(len(genes)),
        "median_log2fc": float(np.nanmedian(log2fc.to_numpy())),
        "mean_log2fc": float(np.nanmean(log2fc.to_numpy())),
        "n_de_significant": int(sig.sum()),
        "fraction_de_significant": float(sig.mean()) if len(sig) else float("nan"),
        "n_up": up,
        "n_down": down,
        "up_down_ratio": (up / down) if down > 0 else (float("inf") if up > 0 else float("nan")),
        "padj_threshold": config.DE_PVALUE_THRESHOLD,
        "abs_log2fc_threshold": config.DE_LOG2FC_THRESHOLD,
    }


def _plot_pca(scores: pd.DataFrame, ev: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    for group, color in (("normal", "darkorange"), ("tumor", "steelblue")):
        sub = scores[scores["group"] == group]
        ax.scatter(sub["PC1"], sub["PC2"], s=12, alpha=0.6, color=color, label=f"{group} (n={len(sub)})",
                   linewidths=0)
    ax.set_xlabel(f"PC1 ({100 * ev['explained_variance_ratio'].iloc[0]:.1f}% var)")
    if len(ev) > 1:
        ax.set_ylabel(f"PC2 ({100 * ev['explained_variance_ratio'].iloc[1]:.1f}% var)")
    ax.set_title("Sample PCA (most variable genes)")
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def _plot_medians(per_sample: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5, 5))
    data = [per_sample.loc[per_sample["group"] == g, "median"].to_numpy() for g in ("tumor", "normal")]
    ax.boxplot(data, patch_artist=True)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["tumor", "normal"])
    ax.set_ylabel("Per-sample median log2(x + 1)")
    ax.set_title("Per-sample median expression")
    plt.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def run_qc_report() -> dict:
    log.info("=" * 60)
    log.info("STEP 2b — SAMPLE-LEVEL QC AND COHORT DIAGNOSTICS")
    log.info("=" * 60)

    tumor = _load_matrix(config.TUMOR_EXPR_PROCESSED)
    normal = _load_matrix(config.NORMAL_EXPR_PROCESSED)
    genes = tumor.index.intersection(normal.index)
    tumor, normal = tumor.loc[genes], normal.loc[genes]
    log.info(f"  Genes shared by both groups: {len(genes):,}  "
             f"(tumor samples: {tumor.shape[1]}, normal samples: {normal.shape[1]})")

    res_dir = Path(config.RESULTS_DIR)
    res_dir.mkdir(parents=True, exist_ok=True)

    per_sample = sample_summary(tumor, normal)
    per_group = group_summary(per_sample)
    per_sample.to_csv(res_dir / "qc_sample_summary.csv", index=False)
    per_group.to_csv(res_dir / "qc_group_summary.csv", index=False)
    log.info("\n" + per_group.to_string(index=False))

    scores, ev = sample_pca(tumor, normal, top_genes=config.QC_PCA_TOP_GENES)
    scores.to_csv(res_dir / "qc_sample_pca.csv", index=False)
    ev.to_csv(res_dir / "qc_pca_explained_variance.csv", index=False)
    sep = auroc_score((scores["group"] == "tumor").astype(int), scores["PC1"])
    sep = max(sep, 1.0 - sep) if np.isfinite(sep) else float("nan")
    log.info(f"  PCA: PC1 explains {100 * ev['explained_variance_ratio'].iloc[0]:.1f}% of the variance; "
             f"PC1 separates tumor from normal with AUROC = {sep:.3f}")
    _plot_pca(scores, ev, Path(config.FIGURES_DIR) / "qc_sample_pca.png")
    _plot_medians(per_sample, Path(config.FIGURES_DIR) / "qc_sample_medians.png")

    diag = global_de_diagnostics(tumor, normal)
    diag["pc1_group_separation_auroc"] = sep
    pd.DataFrame([diag]).to_csv(res_dir / "qc_global_de_summary.csv", index=False)
    log.info(f"  Global DE diagnostics: median log2FC = {diag['median_log2fc']:.3f}; "
             f"{100 * diag['fraction_de_significant']:.1f}% DE-significant "
             f"(up={diag['n_up']}, down={diag['n_down']}, up:down={diag['up_down_ratio']:.3g})")

    warn_lines = evaluate_qc_warnings(
        diag["median_log2fc"], diag["fraction_de_significant"], diag["n_up"], diag["n_down"],
        max_abs_median_log2fc=config.QC_MAX_ABS_MEDIAN_LOG2FC,
        max_fraction_de=config.QC_MAX_FRACTION_DE,
        max_up_down_ratio=config.QC_MAX_UP_DOWN_RATIO,
    )
    report_path = res_dir / "qc_cohort_warning.txt"
    if warn_lines:
        lines = warn_lines + [
            "",
            "These patterns indicate that the contrast may be dominated by cohort / platform /",
            "normalisation differences rather than disease biology.  Do not interpret the DE",
            "and network results as LUAD-specific without a uniformly processed design",
            "(docs/tcga_paired_design.md).",
        ]
        banner = "!" * 70
        log.warning("\n" + banner + "\n  SAMPLE QC WARNING\n" + "\n".join("  " + l for l in lines)
                    + "\n" + banner)
        update_report_section(report_path, "SAMPLE_QC (step2b)", lines)
    else:
        log.info("  Sample QC: no global cohort-artefact warning triggered.")
        update_report_section(report_path, "SAMPLE_QC (step2b)", None)

    log.info("STEP 2b COMPLETE")
    return {"per_sample": per_sample, "per_group": per_group, "pca": scores,
            "explained_variance": ev, "global_de": diag, "warnings": warn_lines}


if __name__ == "__main__":
    out = run_qc_report()
    print(pd.DataFrame([out["global_de"]]).T)
