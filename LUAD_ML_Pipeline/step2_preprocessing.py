"""
step2_preprocessing.py
----------------------
Performs preprocessing and quality control on the raw-loaded
tumor and normal expression matrices.

Operations:
  1. Replace non-positive values with NaN
  2. Log2-transform expression values (log2(x + 1))
  3. Filter low-expression genes
  4. Filter low-variance genes (bottom 10% IQR)
  5. Validate sample overlap between expression matrix and metadata
  6. Generate QC summary statistics and plots

Outputs (to data/processed/):
  - tumor_expression_processed.csv
  - normal_expression_processed.csv
  - tumor_metadata_processed.csv
  - normal_metadata_processed.csv
  - qc_summary.csv
"""

import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import skew, kurtosis

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

config.create_output_dirs()

logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(config.LOG_FILE),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)


def replace_nonpositive_with_nan(df, label):
    n_zero = (df == 0).sum().sum()
    n_neg  = (df  < 0).sum().sum()
    df = df.where(df > 0, other=np.nan)
    log.info(f"  [{label}] Replaced {n_zero} zeros and {n_neg} negatives with NaN.")
    return df


def log2_transform(df, label):
    transformed = np.log2(df + 1)
    log.info(f"  [{label}] log2(x+1) applied. Range: [{transformed.min().min():.3f}, {transformed.max().max():.3f}]")
    return transformed


def filter_low_expression_genes(df, label,
                                 min_fraction=config.MIN_EXPRESSION_FRACTION,
                                 min_value=config.MIN_EXPRESSION_VALUE):
    log_threshold  = np.log2(min_value + 1)
    n_before       = df.shape[0]
    expressed_frac = (df > log_threshold).sum(axis=1) / df.shape[1]
    mask           = expressed_frac >= min_fraction
    df_filtered    = df.loc[mask]
    log.info(f"  [{label}] Low-expr filter: {n_before} → {df_filtered.shape[0]} genes.")
    return df_filtered


def filter_low_variance_genes(df, label, bottom_percentile=10.0):
    n_before  = df.shape[0]
    iqr       = df.quantile(0.75, axis=1) - df.quantile(0.25, axis=1)
    threshold = np.percentile(iqr.dropna(), bottom_percentile)
    mask      = iqr >= threshold
    df_f      = df.loc[mask]
    log.info(f"  [{label}] Low-var filter: {n_before} → {df_f.shape[0]} genes.")
    return df_f


def validate_sample_overlap(expr_df, meta_df, label):
    expr_samples = set(expr_df.columns)
    meta_samples = set(meta_df.index)
    common       = expr_samples & meta_samples
    if len(common) == 0:
        log.warning(f"  [{label}] No common samples — trying TCGA ID truncation.")
        short_map    = {s: s[:12] for s in expr_df.columns}
        expr_df_short = expr_df.rename(columns=short_map)
        common2 = set(expr_df_short.columns) & meta_samples
        if len(common2) > 0:
            expr_df = expr_df_short
            common  = common2
        else:
            raise ValueError(f"[{label}] Zero samples overlap even after truncation.")
    common_sorted = sorted(common)
    expr_aligned  = expr_df[common_sorted]
    meta_aligned  = meta_df.loc[common_sorted]
    log.info(f"  [{label}] Retained {len(common_sorted)} common samples.")
    return expr_aligned, meta_aligned


def check_minimum_samples(df, label, min_n):
    n = df.shape[1]
    if n < min_n:
        raise ValueError(f"[{label}] Only {n} samples — minimum is {min_n}.")
    log.info(f"  [{label}] Sample count OK: {n} >= {min_n}.")


def plot_sample_distributions(df_before, df_after, label, out_path, n_samples=50):
    rng = np.random.default_rng(config.RANDOM_STATE)
    def _sample_medians(df):
        cols = rng.choice(df.columns, size=min(n_samples, df.shape[1]), replace=False)
        return df[cols].median(axis=0)
    med_before = _sample_medians(df_before)
    med_after  = _sample_medians(df_after)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    fig.suptitle(f"{label.capitalize()} — Sample Median Expression", fontsize=12)
    for ax, medians, title in zip(axes,
        [med_before, med_after],
        ["Before QC (raw)", "After QC (log2-transformed)"]):
        ax.violinplot(medians, positions=[0], showmedians=True)
        ax.boxplot(medians, positions=[1], patch_artist=True,
                   boxprops=dict(facecolor="steelblue", alpha=0.6))
        ax.set_title(title, fontsize=10)
        ax.set_ylabel("Median expression value")
        ax.set_xticks([])
    plt.tight_layout()
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def plot_gene_filter_summary(stats, out_path):
    stages = ["Raw loaded", "After low-expr filter", "After low-var filter", "After sample alignment"]
    tumor_counts  = [stats["tumor"].get(s, 0)  for s in stages]
    normal_counts = [stats["normal"].get(s, 0) for s in stages]
    x = np.arange(len(stages))
    width = 0.35
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(x - width/2, tumor_counts,  width, label="Tumor",  color="steelblue")
    ax.bar(x + width/2, normal_counts, width, label="Normal", color="darkorange")
    ax.set_xticks(x)
    ax.set_xticklabels(stages, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("Number of genes")
    ax.set_title("Gene counts at each QC filtering stage")
    ax.legend()
    plt.tight_layout()
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def run_preprocessing() -> dict:
    log.info("=" * 60)
    log.info("STEP 2 — DATA PREPROCESSING & QUALITY CONTROL")
    log.info("=" * 60)

    tumor_expr_raw  = pd.read_csv(config.PROCESSED_DIR / "tumor_expression_raw.csv",  index_col=0)
    normal_expr_raw = pd.read_csv(config.PROCESSED_DIR / "normal_expression_raw.csv", index_col=0)
    tumor_meta      = pd.read_csv(config.PROCESSED_DIR / "tumor_metadata_raw.csv",    index_col=0)
    normal_meta     = pd.read_csv(config.PROCESSED_DIR / "normal_metadata_raw.csv",   index_col=0)

    gene_stats = {"tumor": {"Raw loaded": tumor_expr_raw.shape[0]},
                  "normal": {"Raw loaded": normal_expr_raw.shape[0]}}

    log.info("\n--- Tumor preprocessing ---")
    tumor_expr = replace_nonpositive_with_nan(tumor_expr_raw.copy(), "tumor")
    tumor_expr = log2_transform(tumor_expr, "tumor")
    tumor_expr = filter_low_expression_genes(tumor_expr, "tumor")
    gene_stats["tumor"]["After low-expr filter"] = tumor_expr.shape[0]
    tumor_expr = filter_low_variance_genes(tumor_expr, "tumor")
    gene_stats["tumor"]["After low-var filter"] = tumor_expr.shape[0]
    tumor_expr, tumor_meta = validate_sample_overlap(tumor_expr, tumor_meta, "tumor")
    gene_stats["tumor"]["After sample alignment"] = tumor_expr.shape[0]
    check_minimum_samples(tumor_expr, "tumor", config.MIN_SAMPLES_TUMOR)

    log.info("\n--- Normal preprocessing ---")
    normal_expr = replace_nonpositive_with_nan(normal_expr_raw.copy(), "normal")
    normal_expr = log2_transform(normal_expr, "normal")
    normal_expr = filter_low_expression_genes(normal_expr, "normal")
    gene_stats["normal"]["After low-expr filter"] = normal_expr.shape[0]
    normal_expr = filter_low_variance_genes(normal_expr, "normal")
    gene_stats["normal"]["After low-var filter"] = normal_expr.shape[0]
    gene_stats["normal"]["After sample alignment"] = normal_expr.shape[0]
    check_minimum_samples(normal_expr, "normal", config.MIN_SAMPLES_NORMAL)

    log.info("\nGenerating QC plots …")
    plot_sample_distributions(tumor_expr_raw, tumor_expr, "tumor",
                               config.FIGURES_DIR / "qc_tumor_sample_distributions.png")
    plot_sample_distributions(normal_expr_raw, normal_expr, "normal",
                               config.FIGURES_DIR / "qc_normal_sample_distributions.png")
    plot_gene_filter_summary(gene_stats, config.FIGURES_DIR / "qc_gene_filter_summary.png")

    rows = []
    for dataset in ("tumor", "normal"):
        for stage, count in gene_stats[dataset].items():
            rows.append({"dataset": dataset, "stage": stage, "n_genes": count})
    rows.append({"dataset": "tumor",  "stage": "Final samples", "n_genes": tumor_expr.shape[1]})
    rows.append({"dataset": "normal", "stage": "Final samples", "n_genes": normal_expr.shape[1]})
    qc_summary = pd.DataFrame(rows)

    tumor_expr.to_csv(config.TUMOR_EXPR_PROCESSED)
    normal_expr.to_csv(config.NORMAL_EXPR_PROCESSED)
    tumor_meta.to_csv(config.PROCESSED_DIR / "tumor_metadata_processed.csv")
    normal_meta.to_csv(config.PROCESSED_DIR / "normal_metadata_processed.csv")
    qc_summary.to_csv(config.PROCESSED_DIR / "qc_summary.csv", index=False)

    log.info("\n" + "=" * 60)
    log.info("STEP 2 SUMMARY")
    log.info("=" * 60)
    log.info(f"  Tumor  expr : {tumor_expr.shape[0]:>6} genes x {tumor_expr.shape[1]:>4} samples")
    log.info(f"  Normal expr : {normal_expr.shape[0]:>6} genes x {normal_expr.shape[1]:>4} samples")
    log.info("STEP 2 COMPLETE")
    return {"tumor_expr": tumor_expr, "tumor_meta": tumor_meta,
            "normal_expr": normal_expr, "normal_meta": normal_meta,
            "qc_summary": qc_summary}


if __name__ == "__main__":
    result = run_preprocessing()
    print("\n--- QC Summary ---")
    print(result["qc_summary"].to_string(index=False))
