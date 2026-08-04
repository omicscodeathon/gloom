"""
step5_expression_features.py
-----------------------------
Constructs per-gene expression feature matrix from:
  (A) Tumor expression statistics
  (B) Normal expression statistics
  (C) DE metrics
  (D) Contrast/ratio features
  (E) Percentile-rank features
Output: expression_features.csv
"""
import logging, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import skew, kurtosis

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)
EPSILON = 1e-6

def compute_expression_stats(expr_df, prefix):
    log.info(f"  Computing stats for [{prefix}] ({expr_df.shape}) …")
    vals = expr_df.values.astype(float)
    mean_v   = np.nanmean(vals, axis=1)
    median_v = np.nanmedian(vals, axis=1)
    std_v    = np.nanstd(vals, axis=1, ddof=1)
    q75      = np.nanpercentile(vals, 75, axis=1)
    q25      = np.nanpercentile(vals, 25, axis=1)
    iqr_v    = q75 - q25
    cv_v     = np.where(np.abs(mean_v) > EPSILON, std_v / np.abs(mean_v), 0.0)
    skew_v   = np.array([skew(r[~np.isnan(r)]) if np.sum(~np.isnan(r)) > 2 else 0.0 for r in vals])
    kurt_v   = np.array([kurtosis(r[~np.isnan(r)], fisher=True) if np.sum(~np.isnan(r)) > 3 else 0.0 for r in vals])
    log_thresh = np.log2(config.MIN_EXPRESSION_VALUE + 1)
    pct_expr   = np.nanmean(vals > log_thresh, axis=1)
    return pd.DataFrame({
        f"{prefix}_mean": mean_v, f"{prefix}_median": median_v,
        f"{prefix}_std": std_v,   f"{prefix}_iqr": iqr_v,
        f"{prefix}_cv": cv_v,     f"{prefix}_skewness": skew_v,
        f"{prefix}_kurtosis": kurt_v, f"{prefix}_pct_expressed": pct_expr,
    }, index=expr_df.index)

def extract_de_features(de_df):
    required = ["log2fc","neg_log10_padj","cohens_d","t_stat"]
    df = de_df[required].copy()
    df["abs_log2fc"] = de_df["log2fc"].abs()
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.fillna(0.0, inplace=True)
    return df

def compute_contrast_features(t_stats, n_stats):
    c = pd.DataFrame(index=t_stats.index)
    c["tumor_normal_mean_ratio"] = t_stats["tumor_mean"] / (n_stats["normal_mean"].abs() + EPSILON)
    c["std_ratio"] = t_stats["tumor_std"] / (n_stats["normal_std"] + EPSILON)
    c["iqr_ratio"] = t_stats["tumor_iqr"] / (n_stats["normal_iqr"] + EPSILON)
    for col in c.columns:
        c[col] = c[col].clip(lower=c[col].quantile(0.001), upper=c[col].quantile(0.999))
    return c

def compute_rank_features(feature_df, cols):
    rank_df = pd.DataFrame(index=feature_df.index)
    n = len(feature_df)
    for col in cols:
        if col in feature_df.columns:
            rank_df[f"{col}_rank"] = feature_df[col].rank(method="average", na_option="bottom") / n
    return rank_df

def impute_and_validate(features, label="features"):
    features.replace([np.inf, -np.inf], np.nan, inplace=True)
    for col in features.columns:
        if features[col].isna().any():
            features[col].fillna(features[col].median(), inplace=True)
    assert not features.isna().any().any()
    log.info(f"  [{label}] Validated. Shape: {features.shape}")
    return features

def run_expression_feature_construction():
    log.info("="*60); log.info("STEP 5 — EXPRESSION FEATURE CONSTRUCTION"); log.info("="*60)
    tumor_expr  = pd.read_csv(config.TUMOR_EXPR_HARMONIZED,  index_col=0)
    normal_expr = pd.read_csv(config.NORMAL_EXPR_HARMONIZED, index_col=0)
    de_df       = pd.read_csv(config.DE_RESULTS_FILE,         index_col=0)
    de_df       = de_df.reindex(tumor_expr.index)

    t_stats    = compute_expression_stats(tumor_expr,  "tumor")
    n_stats    = compute_expression_stats(normal_expr, "normal")
    de_feats   = extract_de_features(de_df)
    contrast   = compute_contrast_features(t_stats, n_stats)
    tmp        = pd.concat([t_stats, de_feats], axis=1)
    rank_feats = compute_rank_features(tmp, ["tumor_mean","abs_log2fc","neg_log10_padj","cohens_d","tumor_iqr"])

    features = pd.concat([t_stats, n_stats, de_feats, contrast, rank_feats], axis=1)
    features = impute_and_validate(features, "expression_features")

    fig, ax = plt.subplots(figsize=(14,12))
    corr = features.iloc[:, :30].corr()
    im = ax.imshow(corr.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(corr.columns))); ax.set_yticks(range(len(corr.columns)))
    ax.set_xticklabels(corr.columns, rotation=90, fontsize=7)
    ax.set_yticklabels(corr.columns, fontsize=7)
    ax.set_title("Expression Feature Correlation Matrix", fontsize=12)
    plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    plt.tight_layout(); fig.savefig(config.FIGURES_DIR/"expr_feature_correlation_heatmap.png", dpi=130, bbox_inches="tight"); plt.close(fig)

    features.to_csv(config.EXPR_FEATURES_FILE)
    log.info(f"  Expression features: {features.shape}")
    log.info("STEP 5 COMPLETE")
    return features

if __name__ == "__main__":
    f = run_expression_feature_construction()
    print(f"Shape: {f.shape}")
    print(list(f.columns))
