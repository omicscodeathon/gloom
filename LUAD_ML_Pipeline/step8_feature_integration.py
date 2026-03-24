"""
step8_feature_integration.py - Feature Integration
Concatenates expression (Step 5) and network (Step 7) features.
Removes zero-variance/high-missing features, applies RobustScaler.
Outputs: integrated_features.csv, integrated_features_scaled.csv
"""
import logging, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

def run_feature_integration():
    log.info("="*60); log.info("STEP 8 — FEATURE INTEGRATION"); log.info("="*60)
    expr_features = pd.read_csv(config.EXPR_FEATURES_FILE,    index_col=0)
    net_features  = pd.read_csv(config.NETWORK_FEATURES_FILE, index_col=0)
    log.info(f"  Expression: {expr_features.shape}  Network: {net_features.shape}")

    # Align to common genes
    common = expr_features.index.intersection(net_features.index)
    expr_features = expr_features.loc[common]
    net_features  = net_features.loc[common]

    # Concatenate
    features_raw = pd.concat([expr_features, net_features], axis=1)
    log.info(f"  Concatenated: {features_raw.shape}")

    # Remove duplicate columns
    if features_raw.columns.duplicated().any():
        features_raw = features_raw.loc[:, ~features_raw.columns.duplicated(keep="first")]

    # Remove zero-variance
    variances = features_raw.var(axis=0)
    low_var   = variances[variances < 1e-8].index.tolist()
    if low_var:
        log.warning(f"  Removing {len(low_var)} zero-variance features.")
        features_raw = features_raw.drop(columns=low_var)

    # Fill NaN
    features_raw.replace([np.inf, -np.inf], np.nan, inplace=True)
    for col in features_raw.columns:
        if features_raw[col].isna().any():
            features_raw[col].fillna(features_raw[col].median(), inplace=True)

    # Feature correlations (flag highly correlated pairs)
    corr_matrix = features_raw.corr(method="pearson")
    corr_pairs  = []
    cols = corr_matrix.columns.tolist()
    for i in range(len(cols)):
        for j in range(i+1, len(cols)):
            r = corr_matrix.iloc[i,j]
            if abs(r) >= 0.95:
                corr_pairs.append({"feature_a": cols[i], "feature_b": cols[j], "pearson_r": round(r,4)})
    pairs_df = pd.DataFrame(corr_pairs).sort_values("pearson_r", key=abs, ascending=False) if corr_pairs else pd.DataFrame()
    log.info(f"  Highly correlated pairs (|r|>=0.95): {len(pairs_df)}")

    # RobustScaler
    scaler = RobustScaler()
    scaled_arr = scaler.fit_transform(features_raw.values)
    features_scaled = pd.DataFrame(scaled_arr, index=features_raw.index, columns=features_raw.columns)

    # PCA plot
    pca    = PCA(n_components=3, random_state=config.RANDOM_STATE)
    pcs    = pca.fit_transform(scaled_arr)
    color_vals = features_raw["abs_log2fc_rank"].values if "abs_log2fc_rank" in features_raw.columns else np.zeros(len(features_raw))
    fig, axes = plt.subplots(1,2, figsize=(14,5))
    sc = axes[0].scatter(pcs[:,0], pcs[:,1], c=color_vals, cmap="RdYlBu_r", s=4, alpha=0.5, linewidths=0)
    axes[0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)")
    axes[0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)")
    axes[0].set_title("PCA — PC1 vs PC2")
    plt.colorbar(sc, ax=axes[0], fraction=0.03)
    sc2 = axes[1].scatter(pcs[:,0], pcs[:,2], c=color_vals, cmap="RdYlBu_r", s=4, alpha=0.5, linewidths=0)
    axes[1].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)")
    axes[1].set_ylabel(f"PC3 ({pca.explained_variance_ratio_[2]*100:.1f}% var)")
    axes[1].set_title("PCA — PC1 vs PC3")
    plt.colorbar(sc2, ax=axes[1], fraction=0.03)
    fig.suptitle("PCA of Integrated Feature Matrix", fontsize=11)
    plt.tight_layout(); fig.savefig(config.FIGURES_DIR/"integrated_feature_pca.png", dpi=130, bbox_inches="tight"); plt.close(fig)

    # Save
    features_raw.to_csv(config.INTEGRATED_FEATURES_FILE)
    features_scaled.to_csv(config.PROCESSED_DIR/"integrated_features_scaled.csv")
    if not pairs_df.empty:
        pairs_df.to_csv(config.PROCESSED_DIR/"highly_correlated_features.csv", index=False)
    log.info(f"  Integrated features: {features_raw.shape}")
    log.info("STEP 8 COMPLETE")
    return {"features": features_raw, "features_scaled": features_scaled}

if __name__ == "__main__":
    r = run_feature_integration()
    print(f"Shape: {r['features'].shape}")
    print(list(r['features'].columns))
