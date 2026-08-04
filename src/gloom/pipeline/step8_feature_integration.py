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

    # P1.3 FIX: Remove confirmed perfectly collinear features (|r| >= 0.999).
    # These pairs carry zero marginal information and inflate Logistic Regression
    # coefficient variance via multicollinearity. One from each pair is dropped;
    # the retained feature is noted in the comment.
    REDUNDANT_FEATURES = [
        # P1.3: perfectly collinear pairs
        "weighted_degree",   # r=0.9996 with degree          → keep degree
        "component_size",    # r=1.0000 with in_largest_component → keep in_largest_component
        "min_edge_weight",   # r=0.9991 with mean_edge_weight → keep mean_edge_weight
        "tumor_median",      # r=0.9985 with tumor_mean      → keep tumor_mean
        "normal_median",     # r=0.9988 with normal_mean     → keep normal_mean
        "t_stat",            # r=0.9994 with cohens_d        → keep cohens_d
        # Rank features: computed on full gene universe before train/val split
        # → data leakage; also near-duplicate of their raw counterparts (r>0.96)
        "tumor_iqr_rank",       # dup of tumor_iqr
        "tumor_mean_rank",      # dup of tumor_mean
        "cohens_d_rank",        # dup of cohens_d
        "abs_log2fc_rank",      # dup of abs_log2fc
        "neg_log10_padj_rank",  # dup of neg_log10_padj
        # Zero-importance features: set to 0 by density guards (graph too dense)
        # or structurally redundant with higher-ranked features
        "betweenness_centrality",  # zeroed — tumor graph too dense (>500K edges)
        "tumor_betweenness",       # zeroed — same density limit
        "normal_betweenness",      # zeroed — normal graph too dense
        "delta_betweenness",       # zeroed — derived from above two zeros
        "normal_clustering",       # zeroed — normal graph too dense (>300K edges)
        "delta_clustering",        # zeroed — derived from above
        "clustering_coefficient",  # 0 importance; tumor_clustering covers same concept
        "max_edge_weight",         # r=0.993 with mean_edge_weight; 0 importance
        "std_edge_weight",         # 0 importance; mean_edge_weight sufficient
        "tumor_degree",            # 0 importance; redundant with degree
        "in_largest_component",    # 0 importance
        "log2fc",                  # 0 importance; abs_log2fc already present
        "iqr_ratio",               # 0 importance
        "std_ratio",               # 0 importance
        "tumor_kurtosis",          # 0 importance
        "normal_kurtosis",         # 0 importance
    ]
    to_drop = [c for c in REDUNDANT_FEATURES if c in features_raw.columns]
    if to_drop:
        features_raw = features_raw.drop(columns=to_drop)
        log.info(f"  Removed {len(to_drop)} collinear features: {to_drop}")
    log.info(f"  Features after collinearity removal: {features_raw.shape[1]}")

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
    color_vals = features_raw["abs_log2fc"].values if "abs_log2fc" in features_raw.columns else np.zeros(len(features_raw))
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

    # P3.2: Optionally include differential co-expression features (step7b output)
    diff_feat_path = getattr(config, "DIFFERENTIAL_NETWORK_FEATURES_FILE", None)
    if diff_feat_path and Path(diff_feat_path).exists():
        diff_features = pd.read_csv(diff_feat_path, index_col=0)
        common_diff   = features_raw.index.intersection(diff_features.index)
        diff_features = diff_features.reindex(features_raw.index).fillna(0.0)
        # Remove any collinear differential columns before concatenating
        # Exclude any diff feature already in REDUNDANT_FEATURES (they were zeroed
        # out by density guards and would bypass the earlier collinearity drop).
        new_cols = [c for c in diff_features.columns
                    if c not in features_raw.columns and c not in REDUNDANT_FEATURES]
        features_raw = pd.concat([features_raw, diff_features[new_cols]], axis=1)
        log.info(f"  Differential network features added: {len(new_cols)} columns "
                 f"(coverage: {len(common_diff):,}/{len(features_raw):,} genes)")
    else:
        log.info("  Differential network features not found — run step6b + step7b to include them.")

    # Save
    features_raw.to_csv(config.INTEGRATED_FEATURES_FILE)
    features_scaled_final = pd.DataFrame(
        RobustScaler().fit_transform(features_raw.values),
        index=features_raw.index, columns=features_raw.columns,
    )
    features_scaled_final.to_csv(config.PROCESSED_DIR/"integrated_features_scaled.csv")
    if not pairs_df.empty:
        pairs_df.to_csv(config.PROCESSED_DIR/"highly_correlated_features.csv", index=False)
    log.info(f"  Integrated features: {features_raw.shape}")
    log.info("STEP 8 COMPLETE")
    return {"features": features_raw, "features_scaled": features_scaled_final}

if __name__ == "__main__":
    r = run_feature_integration()
    print(f"Shape: {r['features'].shape}")
    print(list(r['features'].columns))
