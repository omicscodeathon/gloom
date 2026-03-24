"""
step13_feature_importance.py - Feature Importance
Methods: MDI (tree models), |coefficients| (LR), permutation importance.
Consolidates into ranked importance table with feature groups.
Output: feature_importance.csv + importance plots.
"""
import logging, sys, warnings
from pathlib import Path
import joblib
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.inspection import permutation_importance

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

FEATURE_GROUPS = {
    "Tumor Expression Stats"  : ["tumor_mean","tumor_median","tumor_std","tumor_iqr","tumor_cv","tumor_skewness","tumor_kurtosis","tumor_pct_expressed"],
    "Normal Expression Stats" : ["normal_mean","normal_median","normal_std","normal_iqr","normal_cv","normal_skewness","normal_kurtosis","normal_pct_expressed"],
    "Differential Expression" : ["log2fc","abs_log2fc","neg_log10_padj","cohens_d","t_stat"],
    "Contrast Features"       : ["tumor_normal_mean_ratio","std_ratio","iqr_ratio"],
    "Rank Features"           : ["tumor_mean_rank","abs_log2fc_rank","neg_log10_padj_rank","cohens_d_rank","tumor_iqr_rank"],
    "Network Topology"        : ["degree","weighted_degree","avg_neighbor_degree","betweenness_centrality","closeness_centrality","eigenvector_centrality","clustering_coefficient"],
    "Network Edge Weights"    : ["mean_edge_weight","max_edge_weight","min_edge_weight","std_edge_weight"],
    "Network Component"       : ["in_largest_component","component_size"],
}

def normalise_importance(values):
    values = np.clip(values, 0, None)
    v_min = values.min(); v_max = values.max(); denom = v_max - v_min + 1e-12
    return (values - v_min) / denom

def get_native_importance(model, model_name, feature_names):
    if hasattr(model, "feature_importances_"):
        raw = model.feature_importances_; label = f"{model_name} (MDI)"
    elif hasattr(model, "coef_"):
        raw = np.abs(model.coef_[0]); label = f"{model_name} (|coef|)"
    else:
        return None
    s = pd.Series(normalise_importance(raw), index=feature_names, name=label)
    return s

def get_permutation_importance(model, model_name, use_scaled, X_val, X_val_scaled, y_val, n_repeats=10):
    X = X_val_scaled if use_scaled else X_val
    log.info(f"  [{model_name}] Permutation importance (n_repeats={n_repeats}) …")
    result = permutation_importance(model, X.values, y_val, n_repeats=n_repeats,
                                    random_state=config.RANDOM_STATE, scoring="roc_auc", n_jobs=-1)
    s = pd.Series(normalise_importance(result.importances_mean), index=X.columns, name=f"{model_name} (permutation)")
    log.info(f"  Top 5: {s.nlargest(5).to_dict()}")
    return s

def run_feature_importance():
    log.info("="*60); log.info("STEP 13 — FEATURE IMPORTANCE"); log.info("="*60)
    X_val        = pd.read_csv(config.VAL_FEATURES_FILE, index_col=0)
    X_val_scaled = pd.read_csv(config.PROCESSED_DIR/"val_features_scaled.csv", index_col=0)
    y_val        = pd.read_csv(config.VAL_LABELS_FILE).set_index("gene")["label"].values
    feature_names = X_val.columns.tolist()
    best_name = (config.MODELS_DIR/"best_model_name.txt").read_text().strip() if (config.MODELS_DIR/"best_model_name.txt").exists() else "random_forest"
    model_registry = {
        "random_forest":(False,"model_random_forest.joblib"),
        "gradient_boosting":(False,"model_gradient_boosting.joblib"),
        "extra_trees":(False,"model_extra_trees.joblib"),
        "logistic_regression":(True,"model_logistic_regression.joblib"),
    }
    loaded = {}
    for name,(use_scaled,fname) in model_registry.items():
        fpath = config.MODELS_DIR/fname
        if fpath.exists():
            loaded[name] = (joblib.load(fpath), use_scaled)
    all_series = []
    for name,(model,_) in loaded.items():
        s = get_native_importance(model, name, feature_names)
        if s is not None: all_series.append(s)
    if best_name in loaded:
        best_model, best_use_scaled = loaded[best_name]
        perm = get_permutation_importance(best_model, best_name, best_use_scaled, X_val, X_val_scaled, y_val)
        all_series.append(perm)
    # Build table
    df = pd.DataFrame(index=feature_names); df.index.name = "feature"
    for s in all_series:
        if s is not None:
            df[s.name] = s.reindex(feature_names).fillna(0.0)
    method_cols = df.columns.tolist()
    df["mean_importance"] = df[method_cols].mean(axis=1)
    df["std_importance"]  = df[method_cols].std(axis=1)
    feat_to_group = {f: g for g, feats in FEATURE_GROUPS.items() for f in feats}
    df["feature_group"] = [feat_to_group.get(f,"Other") for f in df.index]
    df = df.sort_values("mean_importance", ascending=False)
    df.insert(0, "rank", range(1, len(df)+1))

    # Grouped importance plot
    group_means = df.groupby("feature_group")["mean_importance"].mean().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(10,5))
    colors = plt.cm.Set2(np.linspace(0,1,len(group_means)))
    bars = ax.bar(range(len(group_means)), group_means.values, color=colors, alpha=0.85)
    ax.set_xticks(range(len(group_means))); ax.set_xticklabels(group_means.index, rotation=25, ha="right", fontsize=9)
    ax.set_ylabel("Mean Normalised Importance"); ax.set_title("Feature Importance by Group")
    for bar,val in zip(bars,group_means.values):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.005, f"{val:.4f}", ha="center", va="bottom", fontsize=8)
    plt.tight_layout(); fig.savefig(config.FIGURES_DIR/"feature_importance_grouped.png",dpi=130,bbox_inches="tight"); plt.close(fig)

    df.to_csv(config.FEATURE_IMPORTANCE_FILE)
    log.info(f"  Feature importance shape: {df.shape}")
    log.info(f"  Top feature: {df.index[0]} ({df['mean_importance'].iloc[0]:.4f})")
    log.info("STEP 13 COMPLETE")
    return df

if __name__ == "__main__":
    imp = run_feature_importance()
    print(imp[["rank","mean_importance","std_importance","feature_group"]].head(15).to_string())
