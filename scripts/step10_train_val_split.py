"""
step10_train_val_split.py - Train/Validation Split
Stratified 80/20 split + 5-fold CV. RobustScaler fit on train only.
Outputs: train/val features (scaled+unscaled), labels, fold assignments.
"""
import logging, sys
from pathlib import Path
import joblib
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import RobustScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

def run_train_val_split():
    log.info("="*60); log.info("STEP 10 - TRAIN/VALIDATION SPLIT"); log.info("="*60)
    features  = pd.read_csv(config.INTEGRATED_FEATURES_FILE, index_col=0)
    labels_df = pd.read_csv(config.LABELS_FILE)
    labels_df = labels_df.set_index("gene")["label"]
    common    = features.index.intersection(labels_df.index)
    features  = features.loc[common]; labels_df = labels_df.loc[common]
    log.info(f"  Features: {features.shape}  Labels: {len(labels_df)}")

    # Split
    X_train, X_val, y_train, y_val = train_test_split(
        features, labels_df, test_size=config.TEST_SIZE,
        random_state=config.RANDOM_STATE, stratify=labels_df)
    log.info(f"  Train: {len(y_train)} (pos={y_train.sum()})  Val: {len(y_val)} (pos={y_val.sum()})")

    # CV folds
    skf = StratifiedKFold(n_splits=config.CV_FOLDS, shuffle=True, random_state=config.RANDOM_STATE)
    fold_assignments = np.zeros(len(X_train), dtype=int)
    for fold_idx, (_, val_idx) in enumerate(skf.split(X_train, y_train)):
        fold_assignments[val_idx] = fold_idx + 1
    fold_df = pd.DataFrame({"gene": X_train.index, "cv_fold": fold_assignments, "label": y_train.values})

    # Scale (fit on train only)
    scaler = RobustScaler()
    X_train_scaled = pd.DataFrame(scaler.fit_transform(X_train.values), index=X_train.index, columns=X_train.columns)
    X_val_scaled   = pd.DataFrame(scaler.transform(X_val.values),       index=X_val.index,   columns=X_val.columns)

    # PCA plot
    pca = PCA(n_components=2, random_state=config.RANDOM_STATE)
    pca.fit(X_train_scaled.values)
    train_pcs = pca.transform(X_train_scaled.values)
    val_pcs   = pca.transform(X_val_scaled.values)
    fig, axes = plt.subplots(1,2, figsize=(14,5))
    axes[0].scatter(train_pcs[:,0], train_pcs[:,1], s=4, alpha=0.4, color="steelblue", label=f"Train (n={len(y_train):,})", linewidths=0)
    axes[0].scatter(val_pcs[:,0],   val_pcs[:,1],   s=8, alpha=0.7, color="tomato",    label=f"Val (n={len(y_val):,})", linewidths=0)
    axes[0].set_title("Train vs Validation Split"); axes[0].legend(fontsize=8, markerscale=3)
    neg_mask = y_train.values==0; pos_mask = y_train.values==1
    axes[1].scatter(train_pcs[neg_mask,0], train_pcs[neg_mask,1], s=4, alpha=0.3, color="steelblue", label=f"Neg (n={neg_mask.sum():,})", linewidths=0)
    axes[1].scatter(train_pcs[pos_mask,0], train_pcs[pos_mask,1], s=14, alpha=0.8, color="tomato", label=f"Pos (n={pos_mask.sum()})", linewidths=0)
    axes[1].set_title("Train: Positive vs Negative"); axes[1].legend(fontsize=8, markerscale=2)
    fig.suptitle("PCA Overview of Train/Val Split", fontsize=12)
    plt.tight_layout(); fig.savefig(config.FIGURES_DIR/"split_feature_pca.png", dpi=130, bbox_inches="tight"); plt.close(fig)

    # Save
    X_train.to_csv(config.TRAIN_FEATURES_FILE)
    X_val.to_csv(config.VAL_FEATURES_FILE)
    X_train_scaled.to_csv(config.PROCESSED_DIR/"train_features_scaled.csv")
    X_val_scaled.to_csv(config.PROCESSED_DIR/"val_features_scaled.csv")
    pd.DataFrame({"gene":y_train.index,"label":y_train.values}).to_csv(config.TRAIN_LABELS_FILE, index=False)
    pd.DataFrame({"gene":y_val.index,"label":y_val.values}).to_csv(config.VAL_LABELS_FILE, index=False)
    fold_df.to_csv(config.PROCESSED_DIR/"cv_fold_assignments.csv", index=False)
    joblib.dump(scaler, config.MODELS_DIR/"robust_scaler.joblib")
    log.info("STEP 10 COMPLETE")
    return {"X_train":X_train,"X_val":X_val,"X_train_scaled":X_train_scaled,
            "X_val_scaled":X_val_scaled,"y_train":y_train,"y_val":y_val,"skf":skf,"scaler":scaler}

if __name__ == "__main__":
    r = run_train_val_split()
    print(f"Train: {r['X_train'].shape}  Val: {r['X_val'].shape}")
