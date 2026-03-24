"""
step11_model_training.py - Model Training
Trains: RandomForest, GradientBoosting, LogisticRegression, SVM, ExtraTrees.
Uses stratified k-fold CV. Selects best model by AUROC.
Outputs: model_*.joblib, best_model.joblib, cv_results.csv
"""
import logging, sys, time, warnings
from pathlib import Path
import joblib
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier, GradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

def get_model_definitions():
    return {
        "random_forest":       (RandomForestClassifier(n_estimators=500, max_features="sqrt", class_weight="balanced", n_jobs=-1, random_state=config.RANDOM_STATE, oob_score=True), False),
        "gradient_boosting":   (GradientBoostingClassifier(n_estimators=300, learning_rate=0.05, max_depth=4, subsample=0.8, max_features="sqrt", random_state=config.RANDOM_STATE), False),
        "logistic_regression": (LogisticRegression(penalty="l2", C=1.0, class_weight="balanced", solver="lbfgs", max_iter=2000, random_state=config.RANDOM_STATE, n_jobs=-1), True),
        "svm":                 (SVC(kernel="rbf", C=1.0, gamma="scale", class_weight="balanced", probability=True, random_state=config.RANDOM_STATE), True),
        "extra_trees":         (ExtraTreesClassifier(n_estimators=500, max_features="sqrt", class_weight="balanced", n_jobs=-1, random_state=config.RANDOM_STATE), False),
    }

def compute_sample_weights(y):
    classes, counts = np.unique(y, return_counts=True)
    n_samples = len(y); n_classes = len(classes)
    weight_map = {cls: n_samples/(n_classes*cnt) for cls, cnt in zip(classes, counts)}
    return np.array([weight_map[yi] for yi in y], dtype=np.float64)

def run_model_training():
    log.info("="*60); log.info("STEP 11 — MODEL TRAINING"); log.info("="*60)
    X_train        = pd.read_csv(config.TRAIN_FEATURES_FILE, index_col=0)
    X_train_scaled = pd.read_csv(config.PROCESSED_DIR/"train_features_scaled.csv", index_col=0)
    y_train        = pd.read_csv(config.TRAIN_LABELS_FILE).set_index("gene")["label"]
    log.info(f"  X_train: {X_train.shape}  pos={y_train.sum()}")
    skf = StratifiedKFold(n_splits=config.CV_FOLDS, shuffle=True, random_state=config.RANDOM_STATE)
    model_defs = get_model_definitions()
    cv_results = []; trained_models = {}
    total_start = time.time()
    for model_name, (model, use_scaled) in model_defs.items():
        X = X_train_scaled if use_scaled else X_train
        y = y_train.values
        fold_aurocs = []; fold_auprcs = []
        log.info(f"\n  [{model_name}] CV ...")
        for fold_idx, (tr_idx, va_idx) in enumerate(skf.split(X, y)):
            Xf, Xv, yf, yv = X.iloc[tr_idx], X.iloc[va_idx], y[tr_idx], y[va_idx]
            if model_name == "gradient_boosting":
                model.fit(Xf, yf, sample_weight=compute_sample_weights(yf))
            else:
                model.fit(Xf, yf)
            y_prob = model.predict_proba(Xv)[:,1] if hasattr(model,"predict_proba") else model.decision_function(Xv)
            auroc = roc_auc_score(yv, y_prob); auprc = average_precision_score(yv, y_prob)
            fold_aurocs.append(auroc); fold_auprcs.append(auprc)
            log.info(f"    Fold {fold_idx+1} AUROC={auroc:.4f} AUPRC={auprc:.4f}")
        mean_auroc = np.mean(fold_aurocs); std_auroc = np.std(fold_aurocs)
        mean_auprc = np.mean(fold_auprcs); std_auprc = np.std(fold_auprcs)
        log.info(f"  [{model_name}] AUROC={mean_auroc:.4f}+-{std_auroc:.4f}  AUPRC={mean_auprc:.4f}+-{std_auprc:.4f}")
        cv_results.append({"model_name":model_name,"fold_aurocs":fold_aurocs,"fold_auprcs":fold_auprcs,
                            "mean_auroc":mean_auroc,"std_auroc":std_auroc,"mean_auprc":mean_auprc,"std_auprc":std_auprc})
        # Retrain on full set
        if model_name == "gradient_boosting":
            model.fit(X, y, sample_weight=compute_sample_weights(y))
        else:
            model.fit(X, y)
        trained_models[model_name] = (model, use_scaled)
        joblib.dump(model, config.MODELS_DIR/f"model_{model_name}.joblib")
    total_elapsed = time.time() - total_start
    log.info(f"\nAll models trained in {total_elapsed:.1f}s")
    mean_aurocs = {r["model_name"]: r["mean_auroc"] for r in cv_results}
    best_name   = max(mean_aurocs, key=mean_aurocs.get)
    best_model, _ = trained_models[best_name]
    joblib.dump(best_model, config.MODELS_DIR/"best_model.joblib")
    (config.MODELS_DIR/"best_model_name.txt").write_text(best_name)
    log.info(f"  Best: {best_name} (AUROC={mean_aurocs[best_name]:.4f})")

    cv_rows = []
    for r in cv_results:
        row = {"model":r["model_name"],"mean_auroc":round(r["mean_auroc"],6),"std_auroc":round(r["std_auroc"],6),
               "mean_auprc":round(r["mean_auprc"],6),"std_auprc":round(r["std_auprc"],6)}
        for i,(au,ap) in enumerate(zip(r["fold_aurocs"],r["fold_auprcs"]),1):
            row[f"fold{i}_auroc"] = round(au,6); row[f"fold{i}_auprc"] = round(ap,6)
        cv_rows.append(row)
    cv_results_df = pd.DataFrame(cv_rows).sort_values("mean_auroc", ascending=False).reset_index(drop=True)
    cv_results_df.to_csv(config.MODELS_DIR/"cv_results.csv", index=False)

    # Plot AUROC comparison
    fig, ax = plt.subplots(figsize=(10,5))
    models_cv = cv_results_df["model"].tolist(); means_cv = cv_results_df["mean_auroc"].tolist(); stds_cv = cv_results_df["std_auroc"].tolist()
    colors_cv = ["tomato" if m==best_name else "steelblue" for m in models_cv]
    ax.bar(range(len(models_cv)), means_cv, yerr=stds_cv, capsize=4, color=colors_cv, alpha=0.85, error_kw={"elinewidth":1.2})
    ax.set_xticks(range(len(models_cv))); ax.set_xticklabels(models_cv, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("CV AUROC"); ax.set_title("Model Comparison — CV AUROC")
    plt.tight_layout(); fig.savefig(config.FIGURES_DIR/"cv_auroc_comparison.png", dpi=130, bbox_inches="tight"); plt.close(fig)

    log.info("STEP 11 COMPLETE")
    return {"trained_models":trained_models,"cv_results":cv_results,"cv_results_df":cv_results_df,
            "best_model_name":best_name,"best_model":best_model}

if __name__ == "__main__":
    r = run_model_training()
    print(r["cv_results_df"][["model","mean_auroc","std_auroc","mean_auprc","std_auprc"]].to_string(index=False))
