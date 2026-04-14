"""
step11_model_training.py - Model Training
Trains: RandomForest, GradientBoosting, LogisticRegression, SVM, ExtraTrees.
Uses stratified k-fold CV.
Best model selected by config.CV_METRIC_PRIMARY (default: auprc).
Optional SMOTE resampling inside each CV fold (config.USE_SMOTE).
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
from sklearn.metrics import (average_precision_score, f1_score,
                              matthews_corrcoef, roc_auc_score)
from sklearn.model_selection import StratifiedKFold

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
        "extra_trees":         (ExtraTreesClassifier(n_estimators=500, max_features="sqrt", class_weight="balanced", n_jobs=-1, random_state=config.RANDOM_STATE), False),
    }

def compute_sample_weights(y):
    classes, counts = np.unique(y, return_counts=True)
    n_samples = len(y); n_classes = len(classes)
    weight_map = {cls: n_samples/(n_classes*cnt) for cls, cnt in zip(classes, counts)}
    return np.array([weight_map[yi] for yi in y], dtype=np.float64)

def run_model_training():
    log.info("="*60); log.info("STEP 11 - MODEL TRAINING"); log.info("="*60)
    X_train        = pd.read_csv(config.TRAIN_FEATURES_FILE, index_col=0)
    X_train_scaled = pd.read_csv(config.PROCESSED_DIR/"train_features_scaled.csv", index_col=0)
    y_train        = pd.read_csv(config.TRAIN_LABELS_FILE).set_index("gene")["label"]
    log.info(f"  X_train: {X_train.shape}  pos={y_train.sum()}")
    skf = StratifiedKFold(n_splits=config.CV_FOLDS, shuffle=True, random_state=config.RANDOM_STATE)
    model_defs = get_model_definitions()
    cv_results = []; trained_models = {}
    total_start = time.time()
    # Try to load resampling classes once; warn early if missing
    _smote_cls = None
    if config.USE_SMOTE:
        try:
            from imblearn.over_sampling import SMOTE as _SMOTE
            _smote_cls = _SMOTE
            log.info("  SMOTE enabled (imbalanced-learn found).")
        except ImportError:
            log.warning("  USE_SMOTE=True but imbalanced-learn is not installed - SMOTE skipped. "
                        "Install with: pip install imbalanced-learn")

    _undersample_cls = None
    if getattr(config, "USE_UNDERSAMPLING", False):
        try:
            from imblearn.under_sampling import RandomUnderSampler as _RUS
            _undersample_cls = _RUS
            log.info("  Random undersampling enabled (imbalanced-learn found).")
        except ImportError:
            log.warning("  USE_UNDERSAMPLING=True but imbalanced-learn is not installed - "
                        "undersampling skipped.  Install with: pip install imbalanced-learn")

    for model_name, (model, use_scaled) in model_defs.items():
        X = X_train_scaled if use_scaled else X_train
        y = y_train.values
        fold_aurocs = []; fold_auprcs = []; fold_f1s = []; fold_mccs = []
        log.info(f"\n  [{model_name}] CV ...")
        for fold_idx, (tr_idx, va_idx) in enumerate(skf.split(X, y)):
            Xf, Xv, yf, yv = X.iloc[tr_idx], X.iloc[va_idx], y[tr_idx], y[va_idx]

            # Step 5 - resampling inside training fold only
            # First: SMOTE oversampling (if enabled)
            if _smote_cls is not None:
                try:
                    sm = _smote_cls(random_state=config.RANDOM_STATE)
                    Xf_arr, yf = sm.fit_resample(Xf.values, yf)
                    Xf = pd.DataFrame(Xf_arr, columns=Xf.columns)
                except Exception as smote_err:
                    log.warning(f"    SMOTE failed fold {fold_idx+1}: {smote_err}; skipping.")
            # Second: random undersampling (applied after SMOTE when both enabled)
            if _undersample_cls is not None:
                try:
                    rus = _undersample_cls(random_state=config.RANDOM_STATE)
                    Xf_arr, yf = rus.fit_resample(Xf.values, yf)
                    Xf = pd.DataFrame(Xf_arr, columns=Xf.columns)
                except Exception as rus_err:
                    log.warning(f"    Undersampling failed fold {fold_idx+1}: {rus_err}; skipping.")

            if model_name == "gradient_boosting":
                model.fit(Xf, yf, sample_weight=compute_sample_weights(yf))
            else:
                model.fit(Xf, yf)
            y_prob = model.predict_proba(Xv)[:,1] if hasattr(model,"predict_proba") else model.decision_function(Xv)
            y_pred = (y_prob >= 0.5).astype(int)
            auroc = roc_auc_score(yv, y_prob)
            auprc = average_precision_score(yv, y_prob)
            f1    = f1_score(yv, y_pred, zero_division=0)
            mcc   = matthews_corrcoef(yv, y_pred)
            fold_aurocs.append(auroc); fold_auprcs.append(auprc)
            fold_f1s.append(f1);      fold_mccs.append(mcc)
            log.info(f"    Fold {fold_idx+1} AUROC={auroc:.4f} AUPRC={auprc:.4f} F1={f1:.4f} MCC={mcc:.4f}")
        mean_auroc = np.mean(fold_aurocs); std_auroc = np.std(fold_aurocs)
        mean_auprc = np.mean(fold_auprcs); std_auprc = np.std(fold_auprcs)
        mean_f1    = np.mean(fold_f1s);   std_f1    = np.std(fold_f1s)
        mean_mcc   = np.mean(fold_mccs);  std_mcc   = np.std(fold_mccs)
        log.info(f"  [{model_name}] AUROC={mean_auroc:.4f}±{std_auroc:.4f}  "
                 f"AUPRC={mean_auprc:.4f}±{std_auprc:.4f}  "
                 f"F1={mean_f1:.4f}±{std_f1:.4f}  MCC={mean_mcc:.4f}±{std_mcc:.4f}")
        cv_results.append({"model_name":model_name,
                            "fold_aurocs":fold_aurocs,"fold_auprcs":fold_auprcs,
                            "fold_f1s":fold_f1s,"fold_mccs":fold_mccs,
                            "mean_auroc":mean_auroc,"std_auroc":std_auroc,
                            "mean_auprc":mean_auprc,"std_auprc":std_auprc,
                            "mean_f1":mean_f1,"std_f1":std_f1,
                            "mean_mcc":mean_mcc,"std_mcc":std_mcc})
        # Retrain on full set
        if model_name == "gradient_boosting":
            model.fit(X, y, sample_weight=compute_sample_weights(y))
        else:
            model.fit(X, y)
        trained_models[model_name] = (model, use_scaled)
        joblib.dump(model, config.MODELS_DIR/f"model_{model_name}.joblib")
    total_elapsed = time.time() - total_start
    log.info(f"\nAll models trained in {total_elapsed:.1f}s")
    # Step 3 - select best model by primary metric (config.CV_METRIC_PRIMARY)
    primary = getattr(config, "CV_METRIC_PRIMARY", "auprc")
    primary_scores = {r["model_name"]: r[f"mean_{primary}"] for r in cv_results}
    best_name      = max(primary_scores, key=primary_scores.get)
    best_model, _  = trained_models[best_name]
    joblib.dump(best_model, config.MODELS_DIR/"best_model.joblib")
    (config.MODELS_DIR/"best_model_name.txt").write_text(best_name)
    log.info(f"  Best: {best_name} ({primary.upper()}={primary_scores[best_name]:.4f})")

    cv_rows = []
    for r in cv_results:
        row = {"model":r["model_name"],
               "mean_auroc":round(r["mean_auroc"],6),"std_auroc":round(r["std_auroc"],6),
               "mean_auprc":round(r["mean_auprc"],6),"std_auprc":round(r["std_auprc"],6),
               "mean_f1":round(r["mean_f1"],6),      "std_f1":round(r["std_f1"],6),
               "mean_mcc":round(r["mean_mcc"],6),    "std_mcc":round(r["std_mcc"],6)}
        for i,(au,ap,f1,mc) in enumerate(zip(r["fold_aurocs"],r["fold_auprcs"],r["fold_f1s"],r["fold_mccs"]),1):
            row[f"fold{i}_auroc"] = round(au,6); row[f"fold{i}_auprc"] = round(ap,6)
            row[f"fold{i}_f1"]    = round(f1,6); row[f"fold{i}_mcc"]   = round(mc,6)
        cv_rows.append(row)
    cv_results_df = pd.DataFrame(cv_rows).sort_values(f"mean_{primary}", ascending=False).reset_index(drop=True)
    cv_results_df.to_csv(config.MODELS_DIR/"cv_results.csv", index=False)

    # Plot AUROC and AUPRC comparison side-by-side
    models_cv = cv_results_df["model"].tolist()
    colors_cv = ["tomato" if m==best_name else "steelblue" for m in models_cv]
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    for ax, metric, label in [(axes[0],"auroc","CV AUROC"), (axes[1],"auprc","CV AUPRC")]:
        means = cv_results_df[f"mean_{metric}"].tolist()
        stds  = cv_results_df[f"std_{metric}"].tolist()
        ax.bar(range(len(models_cv)), means, yerr=stds, capsize=4,
               color=colors_cv, alpha=0.85, error_kw={"elinewidth":1.2})
        ax.set_xticks(range(len(models_cv)))
        ax.set_xticklabels(models_cv, rotation=20, ha="right", fontsize=9)
        ax.set_ylabel(label); ax.set_title(f"Model Comparison - {label}")
        if metric == primary:
            ax.set_title(f"Model Comparison - {label} ★ (primary)", fontsize=10)
    plt.tight_layout()
    fig.savefig(config.FIGURES_DIR/"cv_auroc_auprc_comparison.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    log.info("STEP 11 COMPLETE")
    return {"trained_models":trained_models,"cv_results":cv_results,"cv_results_df":cv_results_df,
            "best_model_name":best_name,"best_model":best_model}

if __name__ == "__main__":
    r = run_model_training()
    print(r["cv_results_df"][["model","mean_auroc","std_auroc","mean_auprc","std_auprc"]].to_string(index=False))
