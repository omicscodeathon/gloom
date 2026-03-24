"""
step12_model_evaluation.py - Model Evaluation
Evaluates all models on held-out validation set.
Metrics: AUROC, AUPRC, F1, MCC, Brier Score.
Outputs: model_metrics.csv + ROC/PR/confusion plots.
"""
import logging, sys
from pathlib import Path
import joblib
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (average_precision_score, brier_score_loss,
    confusion_matrix, f1_score, matthews_corrcoef, precision_recall_curve,
    precision_score, recall_score, roc_auc_score, roc_curve)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

def find_optimal_threshold(y_true, y_prob):
    precisions, recalls, thresholds = precision_recall_curve(y_true, y_prob)
    f1_scores = np.where((precisions+recalls)>0, 2*precisions*recalls/(precisions+recalls), 0.0)
    best_idx  = np.argmax(f1_scores[:-1])
    return float(thresholds[best_idx]), float(f1_scores[best_idx])

def evaluate_single_model(model_name, model, use_scaled, X_val, X_val_scaled, y_val):
    X = X_val_scaled if use_scaled else X_val
    if hasattr(model,"predict_proba"):
        y_prob = model.predict_proba(X)[:,1]
    else:
        s = model.decision_function(X)
        y_prob = (s-s.min())/(s.max()-s.min()+1e-12)
    opt_thresh, opt_f1 = find_optimal_threshold(y_val, y_prob)
    y_pred_opt = (y_prob>=opt_thresh).astype(int)
    y_pred_05  = (y_prob>=0.5).astype(int)
    auroc = roc_auc_score(y_val, y_prob)
    auprc = average_precision_score(y_val, y_prob)
    brier = brier_score_loss(y_val, y_prob)
    mcc   = matthews_corrcoef(y_val, y_pred_opt)
    prec  = precision_score(y_val, y_pred_opt, zero_division=0)
    rec   = recall_score(y_val, y_pred_opt, zero_division=0)
    f1_05 = f1_score(y_val, y_pred_05, zero_division=0)
    fpr, tpr, _ = roc_curve(y_val, y_prob)
    prec_c, rec_c, _ = precision_recall_curve(y_val, y_prob)
    cm = confusion_matrix(y_val, y_pred_opt)
    log.info(f"  [{model_name:<22}] AUROC={auroc:.4f} AUPRC={auprc:.4f} MCC={mcc:.4f} F1={opt_f1:.4f}")
    return {"model_name":model_name,"auroc":round(auroc,4),"auprc":round(auprc,4),
            "brier_score":round(brier,4),"mcc":round(mcc,4),"precision_opt":round(prec,4),
            "recall_opt":round(rec,4),"f1_opt":round(opt_f1,4),"f1_at_05":round(f1_05,4),
            "optimal_threshold":round(opt_thresh,4),"_fpr":fpr,"_tpr":tpr,
            "_prec_curve":prec_c,"_rec_curve":rec_c,"_cm":cm,"_y_prob":y_prob,"_y_pred_opt":y_pred_opt}

def run_model_evaluation():
    log.info("="*60); log.info("STEP 12 — MODEL EVALUATION"); log.info("="*60)
    X_val        = pd.read_csv(config.VAL_FEATURES_FILE, index_col=0)
    X_val_scaled = pd.read_csv(config.PROCESSED_DIR/"val_features_scaled.csv", index_col=0)
    y_val        = pd.read_csv(config.VAL_LABELS_FILE).set_index("gene")["label"].values
    best_name    = (config.MODELS_DIR/"best_model_name.txt").read_text().strip() if (config.MODELS_DIR/"best_model_name.txt").exists() else "random_forest"
    log.info(f"  Val: {len(y_val)} genes  pos={y_val.sum()}")
    model_files  = {
        "random_forest":(False,"model_random_forest.joblib"),
        "gradient_boosting":(False,"model_gradient_boosting.joblib"),
        "logistic_regression":(True,"model_logistic_regression.joblib"),
        "svm":(True,"model_svm.joblib"),
        "extra_trees":(False,"model_extra_trees.joblib"),
    }
    loaded = {}
    for name,(use_scaled,fname) in model_files.items():
        fpath = config.MODELS_DIR/fname
        if fpath.exists():
            loaded[name] = (joblib.load(fpath), use_scaled)
    eval_results = []
    for model_name,(model,use_scaled) in loaded.items():
        eval_results.append(evaluate_single_model(model_name,model,use_scaled,X_val,X_val_scaled,y_val))
    scalar_keys = ["model_name","auroc","auprc","brier_score","mcc","precision_opt","recall_opt","f1_opt","f1_at_05","optimal_threshold"]
    metrics_df  = pd.DataFrame([{k:r[k] for k in scalar_keys} for r in eval_results]).sort_values("auroc",ascending=False).reset_index(drop=True)
    best_metrics = next((r for r in eval_results if r["model_name"]==best_name), eval_results[0])

    # ROC plot
    fig, ax = plt.subplots(figsize=(8,7))
    ax.plot([0,1],[0,1],"k--",lw=0.8,alpha=0.5,label="Random")
    colors = plt.cm.tab10(np.linspace(0,0.8,len(eval_results)))
    for result, color in zip(eval_results, colors):
        name = result["model_name"]; is_best = name==best_name
        ax.plot(result["_fpr"],result["_tpr"],color="tomato" if is_best else color,
                lw=2.5 if is_best else 1.5,alpha=0.95 if is_best else 0.75,
                label=f"{'* ' if is_best else ''}{name} (AUROC={result['auroc']:.4f})")
    ax.set_xlabel("False Positive Rate"); ax.set_ylabel("True Positive Rate")
    ax.set_title(f"ROC Curves — Validation Set"); ax.legend(fontsize=8,loc="lower right")
    plt.tight_layout(); fig.savefig(config.FIGURES_DIR/"eval_roc_curves.png",dpi=150,bbox_inches="tight"); plt.close(fig)

    # PR plot
    fig2, ax2 = plt.subplots(figsize=(8,7))
    baseline = y_val.mean()
    ax2.axhline(baseline,color="black",linestyle="--",lw=0.8,alpha=0.5,label=f"No-skill ({baseline:.3f})")
    for result, color in zip(eval_results, colors):
        name = result["model_name"]; is_best = name==best_name
        ax2.plot(result["_rec_curve"],result["_prec_curve"],color="tomato" if is_best else color,
                 lw=2.5 if is_best else 1.5,
                 label=f"{'* ' if is_best else ''}{name} (AUPRC={result['auprc']:.4f})")
    ax2.set_xlabel("Recall"); ax2.set_ylabel("Precision"); ax2.set_title("PR Curves — Validation Set")
    ax2.legend(fontsize=8); plt.tight_layout(); fig2.savefig(config.FIGURES_DIR/"eval_pr_curves.png",dpi=150,bbox_inches="tight"); plt.close(fig2)

    metrics_df.to_csv(config.MODEL_METRICS_FILE, index=False)
    log.info("STEP 12 COMPLETE")
    return {"eval_results":eval_results,"metrics_df":metrics_df,"best_model_name":best_name,"best_metrics":best_metrics}

if __name__ == "__main__":
    r = run_model_evaluation()
    print(r["metrics_df"][["model_name","auroc","auprc","mcc","f1_opt","brier_score"]].to_string(index=False))
