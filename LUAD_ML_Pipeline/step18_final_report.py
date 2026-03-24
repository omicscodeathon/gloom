"""
step18_final_report.py - Final Outputs & Reporting
Generates: pipeline_summary_table.csv, pipeline_report.txt,
pipeline_summary_figure.png/.pdf, and run_pipeline.py master script.
"""
import logging, sys, time
from datetime import datetime
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

def safe_read_csv(path, **kwargs):
    if not path.exists(): log.warning(f"  Not found: {path}"); return None
    return pd.read_csv(path, **kwargs)

def build_pipeline_summary_table():
    ranking = safe_read_csv(config.GENE_RANKINGS_FILE, index_col=0)
    if ranking is None: raise FileNotFoundError("gene_rankings.csv not found.")
    summary = ranking.copy()
    net_feats  = safe_read_csv(config.NETWORK_FEATURES_FILE, index_col=0)
    expr_feats = safe_read_csv(config.EXPR_FEATURES_FILE, index_col=0)
    if net_feats is not None:
        for col in ["degree","weighted_degree","betweenness_centrality","eigenvector_centrality","clustering_coefficient","in_largest_component","mean_edge_weight"]:
            if col in net_feats.columns and col not in summary.columns:
                summary[col] = net_feats[col].reindex(summary.index)
    if expr_feats is not None:
        for col in ["tumor_mean","tumor_std","normal_mean","normal_std","abs_log2fc_rank","neg_log10_padj_rank"]:
            if col in expr_feats.columns and col not in summary.columns:
                summary[col] = expr_feats[col].reindex(summary.index)
    if "rank" in summary.columns: summary = summary.sort_values("rank")
    float_cols = summary.select_dtypes(include=[float]).columns
    summary[float_cols] = summary[float_cols].round(6)
    return summary

def build_pipeline_report(summary_table):
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    sep = "="*70; sep2 = "-"*70
    lines = [sep,"LUAD ML PIPELINE — FINAL REPORT",f"Generated: {now}",sep,""]
    lines += ["SECTION 1 — DATA PROVENANCE",sep2,
              f"  Tumor expression   : {config.TUMOR_EXPR_FILE.name}",
              f"  Normal expression  : {config.NORMAL_EXPR_FILE.name}",
              f"  Cancer gene list   : {config.CANCER_GENE_FILE.name}",""]
    qc = safe_read_csv(config.PROCESSED_DIR/"qc_summary.csv")
    lines += ["SECTION 2 — QC SUMMARY",sep2]
    if qc is not None:
        for _,row in qc.iterrows(): lines.append(f"  [{row['dataset']:<6}] {row['stage']:<35}: {int(row['n_genes']):>7,}")
    lines.append("")
    de_df = safe_read_csv(config.DE_RESULTS_FILE, index_col=0)
    lines += ["SECTION 3 — DIFFERENTIAL EXPRESSION",sep2]
    if de_df is not None:
        n_up = (de_df.get("direction",pd.Series())=="up").sum()
        n_down = (de_df.get("direction",pd.Series())=="down").sum()
        n_sig = de_df.get("significant",pd.Series(False)).sum()
        lines += [f"  Genes tested: {len(de_df):,}  Significant: {int(n_sig):,}  Up: {int(n_up):,}  Down: {int(n_down):,}",""]
    metrics = safe_read_csv(config.MODEL_METRICS_FILE)
    best_name = (config.MODELS_DIR/"best_model_name.txt").read_text().strip() if (config.MODELS_DIR/"best_model_name.txt").exists() else "N/A"
    lines += ["SECTION 4 — MODEL PERFORMANCE",sep2]
    if metrics is not None:
        lines.append(f"  Best model: {best_name}")
        lines.append(f"  {'Model':<25} {'AUROC':>7} {'AUPRC':>7} {'MCC':>7} {'F1':>7}")
        lines.append(f"  {'-'*55}")
        for _,row in metrics.iterrows():
            m = " *" if row["model_name"]==best_name else ""
            lines.append(f"  {row['model_name']:<25} {row['auroc']:>7.4f} {row['auprc']:>7.4f} {row['mcc']:>7.4f} {row['f1_opt']:>7.4f}{m}")
    lines.append("")
    imp = safe_read_csv(config.FEATURE_IMPORTANCE_FILE, index_col=0)
    lines += ["SECTION 5 — TOP 15 FEATURES",sep2]
    if imp is not None:
        lines.append(f"  {'Rank':<5} {'Feature':<40} {'Mean Imp':>9} {'Group'}")
        for _,row in imp.head(15).iterrows():
            lines.append(f"  {int(row['rank']):<5} {row.name:<40} {row['mean_importance']:>9.4f}  {row['feature_group']}")
    lines.append("")
    ranking = safe_read_csv(config.GENE_RANKINGS_FILE, index_col=0)
    novel   = safe_read_csv(config.RESULTS_DIR/"novel_candidates.csv", index_col=0)
    lines += ["SECTION 6 — TOP 20 RANKED GENES",sep2]
    if ranking is not None:
        lines.append(f"  {'Rank':<6} {'Gene':<15} {'Prob':>7} {'Log2FC':>8} {'Dir':>5} {'CGC':>5} {'Novel':>6}")
        for gene,row in ranking.head(20).iterrows():
            lines.append(f"  {int(row['rank']):<6} {gene:<15} {row['predicted_prob']:>7.4f} {row.get('log2fc',0):>8.3f} {row.get('direction','ns'):>5} {'Yes' if row.get('is_cgc_gene',False) else 'No':>5} {'Yes' if row.get('novel_candidate',False) else 'No':>6}")
    lines.append("")
    lines += ["SECTION 7 — NOVEL CANDIDATES (top 15)",sep2]
    if novel is not None:
        lines.append(f"  Total: {len(novel):,}")
        for gene,row in novel.head(15).iterrows():
            lines.append(f"  {int(row['rank']):<6} {gene:<15} prob={row['predicted_prob']:.4f}  log2fc={row.get('log2fc',0):.3f}  dir={row.get('direction','ns')}")
    lines += ["",sep,"END OF REPORT",sep]
    return lines

def build_summary_figure(out_path_prefix):
    fig = plt.figure(figsize=(20,24))
    gs  = gridspec.GridSpec(3,2,figure=fig,hspace=0.42,wspace=0.35)
    ax_A=fig.add_subplot(gs[0,0]); ax_B=fig.add_subplot(gs[0,1])
    ax_C=fig.add_subplot(gs[1,0]); ax_D=fig.add_subplot(gs[1,1])
    ax_E=fig.add_subplot(gs[2,0]); ax_F=fig.add_subplot(gs[2,1])
    # Panel A: QC
    qc = safe_read_csv(config.PROCESSED_DIR/"qc_summary.csv")
    if qc is not None:
        for dataset,color,offset in [("tumor","steelblue",-0.2),("normal","darkorange",0.2)]:
            sub = qc[qc["dataset"]==dataset]
            if len(sub): ax_A.bar(np.arange(len(sub))+offset, sub["n_genes"].values, 0.35, label=dataset.capitalize(), color=color, alpha=0.85)
        ax_A.set_title("(A) Gene Retention Through Pipeline",fontsize=10); ax_A.set_ylabel("Number of genes",fontsize=9); ax_A.legend(fontsize=8)
    # Panel B: Volcano
    de_df = safe_read_csv(config.DE_RESULTS_FILE, index_col=0)
    if de_df is not None and "log2fc" in de_df.columns and "neg_log10_padj" in de_df.columns:
        direction = de_df.get("direction",pd.Series("ns",index=de_df.index)).fillna("ns")
        cmap = {"up":"steelblue","down":"tomato","ns":"lightgrey"}
        for d in ["ns","down","up"]:
            mask = direction==d
            ax_B.scatter(de_df.loc[mask,"log2fc"].fillna(0), de_df.loc[mask,"neg_log10_padj"].fillna(0),
                         c=cmap[d], s=4 if d=="ns" else 8, alpha=0.4 if d=="ns" else 0.7, linewidths=0, rasterized=True)
        ax_B.axvline(config.DE_LOG2FC_THRESHOLD, color="black",linestyle="--",lw=0.8,alpha=0.5)
        ax_B.axvline(-config.DE_LOG2FC_THRESHOLD,color="black",linestyle="--",lw=0.8,alpha=0.5)
        ax_B.axhline(-np.log10(config.DE_PVALUE_THRESHOLD),color="black",linestyle=":",lw=0.8,alpha=0.5)
        ax_B.set_xlabel("Log2 Fold-Change",fontsize=9); ax_B.set_ylabel("-log10(adj P)",fontsize=9); ax_B.set_title("(B) Volcano Plot",fontsize=10)
    # Panel C: CV AUROC
    cv_res = safe_read_csv(config.MODELS_DIR/"cv_results.csv")
    best_name = (config.MODELS_DIR/"best_model_name.txt").read_text().strip() if (config.MODELS_DIR/"best_model_name.txt").exists() else ""
    if cv_res is not None:
        models_cv=cv_res["model"].tolist(); means_cv=cv_res["mean_auroc"].tolist(); stds_cv=cv_res["std_auroc"].tolist()
        colors_cv=["tomato" if m==best_name else "steelblue" for m in models_cv]
        ax_C.bar(range(len(models_cv)),means_cv,yerr=stds_cv,capsize=4,color=colors_cv,alpha=0.85,error_kw={"elinewidth":1.2})
        ax_C.set_xticks(range(len(models_cv))); ax_C.set_xticklabels(models_cv,rotation=20,ha="right",fontsize=8)
        ax_C.set_ylabel("CV AUROC",fontsize=9); ax_C.set_title("(C) Cross-Validation AUROC",fontsize=10)
    # Panel D: ROC curve
    try:
        import joblib
        model_path = config.MODELS_DIR/"best_model.joblib"
        if model_path.exists():
            model = joblib.load(model_path)
            scaled = best_name in {"logistic_regression","svm"}
            X_val_use = safe_read_csv(config.PROCESSED_DIR/"val_features_scaled.csv" if scaled else config.VAL_FEATURES_FILE, index_col=0)
            y_val = pd.read_csv(config.VAL_LABELS_FILE).set_index("gene")["label"].values
            if X_val_use is not None:
                y_prob = model.predict_proba(X_val_use.values)[:,1] if hasattr(model,"predict_proba") else model.decision_function(X_val_use.values)
                fpr,tpr,_ = roc_curve(y_val,y_prob); roc_auc_v = auc(fpr,tpr)
                ax_D.plot(fpr,tpr,color="tomato",lw=2,label=f"{best_name}\n(AUROC={roc_auc_v:.4f})")
                ax_D.plot([0,1],[0,1],"k--",lw=0.8,alpha=0.5,label="Random")
                ax_D.set_xlabel("FPR",fontsize=9); ax_D.set_ylabel("TPR",fontsize=9); ax_D.set_title("(D) ROC Curve",fontsize=10); ax_D.legend(fontsize=8)
    except Exception as e:
        ax_D.text(0.5,0.5,"ROC unavailable",ha="center",va="center",transform=ax_D.transAxes); ax_D.set_title("(D) ROC Curve",fontsize=10)
    # Panel E: Top-20 genes
    ranking = safe_read_csv(config.GENE_RANKINGS_FILE, index_col=0)
    if ranking is not None:
        top20 = ranking.head(20); genes_e = top20.index.tolist(); probs_e = top20["predicted_prob"].values
        colors_e = ["tomato" if top20.loc[g].get("is_cgc_gene",False) else "steelblue" for g in genes_e]
        ax_E.barh(range(len(genes_e)),probs_e,color=colors_e,alpha=0.85)
        ax_E.set_yticks(range(len(genes_e))); ax_E.set_yticklabels(genes_e,fontsize=8); ax_E.invert_yaxis()
        ax_E.set_xlabel("Predicted Probability",fontsize=9); ax_E.set_title("(E) Top-20 Ranked Genes",fontsize=10); ax_E.set_xlim(0,1.1)
        from matplotlib.patches import Patch
        ax_E.legend(handles=[Patch(color="tomato",alpha=0.85,label="CGC"),Patch(color="steelblue",alpha=0.85,label="Novel")],fontsize=7)
    # Panel F: Feature importance
    imp = safe_read_csv(config.FEATURE_IMPORTANCE_FILE, index_col=0)
    if imp is not None:
        top15=imp.head(15); feats_f=top15.index.tolist(); imps_f=top15["mean_importance"].values; stds_f=top15["std_importance"].values
        groups_f=top15["feature_group"].tolist(); unique_groups=list(dict.fromkeys(groups_f))
        palette=plt.cm.Set2(np.linspace(0,1,len(unique_groups))); gcolor={g:palette[i] for i,g in enumerate(unique_groups)}
        colors_f=[gcolor[g] for g in groups_f]
        ax_F.barh(range(len(feats_f)),imps_f,xerr=stds_f,capsize=3,color=colors_f,alpha=0.85,error_kw={"elinewidth":1.0})
        ax_F.set_yticks(range(len(feats_f))); ax_F.set_yticklabels(feats_f,fontsize=7); ax_F.invert_yaxis()
        ax_F.set_xlabel("Normalised Importance",fontsize=9); ax_F.set_title("(F) Top-15 Feature Importances",fontsize=10)
    fig.suptitle("LUAD ML Pipeline — Summary Figure",fontsize=13,fontweight="bold",y=1.01)
    for ext in ("png","pdf"):
        out = Path(str(out_path_prefix)+f".{ext}")
        fig.savefig(out,dpi=150,bbox_inches="tight",format=ext)
        log.info(f"  Summary figure ({ext.upper()}) saved: {out}")
    plt.close(fig)

def write_master_runner(out_path):
    script = '''"""
run_pipeline.py - Master Pipeline Runner
Executes all 18 steps of the LUAD ML Pipeline in sequence.
Usage: cd LUAD_ML_Pipeline && python run_pipeline.py
"""
import sys, time, logging
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()

logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger("run_pipeline")

STEPS = [
    ("Step 0  — Config",              "config",                          "create_output_dirs"),
    ("Step 1  — Data loading",        "step1_data_loading",              "run_data_loading"),
    ("Step 2  — Preprocessing",       "step2_preprocessing",             "run_preprocessing"),
    ("Step 3  — Harmonization",       "step3_harmonization",             "run_harmonization"),
    ("Step 4  — Differential expr",   "step4_differential_expression",   "run_differential_expression"),
    ("Step 5  — Expr features",       "step5_expression_features",       "run_expression_feature_construction"),
    ("Step 6  — Co-expr network",     "step6_coexpression_network",      "run_coexpression_network"),
    ("Step 7  — Network features",    "step7_network_features",          "run_network_feature_extraction"),
    ("Step 8  — Feature integration", "step8_feature_integration",       "run_feature_integration"),
    ("Step 9  — Label construction",  "step9_label_construction",        "run_label_construction"),
    ("Step 10 — Train/val split",     "step10_train_val_split",          "run_train_val_split"),
    ("Step 11 — Model training",      "step11_model_training",           "run_model_training"),
    ("Step 12 — Model evaluation",    "step12_model_evaluation",         "run_model_evaluation"),
    ("Step 13 — Feature importance",  "step13_feature_importance",       "run_feature_importance"),
    ("Step 14 — Gene ranking",        "step14_gene_ranking",             "run_gene_ranking"),
    ("Step 15 — Network annotation",  "step15_network_annotation",       "run_network_annotation"),
    ("Step 16 — Network export",      "step16_network_export",           "run_network_export"),
    ("Step 17 — Visualization",       "step17_interactive_visualization","run_interactive_visualization"),
    ("Step 18 — Final report",        "step18_final_report",             "run_final_report"),
]

def run_pipeline():
    log.info("="*70); log.info("LUAD ML PIPELINE — FULL RUN"); log.info("="*70)
    pipeline_start = time.time(); completed = []; failed_step = None
    for step_label, module_name, func_name in STEPS:
        log.info(f">>> {step_label}"); t0 = time.time()
        try:
            module = __import__(module_name)
            func   = getattr(module, func_name)
            func()
            elapsed = time.time() - t0
            log.info(f"<<< {step_label} DONE ({elapsed:.1f}s)")
            completed.append((step_label, elapsed, "OK"))
        except Exception as e:
            elapsed = time.time() - t0
            log.error(f"<<< {step_label} FAILED: {type(e).__name__}: {e}")
            completed.append((step_label, elapsed, f"FAILED: {e}"))
            failed_step = step_label; break
    total = time.time() - pipeline_start
    log.info("="*70); log.info("EXECUTION SUMMARY"); log.info("="*70)
    for label, elapsed, status in completed:
        icon = "OK" if status=="OK" else "FAIL"
        log.info(f"  [{icon}] {label:<45} {elapsed:>7.1f}s")
    log.info(f"  Total: {total:.1f}s ({total/60:.1f} min)")
    if failed_step: log.error(f"  STOPPED at: {failed_step}"); sys.exit(1)
    else: log.info("  Pipeline COMPLETE.")
    log.info("="*70)

if __name__ == "__main__":
    run_pipeline()
'''
    out_path.write_text(script, encoding="utf-8")
    log.info(f"  Master runner saved: {out_path}")

def run_final_report():
    log.info("="*60); log.info("STEP 18 — FINAL OUTPUTS & REPORTING"); log.info("="*60)
    reports_dir = config.REPORTS_DIR; reports_dir.mkdir(parents=True, exist_ok=True)
    log.info("[1] Building summary table …")
    try:
        summary_table = build_pipeline_summary_table()
        summary_table.to_csv(reports_dir/"pipeline_summary_table.csv")
        log.info(f"  Summary table: {summary_table.shape}")
    except Exception as e:
        log.error(f"  Summary table failed: {e}"); summary_table = pd.DataFrame()
    log.info("[2] Building text report …")
    try:
        report_lines = build_pipeline_report(summary_table)
        (reports_dir/"pipeline_report.txt").write_text("\n".join(report_lines), encoding="utf-8")
        log.info(f"  Report: {len(report_lines)} lines")
    except Exception as e:
        log.error(f"  Report failed: {e}"); report_lines = []
    log.info("[3] Building summary figure …")
    try: build_summary_figure(reports_dir/"pipeline_summary_figure")
    except Exception as e: log.error(f"  Figure failed: {e}")
    log.info("[4] Writing master runner …")
    runner_path = Path(__file__).resolve().parent/"run_pipeline.py"
    write_master_runner(runner_path)
    log.info("="*60); log.info("STEP 18 SUMMARY"); log.info("="*60)
    log.info(f"  Summary table: {reports_dir}/pipeline_summary_table.csv")
    log.info(f"  Text report  : {reports_dir}/pipeline_report.txt")
    log.info(f"  Summary fig  : {reports_dir}/pipeline_summary_figure.png")
    log.info(f"  Master runner: {runner_path}")
    log.info("")
    log.info("  Key deliverables:")
    log.info(f"    Gene rankings    : {config.GENE_RANKINGS_FILE}")
    log.info(f"    Novel candidates : {config.RESULTS_DIR}/novel_candidates.csv")
    log.info(f"    Dashboard        : {config.FIGURES_DIR}/interactive_dashboard.html")
    log.info("STEP 18 COMPLETE")
    return {"summary_table":summary_table,"report_lines":report_lines}

if __name__ == "__main__":
    r = run_final_report()
    if not r["summary_table"].empty:
        cols = [c for c in ["rank","predicted_prob","is_cgc_gene","novel_candidate","log2fc","direction"] if c in r["summary_table"].columns]
        print(r["summary_table"].head(10)[cols].round(4).to_string())
