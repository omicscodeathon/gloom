"""
step14_gene_ranking.py - Gene Prediction & Ranking
Scores ALL genes using best model. Ranks by predicted probability.
Flags novel candidates (high-prob + not in CGC).
Outputs: gene_rankings.csv, novel_candidates.csv + ranking plots.
"""
import logging, sys
from pathlib import Path
import joblib
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)
NOVEL_PROB_THRESHOLD = 0.50

def run_gene_ranking():
    log.info("="*60); log.info("STEP 14 - GENE PREDICTION & RANKING"); log.info("="*60)
    features        = pd.read_csv(config.INTEGRATED_FEATURES_FILE, index_col=0)
    features_scaled = pd.read_csv(config.PROCESSED_DIR/"integrated_features_scaled.csv", index_col=0)
    labels          = pd.read_csv(config.LABELS_FILE).set_index("gene")["label"]
    annotation      = pd.read_csv(config.PROCESSED_DIR/"gene_annotation_table.csv", index_col=0)
    model_path      = config.MODELS_DIR/"best_model.joblib"
    if not model_path.exists():
        raise FileNotFoundError(f"Best model not found: {model_path}")
    model     = joblib.load(model_path)
    best_name = (config.MODELS_DIR/"best_model_name.txt").read_text().strip() if (config.MODELS_DIR/"best_model_name.txt").exists() else "random_forest"
    use_scaled = best_name in {"logistic_regression","svm"}
    opt_thresh = 0.50
    if config.MODEL_METRICS_FILE.exists():
        m = pd.read_csv(config.MODEL_METRICS_FILE)
        match = m[m["model_name"]==best_name]
        if not match.empty: opt_thresh = float(match["optimal_threshold"].iloc[0])
    X = features_scaled if use_scaled else features
    log.info(f"  Scoring {X.shape[0]:,} genes …")
    probs = model.predict_proba(X.values)[:,1] if hasattr(model,"predict_proba") else model.decision_function(X.values)
    probs_s = pd.Series(probs, index=X.index, name="predicted_prob")

    # Build ranking
    ranking = pd.DataFrame(index=probs_s.index)
    ranking["predicted_prob"]  = probs_s
    ranking["rank"]            = probs_s.rank(ascending=False, method="min").astype(int)
    ranking["percentile"]      = (probs_s.rank(ascending=True, pct=True)*100).round(2)
    ranking["predicted_label"] = (probs_s >= opt_thresh).astype(int)
    ranking["label"]           = labels.reindex(probs_s.index).fillna(0).astype(int)
    for col in ["is_cgc_gene","log2fc","pvalue_adj","neg_log10_padj","direction","is_de_significant","mean_tumor","mean_normal"]:
        if col in annotation.columns:
            ranking[col] = annotation[col].reindex(probs_s.index)
    ranking["is_cgc_gene"]       = ranking.get("is_cgc_gene", pd.Series(False,index=probs_s.index)).fillna(False).astype(bool)
    ranking["is_de_significant"] = ranking.get("is_de_significant", pd.Series(False,index=probs_s.index)).fillna(False).astype(bool)
    ranking["log2fc"]     = ranking["log2fc"].fillna(0.0)
    ranking["direction"]  = ranking["direction"].fillna("ns")
    ranking["pvalue_adj"] = ranking["pvalue_adj"].fillna(1.0)
    ranking["novel_candidate"] = (ranking["predicted_prob"] >= NOVEL_PROB_THRESHOLD) & (~ranking["is_cgc_gene"])
    ranking = ranking.sort_values("predicted_prob", ascending=False)

    # Summary
    n_novel = ranking["novel_candidate"].sum()
    cgc_ranks = ranking.loc[ranking["is_cgc_gene"],"rank"]
    median_cgc_rank = cgc_ranks.median() if len(cgc_ranks)>0 else np.nan
    cgc_top100 = ranking.head(100)["is_cgc_gene"].sum()
    log.info(f"  Total: {len(ranking):,}  Novel: {n_novel:,}")
    log.info(f"  Median CGC rank: {median_cgc_rank:.0f}  CGC in top-100: {cgc_top100}")

    # Score distribution plot
    cgc_probs    = ranking.loc[ranking["is_cgc_gene"],"predicted_prob"]
    noncgc_probs = ranking.loc[~ranking["is_cgc_gene"],"predicted_prob"]
    fig, axes = plt.subplots(1,2, figsize=(13,5))
    bins = np.linspace(0,1,50)
    axes[0].hist(noncgc_probs, bins=bins, color="steelblue", alpha=0.6, label=f"Non-CGC (n={len(noncgc_probs):,})", density=True)
    axes[0].hist(cgc_probs,    bins=bins, color="tomato",    alpha=0.75,label=f"CGC (n={len(cgc_probs):,})",         density=True)
    axes[0].axvline(NOVEL_PROB_THRESHOLD, color="black", linestyle="--", lw=1.2, label=f"Novel threshold")
    axes[0].set_xlabel("Predicted Probability"); axes[0].set_ylabel("Density")
    axes[0].set_title("Score Distribution: CGC vs Non-CGC"); axes[0].legend(fontsize=8)
    sorted_r = ranking.sort_values("rank")
    cum_cgc = sorted_r["is_cgc_gene"].cumsum(); total_cgc = sorted_r["is_cgc_gene"].sum()
    frac_captured = cum_cgc/total_cgc; frac_genes = np.arange(1,len(sorted_r)+1)/len(sorted_r)
    axes[1].plot(frac_genes, frac_captured, color="tomato", lw=2, label="Model ranking")
    axes[1].plot([0,1],[0,1],"k--",lw=0.8,alpha=0.5,label="Random")
    axes[1].fill_between(frac_genes,frac_captured,frac_genes,alpha=0.15,color="tomato")
    axes[1].set_xlabel("Fraction of genes screened"); axes[1].set_ylabel("Fraction of CGC captured")
    axes[1].set_title("CGC Gene Enrichment Curve"); axes[1].legend(fontsize=8)
    plt.tight_layout(); fig.savefig(config.FIGURES_DIR/"ranking_score_distribution.png",dpi=150,bbox_inches="tight"); plt.close(fig)

    novel_candidates = ranking[ranking["novel_candidate"]].copy()
    ranking.to_csv(config.GENE_RANKINGS_FILE)
    novel_candidates.to_csv(config.RESULTS_DIR/"novel_candidates.csv")
    log.info("STEP 14 COMPLETE")
    return {"ranking":ranking,"novel_candidates":novel_candidates}

if __name__ == "__main__":
    r = run_gene_ranking()
    print("\n--- Top 20 Ranked Genes ---")
    cols = ["rank","predicted_prob","log2fc","direction","is_cgc_gene","novel_candidate"]
    print(r["ranking"].head(20)[[c for c in cols if c in r["ranking"].columns]].to_string())
