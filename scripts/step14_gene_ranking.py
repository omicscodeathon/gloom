"""
step14_gene_ranking.py - Gene Prediction & Ranking
Scores ALL genes using best model. Ranks by predicted probability.
Flags novel candidates (high-prob + not in LCGene known-positive set).
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
    use_scaled = best_name in {"logistic_regression"}
    opt_thresh = 0.50
    if config.MODEL_METRICS_FILE.exists():
        m = pd.read_csv(config.MODEL_METRICS_FILE)
        match = m[m["model_name"]==best_name]
        if not match.empty: opt_thresh = float(match["optimal_threshold"].iloc[0])
    X = features_scaled if use_scaled else features
    log.info(f"  Scoring {X.shape[0]:,} genes …")
    probs = model.predict_proba(X)[:,1] if hasattr(model,"predict_proba") else model.decision_function(X)
    probs_s = pd.Series(probs, index=X.index, name="predicted_prob")

    # Build ranking
    ranking = pd.DataFrame(index=probs_s.index)
    ranking["predicted_prob"]  = probs_s
    ranking["rank"]            = probs_s.rank(ascending=False, method="min").astype(int)
    ranking["percentile"]      = (probs_s.rank(ascending=True, pct=True)*100).round(2)
    ranking["predicted_label"] = (probs_s >= opt_thresh).astype(int)
    ranking["label"]           = labels.reindex(probs_s.index).fillna(0).astype(int)
    for col in ["is_lcgene_gene","log2fc","pvalue_adj","neg_log10_padj","direction","is_de_significant","mean_tumor","mean_normal"]:
        if col in annotation.columns:
            ranking[col] = annotation[col].reindex(probs_s.index)
    ranking["is_lcgene_gene"]    = ranking.get("is_lcgene_gene", pd.Series(False,index=probs_s.index)).fillna(False).astype(bool)
    ranking["is_de_significant"] = ranking.get("is_de_significant", pd.Series(False,index=probs_s.index)).fillna(False).astype(bool)
    ranking["log2fc"]     = ranking["log2fc"].fillna(0.0)
    ranking["direction"]  = ranking["direction"].fillna("ns")
    ranking["pvalue_adj"] = ranking["pvalue_adj"].fillna(1.0)
    ranking["novel_candidate"] = (ranking["predicted_prob"] >= NOVEL_PROB_THRESHOLD) & (~ranking["is_lcgene_gene"])
    ranking = ranking.sort_values("predicted_prob", ascending=False)

    # PU framing - choose labels for logging/plots (Step 6)
    _pu = getattr(config, "PU_FRAMING", False)
    pos_label_str = "positive (known)" if _pu else "LCGene"
    neg_label_str = "unlabeled"        if _pu else "non-LCGene"

    # Summary + ranking metrics (Step 3)
    n_novel = ranking["novel_candidate"].sum()
    lcgene_ranks     = ranking.loc[ranking["is_lcgene_gene"],"rank"]
    median_lcgene_rank = lcgene_ranks.median() if len(lcgene_ranks)>0 else np.nan
    total_pos        = ranking["is_lcgene_gene"].sum()

    def _recall_at(k):
        return ranking.head(k)["is_lcgene_gene"].sum() / total_pos if total_pos > 0 else 0.0
    def _prec_at(k):
        return ranking.head(k)["is_lcgene_gene"].mean()

    lcgene_top50  = ranking.head(50)["is_lcgene_gene"].sum()
    lcgene_top100 = ranking.head(100)["is_lcgene_gene"].sum()
    lcgene_top200 = ranking.head(200)["is_lcgene_gene"].sum()
    lcgene_top500 = ranking.head(500)["is_lcgene_gene"].sum()
    recall_50  = _recall_at(50);  recall_100 = _recall_at(100)
    recall_200 = _recall_at(200); recall_500 = _recall_at(500)
    prec_50    = _prec_at(50);    prec_100   = _prec_at(100)

    log.info(f"  Total: {len(ranking):,}  Novel candidates: {n_novel:,}")
    log.info(f"  Median {pos_label_str} rank: {median_lcgene_rank:.0f}")
    log.info(f"  {pos_label_str.capitalize()} in top-50:  {lcgene_top50}  "
             f"top-100: {lcgene_top100}  top-200: {lcgene_top200}  top-500: {lcgene_top500}")
    log.info(f"  Recall@50={recall_50:.4f}  Recall@100={recall_100:.4f}  "
             f"Recall@200={recall_200:.4f}  Recall@500={recall_500:.4f}")
    log.info(f"  Precision@50={prec_50:.4f}  Precision@100={prec_100:.4f}")

    # Score distribution plot (Step 3 / Step 6 - PU-aware labels)
    lcgene_probs    = ranking.loc[ranking["is_lcgene_gene"],"predicted_prob"]
    nonlcgene_probs = ranking.loc[~ranking["is_lcgene_gene"],"predicted_prob"]
    fig, axes = plt.subplots(1,2, figsize=(13,5))
    bins = np.linspace(0,1,50)
    axes[0].hist(nonlcgene_probs, bins=bins, color="steelblue", alpha=0.6,
                 label=f"{neg_label_str.capitalize()} (n={len(nonlcgene_probs):,})", density=True)
    axes[0].hist(lcgene_probs,   bins=bins, color="tomato",    alpha=0.75,
                 label=f"{pos_label_str.capitalize()} (n={len(lcgene_probs):,})", density=True)
    axes[0].axvline(NOVEL_PROB_THRESHOLD, color="black", linestyle="--", lw=1.2, label="Novel threshold")
    axes[0].set_xlabel("Predicted Probability"); axes[0].set_ylabel("Density")
    axes[0].set_title(f"Score Distribution: {pos_label_str.capitalize()} vs {neg_label_str.capitalize()}")
    axes[0].legend(fontsize=8)
    sorted_r = ranking.sort_values("rank")
    cum_lcgene = sorted_r["is_lcgene_gene"].cumsum(); total_lcgene = sorted_r["is_lcgene_gene"].sum()
    frac_captured = cum_lcgene/total_lcgene; frac_genes = np.arange(1,len(sorted_r)+1)/len(sorted_r)
    axes[1].plot(frac_genes, frac_captured, color="tomato", lw=2, label="Model ranking")
    axes[1].plot([0,1],[0,1],"k--",lw=0.8,alpha=0.5,label="Random")
    axes[1].fill_between(frac_genes,frac_captured,frac_genes,alpha=0.15,color="tomato")
    axes[1].set_xlabel("Fraction of genes screened")
    axes[1].set_ylabel(f"Fraction of {pos_label_str} captured")
    axes[1].set_title(f"{pos_label_str.capitalize()} Gene Enrichment Curve")
    axes[1].legend(fontsize=8)
    plt.tight_layout(); fig.savefig(config.FIGURES_DIR/"ranking_score_distribution.png",dpi=150,bbox_inches="tight"); plt.close(fig)

    novel_candidates = ranking[ranking["novel_candidate"]].copy()
    ranking.to_csv(config.GENE_RANKINGS_FILE)
    novel_candidates.to_csv(config.RESULTS_DIR/"novel_candidates.csv")

    # Save ranking-level metrics for downstream steps (Step 3)
    ranking_metrics = pd.DataFrame([{
        "median_lcgene_rank": round(median_lcgene_rank, 1) if not np.isnan(median_lcgene_rank) else np.nan,
        "lcgene_top50": int(lcgene_top50), "lcgene_top100": int(lcgene_top100),
        "lcgene_top200": int(lcgene_top200), "lcgene_top500": int(lcgene_top500),
        "recall_at_50": round(recall_50, 4), "recall_at_100": round(recall_100, 4),
        "recall_at_200": round(recall_200, 4), "recall_at_500": round(recall_500, 4),
        "precision_at_50": round(prec_50, 4), "precision_at_100": round(prec_100, 4),
        "n_novel": int(n_novel), "pu_framing": _pu,
    }])
    ranking_metrics.to_csv(config.RESULTS_DIR/"ranking_metrics.csv", index=False)

    log.info("STEP 14 COMPLETE")
    return {"ranking":ranking,"novel_candidates":novel_candidates,"ranking_metrics":ranking_metrics}

if __name__ == "__main__":
    r = run_gene_ranking()
    print("\n--- Top 20 Ranked Genes ---")
    cols = ["rank","predicted_prob","log2fc","direction","is_lcgene_gene","novel_candidate"]
    print(r["ranking"].head(20)[[c for c in cols if c in r["ranking"].columns]].to_string())
