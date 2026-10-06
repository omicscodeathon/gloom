"""
step14_gene_ranking.py - Gene Prediction & Ranking
Ranks ALL genes of the analysis universe and flags non-LCGene candidates
(high score + not in the LCGene reference set).

PRIMARY ranking (v0.2.0): out-of-fold (cross-fitted) PU-bagging scores of step11c
(config.USE_CROSSFIT = True, results/oof_scores.csv).  Every gene is scored by models that
never saw it during training, so the top of the ranking and all top-K statistics are
independent of the training labels.  The score of the model fitted on ALL data (step11b, which
has seen the training positives) is kept in a separate column ``full_model_score``;
``is_training_positive`` marks the LCGene positives that were used to fit that full model.
If cross-fitting is disabled or missing the ranking falls back to the full-data model and the
reported top-K metrics are restricted to genes the model did not train on.

"Candidate" terminology: absence from LCGene does not establish biological novelty, so genes
are reported as ``non_lcgene_candidate`` (see step20 for the external-evidence check).

Outputs: gene_rankings.csv, non_lcgene_candidates.csv (+ _sensitivity), ranking_metrics.csv
         + ranking plots.
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
    handlers=[logging.FileHandler(config.LOG_FILE, encoding="utf-8"), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)
CANDIDATE_PROB_THRESHOLD      = getattr(config, "CANDIDATE_PROB_THRESHOLD",
                                        getattr(config, "NOVEL_PROB_THRESHOLD", 0.70))
CANDIDATE_PROB_THRESHOLD_SENS = getattr(config, "CANDIDATE_PROB_THRESHOLD_SENS",
                                        getattr(config, "NOVEL_PROB_THRESHOLD_SENS", 0.50))
CANDIDATE_MIN_ABS_LOG2FC      = getattr(config, "CANDIDATE_MIN_ABS_LOG2FC", 1.0)
TOP_KS = (50, 100, 200, 500)


def _load_full_model_scores(features):
    """
    Scores of the model fitted on ALL training data (step11b PU bagging; fallback: best_model).
    Returns (Series or None, source name).  These scores are resubstitution scores for the
    training positives — they must not be used for independent evaluation.
    """
    pu_scores_path = config.RESULTS_DIR / "pu_bagging_scores.csv"
    if pu_scores_path.exists():
        pu_df = pd.read_csv(pu_scores_path, index_col=0)
        return pu_df["pu_score"], "pu_bagging_full_model"
    model_path = config.MODELS_DIR / "best_model.joblib"
    if not model_path.exists():
        return None, ""
    log.warning("  PU bagging scores not found — falling back to best_model.joblib.")
    model     = joblib.load(model_path)
    best_name = (config.MODELS_DIR/"best_model_name.txt").read_text().strip() \
                if (config.MODELS_DIR/"best_model_name.txt").exists() else "random_forest"
    train_feats_path = config.TRAIN_FEATURES_FILE
    if train_feats_path.exists():
        selected_cols = pd.read_csv(train_feats_path, index_col=0).columns.tolist()
        selected_cols = [c for c in selected_cols if c in features.columns]
        X = features[selected_cols]
    else:
        X = features
    if best_name == "logistic_regression":
        scaler_path = config.MODELS_DIR / "robust_scaler.joblib"
        if scaler_path.exists():
            scaler = joblib.load(scaler_path)
            X = pd.DataFrame(scaler.transform(X.values), index=X.index, columns=X.columns)
    probs = model.predict_proba(X)[:, 1] if hasattr(model, "predict_proba") \
            else model.decision_function(X)
    return pd.Series(probs, index=X.index, name="predicted_prob"), best_name


def run_gene_ranking():
    log.info("="*60); log.info("STEP 14 — GENE PREDICTION & RANKING"); log.info("="*60)
    features   = pd.read_csv(config.INTEGRATED_FEATURES_FILE, index_col=0)
    labels     = pd.read_csv(config.LABELS_FILE).set_index("gene")["label"]
    annotation = pd.read_csv(config.PROCESSED_DIR/"gene_annotation_table.csv", index_col=0)
    # DE results: fallback annotation for genes that are not in the annotation table.
    de_df = pd.read_csv(config.DE_RESULTS_FILE, index_col=0) \
            if config.DE_RESULTS_FILE.exists() else pd.DataFrame()
    universe = features.index.intersection(labels.index)

    # ── Score sources ────────────────────────────────────────────────────────────────────────
    full_scores, full_source = _load_full_model_scores(features)
    oof_path = config.RESULTS_DIR / "oof_scores.csv"
    oof_scores = None
    if getattr(config, "USE_CROSSFIT", True):
        if oof_path.exists():
            oof_scores = pd.read_csv(oof_path, index_col=0)["oof_score"]
        else:
            log.warning("  USE_CROSSFIT = True but results/oof_scores.csv is missing — run step11c. "
                        "Falling back to the full-data model: its top-K includes training positives "
                        "and is NOT an independent evaluation.")

    if oof_scores is not None:
        probs_s      = oof_scores.reindex(universe)
        score_source = "crossfit_oof_pu_bagging"
        metric_basis = "out_of_fold"
        log.info(f"  PRIMARY ranking: out-of-fold (cross-fitted) scores — {probs_s.notna().sum():,} genes.")
    elif full_scores is not None:
        probs_s      = full_scores.reindex(universe)
        score_source = full_source
        metric_basis = "held_out_genes_only"
        log.warning("  Ranking with the FULL-DATA model: top-K metrics below are computed on genes "
                    "outside the training split only.")
    else:
        raise FileNotFoundError("No scoring source found. Run step11c (cross-fitting), step11b and/or step11 first.")
    probs_s = probs_s.fillna(0.0)
    probs_s.name = "predicted_prob"
    log.info(f"  Score source: {score_source}  |  candidate threshold: {CANDIDATE_PROB_THRESHOLD}")

    # ── Training split (to flag positives the full model was fitted on) ───────────────────────
    if config.TRAIN_FEATURES_FILE.exists():
        train_idx = pd.read_csv(config.TRAIN_FEATURES_FILE, index_col=0).index
    else:
        log.warning("  Training split file not found — is_training_positive cannot be determined.")
        train_idx = pd.Index([])

    # Build ranking
    ranking = pd.DataFrame(index=probs_s.index)
    ranking["predicted_prob"]  = probs_s
    ranking["rank"]            = probs_s.rank(ascending=False, method="min").astype(int)
    ranking["percentile"]      = (probs_s.rank(ascending=True, pct=True)*100).round(2)
    LABEL_THRESHOLD = 0.40  # for PU bagging scores: fraction of base classifiers agreeing
    ranking["predicted_label"] = (probs_s >= LABEL_THRESHOLD).astype(int)
    ranking["label"]           = labels.reindex(probs_s.index).fillna(0).astype(int)
    ranking["score_source"]    = score_source
    if oof_scores is not None:
        ranking["oof_score"]   = oof_scores.reindex(ranking.index)
    ranking["full_model_score"] = (full_scores.reindex(ranking.index)
                                   if full_scores is not None else probs_s)
    # LCGene positives that were used to fit the full-data model (resubstitution scores)
    ranking["is_training_positive"] = (ranking["label"] == 1) & ranking.index.isin(train_idx)
    for col in ["is_lcgene_gene","log2fc","pvalue_adj","neg_log10_padj","direction","is_de_significant","mean_tumor","mean_normal"]:
        if col in annotation.columns:
            ranking[col] = annotation[col].reindex(probs_s.index)
    ranking["is_lcgene_gene"]    = ranking.get("is_lcgene_gene", pd.Series(False,index=probs_s.index)).fillna(False).astype(bool)
    ranking["is_de_significant"] = ranking.get("is_de_significant", pd.Series(False,index=probs_s.index)).fillna(False).astype(bool)
    # Fill missing DE columns from DE results file
    if not de_df.empty:
        de_fallback = {"log2fc": "log2fc", "pvalue_adj": "pvalue_adj",
                       "neg_log10_padj": "neg_log10_padj", "direction": "direction",
                       "mean_tumor": "mean_tumor", "mean_normal": "mean_normal"}
        for rcol, dcol in de_fallback.items():
            if dcol in de_df.columns:
                if rcol not in ranking.columns:
                    ranking[rcol] = np.nan
                missing = ranking[rcol].isna()
                if missing.any():
                    ranking.loc[missing, rcol] = de_df[dcol].reindex(ranking.index[missing]).values
        if "significant" in de_df.columns:
            not_in_annot = ~ranking.index.isin(annotation.index)
            if not_in_annot.any():
                ranking.loc[not_in_annot, "is_de_significant"] = \
                    de_df["significant"].reindex(ranking.index[not_in_annot]).fillna(False).values
    ranking["log2fc"]     = ranking["log2fc"].fillna(0.0)
    ranking["direction"]  = ranking["direction"].fillna("ns")
    ranking["pvalue_adj"] = ranking["pvalue_adj"].fillna(1.0)
    # Mark query genes (user's candidate genes with unknown LUAD relationship)
    query_genes_file = getattr(config, "QUERY_GENES_FILE", None)
    if query_genes_file and Path(query_genes_file).exists():
        q_df = pd.read_csv(Path(query_genes_file), sep="\t")
        query_set = set(q_df["GeneSymbol"].dropna().astype(str).str.strip().str.upper())
        ranking["is_query_gene"] = ranking.index.str.strip().str.upper().isin(query_set)
        n_found = ranking["is_query_gene"].sum()
        log.info(f"  Query genes found in ranking: {n_found} / {len(query_set)}")
        if len(query_set) - n_found > 0:
            missing_q = query_set - set(ranking.index.str.strip().str.upper())
            log.warning(f"  Query genes not found in transcriptome: {sorted(missing_q)}")
    else:
        ranking["is_query_gene"] = False

    # Non-LCGene candidates = genes absent from the LCGene reference set with a high score AND a
    # minimum DE signal (|log2FC| >= CANDIDATE_MIN_ABS_LOG2FC; reporting filter only — it does
    # not influence training or the ranking).  Absence from LCGene does NOT imply novelty:
    # see step20 for the external-evidence classification.
    # Primary list: prob >= CANDIDATE_PROB_THRESHOLD — high-confidence.
    # Sensitivity list: prob >= CANDIDATE_PROB_THRESHOLD_SENS — broader.
    _de_filter = ranking["log2fc"].abs() >= CANDIDATE_MIN_ABS_LOG2FC
    if ranking["is_query_gene"].any():
        scope = ranking["is_query_gene"] & (~ranking["is_lcgene_gene"])
    else:
        scope = ~ranking["is_lcgene_gene"]
    ranking["non_lcgene_candidate"]      = scope & (ranking["predicted_prob"] >= CANDIDATE_PROB_THRESHOLD) & _de_filter
    ranking["non_lcgene_candidate_sens"] = scope & (ranking["predicted_prob"] >= CANDIDATE_PROB_THRESHOLD_SENS) & _de_filter
    ranking = ranking.sort_values("predicted_prob", ascending=False)

    # PU framing — choose labels for logging/plots
    _pu = getattr(config, "PU_FRAMING", False)
    pos_label_str = "positive (known)" if _pu else "LCGene"
    neg_label_str = "unlabeled"        if _pu else "non-LCGene"

    # ── Evaluation view: genes whose score is independent of their own label ──────────────────
    # out-of-fold ranking -> every gene; full-model ranking -> genes outside the training split.
    if metric_basis == "out_of_fold":
        eval_view = ranking
    else:
        eval_view = ranking[~ranking.index.isin(train_idx)]
    eval_view = eval_view.copy()
    eval_view["rank_in_eval"] = eval_view["predicted_prob"].rank(ascending=False, method="min").astype(int)
    total_pos = int(eval_view["is_lcgene_gene"].sum())
    n_novel      = int(ranking["non_lcgene_candidate"].sum())
    n_novel_sens = int(ranking["non_lcgene_candidate_sens"].sum())
    lcgene_ranks = eval_view.loc[eval_view["is_lcgene_gene"], "rank_in_eval"]
    median_lcgene_rank = lcgene_ranks.median() if len(lcgene_ranks) > 0 else np.nan

    metrics_row = {
        "score_source": score_source,
        "metric_basis": metric_basis,
        "n_eval_genes": len(eval_view),
        "n_eval_positives": total_pos,
        "median_lcgene_rank": round(median_lcgene_rank, 1) if not np.isnan(median_lcgene_rank) else np.nan,
    }
    for k in TOP_KS:
        top_eval = eval_view.head(k)
        hits = int(top_eval["is_lcgene_gene"].sum())
        # independent counts (out-of-fold, or held-out genes only)
        metrics_row[f"lcgene_top{k}"]      = hits
        metrics_row[f"recall_at_{k}"]      = round(hits / total_pos, 4) if total_pos else 0.0
        metrics_row[f"precision_at_{k}"]   = round(float(top_eval["is_lcgene_gene"].mean()), 4) if len(top_eval) else np.nan
        if metric_basis == "out_of_fold":
            metrics_row[f"lcgene_top{k}_training"] = np.nan   # not applicable: no gene is scored by a model that saw it
        else:
            # resubstitution counts for transparency (NOT independent)
            metrics_row[f"lcgene_top{k}_training"] = int(ranking.head(k)["is_training_positive"].sum())
    metrics_row["n_non_lcgene_candidates"] = n_novel
    metrics_row["pu_framing"] = _pu

    log.info(f"  Total: {len(ranking):,}  Non-LCGene candidates (prob>={CANDIDATE_PROB_THRESHOLD}): {n_novel:,}  "
             f"Sensitivity list (prob>={CANDIDATE_PROB_THRESHOLD_SENS}): {n_novel_sens:,}")
    log.info(f"  Metric basis: {metric_basis}  ({len(eval_view):,} genes, {total_pos} {pos_label_str})")
    log.info(f"  Median {pos_label_str} rank (evaluation view): {median_lcgene_rank:.0f}")
    log.info("  " + "  ".join(f"{pos_label_str.capitalize()} in top-{k}: {metrics_row[f'lcgene_top{k}']}" for k in TOP_KS))
    log.info("  " + "  ".join(f"Recall@{k}={metrics_row[f'recall_at_{k}']:.4f}" for k in TOP_KS))
    log.info("  " + "  ".join(f"Precision@{k}={metrics_row[f'precision_at_{k}']:.4f}" for k in (50, 100)))
    if metric_basis != "out_of_fold":
        log.warning("  Training positives are NOT counted above; "
                    f"in the all-gene top-100, {metrics_row['lcgene_top100_training']} are training positives "
                    "(resubstitution).")

    # Score distribution plot (evaluation view only)
    lcgene_probs    = eval_view.loc[eval_view["is_lcgene_gene"],"predicted_prob"]
    nonlcgene_probs = eval_view.loc[~eval_view["is_lcgene_gene"],"predicted_prob"]
    fig, axes = plt.subplots(1,2, figsize=(13,5))
    bins = np.linspace(0,1,50)
    basis_txt = "out-of-fold" if metric_basis == "out_of_fold" else "held-out genes"
    axes[0].hist(nonlcgene_probs, bins=bins, color="steelblue", alpha=0.6,
                 label=f"{neg_label_str.capitalize()} (n={len(nonlcgene_probs):,})", density=True)
    axes[0].hist(lcgene_probs,   bins=bins, color="tomato",    alpha=0.75,
                 label=f"{pos_label_str.capitalize()} (n={len(lcgene_probs):,})", density=True)
    axes[0].axvline(CANDIDATE_PROB_THRESHOLD, color="black", linestyle="--", lw=1.2, label="Candidate threshold")
    axes[0].set_xlabel("Predicted Probability"); axes[0].set_ylabel("Density")
    axes[0].set_title(f"Score Distribution ({basis_txt}): {pos_label_str.capitalize()} vs {neg_label_str.capitalize()}")
    axes[0].legend(fontsize=8)
    sorted_r = eval_view.sort_values("rank_in_eval")
    cum_lcgene = sorted_r["is_lcgene_gene"].cumsum(); total_lcgene = sorted_r["is_lcgene_gene"].sum()
    frac_captured = cum_lcgene/total_lcgene if total_lcgene > 0 else cum_lcgene*0
    frac_genes = np.arange(1,len(sorted_r)+1)/max(len(sorted_r), 1)
    axes[1].plot(frac_genes, frac_captured, color="tomato", lw=2, label=f"Model ranking ({basis_txt})")
    axes[1].plot([0,1],[0,1],"k--",lw=0.8,alpha=0.5,label="Random")
    axes[1].fill_between(frac_genes,frac_captured,frac_genes,alpha=0.15,color="tomato")
    axes[1].set_xlabel("Fraction of genes screened")
    axes[1].set_ylabel(f"Fraction of {pos_label_str} captured")
    axes[1].set_title(f"{pos_label_str.capitalize()} Gene Enrichment Curve")
    axes[1].legend(fontsize=8)
    plt.tight_layout(); fig.savefig(config.FIGURES_DIR/"ranking_score_distribution.png",dpi=150,bbox_inches="tight"); plt.close(fig)

    non_lcgene_candidates      = ranking[ranking["non_lcgene_candidate"]].copy()
    non_lcgene_candidates_sens = ranking[ranking["non_lcgene_candidate_sens"]].copy()
    ranking.to_csv(config.GENE_RANKINGS_FILE)
    non_lcgene_candidates.to_csv(config.RESULTS_DIR/"non_lcgene_candidates.csv")
    non_lcgene_candidates_sens.to_csv(config.RESULTS_DIR/"non_lcgene_candidates_sensitivity.csv")
    log.info(f"  High-confidence non-LCGene candidates saved: {len(non_lcgene_candidates)} genes (prob>={CANDIDATE_PROB_THRESHOLD})")
    log.info(f"  Sensitivity non-LCGene candidates saved:     {len(non_lcgene_candidates_sens)} genes (prob>={CANDIDATE_PROB_THRESHOLD_SENS})")

    # Save query-gene-only ranking (the user's candidates scored by the model)
    query_ranking = ranking[ranking["is_query_gene"]].copy() if ranking["is_query_gene"].any() else pd.DataFrame()
    query_ranking.to_csv(config.RESULTS_DIR/"query_gene_rankings.csv")
    log.info(f"  Query gene rankings saved: {len(query_ranking)} genes")

    # Save ranking-level metrics for downstream steps
    ranking_metrics = pd.DataFrame([metrics_row])
    ranking_metrics.to_csv(config.RESULTS_DIR/"ranking_metrics.csv", index=False)

    log.info("STEP 14 COMPLETE")
    return {"ranking":ranking,"non_lcgene_candidates":non_lcgene_candidates,"ranking_metrics":ranking_metrics}

if __name__ == "__main__":
    r = run_gene_ranking()
    print("\n--- Top 20 Ranked Genes ---")
    cols = ["rank","predicted_prob","full_model_score","log2fc","direction","is_lcgene_gene",
            "is_training_positive","non_lcgene_candidate"]
    print(r["ranking"].head(20)[[c for c in cols if c in r["ranking"].columns]].to_string())
