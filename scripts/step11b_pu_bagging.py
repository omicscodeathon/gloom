"""
step11b_pu_bagging.py — Mordelet-Vert PU Bagging
=================================================
Implements true Positive-Unlabeled (PU) learning using the bagging approach
from Mordelet & Vert (2014) "A bagging SVM to learn from positive and
unlabeled examples" (Neurocomputing, 163:73-83).

Algorithm (Option B from the technical review):
  P   = confirmed positive genes (LCGene in training set,  label = 1)
  U   = unlabeled genes          (non-LCGene in training set, label = 0)
  For i = 1 … B  (B = PU_N_ESTIMATORS):
      U_i  ← random subsample of U,  |U_i| = |P| × PU_SUBSAMPLE_RATIO
      T_i  ← P ∪ U_i    (P labelled 1, U_i labelled 0)
      f_i  ← train RandomForest on T_i
      s_i(x) ← f_i.predict_proba(x)[1]   for ALL genes in full universe
  s(x) ← (1/B) Σ_i s_i(x)               (ensemble PU score)

The model-building function lives in gloom_utils.pu_bagging_fit_predict (re-exported here as
``pu_bagging_fit_predict``) so that step11c (cross-fitting) and step13b (ablation) train
EXACTLY the same model.

IMPORTANT — what this step reports (v0.2.0):
  The score of a gene in the *training split* is NOT an independent evaluation: the full-data
  model has seen the training positives.  The all-gene score is therefore stored as the
  "full model" score (see also ``is_training_positive``) and the top-K reporting separates
    * held-out positives (validation split, never seen while fitting), from
    * training positives (resubstitution — shown only for transparency).
  The PRIMARY, fully out-of-sample ranking is produced by step11c (cross-fitting) and used by
  step14.

Key properties vs. standard classifier:
  — Never treats all unlabeled genes as confirmed negatives globally.
  — Each iteration sees a fresh random negative subsample → label noise
    averages out across B iterations.

Outputs:
  results/pu_bagging_scores.csv   — per-gene full-model PU score, rank, flags
  results/pu_bagging_metrics.csv  — validation-set AUROC/AUPRC + held-out vs training top-K
  figures/pu_bagging_score_distribution.png
"""
import logging, sys, time, warnings
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
from gloom_utils import pu_bagging_fit_predict  # noqa: F401  (re-exported for step11c / step13b)
config.create_output_dirs()
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

PU_N_ESTIMATORS    = getattr(config, "PU_N_ESTIMATORS",    100)
PU_SUBSAMPLE_RATIO = getattr(config, "PU_SUBSAMPLE_RATIO", 1.0)
PU_BASE_N_TREES    = getattr(config, "PU_BASE_N_TREES",    100)
CANDIDATE_THRESHOLD = getattr(config, "CANDIDATE_PROB_THRESHOLD",
                              getattr(config, "NOVEL_PROB_THRESHOLD", 0.70))


def run_pu_bagging():
    log.info("=" * 60)
    log.info("STEP 11b — MORDELET-VERT PU BAGGING")
    log.info(f"  B={PU_N_ESTIMATORS}  subsample_ratio={PU_SUBSAMPLE_RATIO}  "
             f"base_trees={PU_BASE_N_TREES}")
    log.info("=" * 60)

    # ── Load data ──────────────────────────────────────────────────────────────
    features = pd.read_csv(config.INTEGRATED_FEATURES_FILE, index_col=0)
    labels   = pd.read_csv(config.LABELS_FILE).set_index("gene")["label"]

    # Restrict scoring universe to the analysis universe (genes in gene_labels.csv).
    features = features.loc[features.index.isin(labels.index)]
    log.info(f"  Scoring universe: {len(features):,} genes")

    # Training gene index: only these genes are used to FIT the classifiers.
    train_idx = pd.read_csv(config.TRAIN_FEATURES_FILE, index_col=0).index
    val_idx   = pd.read_csv(config.VAL_FEATURES_FILE,   index_col=0).index

    X_tr = features.loc[features.index.isin(train_idx)]
    y_tr = labels.reindex(X_tr.index).fillna(0).astype(int)

    X_pos = X_tr[y_tr == 1]
    X_unl = X_tr[y_tr == 0]
    n_pos = len(X_pos)
    n_sub = max(1, int(round(n_pos * PU_SUBSAMPLE_RATIO)))

    log.info(f"  Training positives (P): {n_pos}")
    log.info(f"  Training unlabeled  (U): {len(X_unl):,}")
    log.info(f"  Subsample per iteration: {n_sub}")

    # ── PU Bagging loop (shared implementation) ────────────────────────────────
    t0 = time.time()

    def _progress(i, total):
        if i % 10 == 0 or i == total:
            log.info(f"  Iter {i:>4}/{total}  elapsed: {time.time() - t0:.1f}s")

    scores = pu_bagging_fit_predict(
        X_pos.values, X_unl.values, features.values,
        n_estimators=PU_N_ESTIMATORS,
        subsample_ratio=PU_SUBSAMPLE_RATIO,
        base_n_trees=PU_BASE_N_TREES,
        seed=int(config.SEED),
        progress=_progress,
    )
    pu_scores = pd.Series(scores, index=features.index, name="pu_score")
    log.info(f"  PU bagging complete in {time.time() - t0:.1f}s")

    # ── Build output DataFrame ─────────────────────────────────────────────────
    out = pd.DataFrame(index=pu_scores.index)
    out["pu_score"]       = pu_scores
    out["pu_rank"]        = pu_scores.rank(ascending=False, method="min").astype(int)
    out["label"]          = labels.reindex(pu_scores.index).fillna(0).astype(int)
    out["is_lcgene_gene"] = out["label"] == 1
    out["is_val_gene"]    = out.index.isin(val_idx)
    # Positives that were used to FIT the full model (their score is resubstitution, not
    # independent evidence) vs positives that were held out in the validation split.
    out["is_training_positive"] = out["is_lcgene_gene"] & out.index.isin(train_idx)
    out["is_heldout_positive"]  = out["is_lcgene_gene"] & out["is_val_gene"]

    # Query genes (optional)
    query_genes_file = getattr(config, "QUERY_GENES_FILE", None)
    if query_genes_file and Path(query_genes_file).exists():
        q_df      = pd.read_csv(Path(query_genes_file), sep="\t")
        query_set = set(q_df["GeneSymbol"].dropna().astype(str).str.strip().str.upper())
        out["is_query_gene"] = out.index.str.strip().str.upper().isin(query_set)
    else:
        out["is_query_gene"] = False

    # Non-LCGene candidate flag (mirrors step14 logic)
    if out["is_query_gene"].any():
        out["pu_non_lcgene_candidate"] = (
            out["is_query_gene"] & ~out["is_lcgene_gene"] & (out["pu_score"] >= CANDIDATE_THRESHOLD)
        )
    else:
        out["pu_non_lcgene_candidate"] = (
            ~out["is_lcgene_gene"] & (out["pu_score"] >= CANDIDATE_THRESHOLD)
        )

    out = out.sort_values("pu_score", ascending=False)

    # ── Validation-set evaluation (genes not seen during fitting) ──────────────
    val_mask = out["is_val_gene"]
    auroc = auprc = np.nan
    if val_mask.sum() > 0:
        y_v = out.loc[val_mask, "label"].values
        s_v = out.loc[val_mask, "pu_score"].values
        if y_v.sum() > 0 and y_v.sum() < len(y_v):
            auroc = roc_auc_score(y_v, s_v)
            auprc = average_precision_score(y_v, s_v)
            log.info(f"  Validation AUROC={auroc:.4f}  AUPRC={auprc:.4f}")
        else:
            log.warning("  Validation set has no positive genes — skipping AUROC/AUPRC.")

    # ── Top-K reporting: held-out vs training positives, kept SEPARATE ─────────
    # (a) Held-out ranking: genes of the validation split ranked by the full-model score.
    #     These genes were never seen while fitting, so this is the independent evaluation.
    # (b) Resubstitution: where the TRAINING positives land in the all-gene ranking.  This only
    #     shows that the model re-scores what it has already seen — NOT independent performance.
    held = out[out["is_val_gene"]].sort_values("pu_score", ascending=False)
    n_heldout_pos = int(held["is_lcgene_gene"].sum())
    n_train_pos   = int(out["is_training_positive"].sum())
    ks = (50, 100, 200)

    metrics_row = {
        "n_estimators":      PU_N_ESTIMATORS,
        "subsample_ratio":   PU_SUBSAMPLE_RATIO,
        "base_n_trees":      PU_BASE_N_TREES,
        "n_pos":             n_pos,
        "n_unlabeled_pool":  len(X_unl),
        "n_sub_per_iter":    n_sub,
        "n_heldout_genes":   len(held),
        "n_heldout_positives": n_heldout_pos,
        "n_training_positives": n_train_pos,
        "val_auroc":         round(float(auroc), 4) if not np.isnan(auroc) else np.nan,
        "val_auprc":         round(float(auprc), 4) if not np.isnan(auprc) else np.nan,
    }
    for k in ks:
        top_all  = out.head(k)
        top_held = held.head(k)
        hp = int(top_held["is_lcgene_gene"].sum())
        # independent (held-out genes only)
        metrics_row[f"lcgene_top{k}_heldout"]        = hp
        metrics_row[f"precision_at_{k}_heldout"]     = round(hp / max(len(top_held), 1), 4)
        metrics_row[f"recall_at_{k}_heldout"]        = round(hp / n_heldout_pos, 4) if n_heldout_pos else np.nan
        # resubstitution (training positives in the all-gene ranking) — NOT independent
        metrics_row[f"lcgene_top{k}_training"]       = int(top_all["is_training_positive"].sum())
        metrics_row[f"unlabeled_top{k}_allgenes"]    = int((~top_all["is_lcgene_gene"]).sum())
        log.info(f"  K={k:<4} held-out positives in held-out top-K: {hp}/{len(top_held)}   |   "
                 f"training positives in all-gene top-K: {metrics_row[f'lcgene_top{k}_training']}"
                 f" (resubstitution, not independent)")
    metrics_row["n_non_lcgene_candidates"] = int(out["pu_non_lcgene_candidate"].sum())
    metrics_row["candidate_threshold"]     = CANDIDATE_THRESHOLD
    log.info(f"  Non-LCGene candidates (score>={CANDIDATE_THRESHOLD}): "
             f"{metrics_row['n_non_lcgene_candidates']}")

    # ── Save outputs ───────────────────────────────────────────────────────────
    pu_scores_file = config.RESULTS_DIR / "pu_bagging_scores.csv"
    out.to_csv(pu_scores_file)
    metrics = pd.DataFrame([metrics_row])
    metrics.to_csv(config.RESULTS_DIR / "pu_bagging_metrics.csv", index=False)

    # ── Score distribution + enrichment curve plot ─────────────────────────────
    lcgene_scores    = out.loc[out["is_lcgene_gene"], "pu_score"]
    nonlcgene_scores = out.loc[~out["is_lcgene_gene"], "pu_score"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    bins = np.linspace(0, 1, 50)
    axes[0].hist(nonlcgene_scores, bins=bins, color="steelblue", alpha=0.6,
                 label=f"Unlabeled (n={len(nonlcgene_scores):,})", density=True)
    axes[0].hist(lcgene_scores, bins=bins, color="tomato", alpha=0.75,
                 label=f"LCGene (n={len(lcgene_scores):,})", density=True)
    axes[0].axvline(CANDIDATE_THRESHOLD, color="black", linestyle="--", lw=1.2,
                    label=f"Candidate threshold ({CANDIDATE_THRESHOLD})")
    axes[0].set_xlabel("PU Bagging Score (full model)"); axes[0].set_ylabel("Density")
    axes[0].set_title("PU Score Distribution: LCGene vs Unlabeled")
    axes[0].legend(fontsize=8)

    # Enrichment curve on HELD-OUT genes only (training positives would inflate it)
    sorted_h      = held.sort_values("pu_score", ascending=False)
    cum_lcgene    = sorted_h["is_lcgene_gene"].cumsum()
    total_lcgene  = sorted_h["is_lcgene_gene"].sum()
    frac_captured = cum_lcgene / total_lcgene if total_lcgene > 0 else cum_lcgene * 0
    frac_genes    = np.arange(1, len(sorted_h) + 1) / max(len(sorted_h), 1)
    axes[1].plot(frac_genes, frac_captured, color="tomato", lw=2, label="PU Bagging (held-out genes)")
    axes[1].plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5, label="Random")
    axes[1].fill_between(frac_genes, frac_captured, frac_genes, alpha=0.15, color="tomato")
    axes[1].set_xlabel("Fraction of held-out genes screened")
    axes[1].set_ylabel("Fraction of held-out LCGene genes captured")
    axes[1].set_title("Held-out LCGene Enrichment Curve — PU Bagging")
    axes[1].legend(fontsize=8)

    plt.tight_layout()
    fig.savefig(config.FIGURES_DIR / "pu_bagging_score_distribution.png",
                dpi=150, bbox_inches="tight")
    plt.close(fig)

    log.info(f"  Scores saved → {pu_scores_file}")
    log.info("STEP 11b COMPLETE")
    return {"pu_scores": out, "metrics": metrics}


if __name__ == "__main__":
    r = run_pu_bagging()
    cols = ["pu_score", "pu_rank", "is_lcgene_gene", "is_training_positive", "pu_non_lcgene_candidate"]
    print("\n--- Top 20 PU Bagging Genes (full model) ---")
    print(r["pu_scores"].head(20)[[c for c in cols if c in r["pu_scores"].columns]].to_string())
    print("\n--- Metrics ---")
    print(r["metrics"].T.to_string(header=False))
