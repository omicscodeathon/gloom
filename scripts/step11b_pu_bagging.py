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

Key properties vs. standard classifier:
  — Never treats all unlabeled genes as confirmed negatives globally.
  — Each iteration sees a fresh random negative subsample → label noise
    averages out across B iterations.
  — Genes consistently scored high across many iterations are robust
    novel candidates, not artefacts of a single training configuration.

Outputs:
  results/pu_bagging_scores.csv   — per-gene PU score, rank, novel flag
  results/pu_bagging_metrics.csv  — validation-set AUROC/AUPRC + ranking metrics
  figures/pu_bagging_score_distribution.png
"""
import logging, sys, time, warnings
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

PU_N_ESTIMATORS    = getattr(config, "PU_N_ESTIMATORS",    100)
PU_SUBSAMPLE_RATIO = getattr(config, "PU_SUBSAMPLE_RATIO", 1.0)
PU_BASE_N_TREES    = getattr(config, "PU_BASE_N_TREES",    100)
NOVEL_THRESHOLD    = getattr(config, "NOVEL_PROB_THRESHOLD", 0.70)


def run_pu_bagging():
    log.info("=" * 60)
    log.info("STEP 11b — MORDELET-VERT PU BAGGING")
    log.info(f"  B={PU_N_ESTIMATORS}  subsample_ratio={PU_SUBSAMPLE_RATIO}  "
             f"base_trees={PU_BASE_N_TREES}")
    log.info("=" * 60)

    # ── Load data ──────────────────────────────────────────────────────────────
    features = pd.read_csv(config.INTEGRATED_FEATURES_FILE, index_col=0)
    labels   = pd.read_csv(config.LABELS_FILE).set_index("gene")["label"]

    # Restrict scoring universe to the trimmed gene universe (genes in gene_labels.csv).
    # Without this, genes removed by step9 trimming (abs_log2fc < 2.0) get scored and
    # dominate top rankings despite having no meaningful LUAD DE signal.
    features = features.loc[features.index.isin(labels.index)]
    log.info(f"  Scoring universe: {len(features):,} genes (trimmed universe)")

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

    # ── PU Bagging loop ────────────────────────────────────────────────────────
    scores_accum = np.zeros(len(features), dtype=np.float64)
    X_all_arr    = features.values.astype(np.float32)
    X_pos_arr    = X_pos.values.astype(np.float32)
    X_unl_arr    = X_unl.values.astype(np.float32)

    rng = np.random.default_rng(config.RANDOM_STATE)
    t0  = time.time()

    for i in range(PU_N_ESTIMATORS):
        # Sample unlabeled subset (without replacement when possible)
        replace = n_sub > len(X_unl_arr)
        sel_idx = rng.choice(len(X_unl_arr), size=n_sub, replace=replace)
        X_u_sub = X_unl_arr[sel_idx]

        # Build balanced training set for this iteration
        X_iter = np.vstack([X_pos_arr, X_u_sub])
        y_iter = np.concatenate([np.ones(n_pos, dtype=np.int8),
                                  np.zeros(n_sub, dtype=np.int8)])

        # Train a lightweight RF base classifier
        clf = RandomForestClassifier(
            n_estimators=PU_BASE_N_TREES,
            max_features="sqrt",
            class_weight="balanced",
            n_jobs=-1,
            random_state=int(config.RANDOM_STATE) + i,
        )
        clf.fit(X_iter, y_iter)

        # Accumulate scores over the full gene universe
        scores_accum += clf.predict_proba(X_all_arr)[:, 1]

        if (i + 1) % 10 == 0 or i == PU_N_ESTIMATORS - 1:
            log.info(f"  Iter {i + 1:>4}/{PU_N_ESTIMATORS}  elapsed: {time.time() - t0:.1f}s")

    pu_scores = pd.Series(
        scores_accum / PU_N_ESTIMATORS,
        index=features.index,
        name="pu_score",
    )
    log.info(f"  PU bagging complete in {time.time() - t0:.1f}s")

    # ── Build output DataFrame ─────────────────────────────────────────────────
    out = pd.DataFrame(index=pu_scores.index)
    out["pu_score"]       = pu_scores
    out["pu_rank"]        = pu_scores.rank(ascending=False, method="min").astype(int)
    out["label"]          = labels.reindex(pu_scores.index).fillna(0).astype(int)
    out["is_lcgene_gene"] = out["label"] == 1
    out["is_val_gene"]    = out.index.isin(val_idx)

    # Query genes (optional)
    query_genes_file = getattr(config, "QUERY_GENES_FILE", None)
    if query_genes_file and Path(query_genes_file).exists():
        q_df      = pd.read_csv(Path(query_genes_file), sep="\t")
        query_set = set(q_df["GeneSymbol"].dropna().astype(str).str.strip().str.upper())
        out["is_query_gene"] = out.index.str.strip().str.upper().isin(query_set)
    else:
        out["is_query_gene"] = False

    # Novel candidate flag (mirrors step14 logic)
    if out["is_query_gene"].any():
        out["pu_novel_candidate"] = (
            out["is_query_gene"] & ~out["is_lcgene_gene"] & (out["pu_score"] >= NOVEL_THRESHOLD)
        )
    else:
        out["pu_novel_candidate"] = (
            ~out["is_lcgene_gene"] & (out["pu_score"] >= NOVEL_THRESHOLD)
        )

    out = out.sort_values("pu_score", ascending=False)

    # ── Validation-set evaluation (proxy; genes not seen during fitting) ───────
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

    # ── Ranking metrics (full-genome) ──────────────────────────────────────────
    total_pos = out["is_lcgene_gene"].sum()

    def _recall_at(k):
        return out.head(k)["is_lcgene_gene"].sum() / total_pos if total_pos > 0 else 0.0

    def _prec_at(k):
        return float(out.head(k)["is_lcgene_gene"].mean())

    lcgene_top50  = int(out.head(50)["is_lcgene_gene"].sum())
    lcgene_top100 = int(out.head(100)["is_lcgene_gene"].sum())
    lcgene_top200 = int(out.head(200)["is_lcgene_gene"].sum())
    log.info(f"  LCGene in top-50: {lcgene_top50}  top-100: {lcgene_top100}  top-200: {lcgene_top200}")
    log.info(f"  Recall@50={_recall_at(50):.4f}  Recall@100={_recall_at(100):.4f}  "
             f"Recall@200={_recall_at(200):.4f}")
    log.info(f"  Precision@50={_prec_at(50):.4f}  Precision@100={_prec_at(100):.4f}")
    log.info(f"  PU novel candidates (score>={NOVEL_THRESHOLD}): {out['pu_novel_candidate'].sum()}")

    # ── Save outputs ───────────────────────────────────────────────────────────
    pu_scores_file = config.RESULTS_DIR / "pu_bagging_scores.csv"
    out.to_csv(pu_scores_file)

    metrics = pd.DataFrame([{
        "n_estimators":      PU_N_ESTIMATORS,
        "subsample_ratio":   PU_SUBSAMPLE_RATIO,
        "base_n_trees":      PU_BASE_N_TREES,
        "n_pos":             n_pos,
        "n_unlabeled_pool":  len(X_unl),
        "n_sub_per_iter":    n_sub,
        "val_auroc":         round(float(auroc), 4) if not np.isnan(auroc) else np.nan,
        "val_auprc":         round(float(auprc), 4) if not np.isnan(auprc) else np.nan,
        "lcgene_top50":      lcgene_top50,
        "lcgene_top100":     lcgene_top100,
        "lcgene_top200":     lcgene_top200,
        "recall_at_50":      round(_recall_at(50),  4),
        "recall_at_100":     round(_recall_at(100), 4),
        "recall_at_200":     round(_recall_at(200), 4),
        "precision_at_50":   round(_prec_at(50),  4),
        "precision_at_100":  round(_prec_at(100), 4),
        "n_pu_novel":        int(out["pu_novel_candidate"].sum()),
        "novel_threshold":   NOVEL_THRESHOLD,
    }])
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
    axes[0].axvline(NOVEL_THRESHOLD, color="black", linestyle="--", lw=1.2,
                    label=f"Novel threshold ({NOVEL_THRESHOLD})")
    axes[0].set_xlabel("PU Bagging Score"); axes[0].set_ylabel("Density")
    axes[0].set_title("PU Score Distribution: LCGene vs Unlabeled")
    axes[0].legend(fontsize=8)

    sorted_out    = out.sort_values("pu_rank")
    cum_lcgene    = sorted_out["is_lcgene_gene"].cumsum()
    total_lcgene  = sorted_out["is_lcgene_gene"].sum()
    frac_captured = cum_lcgene / total_lcgene if total_lcgene > 0 else cum_lcgene * 0
    frac_genes    = np.arange(1, len(sorted_out) + 1) / len(sorted_out)
    axes[1].plot(frac_genes, frac_captured, color="tomato", lw=2, label="PU Bagging")
    axes[1].plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5, label="Random")
    axes[1].fill_between(frac_genes, frac_captured, frac_genes, alpha=0.15, color="tomato")
    axes[1].set_xlabel("Fraction of genes screened")
    axes[1].set_ylabel("Fraction of LCGene genes captured")
    axes[1].set_title("LCGene Enrichment Curve — PU Bagging")
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
    cols = ["pu_score", "pu_rank", "is_lcgene_gene", "pu_novel_candidate"]
    print("\n--- Top 20 PU Bagging Genes ---")
    print(r["pu_scores"].head(20)[[c for c in cols if c in r["pu_scores"].columns]].to_string())
    print("\n--- Metrics ---")
    print(r["metrics"].to_string(index=False))
