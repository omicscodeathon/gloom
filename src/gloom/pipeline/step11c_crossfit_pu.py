"""
step11c_crossfit_pu.py — Cross-fitted (out-of-fold) PU-bagging scores
=====================================================================
Gives every gene a score produced by models that NEVER saw that gene during training, so that
the final ranking and all top-K statistics are genuinely out-of-sample (reviewer point 4:
"98 of the top 100 genes are LCGene positives" only demonstrated recovery of training positives).

Procedure (repeated R times with different seeds, default R = 5):
  1. Split ALL genes of the analysis universe into K stratified folds (default K = 5),
     stratified by the PU label (positive / unlabeled).
  2. For every fold k: fit the SAME PU-bagging model as step11b (Mordelet-Vert bagging of
     random forests, gloom_utils.pu_bagging_fit_predict) on the positives + unlabeled genes of
     the OTHER K-1 folds, and predict ONLY the held-out fold k.
  3. Concatenate the held-out predictions into one out-of-fold (OOF) score per gene.
  4. Average the OOF scores over the R repeats  ->  results/oof_scores.csv.

Notes
  * All integrated features are used (no importance-based feature selection, which would be
    computed from labels of the whole universe and leak into the held-out folds).
  * The PU-bagging size defaults to config.PU_N_ESTIMATORS (identical to step11b); set
    config.CROSSFIT_PU_N_ESTIMATORS to a smaller value to shorten the run.
  * Cost: K x R model fits of B base forests each (default 5 x 5 x 300 forests).

Outputs:
  results/oof_scores.csv   — gene, oof_score (mean over repeats), oof_score_sd, oof_rank,
                             label, is_lcgene_gene, repeat_1 … repeat_R
  results/oof_fold_assignments.csv — fold id of every gene in every repeat
"""
import logging, sys, time, warnings
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
from gloom_utils import crossfit_oof_scores, pu_bagging_fit_predict

config.create_output_dirs()
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)


def make_pu_score_fn(n_estimators, subsample_ratio, base_n_trees):
    """Adapter: (X_train, y_train, X_test, seed) -> PU-bagging scores for X_test."""
    def _score(X_train, y_train, X_test, seed):
        return pu_bagging_fit_predict(
            X_train[y_train == 1], X_train[y_train == 0], X_test,
            n_estimators=n_estimators, subsample_ratio=subsample_ratio,
            base_n_trees=base_n_trees, seed=seed,
        )
    return _score


def load_universe_features():
    """Integrated features restricted to the analysis universe + the aligned PU labels."""
    features = pd.read_csv(config.INTEGRATED_FEATURES_FILE, index_col=0)
    labels = pd.read_csv(config.LABELS_FILE).set_index("gene")["label"]
    features = features.loc[features.index.isin(labels.index)]
    y = labels.reindex(features.index).fillna(0).astype(int)
    return features, y


def run_crossfit_pu() -> dict:
    log.info("=" * 60)
    log.info("STEP 11c — CROSS-FITTED PU BAGGING (out-of-fold scores)")
    log.info("=" * 60)
    if not getattr(config, "USE_CROSSFIT", True):
        log.info("  USE_CROSSFIT = False — skipped. Step14 will rank with the full-data model, "
                 "whose top-K contains training positives (NOT an independent evaluation).")
        return {}

    features, y = load_universe_features()
    n_splits = int(getattr(config, "CROSSFIT_FOLDS", 5))
    n_repeats = int(getattr(config, "CROSSFIT_REPEATS", 5))
    n_est = getattr(config, "CROSSFIT_PU_N_ESTIMATORS", None) or getattr(config, "PU_N_ESTIMATORS", 100)
    ratio = getattr(config, "PU_SUBSAMPLE_RATIO", 1.0)
    trees = getattr(config, "PU_BASE_N_TREES", 100)
    log.info(f"  Universe: {len(features):,} genes ({int(y.sum())} positives, "
             f"{int((y == 0).sum()):,} unlabeled), {features.shape[1]} features")
    log.info(f"  K={n_splits} folds x R={n_repeats} repeats, PU bagging B={n_est}, base trees={trees}")

    t0 = time.time()
    total = n_splits * n_repeats

    def _progress(rep, fold, done):
        log.info(f"  repeat {rep + 1}/{n_repeats}  fold {fold + 1}/{n_splits}  "
                 f"({done}/{total} fits, {time.time() - t0:.0f}s)")

    oof_mean, oof_rep, folds = crossfit_oof_scores(
        features, y.to_numpy(),
        make_pu_score_fn(n_est, ratio, trees),
        n_splits=n_splits, n_repeats=n_repeats, seed=int(config.SEED),
        index=features.index, progress=_progress,
    )

    out = pd.DataFrame(index=features.index)
    out.index.name = "gene"
    out["oof_score"] = oof_mean
    out["oof_score_sd"] = oof_rep.std(axis=1, ddof=0)
    out["oof_rank"] = out["oof_score"].rank(ascending=False, method="min").astype(int)
    out["label"] = y
    out["is_lcgene_gene"] = out["label"] == 1
    out = pd.concat([out, oof_rep], axis=1).sort_values("oof_score", ascending=False)

    res_dir = Path(config.RESULTS_DIR)
    out.to_csv(res_dir / "oof_scores.csv")
    folds.to_csv(res_dir / "oof_fold_assignments.csv")
    log.info(f"  OOF scores saved -> {res_dir / 'oof_scores.csv'}  ({time.time() - t0:.0f}s)")
    log.info("STEP 11c COMPLETE")
    return {"oof_scores": out, "fold_assignments": folds}


if __name__ == "__main__":
    r = run_crossfit_pu()
    if r:
        print(r["oof_scores"].head(20)[["oof_score", "oof_score_sd", "oof_rank", "is_lcgene_gene"]].to_string())
