"""
step12b_oof_metrics.py — Out-of-fold ranking metrics with bootstrap confidence intervals
========================================================================================
Evaluates the cross-fitted scores of step11c (results/oof_scores.csv).  Because every gene was
scored by models that never saw it, these are genuinely out-of-sample numbers.

Metrics (all computed over the whole analysis universe, LCGene = positive):
  AUROC, AUPRC (trapezoidal), average precision,
  Precision@K, Recall@K, enrichment factor@K   for K in config.EVAL_K_VALUES (default 10/50/100)

Uncertainty: 95 % percentile bootstrap CIs from a STRATIFIED bootstrap over genes (positives
and unlabeled genes are resampled separately), config.BOOTSTRAP_N resamples (default 1000)
with the fixed seed config.BOOTSTRAP_SEED.

Caveat: in a PU setting the "unlabeled" genes include undiscovered positives, so precision
and AUROC are conservative (lower-bound) estimates of the true performance.

Output: results/oof_metrics.csv  (metric, estimate, ci_low, ci_high, n_boot, n_genes, n_positives)
"""
import logging, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
from gloom_utils import bootstrap_ranking_metrics

config.create_output_dirs()
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)


def run_oof_metrics() -> pd.DataFrame:
    log.info("=" * 60)
    log.info("STEP 12b — OUT-OF-FOLD METRICS (bootstrap 95% CIs)")
    log.info("=" * 60)
    oof_path = Path(config.RESULTS_DIR) / "oof_scores.csv"
    if not oof_path.exists():
        if not getattr(config, "USE_CROSSFIT", True):
            log.info("  USE_CROSSFIT = False and no oof_scores.csv — skipped.")
            return pd.DataFrame()
        raise FileNotFoundError(f"{oof_path} not found. Run step11c first.")

    oof = pd.read_csv(oof_path, index_col=0)
    y = oof["label"].to_numpy().astype(int)
    s = oof["oof_score"].to_numpy(dtype=float)
    ks = tuple(getattr(config, "EVAL_K_VALUES", (10, 50, 100)))
    n_boot = int(getattr(config, "BOOTSTRAP_N", 1000))
    seed = int(getattr(config, "BOOTSTRAP_SEED", config.SEED))
    alpha = float(getattr(config, "BOOTSTRAP_ALPHA", 0.05))
    log.info(f"  {len(y):,} genes, {int(y.sum())} positives; K={ks}; bootstrap n={n_boot}, seed={seed}")

    metrics = bootstrap_ranking_metrics(y, s, ks=ks, n_boot=n_boot, seed=seed, alpha=alpha)
    metrics["n_genes"] = len(y)
    metrics["n_positives"] = int(y.sum())

    # Number of LCGene positives among the top-K OOF-ranked genes (no CI; a plain count)
    order = np.argsort(-s, kind="mergesort")
    extra = []
    for k in ks:
        extra.append({"metric": f"n_lcgene_in_top_{k}", "estimate": float(y[order[:k]].sum()),
                      "ci_low": np.nan, "ci_high": np.nan, "n_boot": 0,
                      "n_genes": len(y), "n_positives": int(y.sum())})
    metrics = pd.concat([metrics, pd.DataFrame(extra)], ignore_index=True)

    out_path = Path(config.RESULTS_DIR) / "oof_metrics.csv"
    metrics.to_csv(out_path, index=False)
    log.info("\n" + metrics.round(4).to_string(index=False))
    log.info(f"  Saved -> {out_path}")
    log.info("STEP 12b COMPLETE")
    return metrics


if __name__ == "__main__":
    m = run_oof_metrics()
    print(m.round(4).to_string(index=False))
