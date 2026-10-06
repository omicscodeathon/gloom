"""
step13b_ablation.py — Feature-set ablation under cross-fitting
==============================================================
Answers the reviewer's question "do network features improve prediction beyond expression /
differential-expression information?" by evaluating, with the SAME cross-fitting protocol
(same folds, same PU-bagging model, same bootstrap), the following scorers:

  a_abs_log2fc              rank genes by |log2FC| alone                      (no model)
  b_adj_pvalue              rank genes by adjusted p-value alone               (no model;
                            neg_log10_padj, ties broken by |t statistic| because BH-adjusted
                            p-values saturate)
  c_expression_only         PU bagging on expression-derived features
  d_network_only            PU bagging on network-derived features
  e_expression_plus_network PU bagging on expression + network features
  f_combined_no_de_stats    (e) minus the explicit DE / log2FC statistics
                            NOTE: tumor_mean and normal_mean remain, so fold change is still
                            implicitly recoverable by the model; (f) removes the *explicit*
                            log2FC-derived columns only.
  g_full_model              PU bagging on every integrated feature (the GLOOM model)

Feature groups come from config.FEATURE_GROUP_* (matched against the columns produced by
steps 5, 7, 7b and 8; only existing columns are used).

Outputs (results/):
  ablation_scores.csv            per-gene scores of every variant (+ label)
  ablation_metrics.csv           variant x metric with 95 % bootstrap CIs
                                 (AUROC, AUPRC, average precision, P@K, R@K, EF@K)
  ablation_paired_comparison.csv paired-bootstrap difference (estimate, CI, p-value) for
                                 expression-only vs full model and other key pairs
  ablation_verdict.txt           plain-text verdict: does network integration add a
                                 measurable gain over expression-only?  (yes/no + CI)

Verdict rule: the gain is called "measurable" only if the 95 % paired-bootstrap CI of
Delta-AUPRC (full model minus expression-only) lies entirely above 0.
"""
import logging, sys, time, warnings
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
from gloom_utils import (bootstrap_ranking_metrics, crossfit_oof_scores,
                         paired_bootstrap_difference, resolve_feature_groups,
                         tie_broken_score)
from step11c_crossfit_pu import make_pu_score_fn, load_universe_features

config.create_output_dirs()
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

DESCRIPTIONS = {
    "a_abs_log2fc":              "|log2FC| alone (no model)",
    "b_adj_pvalue":              "adjusted p-value alone (no model)",
    "c_expression_only":         "expression-derived features",
    "d_network_only":            "network-derived features",
    "e_expression_plus_network": "expression + network features",
    "f_combined_no_de_stats":    "expression + network without explicit DE/log2FC statistics",
    "g_full_model":              "all integrated features (GLOOM full model)",
}

# (variant A, variant B) pairs for the paired bootstrap of metric(A) - metric(B)
PAIRS = [
    ("g_full_model", "c_expression_only"),
    ("e_expression_plus_network", "c_expression_only"),
    ("g_full_model", "a_abs_log2fc"),
]


def _baseline_scores(index) -> dict:
    """Model-free baselines from the DE table, aligned to ``index`` (genes)."""
    de = pd.read_csv(config.DE_RESULTS_FILE, index_col=0).reindex(index)
    abs_fc = de["log2fc"].abs().fillna(0.0).to_numpy()
    neg_log_p = de["neg_log10_padj"].fillna(0.0).to_numpy()
    tie = de["t_stat"].abs().fillna(0.0).to_numpy() if "t_stat" in de.columns else np.zeros(len(de))
    return {
        "a_abs_log2fc": abs_fc,
        "b_adj_pvalue": tie_broken_score(neg_log_p, tie),
    }


def _format_ci(est, lo, hi, digits=3):
    return f"{est:+.{digits}f} (95% CI {lo:+.{digits}f} to {hi:+.{digits}f})"


def write_verdict(paired: pd.DataFrame, groups: dict, path: Path) -> str:
    """Create the one-paragraph plain-text verdict from the paired comparison table."""
    if not groups["network_only"]:
        text = ("Ablation verdict: NOT EVALUABLE - no network-derived feature columns were found "
                "(run steps 6/7/7b/8).")
    else:
        sub = paired[(paired["a"] == "g_full_model") & (paired["b"] == "c_expression_only")]
        row_ap = sub[sub["metric"] == "auprc"].iloc[0]
        row_au = sub[sub["metric"] == "auroc"].iloc[0]
        gain = bool(row_ap["ci_low"] > 0)
        text = (
            f"Does network integration add a measurable gain over expression-only features? "
            f"{'YES' if gain else 'NO'}. "
            f"Full model minus expression-only: delta-AUPRC = "
            f"{_format_ci(row_ap['diff'], row_ap['ci_low'], row_ap['ci_high'])}, "
            f"paired-bootstrap p = {row_ap['p_value']:.3g}; "
            f"delta-AUROC = {_format_ci(row_au['diff'], row_au['ci_low'], row_au['ci_high'])}, "
            f"p = {row_au['p_value']:.3g}. "
            f"Rule: 'yes' only if the 95% CI of delta-AUPRC lies entirely above 0. "
            + ("Network features therefore improve ranking beyond expression-derived features in this run."
               if gain else
               "Network-derived features should be described as supporting contextual interpretation "
               "rather than as the source of predictive improvement in this run.")
        )
    path.write_text(text + "\n", encoding="utf-8")
    return text


def run_ablation() -> dict:
    log.info("=" * 60)
    log.info("STEP 13b — FEATURE-SET ABLATION (cross-fitted)")
    log.info("=" * 60)
    if not getattr(config, "USE_ABLATION", True):
        log.info("  USE_ABLATION = False — skipped.")
        return {}

    res_dir = Path(config.RESULTS_DIR)
    features, y_ser = load_universe_features()
    y = y_ser.to_numpy()
    groups = resolve_feature_groups(
        features.columns, config.FEATURE_GROUP_EXPRESSION,
        config.FEATURE_GROUP_NETWORK, config.FEATURE_GROUP_DE_STATS,
    )
    for name in ("expression_only", "network_only", "expression_plus_network",
                 "combined_without_de_stats", "full"):
        log.info(f"  feature set {name:<26}: {len(groups[name])} features")
    if groups["unassigned"]:
        log.info(f"  columns in no group (used only by the full model): {groups['unassigned']}")

    n_splits = int(getattr(config, "CROSSFIT_FOLDS", 5))
    cf_rep = int(getattr(config, "CROSSFIT_REPEATS", 5))
    cf_est = getattr(config, "CROSSFIT_PU_N_ESTIMATORS", None) or getattr(config, "PU_N_ESTIMATORS", 100)
    n_rep = int(getattr(config, "ABLATION_REPEATS", None) or cf_rep)
    n_est = getattr(config, "ABLATION_PU_N_ESTIMATORS", None) or cf_est
    ratio = getattr(config, "PU_SUBSAMPLE_RATIO", 1.0)
    trees = getattr(config, "PU_BASE_N_TREES", 100)
    ks = tuple(getattr(config, "EVAL_K_VALUES", (10, 50, 100)))
    n_boot = int(getattr(config, "BOOTSTRAP_N", 1000))
    seed = int(getattr(config, "BOOTSTRAP_SEED", config.SEED))
    alpha = float(getattr(config, "BOOTSTRAP_ALPHA", 0.05))
    log.info(f"  Cross-fitting: K={n_splits}, repeats={n_rep}, PU bagging B={n_est}, trees={trees}")

    scores = {"label": y}
    scores.update(_baseline_scores(features.index))

    model_sets = [
        ("c_expression_only", groups["expression_only"]),
        ("d_network_only", groups["network_only"]),
        ("e_expression_plus_network", groups["expression_plus_network"]),
        ("f_combined_no_de_stats", groups["combined_without_de_stats"]),
        ("g_full_model", groups["full"]),
    ]
    score_fn = make_pu_score_fn(n_est, ratio, trees)
    oof_main = res_dir / "oof_scores.csv"
    for name, cols in model_sets:
        if not cols:
            log.warning(f"  {name}: no feature columns available — variant skipped.")
            continue
        if name == "g_full_model" and oof_main.exists() and n_rep == cf_rep and n_est == cf_est:
            reuse = pd.read_csv(oof_main, index_col=0)
            if reuse.index.equals(features.index) or set(reuse.index) == set(features.index):
                scores[name] = reuse["oof_score"].reindex(features.index).to_numpy()
                log.info("  g_full_model: re-using step11c out-of-fold scores (identical settings).")
                continue
        t0 = time.time()
        log.info(f"  {name}: {len(cols)} features …")
        oof_mean, _, _ = crossfit_oof_scores(
            features[cols], y, score_fn, n_splits=n_splits, n_repeats=n_rep,
            seed=int(config.SEED), index=features.index,
        )
        scores[name] = oof_mean.to_numpy()
        log.info(f"  {name}: done in {time.time() - t0:.0f}s")

    score_df = pd.DataFrame(scores, index=features.index)
    score_df.index.name = "gene"
    score_df.to_csv(res_dir / "ablation_scores.csv")

    # ── metrics with bootstrap CIs (long format) ────────────────────────────────────────────
    variants = [v for v in DESCRIPTIONS if v in score_df.columns]
    rows = []
    for v in variants:
        m = bootstrap_ranking_metrics(y, score_df[v].to_numpy(), ks=ks, n_boot=n_boot,
                                      seed=seed, alpha=alpha)
        m.insert(0, "variant", v)
        m.insert(1, "description", DESCRIPTIONS[v])
        n_feat = {"c_expression_only": len(groups["expression_only"]),
                  "d_network_only": len(groups["network_only"]),
                  "e_expression_plus_network": len(groups["expression_plus_network"]),
                  "f_combined_no_de_stats": len(groups["combined_without_de_stats"]),
                  "g_full_model": len(groups["full"])}.get(v, 0)
        m.insert(2, "n_features", n_feat)
        rows.append(m)
    metrics = pd.concat(rows, ignore_index=True)
    metrics.to_csv(res_dir / "ablation_metrics.csv", index=False)
    wide = metrics.pivot(index="variant", columns="metric", values="estimate")
    log.info("\n" + wide[[c for c in ("auroc", "auprc") if c in wide.columns]].round(4).to_string())

    # ── paired bootstrap comparisons ─────────────────────────────────────────────────────────
    prow = []
    for a, b in PAIRS:
        if a in score_df.columns and b in score_df.columns:
            d = paired_bootstrap_difference(y, score_df[a].to_numpy(), score_df[b].to_numpy(),
                                            ks=ks, n_boot=n_boot, seed=seed, alpha=alpha)
            d.insert(0, "b", b)
            d.insert(0, "a", a)
            prow.append(d)
    paired = pd.concat(prow, ignore_index=True) if prow else pd.DataFrame(
        columns=["a", "b", "metric", "estimate_a", "estimate_b", "diff", "ci_low", "ci_high",
                 "p_value", "n_boot"])
    paired.to_csv(res_dir / "ablation_paired_comparison.csv", index=False)

    if (paired["a"] == "g_full_model").any() and "c_expression_only" in score_df.columns:
        verdict = write_verdict(paired, groups, res_dir / "ablation_verdict.txt")
    else:
        verdict = "Ablation verdict: NOT EVALUABLE - full-model or expression-only scores are missing."
        (res_dir / "ablation_verdict.txt").write_text(verdict + "\n", encoding="utf-8")
    log.info(f"  VERDICT: {verdict}")
    log.info("STEP 13b COMPLETE")
    return {"scores": score_df, "metrics": metrics, "paired": paired, "verdict": verdict}


if __name__ == "__main__":
    r = run_ablation()
    if r:
        print(r["verdict"])
