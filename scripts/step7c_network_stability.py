"""
step7c_network_stability.py — Co-expression network stability analysis
=======================================================================
Addresses the question: is the large difference in network density between the tumor and the
normal network (e.g. 110,508 vs 1,276,368 edges at |r| >= 0.60) a robust, interpretable
feature, or a consequence of the correlation threshold / sample size / cohort?

The tumor network (and the normal network when available) is rebuilt with the correlation code
of step6 (compute_pearson_correlation_chunked) in two ways:

1. Threshold sweep, |r| >= each value of config.NETWORK_THRESHOLDS (0.50 … 0.70):
   number of edges, density, mean degree, size of the largest connected component, mean
   clustering coefficient, number of isolated nodes, and the Spearman correlation of the
   degree centrality between consecutive thresholds.  When config.NETWORK_EQUAL_N is True the
   sweep is repeated after sub-sampling the larger group to the size of the smaller group
   ("equal_n"), because the number of samples alone changes the amount of spurious correlation.

2. Sample bootstrap (config.NETWORK_BOOTSTRAP_N resamples, default 20) at the working cut-off
   config.COEXPR_CORRELATION_CUTOFF: samples are resampled (with replacement by default; to the
   same size in both groups when NETWORK_EQUAL_N), the network is rebuilt, and we record
     * edge retention frequency of the reference-network edges,
     * edge recall / precision versus the reference network,
     * hub stability: Jaccard index between the reference hubs (top NETWORK_HUB_FRACTION genes
       by degree, default 5 %) and the hubs of each resample, plus per-hub frequencies,
     * Spearman correlation of the degree vector with the reference.

Outputs (results/):
  network_stability_thresholds.csv           one row per group x sampling x threshold
  network_stability_density_ratio.csv        tumor / normal edge-count ratio per threshold
  network_stability_bootstrap_replicates.csv one row per group x resample
  network_stability_bootstrap_summary.csv    summary per group
  network_stability_hubs.csv                 reference hubs with their bootstrap hub frequency

Optional step: it never stops the pipeline.
"""
import logging, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
from gloom_utils import adjacency_from_edges, network_summary
from step6_coexpression_network import compute_pearson_correlation_chunked

config.create_output_dirs()
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)


def _spearman(a, b) -> float:
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(spearmanr(a, b).correlation)


def _edge_indices(edges: pd.DataFrame, gene_index: pd.Index):
    """Integer endpoints (a < b) of an edge table returned by step6."""
    a = gene_index.get_indexer(edges["gene_a"])
    b = gene_index.get_indexer(edges["gene_b"])
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    return lo.astype(np.int64), hi.astype(np.int64), edges["correlation"].to_numpy(dtype=float)


def threshold_sweep(group: str, sampling: str, expr: pd.DataFrame, thresholds) -> list:
    """Network statistics of one expression matrix at several |r| thresholds."""
    thresholds = sorted(float(t) for t in thresholds)
    genes = expr.index
    n = len(genes)
    edges = compute_pearson_correlation_chunked(
        expr.to_numpy(dtype=np.float32), genes.tolist(), cutoff=thresholds[0])
    a, b, r = _edge_indices(edges, genes)
    rows, prev_deg = [], None
    for t in thresholds:
        keep = np.abs(r) >= t
        A = adjacency_from_edges(a[keep], b[keep], n)
        s = network_summary(A)
        rows.append({
            "group": group, "sampling": sampling, "threshold": t, "n_samples": expr.shape[1],
            "n_nodes": s["n_nodes"], "n_edges": s["n_edges"], "density": s["density"],
            "mean_degree": s["mean_degree"], "largest_component": s["largest_component"],
            "clustering_coefficient": s["clustering_coefficient"],
            "isolated_nodes": s["isolated_nodes"],
            "degree_spearman_vs_previous_threshold":
                _spearman(prev_deg, s["degree"]) if prev_deg is not None else np.nan,
        })
        prev_deg = s["degree"]
        log.info(f"  [{group}/{sampling}] |r|>={t:.2f}: {s['n_edges']:,} edges, "
                 f"density={s['density']:.5f}, largest component={s['largest_component']:,}, "
                 f"isolated={s['isolated_nodes']:,}")
    return rows


def _hubs(degree: np.ndarray, fraction: float) -> np.ndarray:
    k = max(1, int(np.ceil(fraction * len(degree))))
    return np.argsort(-degree, kind="mergesort")[:k]


def bootstrap_stability(group: str, expr: pd.DataFrame, cutoff: float, n_boot: int,
                        n_use: int, replace: bool, hub_fraction: float, rng):
    """Sample bootstrap of one group's network at ``cutoff`` (see module docstring)."""
    genes = expr.index
    n = len(genes)
    X = expr.to_numpy(dtype=np.float32)
    ref_edges = compute_pearson_correlation_chunked(X, genes.tolist(), cutoff=cutoff)
    ra, rb, _ = _edge_indices(ref_edges, genes)
    ref_keys = np.unique(ra * n + rb)
    A_ref = adjacency_from_edges(ra, rb, n)
    ref_deg = np.asarray(A_ref.sum(axis=1)).ravel()
    ref_hubs = _hubs(ref_deg, hub_fraction)
    ref_hub_set = set(ref_hubs.tolist())
    log.info(f"  [{group}] reference network at |r|>={cutoff}: {len(ref_keys):,} edges, "
             f"{len(ref_hub_set)} hubs; {n_boot} resamples of {n_use} samples "
             f"({'with' if replace else 'without'} replacement)")

    freq = np.zeros(len(ref_keys), dtype=np.int64)
    hub_count = np.zeros(n, dtype=np.int64)
    deg_sum = np.zeros(n, dtype=np.float64)
    replicates = []
    for rep in range(n_boot):
        t0 = time.time()
        cols = rng.choice(X.shape[1], size=n_use, replace=replace)
        edges_b = compute_pearson_correlation_chunked(X[:, cols], genes.tolist(), cutoff=cutoff)
        ba, bb, _ = _edge_indices(edges_b, genes)
        keys_b = np.unique(ba * n + bb)
        present = np.isin(ref_keys, keys_b, assume_unique=True)
        freq += present
        A_b = adjacency_from_edges(ba, bb, n)
        deg_b = np.asarray(A_b.sum(axis=1)).ravel()
        hubs_b = set(_hubs(deg_b, hub_fraction).tolist())
        hub_count[list(hubs_b)] += 1
        deg_sum += deg_b
        inter = len(ref_hub_set & hubs_b)
        union = len(ref_hub_set | hubs_b)
        n_common = int(present.sum())
        replicates.append({
            "group": group, "replicate": rep + 1, "n_samples_used": n_use,
            "n_edges": len(keys_b),
            "edge_recall_vs_reference": n_common / max(len(ref_keys), 1),
            "edge_precision_vs_reference": n_common / max(len(keys_b), 1),
            "hub_jaccard_vs_reference": inter / union if union else np.nan,
            "degree_spearman_vs_reference": _spearman(ref_deg, deg_b),
        })
        log.info(f"  [{group}] resample {rep + 1}/{n_boot}: {len(keys_b):,} edges, "
                 f"recall={replicates[-1]['edge_recall_vs_reference']:.3f}, "
                 f"hub Jaccard={replicates[-1]['hub_jaccard_vs_reference']:.3f} ({time.time() - t0:.0f}s)")

    retention = freq / max(n_boot, 1)
    rep_df = pd.DataFrame(replicates)
    summary = {
        "group": group, "reference_cutoff": cutoff, "n_boot": n_boot, "n_samples_used": n_use,
        "reference_edges": len(ref_keys), "hub_fraction": hub_fraction,
        "edge_retention_mean": float(retention.mean()) if len(retention) else np.nan,
        "edge_retention_median": float(np.median(retention)) if len(retention) else np.nan,
        "edge_retention_p25": float(np.percentile(retention, 25)) if len(retention) else np.nan,
        "fraction_edges_retained_ge_50pct": float((retention >= 0.5).mean()) if len(retention) else np.nan,
        "fraction_edges_retained_ge_80pct": float((retention >= 0.8).mean()) if len(retention) else np.nan,
        "fraction_edges_retained_ge_95pct": float((retention >= 0.95).mean()) if len(retention) else np.nan,
        "hub_jaccard_mean": float(rep_df["hub_jaccard_vs_reference"].mean()),
        "hub_jaccard_sd": float(rep_df["hub_jaccard_vs_reference"].std(ddof=1)) if n_boot > 1 else np.nan,
        "degree_spearman_mean": float(rep_df["degree_spearman_vs_reference"].mean()),
    }
    hubs_df = pd.DataFrame({
        "group": group,
        "gene": genes[ref_hubs],
        "reference_degree": ref_deg[ref_hubs],
        "hub_frequency": hub_count[ref_hubs] / max(n_boot, 1),
        "mean_bootstrap_degree": deg_sum[ref_hubs] / max(n_boot, 1),
    }).sort_values("reference_degree", ascending=False)
    return rep_df, summary, hubs_df


def run_network_stability() -> dict:
    log.info("=" * 60)
    log.info("STEP 7c — NETWORK STABILITY (thresholds + sample bootstrap)")
    log.info("=" * 60)

    rng = np.random.default_rng(int(config.SEED))
    groups = {"tumor": pd.read_csv(config.TUMOR_EXPR_HARMONIZED, index_col=0)}
    normal_path = Path(config.NORMAL_EXPR_HARMONIZED)
    if normal_path.exists():
        groups["normal"] = pd.read_csv(normal_path, index_col=0)
    else:
        log.warning(f"  Normal matrix not found ({normal_path}) — tumor network only.")
    for g, df in groups.items():
        log.info(f"  {g}: {df.shape[0]:,} genes x {df.shape[1]:,} samples")

    thresholds = tuple(getattr(config, "NETWORK_THRESHOLDS", (0.50, 0.55, 0.60, 0.65, 0.70)))
    equal_n = bool(getattr(config, "NETWORK_EQUAL_N", True)) and len(groups) == 2
    n_equal = min(df.shape[1] for df in groups.values())
    sizes_differ = len({df.shape[1] for df in groups.values()}) > 1

    # Sampling schemes: the full data, and (optionally) the larger group sub-sampled to equal size
    schemes = {"full": dict(groups)}
    if equal_n and sizes_differ:
        eq = {}
        for g, df in groups.items():
            if df.shape[1] > n_equal:
                cols = rng.choice(df.shape[1], size=n_equal, replace=False)
                eq[g] = df.iloc[:, np.sort(cols)]
            else:
                eq[g] = df
        schemes["equal_n"] = eq
        log.info(f"  NETWORK_EQUAL_N: larger group sub-sampled to {n_equal} samples.")

    # ── 1. threshold sweep ─────────────────────────────────────────────────────────────────────
    rows = []
    for sampling, mats in schemes.items():
        for g, df in mats.items():
            rows += threshold_sweep(g, sampling, df, thresholds)
    thr = pd.DataFrame(rows)
    res_dir = Path(config.RESULTS_DIR)
    thr.to_csv(res_dir / "network_stability_thresholds.csv", index=False)

    if {"tumor", "normal"} <= set(thr["group"]):
        piv = thr.pivot_table(index=["sampling", "threshold"], columns="group", values="n_edges").reset_index()
        piv = piv.rename(columns={"tumor": "tumor_edges", "normal": "normal_edges"})
        piv["normal_to_tumor_edge_ratio"] = piv["normal_edges"] / piv["tumor_edges"].replace(0, np.nan)
        piv.to_csv(res_dir / "network_stability_density_ratio.csv", index=False)
        log.info("\n" + piv.round(3).to_string(index=False))

    # ── 2. sample bootstrap ────────────────────────────────────────────────────────────────────
    n_boot = int(getattr(config, "NETWORK_BOOTSTRAP_N", 20))
    replace = bool(getattr(config, "NETWORK_BOOTSTRAP_REPLACE", True))
    hub_fraction = float(getattr(config, "NETWORK_HUB_FRACTION", 0.05))
    cutoff = float(config.COEXPR_CORRELATION_CUTOFF)
    boot_scheme = "equal_n" if "equal_n" in schemes else "full"
    reps, sums, hubs = [], [], []
    if n_boot > 0:
        for g, df in schemes[boot_scheme].items():
            n_use = df.shape[1] if boot_scheme == "equal_n" else (n_equal if equal_n else df.shape[1])
            rep_df, summary, hubs_df = bootstrap_stability(
                g, df, cutoff, n_boot, n_use, replace, hub_fraction, rng)
            summary["sampling"] = boot_scheme
            reps.append(rep_df); sums.append(summary); hubs.append(hubs_df)
        pd.concat(reps, ignore_index=True).to_csv(res_dir / "network_stability_bootstrap_replicates.csv", index=False)
        summ = pd.DataFrame(sums)
        summ.to_csv(res_dir / "network_stability_bootstrap_summary.csv", index=False)
        pd.concat(hubs, ignore_index=True).to_csv(res_dir / "network_stability_hubs.csv", index=False)
        log.info("\n" + summ.round(3).to_string(index=False))
    else:
        log.info("  NETWORK_BOOTSTRAP_N = 0 — bootstrap skipped.")

    log.info("STEP 7c COMPLETE")
    return {"thresholds": thr}


if __name__ == "__main__":
    out = run_network_stability()
    print(out["thresholds"].round(4).to_string(index=False))
