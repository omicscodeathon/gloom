#!/usr/bin/env python
"""
benchmark_runtime.py — measured run time and peak memory of the main GLOOM computations
=======================================================================================
Provides the evidence needed to say anything about scalability: for several problem sizes
(default 2000 / 5000 / 10000 genes  x  100 / 300 / 500 samples) it times and measures the
peak memory of the main computational steps on SUBSAMPLED input:

    differential_expression   Welch t-test + BH over all genes            (step 4)
    coexpression_network      chunked Pearson correlation + edge list     (step 6, |r| >= cutoff)
    network_summary           sparse degree / components / clustering     (step 7c)
    pu_crossfit               PU-bagging cross-fitting, 5 folds x 1 repeat (steps 11b/11c), small B
    enrichment                hypergeometric over-representation, 200 gene sets (step 19)

Input matrices: the harmonized TCGA/GTEx matrices when they exist (genes and samples are randomly
subsampled with config.SEED), otherwise a synthetic log2-expression matrix with latent factors
(so that the correlation network has realistic density).  Use --source synthetic to force it.

Every size is measured in a fresh worker process so that the peak resident memory is not
contaminated by previous sizes.  Peak memory is reported two ways: ``peak_tracemalloc_mb`` (Python /
NumPy allocations during the step) and ``peak_rss_mb`` (peak resident set size of the worker
process so far; requires ``resource`` on Unix or ``psutil`` on Windows, otherwise empty).

Output: results/runtime_benchmark.csv (one row per size x step, with machine information).
Usage:
    python scripts/benchmark_runtime.py
    python scripts/benchmark_runtime.py --genes 2000 5000 --samples 100 300 --pu-estimators 20
    python scripts/benchmark_runtime.py --source synthetic --repeats 3

Numbers from one machine do not generalise: report the machine columns together with the table.
"""
import argparse
import json
import os
import platform
import subprocess
import sys
import time
import tracemalloc
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
# Prefer the packaged copy of the pipeline (src/gloom/pipeline): its config.py resolves the data and
# output directories relative to the repository root.  Fall back to this folder when run standalone.
_PKG = HERE.parent / "src" / "gloom" / "pipeline"
sys.path.insert(0, str(_PKG if (_PKG / "gloom_utils.py").exists() else HERE))

RESULT_PREFIX = "BENCH_JSON:"


# ── machine information ───────────────────────────────────────────────────────────────────────────

def machine_info() -> dict:
    info = {
        "platform": platform.platform(),
        "processor": platform.processor(),
        "cpu_count_logical": os.cpu_count(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
    try:
        import sklearn
        info["scikit_learn"] = sklearn.__version__
    except ImportError:
        info["scikit_learn"] = ""
    try:
        import scipy
        info["scipy"] = scipy.__version__
    except ImportError:
        info["scipy"] = ""
    try:
        import psutil
        info["ram_total_gb"] = round(psutil.virtual_memory().total / 1e9, 1)
        info["cpu_count_physical"] = psutil.cpu_count(logical=False)
    except ImportError:
        info["ram_total_gb"] = ""
        info["cpu_count_physical"] = ""
    return info


def peak_rss_mb():
    """Peak resident set size of this process in MB (None when it cannot be determined)."""
    try:
        import resource
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return peak / 1e6 if sys.platform == "darwin" else peak / 1e3   # bytes on macOS, KB on Linux
    except ImportError:
        pass
    try:
        import psutil
        mi = psutil.Process().memory_info()
        peak = getattr(mi, "peak_wset", None)        # Windows
        return (peak if peak is not None else mi.rss) / 1e6
    except Exception:
        return None


# ── input data ────────────────────────────────────────────────────────────────────────────────────

def make_input(n_genes, n_samples, source, seed):
    """Return a genes x samples log2-expression matrix (DataFrame) and the source actually used."""
    rng = np.random.default_rng(seed)
    if source in ("auto", "real"):
        try:
            import config
            path = Path(config.TUMOR_EXPR_HARMONIZED)
            if path.exists():
                df = pd.read_csv(path, index_col=0)
                genes = rng.choice(df.index, size=min(n_genes, len(df)), replace=False)
                cols = rng.choice(df.columns, size=min(n_samples, df.shape[1]), replace=False)
                return df.loc[genes, cols], "harmonized_tumor_matrix"
        except Exception:
            pass
        if source == "real":
            raise FileNotFoundError("Harmonized expression matrix not found; run steps 1-3 first.")
    # synthetic: low-rank latent factors + noise -> realistic, moderately dense correlation structure
    k = 12
    loadings = rng.normal(0, 1, size=(n_genes, k)) * (rng.random((n_genes, k)) < 0.15)
    factors = rng.normal(0, 1, size=(k, n_samples))
    expr = 5 + loadings @ factors + rng.normal(0, 1, size=(n_genes, n_samples))
    genes = [f"G{i:05d}" for i in range(n_genes)]
    cols = [f"S{j:04d}" for j in range(n_samples)]
    return pd.DataFrame(expr, index=genes, columns=cols), "synthetic"


# ── measured steps ────────────────────────────────────────────────────────────────────────────────

def _measure(fn):
    tracemalloc.start()
    t0 = time.perf_counter()
    out = fn()
    seconds = time.perf_counter() - t0
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return out, seconds, peak / 1e6


def worker(n_genes, n_samples, source, seed, pu_estimators, cutoff):
    from scipy import stats
    from gloom_utils import (adjacency_from_edges, benjamini_hochberg, crossfit_oof_scores,
                             hypergeom_enrichment, network_summary, pu_bagging_fit_predict)
    from step6_coexpression_network import compute_pearson_correlation_chunked

    expr, used = make_input(n_genes, n_samples, source, seed)
    n_genes, n_samples = expr.shape
    rng = np.random.default_rng(seed)
    results = []

    def record(step, seconds, peak_mb, extra=""):
        results.append({
            "step": step, "genes": n_genes, "samples": n_samples, "input": used,
            "seconds": round(seconds, 3), "peak_tracemalloc_mb": round(peak_mb, 1),
            "peak_rss_mb": None if peak_rss_mb() is None else round(peak_rss_mb(), 1),
            "note": extra,
        })

    # 1) differential expression (half / half split of the samples)
    half = n_samples // 2
    A, B = expr.iloc[:, :half].to_numpy(), expr.iloc[:, half:].to_numpy()

    def de():
        _, p = stats.ttest_ind(A, B, axis=1, equal_var=False, nan_policy="omit")
        return benjamini_hochberg(np.where(np.isnan(p), 1.0, p))
    _, s, m = _measure(de)
    record("differential_expression", s, m)

    # 2) co-expression network (step 6 code)
    X = expr.to_numpy(dtype=np.float32)
    genes = expr.index.tolist()
    edges, s, m = _measure(lambda: compute_pearson_correlation_chunked(X, genes, cutoff=cutoff))
    record("coexpression_network", s, m, f"edges={len(edges)} cutoff={cutoff}")

    # 3) sparse network summary
    gi = pd.Index(genes)
    a_idx, b_idx = gi.get_indexer(edges["gene_a"]), gi.get_indexer(edges["gene_b"])

    def summary():
        return network_summary(adjacency_from_edges(a_idx, b_idx, len(genes)))
    summ, s, m = _measure(summary)
    record("network_summary", s, m, f"largest_component={summ['largest_component']}")

    # 4) PU-bagging cross-fitting on a small expression-derived feature matrix
    feats = np.column_stack([X.mean(axis=1), X.std(axis=1), np.percentile(X, 25, axis=1),
                             np.percentile(X, 75, axis=1), summ["degree"]])
    y = (rng.random(n_genes) < 0.05).astype(int)
    feats[:, 0] += 1.5 * y                                   # give the labels a learnable signal

    def pu_fit(Xtr, ytr, Xte, sd):
        return pu_bagging_fit_predict(Xtr[ytr == 1], Xtr[ytr == 0], Xte,
                                      n_estimators=pu_estimators, base_n_trees=50, seed=sd)

    def crossfit():
        return crossfit_oof_scores(feats, y, pu_fit, n_splits=5, n_repeats=1, seed=seed)
    _, s, m = _measure(crossfit)
    record("pu_crossfit", s, m, f"B={pu_estimators} folds=5 repeats=1")

    # 5) enrichment: 200 random gene sets of 100 genes, query of 200 genes
    gene_sets = {f"set{i}": rng.choice(genes, size=min(100, n_genes), replace=False).tolist()
                 for i in range(200)}
    query = rng.choice(genes, size=min(200, n_genes), replace=False).tolist()
    _, s, m = _measure(lambda: hypergeom_enrichment(query, gene_sets, genes))
    record("enrichment", s, m, "200 gene sets x 100 genes, query 200")

    print(RESULT_PREFIX + json.dumps(results))


# ── driver ────────────────────────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", type=int, nargs="+", default=[2000, 5000, 10000])
    ap.add_argument("--samples", type=int, nargs="+", default=[100, 300, 500])
    ap.add_argument("--repeats", type=int, default=1, help="repeat each size (default 1)")
    ap.add_argument("--source", choices=["auto", "real", "synthetic"], default="auto")
    ap.add_argument("--pu-estimators", type=int, default=20, help="PU-bagging B used in the benchmark")
    ap.add_argument("--cutoff", type=float, default=0.6, help="|r| cut-off of the co-expression network")
    ap.add_argument("--seed", type=int, default=None, help="default: config.SEED (or 42)")
    ap.add_argument("--out", default=None, help="CSV path (default: <results dir>/runtime_benchmark.csv)")
    ap.add_argument("--worker", nargs=2, type=int, metavar=("GENES", "SAMPLES"), help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    seed = args.seed
    if seed is None:
        try:
            import config
            seed = int(config.SEED)
        except Exception:
            seed = 42

    if args.worker:
        worker(args.worker[0], args.worker[1], args.source, seed, args.pu_estimators, args.cutoff)
        return 0

    out = args.out
    if out is None:
        try:
            import config
            out = Path(config.RESULTS_DIR) / "runtime_benchmark.csv"
        except Exception:
            out = HERE.parent / "results" / "runtime_benchmark.csv"
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)

    info = machine_info()
    rows = []
    for g in args.genes:
        for s in args.samples:
            for rep in range(args.repeats):
                cmd = [sys.executable, str(Path(__file__).resolve()), "--worker", str(g), str(s),
                       "--source", args.source, "--pu-estimators", str(args.pu_estimators),
                       "--cutoff", str(args.cutoff), "--seed", str(seed + rep)]
                print(f"[benchmark] genes={g} samples={s} repeat={rep + 1}/{args.repeats} …", flush=True)
                proc = subprocess.run(cmd, capture_output=True, text=True)
                payload = [ln for ln in proc.stdout.splitlines() if ln.startswith(RESULT_PREFIX)]
                if proc.returncode != 0 or not payload:
                    tail = (proc.stderr or proc.stdout).strip().splitlines()[-3:]
                    print(f"[benchmark]   FAILED (exit {proc.returncode}): {' | '.join(tail)}", flush=True)
                    rows.append({"step": "FAILED", "genes": g, "samples": s, "repeat": rep + 1,
                                 "note": " | ".join(tail)[:300], **info})
                    continue
                for r in json.loads(payload[-1][len(RESULT_PREFIX):]):
                    r["repeat"] = rep + 1
                    rows.append({**r, **info})
                    print(f"[benchmark]   {r['step']:<26} {r['seconds']:>8.2f} s   "
                          f"peak tracemalloc {r['peak_tracemalloc_mb']:>8.1f} MB", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print(f"[benchmark] wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
