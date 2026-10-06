"""
gloom_utils.py
--------------
Pure helper functions shared by several pipeline steps.

This module deliberately has NO dependency on ``config`` and creates no files,
directories or log handlers at import time.  Every parameter is passed in
explicitly, which makes the functions easy to unit-test with synthetic data
(see ``test/``).

Contents
  1. Expression sanitising        : sanitize_expression
  2. Cohort design checks         : detect_cohort_collinearity, update_report_section
  3. Ranking metrics + bootstrap  : compute_ranking_metrics, bootstrap_ranking_metrics,
                                    paired_bootstrap_difference, tie_broken_score
  4. Cross-fitting                : make_stratified_folds, crossfit_oof_scores,
                                    pu_bagging_fit_predict
  5. Enrichment                   : benjamini_hochberg, hypergeom_enrichment, read_gmt
  6. Network summaries            : adjacency_from_edges, network_summary
  7. Sample-level QC verdicts     : evaluate_qc_warnings
  8. Patient pairing              : infer_patient_ids, pair_samples_by_patient
  9. Ablation feature groups      : resolve_feature_groups
 10. Candidate evidence classes   : classify_candidate_evidence
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


# ======================================================================
# 1. Expression sanitising
# ======================================================================

def sanitize_expression(df: pd.DataFrame, label: str = "") -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Classify every cell of an expression matrix and return a cleaned copy.

    Rules (a true zero is NOT a missing value):
      * finite value >= 0   -> kept unchanged (zeros stay zeros, so that
                               log2(0 + 1) = 0 downstream)
      * NaN                 -> kept as NaN (missing observation); counted
      * negative value      -> invalid for a count/TPM/RSEM matrix; set to NaN
      * +inf / -inf         -> invalid; set to NaN

    Returns
    -------
    (clean_df, report)
        ``report`` holds the cell counts: n_cells, n_zero_kept, n_nan_kept,
        n_negative_invalid, n_inf_invalid, n_invalid.
    """
    values = df.apply(pd.to_numeric, errors="coerce")
    arr = values.to_numpy(dtype=float)
    is_nan = np.isnan(arr)
    is_inf = np.isinf(arr)
    is_neg = (arr < 0) & ~is_inf
    invalid = is_inf | is_neg

    report = {
        "n_cells":            int(arr.size),
        "n_zero_kept":        int((arr == 0).sum()),
        "n_nan_kept":         int(is_nan.sum()),
        "n_negative_invalid": int(is_neg.sum()),
        "n_inf_invalid":      int(is_inf.sum()),
        "n_invalid":          int(invalid.sum()),
    }

    tag = f"[{label}] " if label else ""
    log.info(
        f"  {tag}Value audit: {report['n_zero_kept']:,} true zeros kept (log2(0+1)=0), "
        f"{report['n_nan_kept']:,} missing (NaN) kept as NaN."
    )
    if report["n_invalid"] > 0:
        log.warning(
            f"  {tag}INVALID values detected and set to NaN: "
            f"{report['n_negative_invalid']:,} negative, {report['n_inf_invalid']:,} +/-inf. "
            f"Check the upstream quantification (a count/TPM/RSEM matrix must be >= 0)."
        )

    clean = values.mask(pd.DataFrame(invalid, index=values.index, columns=values.columns))
    return clean, report


# ======================================================================
# 2. Cohort design checks
# ======================================================================

def detect_cohort_collinearity(groups: pd.Series, cohorts: pd.Series) -> Dict[str, object]:
    """
    Detect perfect collinearity between cohort (batch) and biological group.

    The design is confounded when every cohort contains samples of a single
    group only (e.g. all tumors from cohort A, all normals from cohort B);
    in that case no batch-correction method can separate cohort from disease.

    Parameters
    ----------
    groups, cohorts : pd.Series with the same index (one entry per sample).

    Returns
    -------
    dict with keys: confounded (bool), table (cohort x group crosstab),
    n_cohorts, n_groups.
    """
    df = pd.DataFrame({"group": groups.astype(str), "cohort": cohorts.astype(str)}).dropna()
    table = pd.crosstab(df["cohort"], df["group"])
    n_groups = int(table.shape[1])
    n_cohorts = int(table.shape[0])
    groups_per_cohort = (table > 0).sum(axis=1)
    confounded = bool(n_groups > 1 and n_cohorts > 1 and (groups_per_cohort <= 1).all())
    return {
        "confounded": confounded,
        "table": table,
        "n_cohorts": n_cohorts,
        "n_groups": n_groups,
    }


def update_report_section(path, title: str, lines: Optional[Sequence[str]]) -> Optional[Path]:
    """
    Maintain a small multi-section text report (used for qc_cohort_warning.txt).

    Each section starts with a ``## [title]`` line.  The section called
    ``title`` is replaced by ``lines``; when ``lines`` is empty/None the
    section is removed.  If the file ends up empty it is deleted.  Other
    sections (written by other steps) are preserved.
    """
    path = Path(path)
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    marker = f"## [{title}]"
    blocks = [b for b in re.split(r"(?m)^(?=## \[)", existing) if b.strip()]
    kept = [b.rstrip("\n") + "\n" for b in blocks if not b.startswith(marker)]
    if lines:
        kept.append(marker + "\n" + "\n".join(lines).rstrip("\n") + "\n")
    if kept:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(kept), encoding="utf-8")
        return path
    if path.exists():
        path.unlink()
    return None


# ======================================================================
# 3. Ranking metrics and bootstrap
# ======================================================================

def _as_binary(y) -> np.ndarray:
    return np.asarray(y).astype(int)


def auroc_score(y, s) -> float:
    """Tie-aware AUROC (Mann-Whitney U / (n_pos * n_neg))."""
    from scipy.stats import rankdata

    y = _as_binary(y)
    s = np.asarray(s, dtype=float)
    n_pos = int(y.sum())
    n_neg = len(y) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = rankdata(s)
    return float((ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def _pr_points(y: np.ndarray, s: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Precision/recall at every distinct score threshold (descending score)."""
    order = np.argsort(-s, kind="mergesort")
    ys = y[order]
    ss = s[order]
    tp = np.cumsum(ys)
    fp = np.cumsum(1 - ys)
    last_of_tie = np.r_[np.flatnonzero(np.diff(ss) != 0), len(ss) - 1]
    tp = tp[last_of_tie]
    fp = fp[last_of_tie]
    precision = tp / (tp + fp)
    recall = tp / max(int(y.sum()), 1)
    return precision, recall


def average_precision_score_np(y, s) -> float:
    """Average precision = sum_n (R_n - R_{n-1}) * P_n (step-wise, no interpolation)."""
    y = _as_binary(y)
    s = np.asarray(s, dtype=float)
    if y.sum() == 0:
        return float("nan")
    precision, recall = _pr_points(y, s)
    recall_prev = np.r_[0.0, recall[:-1]]
    return float(np.sum((recall - recall_prev) * precision))


def auprc_trapezoid(y, s) -> float:
    """Area under the precision-recall curve by the trapezoidal rule (curve starts at (0, 1))."""
    y = _as_binary(y)
    s = np.asarray(s, dtype=float)
    if y.sum() == 0:
        return float("nan")
    precision, recall = _pr_points(y, s)
    r = np.r_[0.0, recall]
    p = np.r_[1.0, precision]
    return float(np.sum((r[1:] - r[:-1]) * (p[1:] + p[:-1]) / 2.0))


def compute_ranking_metrics(y, s, ks: Sequence[int] = (10, 50, 100)) -> Dict[str, float]:
    """
    Ranking metrics for a binary label vector ``y`` and a score vector ``s``
    (higher score = more likely positive).

    Keys: auroc, auprc (trapezoidal), average_precision, and for every K in
    ``ks``: precision_at_K, recall_at_K, ef_at_K (enrichment factor =
    precision@K divided by the positive base rate).  Ties in the top-K cut are
    broken by input order (stable sort).
    """
    y = _as_binary(y)
    s = np.asarray(s, dtype=float)
    n = len(y)
    n_pos = int(y.sum())

    out: Dict[str, float] = {
        "auroc": auroc_score(y, s),
        "auprc": auprc_trapezoid(y, s),
        "average_precision": average_precision_score_np(y, s),
    }
    order = np.argsort(-s, kind="mergesort")
    cum_pos = np.cumsum(y[order])
    base_rate = n_pos / n if n else float("nan")
    for k in ks:
        kk = int(min(k, n))
        tp = float(cum_pos[kk - 1]) if kk > 0 else 0.0
        prec = tp / kk if kk > 0 else float("nan")
        out[f"precision_at_{k}"] = prec
        out[f"recall_at_{k}"] = tp / n_pos if n_pos > 0 else float("nan")
        out[f"ef_at_{k}"] = prec / base_rate if base_rate and base_rate > 0 else float("nan")
    return out


def tie_broken_score(primary, secondary) -> np.ndarray:
    """
    Convert a (primary, secondary) pair of sort keys into a single score where
    a larger value means a better rank: higher ``primary`` first, ties broken
    by higher ``secondary``.  Useful for baselines such as "adjusted p-value
    alone" whose values saturate (many genes share the minimum p-value).
    """
    primary = np.asarray(primary, dtype=float)
    secondary = np.asarray(secondary, dtype=float)
    order = np.lexsort((secondary, primary))  # ascending: last key = primary
    score = np.empty(len(primary), dtype=float)
    score[order] = np.arange(len(primary), dtype=float)
    return score


def stratified_bootstrap_indices(y, rng: np.random.Generator) -> np.ndarray:
    """Resample positives and negatives separately (with replacement), keeping class counts fixed."""
    y = _as_binary(y)
    pos = np.flatnonzero(y == 1)
    neg = np.flatnonzero(y == 0)
    return np.concatenate([
        rng.choice(pos, size=len(pos), replace=True),
        rng.choice(neg, size=len(neg), replace=True),
    ])


def bootstrap_ranking_metrics(
    y,
    s,
    ks: Sequence[int] = (10, 50, 100),
    n_boot: int = 1000,
    seed: int = 42,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """
    Point estimates plus percentile bootstrap CIs (stratified resampling over genes).

    Returns a DataFrame with columns: metric, estimate, ci_low, ci_high, n_boot.
    """
    y = _as_binary(y)
    s = np.asarray(s, dtype=float)
    point = compute_ranking_metrics(y, s, ks)
    rng = np.random.default_rng(seed)
    boot = {m: np.empty(n_boot) for m in point}
    for b in range(n_boot):
        idx = stratified_bootstrap_indices(y, rng)
        res = compute_ranking_metrics(y[idx], s[idx], ks)
        for m in point:
            boot[m][b] = res[m]
    rows = []
    for m, est in point.items():
        lo, hi = np.nanpercentile(boot[m], [100 * alpha / 2, 100 * (1 - alpha / 2)])
        rows.append({"metric": m, "estimate": est, "ci_low": lo, "ci_high": hi, "n_boot": n_boot})
    return pd.DataFrame(rows)


def paired_bootstrap_difference(
    y,
    s_a,
    s_b,
    ks: Sequence[int] = (10, 50, 100),
    n_boot: int = 1000,
    seed: int = 42,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """
    Paired (same resampled genes for both scorers) bootstrap of metric(A) - metric(B).

    Returns a DataFrame with columns: metric, estimate_a, estimate_b, diff,
    ci_low, ci_high, p_value, n_boot.  The two-sided p-value is
    2 * min(P(diff <= 0), P(diff >= 0)) with +1 smoothing, capped at 1.
    """
    y = _as_binary(y)
    s_a = np.asarray(s_a, dtype=float)
    s_b = np.asarray(s_b, dtype=float)
    pa = compute_ranking_metrics(y, s_a, ks)
    pb = compute_ranking_metrics(y, s_b, ks)
    rng = np.random.default_rng(seed)
    diffs = {m: np.empty(n_boot) for m in pa}
    for b in range(n_boot):
        idx = stratified_bootstrap_indices(y, rng)
        ra = compute_ranking_metrics(y[idx], s_a[idx], ks)
        rb = compute_ranking_metrics(y[idx], s_b[idx], ks)
        for m in pa:
            diffs[m][b] = ra[m] - rb[m]
    rows = []
    for m in pa:
        d = diffs[m][~np.isnan(diffs[m])]
        if len(d) == 0:
            rows.append({"metric": m, "estimate_a": pa[m], "estimate_b": pb[m],
                         "diff": pa[m] - pb[m], "ci_low": np.nan, "ci_high": np.nan,
                         "p_value": np.nan, "n_boot": n_boot})
            continue
        lo, hi = np.percentile(d, [100 * alpha / 2, 100 * (1 - alpha / 2)])
        p_le = (np.sum(d <= 0) + 1.0) / (len(d) + 1.0)
        p_ge = (np.sum(d >= 0) + 1.0) / (len(d) + 1.0)
        rows.append({
            "metric": m, "estimate_a": pa[m], "estimate_b": pb[m], "diff": pa[m] - pb[m],
            "ci_low": lo, "ci_high": hi, "p_value": float(min(1.0, 2 * min(p_le, p_ge))),
            "n_boot": n_boot,
        })
    return pd.DataFrame(rows)


# ======================================================================
# 4. Cross-fitting
# ======================================================================

def make_stratified_folds(y, n_splits: int, seed: int) -> np.ndarray:
    """
    Assign every gene to one of ``n_splits`` folds, stratified by class.

    The returned integer vector partitions the genes: each gene belongs to
    exactly one fold and the class proportions are (almost) equal in all folds.
    """
    y = _as_binary(y)
    if n_splits < 2:
        raise ValueError("n_splits must be >= 2.")
    rng = np.random.default_rng(seed)
    folds = np.full(len(y), -1, dtype=int)
    for cls in np.unique(y):
        idx = np.flatnonzero(y == cls)
        rng.shuffle(idx)
        offset = int(rng.integers(n_splits))   # spread the remainder over random folds
        folds[idx] = (np.arange(len(idx)) + offset) % n_splits
    return folds


def crossfit_oof_scores(
    X,
    y,
    score_fn: Callable[[np.ndarray, np.ndarray, np.ndarray, int], np.ndarray],
    n_splits: int = 5,
    n_repeats: int = 5,
    seed: int = 42,
    index: Optional[Sequence] = None,
    progress: Optional[Callable[[int, int, int], None]] = None,
) -> Tuple[pd.Series, pd.DataFrame, pd.DataFrame]:
    """
    Repeated stratified K-fold cross-fitting.

    For every repeat and every fold, ``score_fn(X_train, y_train, X_test, seed)``
    is fitted on the OTHER folds and predicts ONLY the held-out fold.  The
    held-out predictions are concatenated into one out-of-fold (OOF) score per
    gene per repeat; the final OOF score is the mean across repeats.  A gene is
    therefore always scored by a model that never saw it during training.

    Parameters
    ----------
    X : (n_genes, n_features) DataFrame or ndarray.
    y : binary labels (1 = positive, 0 = unlabeled), length n_genes.
    score_fn : callable returning a 1-D score array for the rows of X_test.
    index : gene identifiers (defaults to X.index when X is a DataFrame).
    progress : optional callback(repeat, fold, n_done_total).

    Returns
    -------
    oof_mean : pd.Series (mean OOF score across repeats)
    oof_repeats : pd.DataFrame (genes x repeats)
    fold_table : pd.DataFrame (genes x repeats; fold id of every gene per repeat)
    """
    if index is None:
        index = X.index if hasattr(X, "index") else np.arange(len(y))
    X_arr = X.to_numpy() if hasattr(X, "to_numpy") else np.asarray(X)
    y = _as_binary(y)
    n = len(y)
    if X_arr.shape[0] != n:
        raise ValueError("X and y must have the same number of rows.")

    oof = np.full((n, n_repeats), np.nan)
    fold_ids = np.zeros((n, n_repeats), dtype=int)
    done = 0
    for rep in range(n_repeats):
        rep_seed = int(seed) + 1000 * rep
        folds = make_stratified_folds(y, n_splits, rep_seed)
        fold_ids[:, rep] = folds
        for f in range(n_splits):
            test_mask = folds == f
            if not test_mask.any():
                continue
            train_mask = ~test_mask
            if y[train_mask].sum() == 0:
                raise ValueError("A training split has no positive gene; reduce n_splits.")
            scores = np.asarray(
                score_fn(X_arr[train_mask], y[train_mask], X_arr[test_mask], rep_seed + f),
                dtype=float,
            )
            if scores.shape[0] != int(test_mask.sum()):
                raise ValueError("score_fn must return one score per held-out row.")
            oof[test_mask, rep] = scores
            done += 1
            if progress is not None:
                progress(rep, f, done)

    if np.isnan(oof).any():
        raise RuntimeError("Cross-fitting left genes without an out-of-fold score.")
    cols = [f"repeat_{r + 1}" for r in range(n_repeats)]
    oof_repeats = pd.DataFrame(oof, index=index, columns=cols)
    fold_table = pd.DataFrame(fold_ids, index=index, columns=cols)
    return oof_repeats.mean(axis=1).rename("oof_score"), oof_repeats, fold_table


def pu_bagging_fit_predict(
    X_pos,
    X_unl,
    X_pred,
    n_estimators: int = 100,
    subsample_ratio: float = 1.0,
    base_n_trees: int = 100,
    seed: int = 42,
    n_jobs: int = -1,
    progress: Optional[Callable[[int, int], None]] = None,
) -> np.ndarray:
    """
    Mordelet-Vert PU bagging (the model of step11b).

    For i in 1..B: draw a random subsample of the unlabeled set of size
    round(|P| * subsample_ratio), train a RandomForest on P (label 1) +
    subsample (label 0) and predict ``X_pred``.  The ensemble score is the mean
    predicted probability over the B base classifiers.
    """
    from sklearn.ensemble import RandomForestClassifier

    X_pos = np.asarray(X_pos, dtype=np.float32)
    X_unl = np.asarray(X_unl, dtype=np.float32)
    X_pred = np.asarray(X_pred, dtype=np.float32)
    n_pos = len(X_pos)
    if n_pos == 0 or len(X_unl) == 0:
        raise ValueError("PU bagging needs at least one positive and one unlabeled gene.")
    n_sub = max(1, int(round(n_pos * subsample_ratio)))

    rng = np.random.default_rng(seed)
    accum = np.zeros(len(X_pred), dtype=np.float64)
    for i in range(n_estimators):
        replace = n_sub > len(X_unl)
        sel = rng.choice(len(X_unl), size=n_sub, replace=replace)
        X_iter = np.vstack([X_pos, X_unl[sel]])
        y_iter = np.concatenate([np.ones(n_pos, dtype=np.int8), np.zeros(n_sub, dtype=np.int8)])
        clf = RandomForestClassifier(
            n_estimators=base_n_trees,
            max_features="sqrt",
            class_weight="balanced",
            n_jobs=n_jobs,
            random_state=int(seed) + i,
        )
        clf.fit(X_iter, y_iter)
        accum += clf.predict_proba(X_pred)[:, 1]
        if progress is not None:
            progress(i + 1, n_estimators)
    return accum / n_estimators


# ======================================================================
# 5. Enrichment
# ======================================================================

def benjamini_hochberg(pvals) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN entries stay NaN and are not counted)."""
    p = np.asarray(pvals, dtype=float)
    adj = np.full(p.shape, np.nan)
    valid = ~np.isnan(p)
    pv = p[valid]
    m = len(pv)
    if m == 0:
        return adj
    order = np.argsort(pv, kind="mergesort")
    scaled = pv[order] * m / np.arange(1, m + 1)
    scaled = np.minimum.accumulate(scaled[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(scaled, 0.0, 1.0)
    adj[valid] = out
    return adj


def read_gmt(path) -> Dict[str, List[str]]:
    """Read a GMT gene-set file into {term: [GENE, ...]} (symbols upper-cased)."""
    sets: Dict[str, List[str]] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets[parts[0]] = [g.strip().upper() for g in parts[2:] if g.strip()]
    return sets


def hypergeom_enrichment(
    query: Iterable[str],
    gene_sets: Mapping[str, Iterable[str]],
    background: Iterable[str],
    min_term_size: int = 1,
) -> pd.DataFrame:
    """
    Over-representation analysis with an EXPLICIT background and BH correction.

    * background  : the genes actually analysed (the analysis universe), NOT all
                    genes of the gene-set library.  Gene-set members outside the
                    background are ignored.
    * query       : genes of interest; only those inside the background are used.
    * Test        : one-sided hypergeometric P(X >= n_overlap).
    * Correction  : Benjamini-Hochberg across ALL tested terms.  No filtering
                    on overlap size or on keywords is applied.

    Returns one row per tested term with the columns: term, n_query,
    n_term_in_background, n_term_total, n_background, n_overlap, fold_enrichment,
    p_raw, p_adj, genes.  Sorted by p_adj then p_raw.
    """
    from scipy.stats import hypergeom

    bg = {str(g).strip().upper() for g in background if str(g).strip()}
    q_all = {str(g).strip().upper() for g in query if str(g).strip()}
    q = q_all & bg
    N = len(bg)
    n = len(q)
    if q_all - bg:
        log.warning(f"  {len(q_all - bg)} query genes are not in the background and were ignored.")

    rows = []
    for term, genes in gene_sets.items():
        members = {str(g).strip().upper() for g in genes}
        in_bg = members & bg
        K = len(in_bg)
        if K < max(1, int(min_term_size)):
            continue
        overlap = q & in_bg
        k = len(overlap)
        p = float(hypergeom.sf(k - 1, N, K, n)) if N > 0 else 1.0
        fold = (k / n) / (K / N) if n > 0 and N > 0 else np.nan
        rows.append({
            "term": term,
            "n_query": n,
            "n_term_in_background": K,
            "n_term_total": len(members),
            "n_background": N,
            "n_overlap": k,
            "fold_enrichment": fold,
            "p_raw": p,
            "genes": ";".join(sorted(overlap)),
        })
    cols = ["term", "n_query", "n_term_in_background", "n_term_total", "n_background",
            "n_overlap", "fold_enrichment", "p_raw", "p_adj", "genes"]
    if not rows:
        return pd.DataFrame(columns=cols)
    df = pd.DataFrame(rows)
    df["p_adj"] = benjamini_hochberg(df["p_raw"].to_numpy())
    return df[cols].sort_values(["p_adj", "p_raw"], kind="mergesort").reset_index(drop=True)


# ======================================================================
# 6. Network summaries (sparse, unweighted)
# ======================================================================

def adjacency_from_edges(idx_a, idx_b, n_nodes: int):
    """Symmetric, binary, loop-free scipy CSR adjacency matrix from integer edge endpoints."""
    from scipy import sparse

    idx_a = np.asarray(idx_a, dtype=np.int64)
    idx_b = np.asarray(idx_b, dtype=np.int64)
    keep = idx_a != idx_b                                  # drop self-loops
    idx_a, idx_b = idx_a[keep], idx_b[keep]
    A = sparse.coo_matrix((np.ones(len(idx_a), dtype=np.float32), (idx_a, idx_b)),
                          shape=(n_nodes, n_nodes)).tocsr()
    A = A + A.T                                            # symmetrise (duplicates are summed)
    A.data[:] = 1.0                                        # binary adjacency
    return A


def local_clustering(A, degree: np.ndarray, chunk: int = 500) -> np.ndarray:
    """Unweighted local clustering coefficient per node (0 for nodes of degree < 2)."""
    n = A.shape[0]
    closed_walks = np.zeros(n, dtype=np.float64)
    for start in range(0, n, chunk):
        end = min(start + chunk, n)
        sub = A[start:end]
        closed_walks[start:end] = np.asarray((sub @ A).multiply(sub).sum(axis=1)).ravel()
    denom = degree * (degree - 1.0)
    out = np.zeros(n, dtype=np.float64)
    ok = denom > 0
    out[ok] = closed_walks[ok] / denom[ok]      # closed walks of length 3 = 2 * triangles
    return out


def network_summary(A) -> Dict[str, object]:
    """
    Summary statistics of an undirected unweighted graph given as a CSR adjacency.

    Keys: n_nodes, n_edges, density, mean_degree, largest_component,
    clustering_coefficient (mean over all nodes), isolated_nodes, degree (ndarray).
    """
    from scipy.sparse.csgraph import connected_components

    n = A.shape[0]
    degree = np.asarray(A.sum(axis=1)).ravel().astype(np.float64)
    m = int(round(degree.sum() / 2.0))
    _, labels = connected_components(A, directed=False)
    largest = int(np.bincount(labels).max()) if n > 0 else 0
    clust = local_clustering(A, degree)
    return {
        "n_nodes": int(n),
        "n_edges": m,
        "density": (2.0 * m / (n * (n - 1))) if n > 1 else 0.0,
        "mean_degree": float(degree.mean()) if n > 0 else 0.0,
        "largest_component": largest,
        "clustering_coefficient": float(clust.mean()) if n > 0 else 0.0,
        "isolated_nodes": int((degree == 0).sum()),
        "degree": degree,
    }


# ======================================================================
# 7. Sample-level QC verdicts
# ======================================================================

def evaluate_qc_warnings(
    median_log2fc: float,
    fraction_de: float,
    n_up: int,
    n_down: int,
    max_abs_median_log2fc: float = 1.0,
    max_fraction_de: float = 0.80,
    max_up_down_ratio: float = 10.0,
) -> List[str]:
    """
    Return a list of warning strings describing global signs of a cohort / normalisation
    artefact in a tumor-vs-normal comparison (empty list = no warning):

      * |median log2FC over all genes| > max_abs_median_log2fc
      * fraction of DE-significant genes > max_fraction_de
      * up:down ratio > max_up_down_ratio or < 1 / max_up_down_ratio
    """
    warnings_out: List[str] = []
    if np.isfinite(median_log2fc) and abs(median_log2fc) > max_abs_median_log2fc:
        warnings_out.append(
            f"|median log2FC| = {abs(median_log2fc):.2f} > {max_abs_median_log2fc:g}: the typical gene "
            f"is shifted between groups (global offset, typical of cohort/normalisation effects)."
        )
    if np.isfinite(fraction_de) and fraction_de > max_fraction_de:
        warnings_out.append(
            f"{100 * fraction_de:.1f}% of genes are DE-significant (> {100 * max_fraction_de:.0f}%): "
            f"implausible for a purely biological contrast."
        )
    if n_up + n_down > 0:
        if n_down == 0:
            ratio = float("inf")
        else:
            ratio = n_up / n_down
        if ratio > max_up_down_ratio or ratio < 1.0 / max_up_down_ratio:
            ratio_txt = "inf" if np.isinf(ratio) else f"{ratio:.3g}"
            warnings_out.append(
                f"up:down ratio = {ratio_txt} (up={n_up}, down={n_down}) is outside "
                f"[{1.0 / max_up_down_ratio:.3g}, {max_up_down_ratio:g}]: extreme asymmetry."
            )
    return warnings_out


# ======================================================================
# 8. Patient pairing (paired tumor / adjacent-normal designs)
# ======================================================================

_TCGA_PATIENT_RE = re.compile(r"^(TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4})")


def infer_patient_ids(sample_ids: Sequence[str], patient_map: Optional[Mapping[str, str]] = None) -> pd.Series:
    """
    Map sample ids to patient ids: use ``patient_map`` when it knows the sample, otherwise
    the first 12 characters of a TCGA barcode (TCGA-XX-XXXX); unknown samples get NaN.
    """
    out = {}
    for s in sample_ids:
        if patient_map is not None and s in patient_map and pd.notna(patient_map[s]):
            out[s] = str(patient_map[s])
            continue
        m = _TCGA_PATIENT_RE.match(str(s))
        out[s] = m.group(1) if m else np.nan
    return pd.Series(out, dtype=object)


def pair_samples_by_patient(tumor_patients: pd.Series, normal_patients: pd.Series) -> Tuple[List[str], List[str]]:
    """
    Pair tumor and normal samples by patient.

    Both inputs are Series indexed by sample id with the patient id as value.  Patients with
    several samples in a group contribute their first sample (sorted by sample id).  Returns
    two equally long, patient-aligned lists of sample ids (tumor, normal), sorted by patient.
    """
    t = tumor_patients.dropna().sort_index()
    n = normal_patients.dropna().sort_index()
    t_first = t[~t.duplicated(keep="first")]
    n_first = n[~n.duplicated(keep="first")]
    t_by_patient = pd.Series(t_first.index, index=t_first.values)
    n_by_patient = pd.Series(n_first.index, index=n_first.values)
    common = sorted(set(t_by_patient.index) & set(n_by_patient.index))
    return [t_by_patient[p] for p in common], [n_by_patient[p] for p in common]


# ======================================================================
# 9. Feature groups for the ablation study
# ======================================================================

def resolve_feature_groups(
    columns: Sequence[str],
    expression: Iterable[str],
    network: Iterable[str],
    de_stats: Iterable[str],
) -> Dict[str, List[str]]:
    """
    Map the feature columns that actually exist onto the ablation feature sets.

    Returns a dict with the keys
      expression_only            : expression-derived features (incl. DE statistics)
      network_only               : network-derived features
      expression_plus_network    : union of both, in the original column order
      combined_without_de_stats  : expression_plus_network minus the explicit DE / log2FC
                                   statistics (``de_stats``)
      full                       : every column (also columns not assigned to any group)
      unassigned                 : columns belonging to no group (informational)
    Column order follows ``columns``.  Only existing columns are returned.
    """
    cols = list(columns)
    expr_set, net_set, de_set = set(expression), set(network), set(de_stats)
    expr = [c for c in cols if c in expr_set]
    net = [c for c in cols if c in net_set and c not in expr_set]
    both = [c for c in cols if c in expr_set or c in net_set]
    return {
        "expression_only": expr,
        "network_only": net,
        "expression_plus_network": both,
        "combined_without_de_stats": [c for c in both if c not in de_set],
        "full": cols,
        "unassigned": [c for c in cols if c not in expr_set and c not in net_set],
    }


# ======================================================================
# 10. Candidate evidence classification (step20)
# ======================================================================

EVIDENCE_CLASSES = (
    "already_LUAD",
    "lung_cancer_unspecified",
    "other_cancer",
    "indirect_only",
    "no_association_found",
    "not_assessed",
)


def classify_candidate_evidence(
    ot_luad_score,
    ot_lung_score,
    ot_cancer_score,
    epmc_luad_hits,
    epmc_lung_hits,
    epmc_cancer_hits,
    min_ot_score: float = 0.05,
    min_hits: int = 3,
) -> str:
    """
    Classify a non-LCGene candidate from external evidence (None / NaN = not retrieved).

      already_LUAD             Open Targets LUAD score >= min_ot_score or Europe PMC LUAD co-mentions >= min_hits
      lung_cancer_unspecified  (not above) lung-carcinoma score / 'lung cancer' co-mentions pass the thresholds
      other_cancer             (not above) generic cancer score / 'cancer' co-mentions pass the thresholds
      indirect_only            some evidence exists (score > 0 or >= 1 hit) but below every threshold
      no_association_found     BOTH sources were queried and returned nothing
      not_assessed             a source could not be queried and no positive evidence was found
                               (never interpreted as novelty)

    Only ``no_association_found`` may be described as "potentially novel" — and only with respect
    to the sources searched.
    """
    def _num(x):
        try:
            if x is None or (isinstance(x, float) and np.isnan(x)):
                return None
            return float(x)
        except (TypeError, ValueError):
            return None

    s = [_num(ot_luad_score), _num(ot_lung_score), _num(ot_cancer_score)]
    h = [_num(epmc_luad_hits), _num(epmc_lung_hits), _num(epmc_cancer_hits)]
    if all(v is None for v in s + h):
        return "not_assessed"

    def _ge(vals, thr):
        return any(v is not None and v >= thr for v in vals)

    if _ge([s[0]], min_ot_score) or _ge([h[0]], min_hits):
        return "already_LUAD"
    if _ge([s[1]], min_ot_score) or _ge([h[1]], min_hits):
        return "lung_cancer_unspecified"
    if _ge([s[2]], min_ot_score) or _ge([h[2]], min_hits):
        return "other_cancer"
    # A negative finding is only meaningful if BOTH sources were queried successfully.
    if all(v is None for v in s) or all(v is None for v in h):
        return "not_assessed"
    if any(v is not None and v > 0 for v in s + h):
        return "indirect_only"
    return "no_association_found"
