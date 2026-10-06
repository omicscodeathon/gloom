"""
step4_differential_expression.py
---------------------------------
Differential expression analysis between tumor and normal samples.
Outputs: differential_expression_results.csv + volcano/heatmap plots.

Test selected by config.DE_METHOD:
  "welch"      Welch t-test on log2(x+1) values + BH FDR + Cohen's d  (default, unpaired).
  "paired"     Paired t-test on patients that have BOTH a tumor and an adjacent-normal sample
               (TCGA-LUAD tumor vs adjacent normal design).  Needs a patient_id per sample
               (config.PAIRING_FILE, the 'patient_id' column of the processed sample metadata,
               or a TCGA barcode).  log2FC = mean paired difference; Cohen's d = d_z.
               Falls back to Welch when fewer than 3 pairs are found.
  "limma_voom" limma-voom on raw counts (config.GDC_COUNTS_FILE) through the OPTIONAL rpy2
               bridge to R (packages limma + edgeR).  Uses a patient blocking factor when
               pairs are available.  rpy2 is NOT a hard dependency: when rpy2/R/limma or the
               counts matrix is missing, the method is skipped with a clear message and the
               Welch test is used instead.

The method actually used is stored in the ``de_method`` column of the result table.

FIX: sort_values key lambda rewritten as an explicit key function that
     is safe across all pandas versions (>=1.1).  The original inline
     lambda worked in newer pandas but was fragile on older releases.
"""
import logging
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
from gloom_utils import infer_patient_ids, pair_samples_by_patient

config.create_output_dirs()
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

MIN_PAIRS = 3


# ----------------------------------------------------------------------
# Welch (unpaired)
# ----------------------------------------------------------------------

def _de_welch(tumor_expr, normal_expr):
    mean_tumor  = tumor_expr.mean(axis=1)
    mean_normal = normal_expr.mean(axis=1)
    log2fc      = mean_tumor - mean_normal

    log.info(f"  Running Welch t-test for {tumor_expr.shape[0]} genes …")
    t_stats, pvalues = stats.ttest_ind(
        tumor_expr.values,
        normal_expr.values,
        axis=1,
        equal_var=False,
        nan_policy="omit",
    )

    n1, n2 = tumor_expr.shape[1], normal_expr.shape[1]
    std1   = tumor_expr.std(axis=1, ddof=1)
    std2   = normal_expr.std(axis=1, ddof=1)
    pooled = np.sqrt(((n1 - 1) * std1 ** 2 + (n2 - 1) * std2 ** 2) / (n1 + n2 - 2))
    cohens_d = np.where(pooled > 0, (mean_tumor - mean_normal) / pooled, 0.0)
    return {
        "mean_tumor": mean_tumor, "mean_normal": mean_normal, "log2fc": log2fc,
        "t_stat": np.asarray(t_stats), "pvalue": np.asarray(pvalues),
        "cohens_d": np.asarray(cohens_d), "method": "welch",
    }


# ----------------------------------------------------------------------
# Pairing helpers
# ----------------------------------------------------------------------

def _patient_map():
    """sample_id -> patient_id from PAIRING_FILE or the processed metadata files (if available)."""
    mapping = {}
    pairing_file = getattr(config, "PAIRING_FILE", None)
    if pairing_file and Path(pairing_file).exists():
        df = pd.read_csv(pairing_file)
        if {"sample_id", "patient_id"} <= set(df.columns):
            mapping.update(dict(zip(df["sample_id"].astype(str), df["patient_id"].astype(str))))
        else:
            log.warning(f"  PAIRING_FILE {pairing_file} lacks sample_id/patient_id columns — ignored.")
    for fname in ("tumor_metadata_processed.csv", "normal_metadata_processed.csv"):
        p = config.PROCESSED_DIR / fname
        if p.exists():
            meta = pd.read_csv(p, index_col=0)
            if "patient_id" in meta.columns:
                for s, pid in meta["patient_id"].items():
                    mapping.setdefault(str(s), str(pid))
    return mapping


def _find_pairs(tumor_expr, normal_expr):
    pmap = _patient_map()
    t_pat = infer_patient_ids(list(tumor_expr.columns), pmap)
    n_pat = infer_patient_ids(list(normal_expr.columns), pmap)
    t_ids, n_ids = pair_samples_by_patient(t_pat, n_pat)
    return t_ids, n_ids


# ----------------------------------------------------------------------
# Paired t-test
# ----------------------------------------------------------------------

def _de_paired(tumor_expr, normal_expr):
    t_ids, n_ids = _find_pairs(tumor_expr, normal_expr)
    if len(t_ids) < MIN_PAIRS:
        log.warning(f"  DE_METHOD='paired' but only {len(t_ids)} tumor/normal pair(s) were found "
                    f"(need >= {MIN_PAIRS}). Provide patient_id per sample (config.PAIRING_FILE or "
                    f"TCGA barcodes). Falling back to the Welch test.")
        return None
    log.info(f"  Paired t-test on {len(t_ids)} patients with both tumor and normal samples …")
    T = tumor_expr[t_ids].to_numpy(dtype=float)
    N = normal_expr[n_ids].to_numpy(dtype=float)
    diff = T - N
    t_stats, pvalues = stats.ttest_rel(T, N, axis=1, nan_policy="omit")
    mean_diff = np.nanmean(diff, axis=1)
    sd_diff = np.nanstd(diff, axis=1, ddof=1)
    dz = np.where(sd_diff > 0, mean_diff / sd_diff, 0.0)
    idx = tumor_expr.index
    return {
        "mean_tumor": pd.Series(np.nanmean(T, axis=1), index=idx),
        "mean_normal": pd.Series(np.nanmean(N, axis=1), index=idx),
        "log2fc": pd.Series(mean_diff, index=idx),
        "t_stat": np.asarray(t_stats), "pvalue": np.asarray(pvalues),
        "cohens_d": np.asarray(dz), "method": "paired",
    }


# ----------------------------------------------------------------------
# limma-voom (optional, via rpy2)
# ----------------------------------------------------------------------

_LIMMA_VOOM_R = r"""
function(counts_csv, group_csv, out_csv) {
  suppressPackageStartupMessages({ library(limma); library(edgeR) })
  counts <- as.matrix(read.csv(counts_csv, row.names = 1, check.names = FALSE))
  info   <- read.csv(group_csv, stringsAsFactors = FALSE)
  info   <- info[match(colnames(counts), info$sample_id), ]
  group  <- factor(info$group, levels = c("normal", "tumor"))
  dge    <- DGEList(counts = counts)
  keep   <- filterByExpr(dge, group = group)
  dge    <- dge[keep, , keep.lib.sizes = FALSE]
  dge    <- calcNormFactors(dge)
  if ("patient_id" %in% colnames(info) && !any(is.na(info$patient_id))) {
    patient <- factor(info$patient_id)
    design  <- model.matrix(~ patient + group)
  } else {
    design  <- model.matrix(~ group)
  }
  v   <- voom(dge, design)
  fit <- eBayes(lmFit(v, design))
  tt  <- topTable(fit, coef = "grouptumor", number = Inf, sort.by = "none")
  write.csv(tt, out_csv)
  invisible(NULL)
}
"""


def _de_limma_voom(tumor_expr, normal_expr):
    """Return a result dict, or None (with a clear log message) when limma-voom cannot run."""
    try:
        import rpy2.robjects as ro  # optional dependency
        from rpy2.robjects.packages import importr
        importr("limma")
        importr("edgeR")
    except Exception as exc:
        log.warning(f"  DE_METHOD='limma_voom' SKIPPED: rpy2 / R packages limma + edgeR are not "
                    f"available ({type(exc).__name__}: {exc}). Install with: pip install rpy2 and "
                    f"BiocManager::install(c('limma','edgeR')). Falling back to the Welch test.")
        return None

    counts_path = Path(getattr(config, "GDC_COUNTS_FILE", ""))
    if not counts_path.exists():
        log.warning(f"  DE_METHOD='limma_voom' SKIPPED: raw counts matrix not found ({counts_path}). "
                    f"limma-voom needs raw counts (run scripts/fetch_gdc_tcga_luad.py). "
                    f"Falling back to the Welch test.")
        return None

    counts = pd.read_csv(counts_path, index_col=0)
    counts.index = counts.index.astype(str).str.strip().str.upper()
    if counts.index.duplicated().any():
        counts["_mean"] = counts.mean(axis=1)
        counts = counts.sort_values("_mean", ascending=False)
        counts = counts[~counts.index.duplicated(keep="first")].drop(columns=["_mean"])
    samples_t = [s for s in tumor_expr.columns if s in counts.columns]
    samples_n = [s for s in normal_expr.columns if s in counts.columns]
    genes = tumor_expr.index.intersection(counts.index)
    if len(samples_t) < 3 or len(samples_n) < 3 or len(genes) == 0:
        log.warning("  DE_METHOD='limma_voom' SKIPPED: counts matrix does not match the processed "
                    "samples/genes. Falling back to the Welch test.")
        return None

    # Patient blocking factor when pairs exist; samples are then restricted to paired patients.
    t_ids, n_ids = _find_pairs(tumor_expr[samples_t], normal_expr[samples_n])
    paired = len(t_ids) >= MIN_PAIRS
    if paired:
        samples_t, samples_n = t_ids, n_ids
        pmap = _patient_map()
        pats = infer_patient_ids(samples_t + samples_n, pmap)
        log.info(f"  limma-voom with patient blocking on {len(t_ids)} pairs …")
    else:
        pats = pd.Series(np.nan, index=samples_t + samples_n, dtype=object)
        log.info("  limma-voom (unpaired design ~ group) …")

    cols = samples_t + samples_n
    info = pd.DataFrame({
        "sample_id": cols,
        "group": ["tumor"] * len(samples_t) + ["normal"] * len(samples_n),
        "patient_id": [pats.get(s, np.nan) for s in cols] if paired else np.nan,
    })
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        counts.loc[genes, cols].to_csv(tmp / "counts.csv")
        info.to_csv(tmp / "info.csv", index=False)
        try:
            fn = ro.r(_LIMMA_VOOM_R)
            fn(str(tmp / "counts.csv"), str(tmp / "info.csv"), str(tmp / "tt.csv"))
            tt = pd.read_csv(tmp / "tt.csv", index_col=0)
        except Exception as exc:
            log.warning(f"  DE_METHOD='limma_voom' FAILED in R ({type(exc).__name__}: {exc}). "
                        f"Falling back to the Welch test.")
            return None

    idx = tumor_expr.index
    welch = _de_welch(tumor_expr, normal_expr)           # fills genes removed by filterByExpr
    log2fc = welch["log2fc"].copy()
    t_stat = pd.Series(welch["t_stat"], index=idx)
    pval = pd.Series(welch["pvalue"], index=idx)
    # genes dropped by filterByExpr are not tested by limma: no evidence -> p = 1, t = 0
    untested = ~idx.isin(tt.index)
    pval[untested] = 1.0
    t_stat[untested] = 0.0
    hit = idx[idx.isin(tt.index)]
    log2fc.loc[hit] = tt.loc[hit, "logFC"].to_numpy()
    t_stat.loc[hit] = tt.loc[hit, "t"].to_numpy()
    pval.loc[hit] = tt.loc[hit, "P.Value"].to_numpy()
    return {
        "mean_tumor": welch["mean_tumor"], "mean_normal": welch["mean_normal"], "log2fc": log2fc,
        "t_stat": t_stat.to_numpy(), "pvalue": pval.to_numpy(),
        "cohens_d": welch["cohens_d"], "method": "limma_voom",
    }


def run_differential_expression():
    log.info("=" * 60)
    log.info("STEP 4 — DIFFERENTIAL EXPRESSION")
    log.info("=" * 60)

    tumor_expr  = pd.read_csv(config.TUMOR_EXPR_HARMONIZED,  index_col=0)
    normal_expr = pd.read_csv(config.NORMAL_EXPR_HARMONIZED, index_col=0)
    log.info(f"  Tumor: {tumor_expr.shape}  Normal: {normal_expr.shape}")

    # ------------------------------------------------------------------
    # Test (welch | paired | limma_voom) -> means, log2FC, statistics, p-values
    # ------------------------------------------------------------------
    requested = str(getattr(config, "DE_METHOD", "welch")).lower()
    if requested not in {"welch", "paired", "limma_voom"}:
        log.warning(f"  Unknown DE_METHOD '{requested}' — using 'welch'.")
        requested = "welch"
    log.info(f"  DE_METHOD requested: {requested}")

    res = None
    if requested == "paired":
        res = _de_paired(tumor_expr, normal_expr)
    elif requested == "limma_voom":
        res = _de_limma_voom(tumor_expr, normal_expr)
    if res is None:
        res = _de_welch(tumor_expr, normal_expr)
    method_used = res["method"]
    log.info(f"  DE method used: {method_used}")

    mean_tumor, mean_normal, log2fc = res["mean_tumor"], res["mean_normal"], res["log2fc"]
    t_stats, pvalues, cohens_d = res["t_stat"], res["pvalue"], res["cohens_d"]
    pvalues = np.where(np.isnan(pvalues), 1.0, pvalues)
    t_stats = np.where(np.isnan(t_stats), 0.0, t_stats)

    # ------------------------------------------------------------------
    # Benjamini-Hochberg FDR correction
    # ------------------------------------------------------------------
    _, pvalues_adj, _, _ = multipletests(
        pvalues, alpha=config.DE_PVALUE_THRESHOLD, method="fdr_bh"
    )
    n_sig = (pvalues_adj <= config.DE_PVALUE_THRESHOLD).sum()
    log.info(f"  Significant at FDR<={config.DE_PVALUE_THRESHOLD}: {n_sig}")

    # ------------------------------------------------------------------
    # Build result DataFrame
    # ------------------------------------------------------------------
    de_df = pd.DataFrame(
        {
            "mean_tumor":      mean_tumor,
            "mean_normal":     mean_normal,
            "log2fc":          log2fc,
            "t_stat":          t_stats,
            "pvalue":          pvalues,
            "pvalue_adj":      pvalues_adj,
            "neg_log10_padj":  -np.log10(np.clip(pvalues_adj, 1e-300, 1.0)),
            "cohens_d":        cohens_d,
        },
        index=tumor_expr.index,
    )
    de_df["de_method"] = method_used

    # ------------------------------------------------------------------
    # Significance labels
    # ------------------------------------------------------------------
    sig_mask = (
        (de_df["pvalue_adj"] <= config.DE_PVALUE_THRESHOLD)
        & (de_df["log2fc"].abs() >= config.DE_LOG2FC_THRESHOLD)
    )
    de_df["significant"] = sig_mask
    de_df["direction"]   = "ns"
    de_df.loc[sig_mask & (de_df["log2fc"] >= config.DE_LOG2FC_THRESHOLD),  "direction"] = "up"
    de_df.loc[sig_mask & (de_df["log2fc"] <= -config.DE_LOG2FC_THRESHOLD), "direction"] = "down"

    # FIX: sort explicitly — first by adjusted p-value ascending,
    #      then by absolute log2FC descending — without relying on a
    #      multi-column key lambda that behaves differently across pandas
    #      versions.
    de_df["_abs_log2fc"] = de_df["log2fc"].abs()
    de_df = de_df.sort_values(
        ["pvalue_adj", "_abs_log2fc"],
        ascending=[True, False],
    ).drop(columns=["_abs_log2fc"])

    # ------------------------------------------------------------------
    # Volcano plot
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 6))
    color_map = {"up": "steelblue", "down": "tomato", "ns": "lightgrey"}
    for direction, group in de_df.groupby("direction"):
        ax.scatter(
            group["log2fc"],
            group["neg_log10_padj"],
            c=color_map[direction],
            s=8 if direction != "ns" else 4,
            alpha=0.7,
            label=direction,
            linewidths=0,
        )
    ax.axvline( config.DE_LOG2FC_THRESHOLD, color="black", linestyle="--", lw=0.8, alpha=0.6)
    ax.axvline(-config.DE_LOG2FC_THRESHOLD, color="black", linestyle="--", lw=0.8, alpha=0.6)
    ax.axhline(-np.log10(config.DE_PVALUE_THRESHOLD), color="black", linestyle=":", lw=0.8, alpha=0.6)
    ax.set_xlabel("Log2 Fold-Change")
    ax.set_ylabel("-log10(adj P)")
    ax.set_title(f"Volcano Plot — LUAD Tumor vs Normal ({method_used})")
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(config.FIGURES_DIR / "de_volcano_plot.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------------------------------------------
    # Log2FC distribution plot
    # ------------------------------------------------------------------
    fig2, ax2 = plt.subplots(figsize=(8, 4))
    ax2.hist(de_df["log2fc"], bins=100, color="steelblue", alpha=0.75, edgecolor="none")
    ax2.axvline( config.DE_LOG2FC_THRESHOLD, color="tomato",     linestyle="--", lw=1.2)
    ax2.axvline(-config.DE_LOG2FC_THRESHOLD, color="darkorange", linestyle="--", lw=1.2)
    ax2.set_xlabel("Log2 Fold-Change")
    ax2.set_ylabel("Number of Genes")
    ax2.set_title("Distribution of Log2 Fold-Change Values")
    plt.tight_layout()
    fig2.savefig(config.FIGURES_DIR / "de_log2fc_distribution.png", dpi=120, bbox_inches="tight")
    plt.close(fig2)

    # ------------------------------------------------------------------
    # Save and report
    # ------------------------------------------------------------------
    de_df.to_csv(config.DE_RESULTS_FILE)
    n_up   = (de_df["direction"] == "up").sum()
    n_down = (de_df["direction"] == "down").sum()
    log.info(f"  Up: {n_up}  Down: {n_down}  NS: {len(de_df) - n_up - n_down}")
    log.info("STEP 4 COMPLETE")
    return de_df


if __name__ == "__main__":
    de = run_differential_expression()
    print(
        de.head(10)[
            ["mean_tumor", "mean_normal", "log2fc", "pvalue_adj", "direction"]
        ].round(4)
    )
