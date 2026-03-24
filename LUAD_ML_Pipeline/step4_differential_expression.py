"""
step4_differential_expression.py
---------------------------------
Differential expression analysis between tumor and normal samples.
Uses Welch t-test + BH FDR correction + Cohen's d effect size.
Outputs: differential_expression_results.csv + volcano/heatmap plots.
"""
import logging, sys
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

config.create_output_dirs()
logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

def run_differential_expression():
    log.info("="*60); log.info("STEP 4 — DIFFERENTIAL EXPRESSION"); log.info("="*60)
    tumor_expr  = pd.read_csv(config.TUMOR_EXPR_HARMONIZED,  index_col=0)
    normal_expr = pd.read_csv(config.NORMAL_EXPR_HARMONIZED, index_col=0)
    log.info(f"  Tumor: {tumor_expr.shape}  Normal: {normal_expr.shape}")

    # Means and log2FC
    mean_tumor  = tumor_expr.mean(axis=1)
    mean_normal = normal_expr.mean(axis=1)
    log2fc      = mean_tumor - mean_normal

    # Welch t-test
    log.info(f"  Running Welch t-test for {tumor_expr.shape[0]} genes …")
    t_stats, pvalues = stats.ttest_ind(
        tumor_expr.values, normal_expr.values,
        axis=1, equal_var=False, nan_policy="omit")
    pvalues = np.where(np.isnan(pvalues), 1.0, pvalues)
    t_stats = np.where(np.isnan(t_stats), 0.0, t_stats)

    # BH correction
    _, pvalues_adj, _, _ = multipletests(pvalues, alpha=config.DE_PVALUE_THRESHOLD, method="fdr_bh")
    n_sig = (pvalues_adj <= config.DE_PVALUE_THRESHOLD).sum()
    log.info(f"  Significant at FDR<={config.DE_PVALUE_THRESHOLD}: {n_sig}")

    # Cohen's d
    n1, n2 = tumor_expr.shape[1], normal_expr.shape[1]
    std1 = tumor_expr.std(axis=1, ddof=1); std2 = normal_expr.std(axis=1, ddof=1)
    pooled = np.sqrt(((n1-1)*std1**2 + (n2-1)*std2**2) / (n1+n2-2))
    cohens_d = np.where(pooled > 0, (mean_tumor - mean_normal) / pooled, 0.0)

    # Build result dataframe
    de_df = pd.DataFrame({
        "mean_tumor": mean_tumor, "mean_normal": mean_normal, "log2fc": log2fc,
        "t_stat": t_stats, "pvalue": pvalues, "pvalue_adj": pvalues_adj,
        "neg_log10_padj": -np.log10(np.clip(pvalues_adj, 1e-300, 1.0)),
        "cohens_d": cohens_d,
    }, index=tumor_expr.index)

    # Significance
    sig_mask = (de_df["pvalue_adj"] <= config.DE_PVALUE_THRESHOLD) &                (de_df["log2fc"].abs() >= config.DE_LOG2FC_THRESHOLD)
    de_df["significant"] = sig_mask
    de_df["direction"] = "ns"
    de_df.loc[sig_mask & (de_df["log2fc"] >= config.DE_LOG2FC_THRESHOLD), "direction"] = "up"
    de_df.loc[sig_mask & (de_df["log2fc"] <= -config.DE_LOG2FC_THRESHOLD), "direction"] = "down"
    de_df = de_df.sort_values(["pvalue_adj", "log2fc"], ascending=[True, False],
                               key=lambda c: c.abs() if c.name == "log2fc" else c)

    # Volcano plot
    fig, ax = plt.subplots(figsize=(9,6))
    color_map = {"up":"steelblue","down":"tomato","ns":"lightgrey"}
    for direction, group in de_df.groupby("direction"):
        ax.scatter(group["log2fc"], group["neg_log10_padj"],
                   c=color_map[direction], s=8 if direction!="ns" else 4,
                   alpha=0.7, label=direction, linewidths=0)
    ax.axvline(config.DE_LOG2FC_THRESHOLD,  color="black", linestyle="--", lw=0.8, alpha=0.6)
    ax.axvline(-config.DE_LOG2FC_THRESHOLD, color="black", linestyle="--", lw=0.8, alpha=0.6)
    ax.axhline(-np.log10(config.DE_PVALUE_THRESHOLD), color="black", linestyle=":", lw=0.8, alpha=0.6)
    ax.set_xlabel("Log2 Fold-Change"); ax.set_ylabel("-log10(adj P)")
    ax.set_title("Volcano Plot — LUAD Tumor vs GTEx Normal"); ax.legend(fontsize=9)
    plt.tight_layout(); fig.savefig(config.FIGURES_DIR/"de_volcano_plot.png", dpi=150, bbox_inches="tight"); plt.close(fig)

    # Log2FC distribution
    fig2, ax2 = plt.subplots(figsize=(8,4))
    ax2.hist(de_df["log2fc"], bins=100, color="steelblue", alpha=0.75, edgecolor="none")
    ax2.axvline(config.DE_LOG2FC_THRESHOLD, color="tomato", linestyle="--", lw=1.2)
    ax2.axvline(-config.DE_LOG2FC_THRESHOLD, color="darkorange", linestyle="--", lw=1.2)
    ax2.set_xlabel("Log2 Fold-Change"); ax2.set_ylabel("Number of Genes")
    ax2.set_title("Distribution of Log2 Fold-Change Values")
    plt.tight_layout(); fig2.savefig(config.FIGURES_DIR/"de_log2fc_distribution.png", dpi=120, bbox_inches="tight"); plt.close(fig2)

    de_df.to_csv(config.DE_RESULTS_FILE)
    n_up = (de_df["direction"]=="up").sum(); n_down = (de_df["direction"]=="down").sum()
    log.info(f"  Up: {n_up}  Down: {n_down}  NS: {len(de_df)-n_up-n_down}")
    log.info("STEP 4 COMPLETE")
    return de_df

if __name__ == "__main__":
    de = run_differential_expression()
    print(de.head(10)[["mean_tumor","mean_normal","log2fc","pvalue_adj","direction"]].round(4))
