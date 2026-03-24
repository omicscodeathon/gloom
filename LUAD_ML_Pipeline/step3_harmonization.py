"""
step3_harmonization.py
----------------------
Harmonizes gene sets between tumor and normal expression matrices.
Steps: normalize symbols, find intersection, check CGC retention,
subset both matrices, validate alignment.
Run: python step3_harmonization.py
"""
import logging, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

config.create_output_dirs()
logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

def normalize_gene_symbols(df, label):
    n_before = len(df)
    df.index = df.index.astype(str).str.strip().str.upper()
    invalid_mask = df.index.isin(["NAN","","NA","NONE"]) | df.index.str.fullmatch(r"\d+")
    if invalid_mask.sum() > 0:
        log.warning(f"  [{label}] Dropping {invalid_mask.sum()} invalid symbols.")
        df = df.loc[~invalid_mask]
    log.info(f"  [{label}] Normalised: {n_before} -> {len(df)} genes.")
    return df

def find_common_genes(tumor_genes, normal_genes):
    common = set(tumor_genes) & set(normal_genes)
    log.info(f"  Common genes: {len(common)}")
    if len(common) == 0:
        raise ValueError("Zero genes overlap between tumor and normal.")
    return pd.Index(sorted(common), name="gene")

def check_cancer_gene_retention(common_genes, cancer_genes, label="harmonized"):
    cancer_set = set(cancer_genes.str.upper())
    retained   = cancer_set & set(common_genes)
    lost       = cancer_set - set(common_genes)
    pct        = 100 * len(retained) / len(cancer_set) if cancer_set else 0.0
    log.info(f"  CGC retained: {len(retained)} / {len(cancer_set)} ({pct:.1f}%)")
    return pd.DataFrame([{"stage": label, "total_cancer_genes": len(cancer_set),
                          "retained": len(retained), "lost": len(lost),
                          "pct_retained": round(pct, 2),
                          "lost_symbols": "|".join(sorted(lost))}])

def plot_venn(tumor_genes, normal_genes, common_genes, out_path):
    n_t = len(tumor_genes); n_n = len(normal_genes); n_c = len(common_genes)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_xlim(0,10); ax.set_ylim(0,6); ax.axis("off")
    ct = mpatches.Circle((3.5,3), 2.4, color="steelblue", alpha=0.4, label="Tumor")
    cn = mpatches.Circle((6.5,3), 2.4, color="darkorange", alpha=0.4, label="Normal")
    ax.add_patch(ct); ax.add_patch(cn)
    ax.text(2.2,3.0,f"{n_t-n_c:,}", ha="center", va="center", fontsize=14, fontweight="bold", color="steelblue")
    ax.text(7.8,3.0,f"{n_n-n_c:,}", ha="center", va="center", fontsize=14, fontweight="bold", color="darkorange")
    ax.text(5.0,3.0,f"{n_c:,}", ha="center", va="center", fontsize=14, fontweight="bold", color="black")
    ax.set_title(f"Gene Set Overlap — Tumor: {n_t:,}  Normal: {n_n:,}  Common: {n_c:,}", fontsize=11)
    ax.legend(handles=[ct,cn], loc="upper right", fontsize=9)
    plt.tight_layout(); fig.savefig(out_path, dpi=120, bbox_inches="tight"); plt.close(fig)

def run_harmonization():
    log.info("="*60); log.info("STEP 3 — GENE ID HARMONIZATION"); log.info("="*60)
    tumor_expr   = pd.read_csv(config.TUMOR_EXPR_PROCESSED,  index_col=0)
    normal_expr  = pd.read_csv(config.NORMAL_EXPR_PROCESSED, index_col=0)
    cancer_genes = pd.read_csv(config.PROCESSED_DIR/"cancer_genes_raw.csv", header=0).squeeze("columns")
    tumor_expr   = normalize_gene_symbols(tumor_expr,  "tumor")
    normal_expr  = normalize_gene_symbols(normal_expr, "normal")
    cancer_genes = cancer_genes.astype(str).str.strip().str.upper().drop_duplicates().reset_index(drop=True)
    cancer_genes.name = "gene"
    common_genes = find_common_genes(tumor_expr.index, normal_expr.index)
    cgc_report   = check_cancer_gene_retention(common_genes, cancer_genes)
    tumor_expr   = tumor_expr.loc[common_genes]
    normal_expr  = normal_expr.loc[common_genes]
    if not tumor_expr.index.equals(normal_expr.index):
        raise ValueError("Gene indexes are NOT aligned after harmonization.")
    plot_venn(tumor_expr.index, normal_expr.index, common_genes,
              config.FIGURES_DIR/"harmonization_venn.png")
    tumor_expr.to_csv(config.TUMOR_EXPR_HARMONIZED)
    normal_expr.to_csv(config.NORMAL_EXPR_HARMONIZED)
    pd.Series(common_genes, name="gene").to_csv(config.PROCESSED_DIR/"common_genes.csv", index=False)
    cgc_report.to_csv(config.PROCESSED_DIR/"cgc_retention_report.csv", index=False)
    log.info(f"  Common genes: {len(common_genes)}")
    log.info("STEP 3 COMPLETE")
    return {"tumor_expr": tumor_expr, "normal_expr": normal_expr,
            "common_genes": common_genes, "cancer_genes": cancer_genes}

if __name__ == "__main__":
    r = run_harmonization()
    print(f"\nCommon genes: {len(r['common_genes'])}")
    print(r['tumor_expr'].iloc[:3,:4].round(4))
