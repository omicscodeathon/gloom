"""
step9_label_construction.py - Label Construction
Assigns binary labels: 1=CGC cancer gene, 0=non-CGC.
Outputs: gene_labels.csv, gene_annotation_table.csv
"""
import logging, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

def run_label_construction():
    log.info("="*60); log.info("STEP 9 — LABEL CONSTRUCTION"); log.info("="*60)
    features  = pd.read_csv(config.INTEGRATED_FEATURES_FILE, index_col=0)
    de_df     = pd.read_csv(config.DE_RESULTS_FILE,           index_col=0)
    cgc_genes = pd.read_csv(config.PROCESSED_DIR/"cancer_genes_raw.csv", header=0).squeeze("columns")
    cgc_set   = set(cgc_genes.dropna().astype(str).str.strip().str.upper())
    log.info(f"  CGC symbols: {len(cgc_set)}")

    gene_index = features.index
    normalised = pd.Index(gene_index.astype(str)).str.strip().str.upper()

    labels = pd.Series(
        normalised.isin(cgc_set).astype(int),
        index=gene_index,
        name="label"
    )

    n_pos = labels.sum(); n_neg = (labels==0).sum(); n_total = len(labels)
    log.info(f"  Total: {n_total}  Positive: {n_pos} ({100*n_pos/n_total:.2f}%)  Negative: {n_neg}")

    # Annotation table
    annotation = pd.DataFrame(index=gene_index)
    annotation["label"]      = labels
    annotation["is_cgc_gene"] = normalised.isin(cgc_set)
    de_cols = ["log2fc","pvalue_adj","neg_log10_padj","cohens_d","direction","significant","mean_tumor","mean_normal"]
    for col in de_cols:
        if col in de_df.columns:
            annotation[col] = de_df[col].reindex(gene_index)
    annotation["is_de_significant"] = annotation.get("significant", pd.Series(False, index=gene_index)).fillna(False).astype(bool)
    annotation["is_cgc_and_de"] = annotation["is_cgc_gene"] & annotation["is_de_significant"]
    annotation["log2fc"]      = annotation["log2fc"].fillna(0.0)
    annotation["pvalue_adj"]  = annotation["pvalue_adj"].fillna(1.0)
    annotation["direction"]   = annotation["direction"].fillna("ns")
    annotation["significant"] = annotation["significant"].fillna(False)

    # Plots
    fig, axes = plt.subplots(1,2, figsize=(10,4))
    bars = axes[0].bar(["Negative","Positive"], [n_neg, n_pos],
                       color=["steelblue","tomato"], alpha=0.85, width=0.5)
    for bar, val in zip(bars, [n_neg, n_pos]):
        axes[0].text(bar.get_x()+bar.get_width()/2, bar.get_height()+max(n_neg,n_pos)*0.01,
                     f"{val:,}", ha="center", va="bottom", fontsize=10)
    axes[0].set_title("Class Counts"); axes[0].set_ylabel("Number of genes")
    axes[1].pie([n_neg, n_pos],
                labels=[f"Negative\n{n_neg:,}",f"Positive\n{n_pos:,}"],
                colors=["steelblue","tomato"], autopct=None, startangle=90, wedgeprops={"alpha":0.85})
    axes[1].set_title("Class Proportions")
    if n_pos > 0:
        ratio_text = f"{n_neg/n_pos:.1f}:1"
    else:
        ratio_text = "undefined (no positive genes)"

    fig.suptitle(f"Label Class Balance  Ratio: {ratio_text}", fontsize=11)
    plt.tight_layout(); fig.savefig(config.FIGURES_DIR/"label_class_balance.png", dpi=130, bbox_inches="tight"); plt.close(fig)

    label_df = labels.reset_index(); label_df.columns = ["gene","label"]
    label_df.to_csv(config.LABELS_FILE, index=False)
    annotation.to_csv(config.PROCESSED_DIR/"gene_annotation_table.csv")
    unmatched = sorted(cgc_set - set(normalised))
    pd.Series(unmatched, name="unmatched_cgc_gene").to_csv(config.PROCESSED_DIR/"unmatched_cgc_genes.csv", index=False)
    log.info(f"  CGC genes unmatched: {len(unmatched)}")
    log.info("STEP 9 COMPLETE")
    return {"labels": labels, "annotation": annotation, "cgc_set": cgc_set}

if __name__ == "__main__":
    r = run_label_construction()
    print(r["labels"].value_counts().rename({0:"Negative",1:"Positive"}))
