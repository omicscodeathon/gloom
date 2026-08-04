"""
config.py
---------
Central configuration for the LUAD ML Pipeline.
All paths, constants, and parameters are defined here.
No hardcoded values should appear in other modules.
"""

import os
from pathlib import Path

# ==================================================
# ROOT PATHS
# ==================================================

# Project root: four levels up from src/gloom/pipeline/config.py
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent

# Data root: data/ at the project root
DATA_ROOT = PROJECT_ROOT / "data"

# ==================================================
# RAW INPUT FILE PATHS
# ==================================================

RAW_DIR = DATA_ROOT / "raw"

# --- Tumor (cBioPortal LUAD) ---
TUMOR_EXPR_FILE   = RAW_DIR / "cBioPortal (RNA Seq Data)" / "data_mrna_seq_v2_rsem.txt"
TUMOR_META_FILE   = RAW_DIR / "cBioPortal (RNA Seq Data)" / "data_clinical_patient.txt"

# --- Normal (GTEx lung) ---
NORMAL_EXPR_FILE  = RAW_DIR / "Gtex (normal samples)" / "gene_tpm_v11_lung.gct.gz"
# Corrected: TSV file, not Excel; filename matches the actual provided file
NORMAL_META_FILE  = RAW_DIR / "Gtex (normal samples)" / \
                    "GTEx_Analysis_v11_Annotations_SampleAttributesDS - LUAD.txt"

# --- Labeled LUAD gene list (LCGene database) ---
# MODIFIED: replaced Cancer Gene Census CSV with LCGene_human_LUAD_filtered.tsv
CANCER_GENE_FILE  = RAW_DIR / "LCGene (Labeled LUAD Data)" / \
                    "LCGene_human_LUAD_filtered.tsv"

# ==================================================
# PROCESSED / INTERMEDIATE OUTPUT PATHS
# ==================================================

PROCESSED_DIR       = DATA_ROOT / "processed"

# All generated outputs live under trial/outputs/ — outside the pipeline code folder
OUTPUTS_ROOT = PROJECT_ROOT.parent / "outputs"
RESULTS_DIR  = OUTPUTS_ROOT / "results"
FIGURES_DIR  = OUTPUTS_ROOT / "figures"
MODELS_DIR   = OUTPUTS_ROOT / "models"
LOGS_DIR     = OUTPUTS_ROOT / "logs"
NETWORK_DIR     = RESULTS_DIR  / "network"
REPORTS_DIR     = RESULTS_DIR  / "reports"
ENRICHMENT_DIR  = RESULTS_DIR  / "enrichment"

# Processed expression matrices (after QC + harmonization)
TUMOR_EXPR_PROCESSED  = PROCESSED_DIR / "tumor_expression_processed.csv"
NORMAL_EXPR_PROCESSED = PROCESSED_DIR / "normal_expression_processed.csv"

# Harmonized (same gene set, same order)
TUMOR_EXPR_HARMONIZED  = PROCESSED_DIR / "tumor_expression_harmonized.csv"
NORMAL_EXPR_HARMONIZED = PROCESSED_DIR / "normal_expression_harmonized.csv"

# Differential expression results
DE_RESULTS_FILE = RESULTS_DIR / "differential_expression_results.csv"

# Expression features per gene
EXPR_FEATURES_FILE = PROCESSED_DIR / "expression_features.csv"

# Co-expression network edge list
NETWORK_EDGES_FILE  = NETWORK_DIR / "coexpression_network_edges.csv"
NETWORK_GRAPH_FILE  = NETWORK_DIR / "coexpression_network.graphml"

# Network features per gene
NETWORK_FEATURES_FILE = PROCESSED_DIR / "network_features.csv"

# Integrated feature matrix (expression + network)
INTEGRATED_FEATURES_FILE = PROCESSED_DIR / "integrated_features.csv"

# Label vector
LABELS_FILE = PROCESSED_DIR / "gene_labels.csv"

# Train / validation splits
TRAIN_FEATURES_FILE = PROCESSED_DIR / "train_features.csv"
VAL_FEATURES_FILE   = PROCESSED_DIR / "val_features.csv"
TRAIN_LABELS_FILE   = PROCESSED_DIR / "train_labels.csv"
VAL_LABELS_FILE     = PROCESSED_DIR / "val_labels.csv"

# Model outputs
GENE_RANKINGS_FILE      = RESULTS_DIR / "gene_rankings.csv"
FEATURE_IMPORTANCE_FILE = RESULTS_DIR / "feature_importance.csv"
MODEL_METRICS_FILE      = RESULTS_DIR / "model_metrics.csv"

# Annotated network
ANNOTATED_NETWORK_FILE  = NETWORK_DIR / "annotated_network.graphml"
ANNOTATED_EDGES_FILE    = NETWORK_DIR / "annotated_edges.csv"
ANNOTATED_NODES_FILE    = NETWORK_DIR / "annotated_nodes.csv"

# Interactive visualization
INTERACTIVE_HTML_FILE   = FIGURES_DIR / "interactive_network.html"

# Final report
FINAL_REPORT_FILE       = REPORTS_DIR / "pipeline_summary_report.csv"

# KEGG enrichment (Step 19)
KEGG_ALL_FILE           = ENRICHMENT_DIR / "kegg_all_candidates.csv"
KEGG_UP_FILE            = ENRICHMENT_DIR / "kegg_upregulated.csv"
KEGG_DOWN_FILE          = ENRICHMENT_DIR / "kegg_downregulated.csv"
KEGG_SUMMARY_FILE       = ENRICHMENT_DIR / "kegg_summary.csv"

# ==================================================
# DATA FORMAT ASSUMPTIONS
# ==================================================

TUMOR_GENE_ID_COL   = "Hugo_Symbol"
TUMOR_SAMPLE_START  = 2

NORMAL_GENE_ID_COL  = "Description"
NORMAL_SAMPLE_START = 2

NORMAL_META_SAMPLE_COL = "SAMPID"
TUMOR_META_SAMPLE_COL  = "PATIENT_ID"

# LCGene TSV uses 'GeneSymbol' as the gene-symbol column
LCGENE_GENE_COL = "GeneSymbol"

# ==================================================
# QC PARAMETERS
# ==================================================

MIN_EXPRESSION_FRACTION = 0.10
MIN_EXPRESSION_VALUE    = 1.0
MIN_SAMPLES_TUMOR  = 50
MIN_SAMPLES_NORMAL = 50

# ==================================================
# DIFFERENTIAL EXPRESSION PARAMETERS
# ==================================================

# P2.2 FIX: Tightened from |log2FC|≥1.0 / FDR≤0.05 to reduce the near-total
# DE significance (~95.3% of genes) caused by cohort/platform batch effects.
# Rerun from step4 through step14 after changing these values.
DE_LOG2FC_THRESHOLD  = 2.0
DE_PVALUE_THRESHOLD  = 0.001

# ==================================================
# CO-EXPRESSION NETWORK PARAMETERS
# ==================================================

COEXPR_CORRELATION_METHOD  = "pearson"
# P2.3 FIX: Lowered from 0.70 to 0.60 to reduce isolated nodes and increase
# network feature variance. More genes will be connected, improving discriminative
# power of network features. Rerun from step6 through step14 after changing.
COEXPR_CORRELATION_CUTOFF  = 0.60
COEXPR_MIN_SAMPLES         = 30

# ==================================================
# ML PARAMETERS
# ==================================================

RANDOM_STATE      = 42
TEST_SIZE         = 0.20
CV_FOLDS          = 5
POSITIVE_LABEL    = 1
NEGATIVE_LABEL    = 0

# ==================================================
# EVALUATION SETTINGS
# ==================================================

# Primary metric for selecting the best model: "auroc" or "auprc"
# AUPRC is preferred for imbalanced datasets (Step 3)
CV_METRIC_PRIMARY = "auprc"

# Threshold strategy for binary classification decisions (Step 4):
#   "f1"               — maximise F1 on the PR curve (default)
#   "target_recall"    — smallest threshold reaching THRESHOLD_TARGET_RECALL
#   "target_precision" — smallest threshold reaching THRESHOLD_TARGET_PRECISION
#   "top_k"            — flag top-THRESHOLD_TOP_K ranked genes as positive
THRESHOLD_STRATEGY         = "f1"
THRESHOLD_TARGET_RECALL    = 0.80
THRESHOLD_TARGET_PRECISION = 0.50
THRESHOLD_TOP_K            = 200

# SMOTE resampling — applied inside each CV training fold only (Step 5)
# Requires: pip install imbalanced-learn
USE_SMOTE = False

# Random undersampling — applied inside each CV training fold only (Step 5)
# Applied after SMOTE when both are True (SMOTE-then-undersample pipeline)
# Requires: pip install imbalanced-learn
USE_UNDERSAMPLING = True
# Only change UNDERSAMPLING_RATIO if USE_SMOTE=False. With SMOTE active, keep at 1.0.
UNDERSAMPLING_RATIO = 3.0

# REC 3: Hyperparameter tuning — set True for a full tuned run (much slower)
USE_HYPERPARAMETER_TUNING = True
HP_N_ITER = 50

USE_FEATURE_SELECTION = True
FEATURE_SELECTION_TOP_N = 50

# Positive-Unlabeled framing (Step 6)
# When True: negatives are treated as "unlabeled" (unannotated) rather than
# confirmed non-cancer genes.  Labels (0/1) are unchanged; only language changes.
PU_FRAMING = True

# ==================================================
# PHASE 3 — PU BAGGING (Step 11b)
# ==================================================

# Mordelet-Vert bagging: number of base classifiers
PU_N_ESTIMATORS    = 300
# Ratio of unlabeled subsample size to positive set size per iteration
PU_SUBSAMPLE_RATIO = 1.0
# Trees per base RF classifier inside each bagging iteration
PU_BASE_N_TREES    = 100

# ==================================================
# PHASE 3 — DIFFERENTIAL CO-EXPRESSION (Steps 6b / 7b)
# ==================================================

# Normal co-expression network outputs
NORMAL_NETWORK_GRAPH_FILE         = NETWORK_DIR / "normal_coexpression_network.graphml"
NORMAL_NETWORK_EDGES_FILE         = NETWORK_DIR / "normal_coexpression_network_edges.csv"
DIFFERENTIAL_NETWORK_FEATURES_FILE = PROCESSED_DIR / "differential_network_features.csv"

# ==================================================
# PHASE 3 — INDEPENDENT COHORT VALIDATION (Step 15)
# ==================================================

# Set these to paths of an external dataset to enable independent validation.
# EXTERNAL_EXPR_FILE:   CSV file, genes × samples (same gene symbols as pipeline)
# EXTERNAL_LABELS_FILE: CSV with columns [gene, label] where label ∈ {0, 1}
EXTERNAL_EXPR_FILE   = None
EXTERNAL_LABELS_FILE = None

# ==================================================
# PHASE 3 — BATCH CORRECTION (Step 1b)
# ==================================================

# Set True and install pyComBat (pip install inmoose) to apply ComBat-seq
# correction to the combined tumor+normal matrix before differential expression.
USE_BATCH_CORRECTION        = False
BATCH_CORRECTED_TUMOR_FILE  = PROCESSED_DIR / "tumor_expression_batch_corrected.csv"
BATCH_CORRECTED_NORMAL_FILE = PROCESSED_DIR / "normal_expression_batch_corrected.csv"

# Novel candidate probability thresholds (Step 14)
# NOVEL_PROB_THRESHOLD      — primary high-confidence list (used in paper)
# NOVEL_PROB_THRESHOLD_SENS — sensitivity/extended list (Supplementary Table)
NOVEL_PROB_THRESHOLD      = 0.85
NOVEL_PROB_THRESHOLD_SENS = 0.50

# ==================================================
# LOGGING
# ==================================================

LOG_FILE   = LOGS_DIR / "pipeline.log"
LOG_LEVEL  = "INFO"


def create_output_dirs() -> None:
    """Create all necessary output directories if they do not exist."""
    dirs = [
        OUTPUTS_ROOT,
        PROCESSED_DIR,
        RESULTS_DIR,
        FIGURES_DIR,
        MODELS_DIR,
        LOGS_DIR,
        NETWORK_DIR,
        REPORTS_DIR,
        ENRICHMENT_DIR,
    ]
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)


def validate_input_files() -> None:
    """
    Check that every raw input file declared in this config actually exists
    on disk.  Raises FileNotFoundError with a clear message if any are missing.
    """
    required = {
        "Tumor expression matrix" : TUMOR_EXPR_FILE,
        "Tumor metadata"           : TUMOR_META_FILE,
        "Normal expression matrix" : NORMAL_EXPR_FILE,
        "Normal metadata"          : NORMAL_META_FILE,
        "Labeled LUAD gene list"   : CANCER_GENE_FILE,
    }
    missing = []
    for label, path in required.items():
        if not path.exists():
            missing.append(f"  [{label}]  ->  {path}")

    if missing:
        msg = "The following required input files were NOT found:\n" + "\n".join(missing)
        raise FileNotFoundError(msg)

    print("[config] All required input files found.")


if __name__ == "__main__":
    print("=== LUAD ML Pipeline — Configuration ===")
    print(f"Project root : {PROJECT_ROOT}")
    print(f"Data root    : {DATA_ROOT}")
    print()
    create_output_dirs()
    print("[config] Output directories created / verified.")
    try:
        validate_input_files()
    except FileNotFoundError as e:
        print(f"[config] WARNING — {e}")
    print()
    print("--- Key paths ---")
    print(f"  Tumor expr        : {TUMOR_EXPR_FILE}")
    print(f"  Normal expr       : {NORMAL_EXPR_FILE}")
    print(f"  Labeled gene list : {CANCER_GENE_FILE}")
    print(f"  Results           : {RESULTS_DIR}")
    print(f"  Figures           : {FIGURES_DIR}")
    print()
    print("[config] Step 0 complete.")
