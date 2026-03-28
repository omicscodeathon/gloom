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

# Project root: directory containing this config.py
PROJECT_ROOT = Path(__file__).resolve().parent

# Data root: sits one level above the project folder, inside LUML/
DATA_ROOT = PROJECT_ROOT.parent / "data"

# ==================================================
# RAW INPUT FILE PATHS
# ==================================================

RAW_DIR = DATA_ROOT / "raw"

# --- Tumor (cBioPortal LUAD) ---
TUMOR_EXPR_FILE   = RAW_DIR / "cBioPortal (RNA Seq Data)" / "data_mrna_seq_v2_rsem.txt"
TUMOR_META_FILE   = RAW_DIR / "cBioPortal (RNA Seq Data)" / "data_clinical_patient.txt"

# --- Normal (GTEx lung) ---
NORMAL_EXPR_FILE  = RAW_DIR / "Gtex (normal samples)" / "gene_tpm_v11_lung.gct.gz"
NORMAL_META_FILE  = RAW_DIR / "Gtex (normal samples)" / "GTEx_Analysis_v11_Annotations_SampleAttributesDD.xlsx"

# --- Known cancer / LUAD genes (Cancer Gene Census) ---
CANCER_GENE_FILE  = RAW_DIR / "Cancer Gene Census (Labeled Data)" / \
                    "Census_allWed Mar 18 05_48_31 2026.csv"

# ==================================================
# PROCESSED / INTERMEDIATE OUTPUT PATHS
# ==================================================

PROCESSED_DIR       = DATA_ROOT / "processed"
RESULTS_DIR         = PROJECT_ROOT / "results"
FIGURES_DIR         = PROJECT_ROOT / "figures"
MODELS_DIR          = PROJECT_ROOT / "models"
LOGS_DIR            = PROJECT_ROOT / "logs"
NETWORK_DIR         = RESULTS_DIR / "network"
REPORTS_DIR         = RESULTS_DIR / "reports"

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

# ==================================================
# DATA FORMAT ASSUMPTIONS
# ==================================================

TUMOR_GENE_ID_COL   = "Hugo_Symbol"
TUMOR_SAMPLE_START  = 2

NORMAL_GENE_ID_COL  = "Description"
NORMAL_SAMPLE_START = 2

NORMAL_META_SAMPLE_COL = "SAMPID"
TUMOR_META_SAMPLE_COL  = "PATIENT_ID"
CGC_GENE_COL = "Gene Symbol"

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

DE_LOG2FC_THRESHOLD  = 1.0
DE_PVALUE_THRESHOLD  = 0.05

# ==================================================
# CO-EXPRESSION NETWORK PARAMETERS
# ==================================================

COEXPR_CORRELATION_METHOD  = "pearson"
COEXPR_CORRELATION_CUTOFF  = 0.70
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
# LOGGING
# ==================================================

LOG_FILE   = LOGS_DIR / "pipeline.log"
LOG_LEVEL  = "INFO"


def create_output_dirs() -> None:
    """Create all necessary output directories if they do not exist."""
    dirs = [
        PROCESSED_DIR,
        RESULTS_DIR,
        FIGURES_DIR,
        MODELS_DIR,
        LOGS_DIR,
        NETWORK_DIR,
        REPORTS_DIR,
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
        "Cancer gene list"         : CANCER_GENE_FILE,
    }
    missing = []
    for label, path in required.items():
        if not path.exists():
            missing.append(f"  [{label}]  →  {path}")

    if missing:
        msg = "The following required input files were NOT found:\n" + "\n".join(missing)
        raise FileNotFoundError(msg)

    print("[config] All required input files found.")


if __name__ == "__main__":
    print("=== LUAD ML Pipeline - Configuration ===")
    print(f"Project root : {PROJECT_ROOT}")
    print(f"Data root    : {DATA_ROOT}")
    print()
    create_output_dirs()
    print("[config] Output directories created / verified.")
    try:
        validate_input_files()
    except FileNotFoundError as e:
        print(f"[config] WARNING - {e}")
    print()
    print("--- Key paths ---")
    print(f"  Tumor expr   : {TUMOR_EXPR_FILE}")
    print(f"  Normal expr  : {NORMAL_EXPR_FILE}")
    print(f"  Cancer genes : {CANCER_GENE_FILE}")
    print(f"  Results      : {RESULTS_DIR}")
    print(f"  Figures      : {FIGURES_DIR}")
    print()
    print("[config] Step 0 complete.")
