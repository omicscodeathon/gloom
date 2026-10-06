"""
config.py
---------
Central configuration for the LUAD ML Pipeline.
All paths, constants, and parameters are defined here.
No hardcoded values should appear in other modules.
"""

import os
import sys
from pathlib import Path


def _make_console_output_safe() -> None:
    """Avoid UnicodeEncodeError on Windows code pages such as cp1256."""
    for stream_name in ("stdout", "stderr"):
        stream = getattr(sys, stream_name, None)
        if hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(errors="backslashreplace")
            except Exception:
                pass


_make_console_output_safe()

# ==================================================
# ROOT PATHS
# ==================================================

# Project root: directory containing this config.py
PROJECT_ROOT = Path(__file__).resolve().parents[3]

# Data root: sits one level above the project folder, inside LUML/
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
OUTPUTS_ROOT = PROJECT_ROOT / "outputs"
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
TUMOR_EXPR_HARMONIZED  = PROCESSED_DIR / "tumor_expression_batch_corrected.csv"
NORMAL_EXPR_HARMONIZED = PROCESSED_DIR / "normal_expression_batch_corrected.csv"

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
# kegg_all_candidates.csv is the complete UNFILTERED table; the lung/cancer keyword subset is
# written separately to ENRICHMENT_DIR / "kegg_lung_cancer_subset.csv" (interpretive only).
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
# DATA SOURCE / COHORT DESIGN  (v0.2.0)
# ==================================================

# "cbioportal_gtex" : legacy design — tumors from TCGA/cBioPortal (RSEM), normals from GTEx (TPM).
#                     Tumor status is PERFECTLY CONFOUNDED with cohort/pipeline in this design
#                     (see COHORT_DESIGN check in step1 and docs/revision/).
# "tcga_gdc"        : uniformly processed TCGA-LUAD primary tumors vs TCGA-LUAD adjacent normals
#                     (GDC STAR counts), produced by scripts/fetch_gdc_tcga_luad.py.
DATA_SOURCE = "cbioportal_gtex"

GDC_DIR               = RAW_DIR / "tcga_gdc"
GDC_SAMPLE_SHEET_FILE = GDC_DIR / "sample_sheet.csv"     # sample_id, patient_id, group
GDC_COUNTS_FILE       = GDC_DIR / "counts_matrix.csv"    # genes x samples, raw STAR unstranded counts
GDC_TPM_FILE          = GDC_DIR / "tpm_matrix.csv"       # genes x samples, STAR tpm_unstranded
# Matrix fed to the pipeline when DATA_SOURCE == "tcga_gdc": "tpm" (default) or "counts".
GDC_EXPRESSION_MATRIX = "tpm"

# Cohort labels used by the COHORT_DESIGN collinearity check in step1.  A sample-level
# column named one of COHORT_COLUMN_CANDIDATES in the metadata takes precedence.
COHORT_COLUMN_CANDIDATES = ("cohort", "batch", "dataset", "source")
# Used for the legacy design; ignored when DATA_SOURCE == "tcga_gdc" (both groups then
# belong to the same cohort, TCGA-LUAD / GDC STAR).
TUMOR_COHORT_NAME  = "TCGA-LUAD (cBioPortal RSEM)"
NORMAL_COHORT_NAME = "GTEx v11 lung (TPM)"

# Thresholds for the sample-level QC warnings of step2b (written to results/qc_cohort_warning.txt)
QC_MAX_ABS_MEDIAN_LOG2FC = 1.0    # warn if |median log2FC over genes| exceeds this
QC_MAX_FRACTION_DE       = 0.80   # warn if more than this fraction of genes is DE-significant
QC_MAX_UP_DOWN_RATIO     = 10.0   # warn if up:down > 10 or < 1/10
QC_PCA_TOP_GENES         = 5000   # most variable genes used for the sample PCA

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

# Differential-expression test used by step4:
#   "welch"      — Welch t-test on log2(x+1) values (default; unpaired)
#   "paired"     — paired t-test on patients that have BOTH a tumor and an adjacent-normal
#                  sample (needs a patient_id for every sample; see PAIRING_FILE below and
#                  docs/tcga_paired_design.md).  Falls back to Welch if no pairs are found.
#   "limma_voom" — limma-voom on raw counts via rpy2 + R (optional dependency; if rpy2/limma
#                  or a counts matrix is unavailable the step logs a clear message and
#                  falls back to Welch).  Needs DATA_SOURCE = "tcga_gdc" (counts matrix).
DE_METHOD = "welch"

# Optional CSV with columns sample_id, patient_id (tumor AND normal samples).  When None,
# the patient_id column of the processed sample metadata is used, and as a last resort the
# first 12 characters of a TCGA barcode (TCGA-XX-XXXX).
PAIRING_FILE = None

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

SEED              = 42      # single source of randomness for the whole pipeline
RANDOM_STATE      = SEED    # backward-compatible alias (do not set independently)
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
# UNDERSAMPLING_RATIO: only change this if USE_SMOTE=False.
# With SMOTE active, SMOTE already creates 1:1 balance, so RUS must stay at 1.0.
# Setting 0.33 here while SMOTE is on causes RUS to fail (can't add majority samples).
UNDERSAMPLING_RATIO = 3.0

# REC 3: Hyperparameter tuning — set True for a full tuned run (much slower)
# HP_N_ITER: RandomizedSearchCV iterations per model (20 = fast, 50 = thorough)
USE_HYPERPARAMETER_TUNING = True
HP_N_ITER = 20

# REC 5: Feature selection — uses top-N features from the last feature_importance.csv run.
# Set True after step13 has been run at least once.
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
# ==================================================


# ==================================================
# PHASE 3 — BATCH CORRECTION (Step 1b)
# ==================================================

# Set True and install pyComBat (pip install inmoose) to apply ComBat-seq
# correction to the combined tumor+normal matrix before differential expression.
#
# DEFAULT IS False — and for the legacy TCGA(cBioPortal) vs GTEx design it SHOULD stay False:
# when every tumor comes from one cohort and every normal from another, the batch variable is
# perfectly collinear with disease status, so ComBat / ComBat-seq cannot tell cohort effects
# from disease effects (it would either remove the disease signal or leave the confounding in
# place).  Batch correction is only meaningful when each cohort contains BOTH tumor and normal
# samples.  The recommended remedy is a uniformly processed design
# (DATA_SOURCE = "tcga_gdc": TCGA-LUAD tumor vs TCGA-LUAD adjacent normal).  See
# step1b_batch_correction.py and docs/tcga_paired_design.md.
USE_BATCH_CORRECTION        = False
BATCH_CORRECTED_TUMOR_FILE  = PROCESSED_DIR / "tumor_expression_batch_corrected.csv"
BATCH_CORRECTED_NORMAL_FILE = PROCESSED_DIR / "normal_expression_batch_corrected.csv"

# Non-LCGene candidate probability thresholds (Step 14)
# Absence from the LCGene reference set does NOT establish biological novelty; these genes are
# reported as "non-LCGene candidates" and are checked against external resources in step20.
# CANDIDATE_PROB_THRESHOLD      — primary high-confidence list
# CANDIDATE_PROB_THRESHOLD_SENS — sensitivity/extended list (Supplementary Table)
CANDIDATE_PROB_THRESHOLD      = 0.85
CANDIDATE_PROB_THRESHOLD_SENS = 0.50
# Minimum |log2FC| for a gene to be listed as a candidate (post-hoc reporting filter only —
# it does not influence training or the ranking itself).
CANDIDATE_MIN_ABS_LOG2FC      = 1.0
# Deprecated aliases (v0.1.x names), kept so that old scripts still run.
NOVEL_PROB_THRESHOLD      = CANDIDATE_PROB_THRESHOLD
NOVEL_PROB_THRESHOLD_SENS = CANDIDATE_PROB_THRESHOLD_SENS

# ==================================================
# v0.2.0 — REVISION: LABEL-INDEPENDENT UNIVERSE, CROSS-FITTING, ABLATION, ...
# ==================================================

# ---- Step 9: analysis universe -------------------------------------------------------------
# True  (default): the analysis universe is fixed BEFORE labels are assigned, with criteria that
#                  are identical for every gene (detectable expression + variance filter of
#                  step2, intersection of tumor/normal genes in step3).  No |log2FC| filter is
#                  applied to unlabeled genes.
# False (DEPRECATED — label-dependent, NOT recommended): reproduces v0.1.x, which removed
#                  unlabeled genes with |log2FC| < 2 while keeping all positives.  Because
#                  log2FC-derived variables are model features, this inflates class separation.
LABEL_INDEPENDENT_UNIVERSE = True

# ---- Steps 11c / 12b: cross-fitting evaluation (PRIMARY ranking) -----------------------------
# Stratified K-fold cross-fitting of the PU-bagging model, repeated with different seeds.
# Every gene is scored only by models that never saw it; the ranking and all top-K metrics are
# computed from these out-of-fold (OOF) scores.
USE_CROSSFIT          = True
CROSSFIT_FOLDS        = 5
CROSSFIT_REPEATS      = 5
# PU-bagging size used inside cross-fitting (None -> PU_N_ESTIMATORS, i.e. identical to step11b).
CROSSFIT_PU_N_ESTIMATORS = None
EVAL_K_VALUES         = (10, 50, 100)   # K for Precision@K / Recall@K / enrichment factor@K
BOOTSTRAP_N           = 1000            # stratified bootstrap resamples over genes
BOOTSTRAP_SEED        = SEED
BOOTSTRAP_ALPHA       = 0.05            # 95 % CI

# ---- Step 13b: feature-set ablation ------------------------------------------------------------
USE_ABLATION = True
# PU-bagging size and repeats used for the ablation variants.  Smaller than the main
# cross-fitting to keep the run time manageable (7 variants); set to None to use the
# CROSSFIT_* values instead.
ABLATION_PU_N_ESTIMATORS = 100
ABLATION_REPEATS         = 3

# Feature groups are matched against the column names produced by steps 5, 7, 7b and 8 (only
# columns that actually exist are used).  Explicit lists are used because the prefix "tumor_"
# appears both in expression features (tumor_mean) and in network features (tumor_degree).
_EXPR_STATS = ("mean", "median", "std", "iqr", "cv", "skewness", "kurtosis", "pct_expressed")
# Differential-expression / fold-change derived features (explicit DE statistics).
FEATURE_GROUP_DE_STATS = [
    "log2fc", "abs_log2fc", "neg_log10_padj", "cohens_d", "t_stat",
    "tumor_normal_mean_ratio",
    "abs_log2fc_rank", "neg_log10_padj_rank", "cohens_d_rank",
]
# Expression-derived features (step5): per-group stats, contrasts, ranks + the DE statistics.
FEATURE_GROUP_EXPRESSION = (
    [f"tumor_{s}" for s in _EXPR_STATS]
    + [f"normal_{s}" for s in _EXPR_STATS]
    + ["std_ratio", "iqr_ratio", "tumor_mean_rank", "tumor_iqr_rank"]
    + FEATURE_GROUP_DE_STATS
)
# Network-derived features (step7 tumor-network topology + step7b differential co-expression).
FEATURE_GROUP_NETWORK = [
    "degree", "weighted_degree", "avg_neighbor_degree", "betweenness_centrality",
    "closeness_centrality", "eigenvector_centrality", "clustering_coefficient",
    "mean_edge_weight", "max_edge_weight", "min_edge_weight", "std_edge_weight",
    "in_largest_component", "component_size",
    "tumor_degree", "normal_degree", "delta_degree",
    "tumor_clustering", "normal_clustering", "delta_clustering",
    "tumor_mean_edge_weight", "normal_mean_edge_weight", "delta_mean_edge_weight",
    "tumor_betweenness", "normal_betweenness", "delta_betweenness",
    "tumor_specific_degree", "normal_specific_degree", "rewiring_ratio",
]

# ---- Step 7c: network stability --------------------------------------------------------------
NETWORK_THRESHOLDS      = (0.50, 0.55, 0.60, 0.65, 0.70)   # |r| cut-offs
NETWORK_EQUAL_N         = True    # sub-sample the larger group to the size of the smaller one
NETWORK_BOOTSTRAP_N     = 20      # sample bootstrap resamples (per group)
NETWORK_BOOTSTRAP_REPLACE = True  # resample samples with replacement
NETWORK_HUB_FRACTION    = 0.05    # hubs = top 5 % of genes by degree

# ---- Step 19: enrichment -----------------------------------------------------------------------
# Primary output: complete, UNFILTERED over-representation table (explicit background = the
# analysis universe, BH correction).  The lung/cancer keyword subset below is an interpretive
# display only and is never used to decide significance.
ENRICHMENT_THEMATIC_SUBSET = True
# Optional local GMT file (e.g. a KEGG GMT).  When None, the library below is downloaded with
# gseapy.get_library (needs internet and gseapy).
ENRICHMENT_GENESETS_FILE   = None
ENRICHMENT_LIBRARY         = "KEGG_2021_Human"
ENRICHMENT_MIN_TERM_SIZE   = 1      # terms with fewer background genes are not tested at all

# ---- Step 20: external evidence for non-LCGene candidates (optional, needs internet) ----
USE_CANDIDATE_EVIDENCE   = True     # step20 never breaks the pipeline; it is skipped offline
EVIDENCE_TOP_N           = 200      # top-N non-LCGene genes of the ranking are checked
EVIDENCE_MIN_EPMC_HITS   = 3        # Europe PMC hit count regarded as "evidence"
EVIDENCE_MIN_OT_SCORE    = 0.05     # Open Targets association score regarded as "evidence"
EVIDENCE_TIMEOUT_S       = 20
EVIDENCE_RETRIES         = 3
EVIDENCE_WORKERS         = 4

# ==================================================
# LOGGING
# ==================================================

LOG_FILE   = LOGS_DIR / "pipeline.log"
LOG_LEVEL  = "INFO"


def set_global_seed(seed: int = None) -> int:
    """
    Seed Python's and NumPy's global generators with ``config.SEED`` (default).

    All steps already pass ``config.SEED`` explicitly to every estimator / generator they use;
    this call additionally pins the global generators for any third-party code that falls back
    on them.  Returns the seed that was applied.
    """
    import random

    seed = SEED if seed is None else int(seed)
    random.seed(seed)
    try:
        import numpy as np
        np.random.seed(seed)
    except ImportError:  # pragma: no cover - numpy is a hard dependency of the pipeline
        pass
    return seed


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
    if DATA_SOURCE == "tcga_gdc":
        expr_file = GDC_COUNTS_FILE if GDC_EXPRESSION_MATRIX == "counts" else GDC_TPM_FILE
        required = {
            "GDC expression matrix"  : expr_file,
            "GDC sample sheet"       : GDC_SAMPLE_SHEET_FILE,
            "Labeled LUAD gene list" : CANCER_GENE_FILE,
        }
    else:
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
