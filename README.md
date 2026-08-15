# 🌌 GLOOM

> **Gene network Learning and Organization through Optimized Machine intelligence**  
> A reproducible package for **gene prioritization, co-expression network construction, and machine-learning-based candidate gene discovery**.

<p align="center">
  <img src="https://github.com/user-attachments/assets/34d8ddad-94ed-49e0-8c3b-fa7006642db5" alt="new logo" width="500" />
</p>

![Python](https://img.shields.io/badge/Python-3.12%2B-3776AB?style=for-the-badge&logo=python&logoColor=white)
![Machine Learning](https://img.shields.io/badge/Machine%20Learning-Gene%20Prioritization-8A2BE2?style=for-the-badge)
![Network Biology](https://img.shields.io/badge/Network%20Biology-Co--expression-00BFA6?style=for-the-badge)
![Cancer Genomics](https://img.shields.io/badge/Cancer%20Genomics-Extensible-FF6B6B?style=for-the-badge)
![License](https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge)

> 🧪 **Starting workflow:** LUAD  
> 🕸️ **Core idea:** combine expression + network topology + machine learning  
> 📈 **Main output:** ranked known and novel candidate genes  
> 🔬 **Goal:** support biological discovery and downstream validation across disease contexts

## 🌐 Project Landing Page

Explore the interactive GLOOM landing page here:

👉 **[Open GLOOM Landing Page](https://beta.sma-it.com/UI/gloom/)**

👉 **Cheat Sheet -> [preview.html](https://github.com/user-attachments/files/27089352/preview.html)**

---

## 🧭 Table of Contents

- 🌌 [Overview](#-overview)
- 💡 [Why GLOOM?](#-why-gloom)
- ⚙️ [What GLOOM Does](#️-what-gloom-does)
- 📊 [Key Results from the LUAD Case Study](#-key-results-from-the-luad-case-study)
- 🔄 [Workflow](#-workflow)
- 🧬 [Data Sources](#-data-sources)
- 🛠️ [Pipeline Stages](#️-pipeline-stages)
- 📦 [Installation](#-installation)
- 💻 [CLI Commands Reference](#-cli-commands-reference)
- 🚀 [Quick Start: gloom prioritize](#-quick-start-gloom-prioritize)
- 📁 [Bring Your Own Data: gloom run](#-bring-your-own-data-gloom-run)
- 🔀 [Understanding the Two Commands](#-understanding-the-two-commands)
- 📝 [Input File Formats](#-input-file-formats)
- 📂 [Output Files](#-output-files)
- 🤖 [Model Performance in the LUAD Case Study](#-model-performance-in-the-luad-case-study)
- 🔬 [Biological Validation](#-biological-validation)
- ⭐ [Feature Importance](#-feature-importance)
- ⚠️ [Limitations](#️-limitations)
- 🔮 [Future Work](#-future-work)
- 🗂️ [Repository Structure](#️-repository-structure)
- 🎯 [Example Use Case](#-example-use-case)
- 👥 [Contributors](#-contributors)
- 📚 [Acknowledgments](#-acknowledgments)
- ✅ [Summary](#-summary)

---

## 🌌 Overview

**GLOOM** is a reproducible Python package for prioritizing disease-associated genes using integrated expression analysis, co-expression network construction, network topology, and machine learning.

The package is designed as a **general framework** that can be adapted when suitable disease/control expression matrices and reference gene labels are provided. The current workflow starts with **lung adenocarcinoma (LUAD)** as the first demonstrated disease case study.

Traditional differential expression analysis identifies genes whose average expression differs between disease and control samples. However, some biologically important genes may not show the strongest fold-change but may occupy important positions in molecular interaction or co-expression networks. GLOOM addresses this limitation by combining:

- 🧪 Disease and control expression statistics
- 📈 Differential expression metrics
- 🕸️ Co-expression network topology
- 🤖 Supervised machine-learning classification (PU learning framework)
- 🧬 Gene-level ranking and biological validation

The final output is a ranked list of known and novel candidate genes supported by model scores, feature importance, pathway enrichment, and annotated network exports.

---

## 💡 Why GLOOM?

Complex diseases are often driven by many interacting genes rather than a single dominant alteration. This makes candidate gene discovery difficult when using only single-gene statistical methods.

GLOOM was designed to provide a unified and reusable solution for **gene prioritization and biological network analysis**.

| Challenge | Traditional Workflow | GLOOM Solution |
|---|---|---|
| Many separate tools are required | Manual preprocessing, DEA, network analysis, ML, visualization | One scriptable end-to-end package |
| Expression-only methods miss network-driven genes | Genes are analyzed independently | Co-expression topology is integrated |
| Candidate lists are difficult to prioritize | Long DEG lists without ranking context | ML-based probability ranking |
| Reproducibility is hard | Many manual intermediate steps | Modular pipeline with cached outputs |
| Biological interpretation is fragmented | Tables, plots, and networks generated separately | Reports, enrichment, and network exports generated together |

---

## ⚙️ What GLOOM Does

GLOOM converts gene expression inputs into biologically interpretable candidate-gene rankings and network outputs.

```mermaid
flowchart LR
    A[Raw Expression Data] --> B[Preprocessing & QC]
    B --> C[Gene Harmonization]
    C --> D[Differential Expression]
    D --> E[Co-expression Network]
    E --> F[Network Features]
    F --> G[Feature Integration]
    G --> H[Machine Learning]
    H --> I[Gene Ranking]
    I --> J[Novel Candidate Genes]
    J --> K[Pathway Enrichment & Reports]
```

### Main capabilities

- 📥 Load disease and control expression matrices
- 🧹 Clean and harmonize gene identifiers
- 📊 Perform differential expression analysis (Welch t-test with BH correction)
- 🕸️ Construct a Pearson-correlation co-expression network
- 🔗 Extract graph-theoretic features (degree, centrality, clustering, components)
- 🧩 Integrate expression and network evidence into a unified feature matrix
- 🤖 Train and compare machine-learning models (Random Forest, Gradient Boosting, Extra Trees, Logistic Regression, and XGBoost when installed)
- 🏆 Build a calibrated soft-voting ensemble
- 📈 Rank genes by predicted disease relevance
- ✨ Identify high-confidence novel candidates
- 🛤️ Run KEGG pathway enrichment analysis
- 📤 Export interactive dashboards and network files (GraphML, Cytoscape)

---

## 📊 Key Results from the LUAD Case Study

The current workflow starts with lung adenocarcinoma as the first disease setting. The latest full local LUAD run in this repository (completed on August 14, 2026) produced the following results:

| Result | Value |
|---|--:|
| Shared genes analyzed | **10,986** |
| LUAD tumor samples (TCGA/cBioPortal) | **510** |
| Normal lung samples (GTEx v11) | **604** |
| Known LUAD genes retained (LCGene) | **438 / 517** |
| Genes tested for differential expression | **10,986** |
| Significant DE genes (FDR <= 0.001) | **9,715** |
| Upregulated genes | **9,700** |
| Downregulated genes | **15** |
| Ranked genes | **10,986** |
| High-confidence novel candidates | **203** |
| Co-expression network nodes (annotated) | **10,986** |
| Co-expression network edges | **110,508** |
| Largest connected component | **4,591 nodes** |
| Integrated features | **30** (21 after selection) |
| Best CV model (XGBoost) Mean AUPRC | **0.4951** |
| Best held-out model (Extra Trees) Val AUPRC | **0.5443** |
| Ensemble Val AUROC | **0.9347** |
| LCGene genes recovered in top 100 ranks | **98** |
| Total pipeline runtime | **~12.2 min** |

---

## 🔄 Workflow

The GLOOM workflow follows a complete analysis path from raw data to biological interpretation.

```mermaid
flowchart TD
    A[Data Sources] --> A1[Disease Expression Data]
    A --> A2[Control / Normal Expression Data]
    A --> A3[Known Disease Gene Labels]

    A1 --> B[Preprocessing]
    A2 --> B
    A3 --> C[Label Construction]

    B --> D[Gene Set Harmonization]
    D --> E[Differential Expression Analysis]
    D --> F[Co-expression Network Construction]

    E --> G[Expression & DE Features]
    F --> H[Network Topology Features]

    G --> I[Integrated Feature Matrix]
    H --> I
    C --> I

    I --> J[Model Training & Validation]
    J --> K[Gene Ranking]
    K --> L[Novel Candidate Detection]
    L --> M[Pathway Enrichment]
    L --> N[Interactive Dashboard]
    L --> O[Annotated Network Export]
```

---

<p align="center">
  <img src="https://github.com/omicscodeathon/gloom/blob/main/workflow/Simplified Workflow.png" alt="workflow overview" width="900"/>
</p>

---

## 🧬 Data Sources

GLOOM can be adapted to any disease context where the user provides compatible disease/control expression matrices and a disease-related reference gene list.

In the starting LUAD workflow, the package used:

| Data Type | Source | Details |
|---|---|---|
| Disease expression | TCGA / cBioPortal LUAD | 20,531 genes, 512 samples (510 after QC) |
| Control expression | GTEx v11 lung tissue | 74,628 transcripts, 604 RNA-seq samples |
| Known disease genes | LCGene (LUAD-filtered) | 517 curated LUAD-associated genes |
| Disease metadata | cBioPortal clinical | 566 patients, 35 clinical variables |
| Normal metadata | GTEx sample attributes | 604 samples, 10 QC columns (RIN, autolysis, map rate) |

After harmonization:

- **10,986 shared genes** across tumor and normal datasets
- **510 LUAD tumor samples**
- **604 normal lung samples**
- **438 LCGene genes** retained in the shared gene universe (84.7%)

### ⚠️ Important note about LCGene

LCGene is a curated database of **expression-based LUAD biomarkers** (297 upregulated, 220 downregulated). It does **not** include many well-known mutation-driven oncogenes/tumor suppressors (e.g., EGFR, KRAS, TP53, BRAF, MET). The model therefore learns expression-signature patterns, not mutation-driven oncogenesis. Users studying mutation-driven genes should provide a custom `--labels` file that includes their genes of interest as known positives.

---

## 🛠️ Pipeline Stages

GLOOM runs 20 core stages (0-19) in sequence. Each stage caches its outputs so the pipeline can be resumed from any point. In addition, optional refinement steps can be enabled between core stages to improve robustness and help make the final rankings more accurate when the required inputs or settings are available.

| Stage | Name | Description |
|--:|---|---|
| 0 | Config | Create output directories and validate configuration |
| 1 | Data loading | Load disease, control, metadata, and label files |
| 2 | Preprocessing | Clean values, log2-transform, remove low-expression/low-variance genes |
| 3 | Gene harmonization | Standardize gene symbols and intersect to shared genes |
| 4 | Differential expression | Welch t-test with Benjamini-Hochberg FDR correction |
| 5 | Expression features | Generate 29 per-gene expression statistics (mean, std, CV, skewness, kurtosis, ranks) |
| 6 | Co-expression network | Build Pearson-correlation network (default cutoff \|r\| >= 0.6) |
| 7 | Network features | Extract 13 topology features (degree, centrality, clustering, components) |
| 8 | Feature integration | Merge expression + network features, remove collinear features (30 final) |
| 9 | Label construction | Assign PU labels using reference gene list, trim noisy unlabeled genes |
| 10 | Train/validation split | Stratified 80/20 split with feature selection (top 21 features) |
| 11 | Model training | Train core models with CV, hyperparameter tuning, and calibration (XGBoost is included when installed) |
| 12 | Model evaluation | Compute AUROC, AUPRC, F1, MCC, Brier score; build calibrated ensemble |
| 13 | Feature importance | Model-based and permutation importance analysis |
| 14 | Gene ranking | Score all genes, flag novel candidates based on probability thresholds |
| 15 | Network annotation | Add ranking and biological metadata to graph nodes |
| 16 | Network export | Export annotated network as GraphML, Cytoscape XML, edge/node tables |
| 17 | Interactive visualization | Generate volcano plots, ROC/PR curves, dashboards (Plotly HTML) |
| 18 | Final report | Generate final reporting artifacts (summary table, text report, summary figure) for the completed run |
| 19 | KEGG enrichment | Pathway enrichment analysis for candidate gene sets |

### Optional refinement steps

Beyond the core 0-19 flow, GLOOM also includes optional supplemental steps that can be enabled when you want additional correction, comparison, or ranking refinement:

| Stage | Name | Purpose |
|--:|---|---|
| 1b | Batch correction | Adjust tumor and control matrices for cross-cohort technical effects before downstream analysis |
| 6b | Normal co-expression network | Build a control-network reference that supports stronger tumor-vs-normal comparison |
| 7b | Differential network features | Add network rewiring features derived from tumor and normal graphs |
| 11b | PU bagging | Add a supplemental positive-unlabeled scoring pass to strengthen candidate ranking robustness |

---

## 📦 Installation

### Prerequisites

| Software | Recommended Version |
|---|---|
| Python | 3.12+ |
| Conda | Latest Anaconda or Miniconda |
| pip | Latest stable version |

### Option 1: Install from Bioconda (recommended)

Install GLOOM directly from Bioconda using the Anaconda Prompt or terminal:

```bash
conda install -c bioconda -c conda-forge gloom
```

This installs GLOOM and all core dependencies (click, pandas, numpy, scipy, statsmodels, scikit-learn, joblib, networkx, matplotlib).

After installation, install the optional features for the full pipeline experience:

```bash
conda install -c conda-forge plotly openpyxl imbalanced-learn
pip install gseapy>=0.10.8
```

If you use `mamba`, the equivalent command is:

```bash
mamba install -c bioconda -c conda-forge gloom
```

### Option 2: Install from source (development)

Clone the repository and install locally:

```bash
git clone https://github.com/omicscodeathon/gloom.git
cd gloom
```

**Using conda environment file (installs all dependencies at once):**

```bash
conda env create -f environment.yml
conda activate gloom
pip install -e .
```

**Or manually:**

```bash
conda create -n gloom python=3.12
conda activate gloom
pip install -e .
```

To install with all optional features:

```bash
pip install -e ".[full]"
```

Available optional groups:

| Group | What it enables | Install command |
|---|---|---|
| `interactive` | Plotly HTML dashboards (Step 17) | `pip install -e ".[interactive]"` |
| `kegg` | KEGG pathway enrichment (Step 19) | `pip install -e ".[kegg]"` |
| `excel` | Excel output (`--format excel`) | `pip install -e ".[excel]"` |
| `resampling` | SMOTE / undersampling (Step 11) | `pip install -e ".[resampling]"` |
| `full` | All optional features | `pip install -e ".[full]"` |
| `dev` | Developer tools (pytest, black, ruff, mypy) | `pip install -e ".[dev]"` |
| `all` | Everything (full + dev) | `pip install -e ".[all]"` |

### Option 3: Build the conda package locally

```bash
conda build conda.recipe/
conda install --use-local gloom
```

### Verify installation

```bash
gloom --version
gloom info
gloom diseases
```

---

## 💻 CLI Commands Reference

GLOOM provides the following CLI commands:

| Command | Purpose |
|---|---|
| `gloom info` | Show version, bundled reference paths, default thresholds, and cached runs |
| `gloom diseases` | List supported `--disease` values and their data sources |
| `gloom validate --data-dir DIR` | Check that all required input files are present before running |
| `gloom cache clear [--output DIR]` | Delete cached intermediate files to force a full re-run |
| `gloom prioritize` | Run the full pipeline using bundled reference data (cBioPortal + GTEx) |
| `gloom run` | Run the full pipeline using your own expression matrices |

### Shared options (`prioritize` and `run`)

| Option | Description | Default |
|---|---|---|
| `--output DIR` | Output directory (created if needed) | *required* |
| `--labels FILE` | Custom positive gene list (CSV/TSV with GeneSymbol column) | bundled LCGene |
| `--fdr FLOAT` | FDR cutoff for differential expression | 0.05 |
| `--log2fc FLOAT` | Log2 fold-change threshold for DE | 1.0 |
| `--prob-threshold FLOAT` | Minimum probability to flag novel candidates | 0.50 |
| `--top-k N` | Keep only the top N candidates in output | all |
| `--format {csv,excel,json}` | Output file format for candidate tables | csv |
| `--from-step STEP` | Resume pipeline from this step key | `0` (`prioritize`) / `2` (`run`) |
| `--to-step STEP` | Stop after this step key | `19` |
| `--skip-optional` | Exclude optional refinement steps from the selected run | off |
| `--skip-step STEP` | Skip one specific step key (may be repeated) | none |
| `--verbose` | Enable DEBUG-level logging | off |

### `gloom prioritize` only

| Option | Description | Default |
|---|---|---|
| `--genes FILE` | Query gene list (.txt, one symbol per line, or TSV with `GeneSymbol`) | *required* |
| `--disease {luad}` | Disease context | `luad` |
| `--data-dir DIR` | Override the bundled raw-data root | bundled data paths |
| `--no-cache` | Force a full re-run by clearing the existing cache in `--output` first | off |
| `--dry-run` | Validate inputs and print the execution plan without running the pipeline | off |

### `gloom run` only

| Option | Description | Default |
|---|---|---|
| `--tumor-expr FILE` | Tumor/disease expression matrix (CSV) | *required* |
| `--normal-expr FILE` | Normal/control expression matrix (CSV) | *required* |
| `--tumor-meta FILE` | Optional tumor sample metadata (CSV) | use all tumor samples |
| `--normal-meta FILE` | Optional normal sample metadata (CSV) | use all normal samples |
| `--genes FILE` | Optional query gene list to tag candidate genes in the final ranking | all genes ranked |

---

## 🚀 Quick Start: gloom prioritize

`gloom prioritize` runs the **full pipeline** (Steps 0-19) using the **bundled LUAD reference data** (cBioPortal tumor + GTEx normal + LCGene labels). The `--genes` file is used to tag query genes of interest in the final ranking.

```bash
gloom prioritize --genes genes.txt \
  --disease luad \
  --output results/
```

### Resume from a specific step

```bash
gloom prioritize --genes genes.txt \
  --output results/ \
  --from-step 11
```

### Change statistical thresholds

```bash
gloom prioritize --genes genes.txt \
  --fdr 0.05 \
  --log2fc 1.0 \
  --prob-threshold 0.5 \
  --output results/
```

### Export top genes as Excel

```bash
gloom prioritize --genes genes.txt \
  --top-k 100 \
  --format excel \
  --output results/
```

### Use custom labels instead of LCGene

```bash
gloom prioritize --genes genes.txt \
  --labels my_known_genes.csv \
  --output results/
```

### Preview what would run (dry run)

```bash
gloom prioritize --genes genes.txt \
  --output results/ \
  --dry-run
```

---

## 📁 Bring Your Own Data: gloom run

`gloom run` trains a **new model from scratch** on **your own expression matrices**. Step 1 (data loading) is bypassed — your files are staged directly into the pipeline cache and preprocessing starts at Step 2.

```bash
gloom run \
  --tumor-expr tumor_expression.csv \
  --normal-expr normal_expression.csv \
  --output results/
```

The flag names `--tumor-expr` and `--normal-expr` reflect the cancer-first design. Conceptually, these correspond to **disease** and **control** expression matrices for any context.

### With metadata files

```bash
gloom run \
  --tumor-expr disease_expr.csv \
  --normal-expr control_expr.csv \
  --tumor-meta disease_meta.csv \
  --normal-meta control_meta.csv \
  --output results/
```

### With custom labels and candidate genes

```bash
gloom run \
  --tumor-expr disease_expr.csv \
  --normal-expr control_expr.csv \
  --genes my_candidates.txt \
  --labels known_genes.tsv \
  --output results/
```

### Full example with all options

```bash
gloom run \
  --tumor-expr tumor.csv \
  --normal-expr normal.csv \
  --genes candidates.txt \
  --labels known_positive_genes.tsv \
  --fdr 0.05 --log2fc 1.0 --prob-threshold 0.5 \
  --top-k 200 --format excel \
  --output my_results/
```

---

## 🔀 Understanding the Two Commands

| | `gloom prioritize` | `gloom run` |
|---|---|---|
| **Data source** | Bundled cBioPortal + GTEx (hardcoded) | Your own `--tumor-expr` / `--normal-expr` |
| **Step 1** | Loads from bundled raw files | Skipped (your files staged directly) |
| **Steps 2-19** | Runs fully | Runs fully |
| **Model training** | Yes, from scratch on bundled data | Yes, from scratch on your data |
| **`--genes` file** | Tags query genes in final ranking | Tags query genes in final ranking |
| **When to use** | Quick analysis using LUAD reference cohort | When you have your own expression data |

**Important:** Both commands train the model from scratch every time. The `--genes` file does **not** influence model training in either command — it only marks which genes are flagged as "query genes" for the novel candidate check at Step 14.

To train on your own data, use `gloom run`. To use the built-in LUAD reference cohort, use `gloom prioritize`.

---

## 📝 Input File Formats

### Query gene file (`--genes`)

One gene symbol per line. Used to tag genes of interest in the ranking output.

```text
EGFR
KRAS
TP53
MET
ALK
MMP11
COL1A1
CXCL9
```

### Expression matrix format (`--tumor-expr` / `--normal-expr`)

CSV with gene symbols as the row index and sample IDs as column headers. Values should be raw expression counts or TPM (the pipeline applies log2 transformation).

```
gene,Sample_1,Sample_2,Sample_3
EGFR,1234.5,2345.6,5678.9
KRAS,567.8,678.9,789.0
TP53,890.1,901.2,123.4
```

### Labels file (`--labels`)

CSV or TSV with a column containing known positive gene symbols. Recognized column names: `GeneSymbol`, `gene_symbol`, `gene`, `symbol`, `hgnc_symbol`.

```
GeneSymbol
EGFR
KRAS
TP53
BRCA1
```

### Metadata files (`--tumor-meta` / `--normal-meta`)

Optional. CSV with sample IDs as the row index. Any additional columns are preserved. If omitted, all samples in the expression matrix are used.

---

## 📂 Output Files

The latest full repository run wrote raw pipeline artifacts under `outputs/`:

```text
outputs/
|-- figures/
|   |-- interactive_dashboard.html
|   |-- interactive_volcano.html
|   |-- interactive_ranking.html
|   |-- interactive_network.html
|   |-- interactive_feature_importance.html
|   |-- interactive_kegg.html
|   |-- eval_roc_curves.png
|   |-- eval_pr_curves.png
|   |-- eval_confusion_matrices.png
|   |-- eval_confusion_multithreshold.png
|   |-- ranking_score_distribution.png
|   |-- pu_bagging_score_distribution.png
|   `-- pipeline-level QC, DE, PCA, feature, and KEGG PNGs
|-- logs/
|   `-- pipeline.log
|-- models/
|   |-- best_model.joblib
|   |-- best_model_calibrated.joblib
|   |-- best_model_name.txt
|   |-- cv_results.csv
|   |-- model_extra_trees.joblib
|   |-- model_gradient_boosting.joblib
|   |-- model_logistic_regression.joblib
|   |-- model_random_forest.joblib
|   |-- model_xgboost.joblib
|   `-- robust_scaler.joblib
`-- results/
    |-- differential_expression_results.csv
    |-- model_metrics.csv
    |-- feature_importance.csv
    |-- gene_rankings.csv
    |-- novel_candidates.csv
    |-- novel_candidates_sensitivity.csv
    |-- query_gene_rankings.csv
    |-- ensemble_soft_vote_probs.csv
    |-- pu_bagging_scores.csv
    |-- pu_bagging_metrics.csv
    |-- ranking_metrics.csv
    |-- classification_reports/
    |   `-- classification_report_*.txt
    |-- enrichment/
    |   |-- kegg_all_candidates.csv
    |   |-- kegg_upregulated.csv
    |   `-- kegg_summary.csv
    |-- network/
    |   |-- coexpression_network.graphml
    |   |-- normal_coexpression_network.graphml
    |   |-- annotated_network.graphml
    |   |-- annotated_nodes.csv
    |   |-- annotated_edges.csv
    |   `-- exports/
    |       |-- network_full.graphml
    |       |-- network_full.gml
    |       |-- network_full_edgelist.tsv
    |       |-- network_cytoscape.json
    |       |-- subnetwork_candidates.graphml
    |       |-- subnetwork_novel.graphml
    |       `-- network_statistics_report.{txt,csv}
    `-- reports/
        |-- pipeline_summary_table.csv
        |-- pipeline_summary_figure.png
        `-- pipeline_summary_figure.pdf
```

If you use the packaged CLI, these raw artifacts can also be reorganized into cleaner user-facing folders such as `candidates/`, `tables/`, `models/`, `plots/`, `network/`, and `reports/`.

### Main output files

| Output | Description |
|---|---|
| `outputs/figures/interactive_dashboard.html` | Consolidated HTML dashboard with ranking, network, enrichment, and evaluation views |
| `outputs/results/gene_rankings.csv` | Genome-wide ranking across all 10,986 scored genes |
| `outputs/results/novel_candidates.csv` | High-confidence non-LCGene candidates (`predicted_prob >= 0.70` with DE support) |
| `outputs/results/novel_candidates_sensitivity.csv` | Broader sensitivity candidate list (`predicted_prob >= 0.50`) |
| `outputs/results/query_gene_rankings.csv` | Query-only ranking when a gene list is supplied; empty in the latest full run |
| `outputs/results/model_metrics.csv` | Held-out AUROC, AUPRC, Brier, fixed-threshold classification metrics, and top-k retrieval metrics |
| `outputs/models/cv_results.csv` | Five-fold CV model comparison used for primary model selection |
| `outputs/results/network/annotated_network.graphml` | Annotated co-expression network for Gephi or Cytoscape |
| `outputs/results/network/exports/network_statistics_report.txt` | Network summary plus top hub genes |
| `outputs/results/enrichment/kegg_summary.csv` | KEGG summary by candidate subset |
| `outputs/results/reports/pipeline_summary_table.csv` | Per-gene summary table used for the final reporting figure |
| `outputs/logs/pipeline.log` | End-to-end execution log with step timings |

### Columns in gene_rankings.csv

| Column | Description |
|---|---|
| `gene` | Gene symbol |
| `predicted_prob` | ML-predicted probability of disease association |
| `rank` | Rank by predicted probability (1 = highest) |
| `percentile` | Percentile rank across all scored genes |
| `predicted_label` | Binary prediction (1 = predicted positive) |
| `label` | Training label (1 = known positive from reference list, 0 = unlabeled) |
| `is_lcgene_gene` | Whether the gene is in the reference label set |
| `log2fc` | Log2 fold-change (disease vs. control) |
| `pvalue_adj` | FDR-adjusted p-value |
| `direction` | DE direction: `up`, `down`, or `ns` (not significant) |
| `is_de_significant` | Whether the gene passed the configured DE significance rule |
| `is_query_gene` | Whether the gene was in the user's `--genes` file |
| `novel_candidate` | High-confidence novel-candidate flag. If no query list is supplied, this is evaluated across the whole ranked universe |
| `novel_candidate_sens` | Sensitivity novel-candidate flag using the broader probability cutoff |

---

## 🤖 Model Performance in the LUAD Case Study

Five machine-learning models plus a calibrated ensemble were compared in the LUAD workflow.

### Cross-Validation Results (5-fold, training set)

| Model | CV AUROC | CV AUPRC | CV F1 | CV MCC |
|---|--:|--:|--:|--:|
| XGBoost | 0.9383 +/- 0.0074 | 0.4951 +/- 0.0147 | 0.3495 +/- 0.0189 | 0.3981 +/- 0.0153 |
| Random Forest | 0.9384 +/- 0.0065 | 0.4831 +/- 0.0052 | 0.4421 +/- 0.0133 | 0.4485 +/- 0.0093 |
| Gradient Boosting | 0.9400 +/- 0.0063 | 0.4796 +/- 0.0319 | 0.4078 +/- 0.0209 | 0.4397 +/- 0.0179 |
| Extra Trees | 0.9368 +/- 0.0073 | 0.4725 +/- 0.0185 | 0.4378 +/- 0.0302 | 0.4406 +/- 0.0369 |
| Logistic Regression | 0.7850 +/- 0.0256 | 0.1498 +/- 0.0310 | 0.1282 +/- 0.0047 | 0.1443 +/- 0.0100 |

### Validation Set Results

| Model | Val AUROC | Val AUPRC | Brier | F1 @ 0.15 | MCC @ 0.15 | Precision@100 | Recall@100 |
|---|--:|--:|--:|--:|--:|--:|--:|
| **Extra Trees** | **0.9330** | **0.5443** | **0.0288** | **0.4955** | **0.4782** | **0.4700** | **0.5341** |
| **Ensemble (soft vote)** | **0.9347** | **0.5392** | **0.0290** | **0.4933** | **0.4761** | **0.5000** | **0.5682** |
| XGBoost | 0.9293 | 0.5288 | 0.0309 | 0.4675 | 0.4508 | 0.4600 | 0.5227 |
| Random Forest | 0.9269 | 0.5113 | 0.0303 | 0.4904 | 0.4688 | 0.4900 | 0.5568 |
| Gradient Boosting | 0.9299 | 0.5077 | 0.0301 | 0.4453 | 0.4319 | 0.4600 | 0.5227 |
| Logistic Regression | 0.8473 | 0.2709 | 0.0374 | 0.3619 | 0.3317 | 0.3200 | 0.3636 |

### Interpretation

- **XGBoost** achieved the best mean cross-validation AUPRC (0.4951) and is the model recorded in `outputs/models/best_model_name.txt`.
- **Extra Trees** achieved the best held-out validation AUPRC (0.5443).
- The **calibrated ensemble** (AUPRC-weighted soft vote) provided the highest held-out AUROC (0.9347) and the highest Precision@100 (0.5000).
- **Step 14 prefers PU bagging scores** from Step 11b for whole-universe ranking when they are available, so the final ranking source is intentionally decoupled from the single-model validation winner.
- **Tree-based models** consistently outperformed Logistic Regression, which is expected given the non-linear feature interactions in co-expression and network topology data.
- The task is inherently difficult due to the **positive-unlabeled (PU) learning** framework: unlabeled genes include both true negatives and hidden positives.
- All models were calibrated using sigmoid calibration on a holdout set.

---

## 🔬 Biological Validation

### LCGene gene recovery

In the starting LUAD workflow, 438 out of 517 LCGene genes were retained after gene harmonization (84.7%). In the latest full ranking:

- Median LCGene rank was **232.5**
- **50** LCGene genes appeared in the top 50 ranks
- **98** LCGene genes appeared in the top 100 ranks
- **190** LCGene genes appeared in the top 200 ranks
- **345** LCGene genes appeared in the top 500 ranks

The top-ranked known genes in the latest run include:

| Gene | Predicted Prob | Rank | Percentile | log2FC | Direction |
|---|--:|--:|--:|--:|---|
| SLCO4A1 | 0.9945 | 1 | 100.0 | 0.32 | ns |
| ID1 | 0.9939 | 2 | 99.99 | 1.80 | ns |
| CLDN18 | 0.9927 | 3 | 99.98 | 0.31 | ns |
| PLA2G2A | 0.9923 | 4 | 99.97 | -0.74 | ns |
| CXCL2 | 0.9920 | 5 | 99.96 | 1.71 | ns |
| MT1G | 0.9919 | 6 | 99.95 | 0.67 | ns |
| TCF21 | 0.9908 | 7 | 99.95 | 0.47 | ns |
| RRAD | 0.9904 | 8 | 99.94 | 1.61 | ns |

### Novel candidate detection

Novel candidates are defined as genes that:
1. **Are not** in the reference label set (LCGene)
2. Have `predicted_prob >= 0.70` for the high-confidence list (`>= 0.50` for the sensitivity list)
3. Show minimum DE signal (`|log2FC| >= 1.0`)
4. Must also be in the query set if a query gene list is supplied

The latest full run was executed **without** a query gene list, so `query_gene_rankings.csv` is empty and novel-candidate detection was applied across the full ranked universe. This produced **203** high-confidence candidates, of which **199** were upregulated and **4** were downregulated.

Top examples from `novel_candidates.csv` include:

| Gene | Predicted Prob | Rank | Percentile | log2FC | Direction |
|---|--:|--:|--:|--:|---|
| CDC20 | 0.9736 | 68 | 99.39 | 6.83 | up |
| MYBL2 | 0.9685 | 93 | 99.16 | 7.39 | up |
| MARCO | 0.9643 | 104 | 99.06 | 2.05 | up |
| RASD1 | 0.9577 | 132 | 98.81 | 4.46 | up |
| CTSE | 0.9475 | 171 | 98.45 | 5.21 | up |
| CHI3L2 | 0.9470 | 172 | 98.44 | 2.31 | up |
| C15ORF48 | 0.9437 | 183 | 98.34 | 6.40 | up |
| TPSAB1 | 0.9426 | 188 | 98.30 | 2.19 | up |

### 🛤️ Enriched KEGG pathways

KEGG pathway enrichment was performed on high-scoring candidate gene sets. The latest run produced the following summary:

| Subset | Input genes | Significant lung/cancer pathways | Top pathway | Top adjusted p-value |
|---|--:|--:|---|--:|
| All candidates | 203 | 4 | ECM-receptor interaction | 0.000500 |
| Upregulated | 199 | 5 | ECM-receptor interaction | 0.000424 |
| Downregulated | 4 | 0 | - | - |

The strongest recurring pathway themes were:

- ECM-receptor interaction
- Focal adhesion
- Complement and coagulation cascades
- Cytokine-cytokine receptor interaction
- PI3K-Akt signaling

No downregulated-subset pathway passed the minimum overlap threshold of 3 genes in the latest run.

---

## ⭐ Feature Importance

The pipeline selects 21 features after collinearity removal. Grouped feature importance shows:

### Top 10 features (normalized importance)

| Rank | Feature | Category | Importance |
|--:|---|---|--:|
| 1 | tumor_std | Tumor expression stats | 0.8097 |
| 2 | tumor_iqr | Tumor expression stats | 0.4568 |
| 3 | normal_std | Normal expression stats | 0.4402 |
| 4 | neg_log10_padj | Differential expression | 0.3814 |
| 5 | tumor_cv | Tumor expression stats | 0.3116 |
| 6 | normal_iqr | Normal expression stats | 0.2341 |
| 7 | cohens_d | Differential expression | 0.2292 |
| 8 | abs_log2fc | Differential expression | 0.1589 |
| 9 | normal_cv | Normal expression stats | 0.1525 |
| 10 | normal_skewness | Normal expression stats | 0.1492 |

### Feature group contributions

| Feature Group | Role |
|---|---|
| Tumor / disease expression statistics | Primary signal (~50% of total importance) |
| Normal / control expression statistics | Strong complementary signal |
| Differential expression metrics | Core discriminative features |
| Network topology (degree, centrality) | Biological refinement signal |
| Network edge weight statistics | Supplementary co-expression context |
| Differential network features | Tumor-vs-normal network rewiring |

The LUAD workflow is mainly expression-driven, but network features help refine interpretation and support biologically coherent ranking.

---

## ⚠️ Limitations

GLOOM provides a prioritized candidate list, not a final set of validated biomarkers.

| Limitation | Explanation |
|---|---|
| Label scope | LCGene captures expression-based biomarkers, not mutation-driven oncogenes. Users should provide custom labels when studying driver genes. |
| Positive-unlabeled setting | Some genes treated as negative may actually be hidden positives |
| Class imbalance | Positive rate is ~4.5% — AUROC can appear high while positive-class recall remains challenging |
| Batch effects | Disease and control datasets from different sources may contain residual technical differences |
| Network contribution | Network features currently act more as refinement signals than dominant predictors |
| Gene coverage | Only genes present in both tumor and normal expression matrices are scored |
| Experimental validation | Novel candidates require laboratory and clinical validation |

---

## 🔮 Future Work

Planned and recommended extensions include:

- 🌍 Apply GLOOM to additional cancer types and non-cancer disease contexts
- 🧬 Integrate methylation, mutation, copy-number, and proteomics data
- 🧠 Add graph neural network embeddings
- 🔄 Test random-walk and network-diffusion features
- 📏 Benchmark against other gene prioritization tools
- 🎯 Improve positive-unlabeled learning strategies
- 🧫 Experimentally validate selected novel candidates
- 🎨 Expand visualization and reporting options

---

## 🗂️ Repository Structure

```text
gloom/
├── data/
│   └── raw/
│       ├── cBioPortal (RNA Seq Data)/             # TCGA LUAD expression + clinical metadata
│       ├── Gtex (normal samples)/                 # GTEx lung expression + sample attributes
│       └── LCGene (Labeled LUAD Data)/            # LCGene curated LUAD gene list
├── src/
│   └── gloom/
│       ├── __init__.py
│       ├── cli.py                                 # CLI entry point (prioritize, run, info, etc.)
│       ├── output.py                              # Output organization and formatting
│       └── pipeline/
│           ├── config.py                          # Configuration and path management
│           ├── run_pipeline.py                    # Step registry and pipeline runner
│           ├── step1_data_loading.py
│           ├── step1b_batch_correction.py         # Optional batch-correction refinement
│           ├── step2_preprocessing.py
│           ├── step3_harmonization.py
│           ├── step4_differential_expression.py
│           ├── step5_expression_features.py
│           ├── step6_coexpression_network.py
│           ├── step6b_normal_network.py           # Optional normal-network reference
│           ├── step7_network_features.py
│           ├── step7b_differential_network_features.py  # Optional tumor-vs-normal rewiring features
│           ├── step8_feature_integration.py
│           ├── step9_label_construction.py
│           ├── step10_train_val_split.py
│           ├── step11_model_training.py
│           ├── step11b_pu_bagging.py              # Optional PU bagging refinement
│           ├── step12_model_evaluation.py
│           ├── step13_feature_importance.py
│           ├── step14_gene_ranking.py
│           ├── step15_network_annotation.py
│           ├── step16_network_export.py
│           ├── step17_interactive_visualization.py
│           ├── step18_final_report.py
│           └── step19_kegg_enrichment.py
├── tests/
│   └── genes.txt                                  # Example query gene list
├── pyproject.toml
├── environment.yml
├── LICENSE
└── README.md
```

---

## 🎯 Example Use Case

### Scenario 1: Using bundled LUAD reference data

A researcher has a list of candidate genes and wants to see how they rank in the LUAD expression landscape:

```bash
gloom prioritize --genes my_candidates.txt \
  --disease luad \
  --output luad_results/
```

### Scenario 2: Using your own expression data

A researcher has disease and control expression matrices and wants to discover biologically relevant candidate genes:

```bash
gloom run \
  --tumor-expr disease_expression.csv \
  --normal-expr control_expression.csv \
  --labels known_disease_genes.csv \
  --output disease_results/
```

### What GLOOM produces

Both commands generate:

1. Differential expression results with fold-change and FDR
2. Co-expression network (GraphML + Cytoscape)
3. Integrated feature matrix (expression + network)
4. Trained and calibrated ML models (core models, optional XGBoost when available, and ensemble outputs)
5. Ranked gene list with predicted probabilities
6. Novel candidate list (high-confidence genes not in reference set)
7. KEGG pathway enrichment
8. Interactive HTML report and dashboard
9. Annotated network files for visualization

---

## 👥 Contributors

| Names | Affiliation(s) | Role(s) |
|---|---|---|
| Rahma Yasser Mahmoud | Faculty of Computers and Information, Assiut University, Assiut, Egypt | Team Leader |
| Khadija Adam Rogo | Department of Bioinformatics, Kalinga University, Raipur, India | |
| Rana Hamed Abu-Zeid | Department of Artificial Intelligence, Badya University, Giza, Egypt | |
| Malick Traore | African Center of Excellence in Bioinformatics and Data Science, USTTB, Mali | |
| Olaitan I. Awe | Institute for Genomic Medicine Research (IGMR); African Society for Bioinformatics and Computational Biology (ASBCB) | Project Advisor |

<br>📧 **Rahma Yasser Mahmoud**: rahmayasserm@gmail.com
<br>📧 **Rana Hamed Abu-Zeid**: ranahamed2111@gmail.com
<br>📧 **Khadija Adam Rogo**: khadijarogo212@gmail.com
<br>📧 **Malick Traore**: malicktra100@gmail.com
<br>📧 **Olaitan I. Awe, Ph.D.**: laitanawe@gmail.com

---

## 📚 Acknowledgments

We thank the contributors to open biomedical datasets (TCGA, GTEx, LCGene) and acknowledge the Python and bioinformatics open-source communities. Special thanks to all collaborators and advisors in LUAD and genomics research.

---

## ✅ Summary

GLOOM provides a transparent and reusable package for moving from gene expression data to biologically meaningful candidate-gene rankings. The current workflow supports optional batch correction, normal-network comparison, PU bagging, KEGG enrichment, and interactive reporting.

In the latest full LUAD run in this repository (August 14, 2026), the pipeline ranked **10,986** genes, recovered **98** LCGene genes within the top 100 positions, and surfaced **203** high-confidence novel candidates in about **12.2 minutes**.

It offers two main modes of operation:

- **`gloom prioritize`** for quick analysis using the bundled LUAD reference cohort
- **`gloom run`** for training on user-provided expression data

By combining expression statistics, differential expression, co-expression network topology, machine learning (with PU learning), enrichment analysis, and interactive reporting, GLOOM supports reproducible gene prioritization and downstream experimental hypothesis generation.
