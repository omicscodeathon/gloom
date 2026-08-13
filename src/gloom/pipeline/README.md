# LUAD ML Pipeline

Machine learning pipeline for LUAD candidate-gene prioritization, co-expression analysis, and downstream reporting.

## Overview

The pipeline currently has 19 numbered steps plus optional branch steps:

- Required path: `0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18`
- Optional steps: `1b, 6b, 7b, 11b, 19`

Important optional-step behavior:

- Optional means the runner continues if the step fails.
- Optional outputs can still affect later steps if their files already exist.
- `step11b` is especially important because `step14` prefers `pu_bagging_scores.csv` when it is present.

## Project Structure

```text
src/gloom/pipeline/
|-- config.py
|-- run_pipeline.py
|-- step1_data_loading.py
|-- step1b_batch_correction.py          # optional
|-- step2_preprocessing.py
|-- step3_harmonization.py
|-- step4_differential_expression.py
|-- step5_expression_features.py
|-- step6_coexpression_network.py
|-- step6b_normal_network.py            # optional
|-- step7_network_features.py
|-- step7b_differential_network_features.py  # optional
|-- step8_feature_integration.py
|-- step9_label_construction.py
|-- step10_train_val_split.py
|-- step11_model_training.py
|-- step11b_pu_bagging.py               # optional
|-- step12_model_evaluation.py
|-- step13_feature_importance.py
|-- step14_gene_ranking.py
|-- step15_network_annotation.py
|-- step16_network_export.py
|-- step17_interactive_visualization.py
|-- step18_final_report.py
|-- step19_kegg_enrichment.py           # optional
`-- requirements.txt
```

## Required Input Data

Place the raw files under `data/raw/` relative to the `trial/` project root:

```text
data/raw/
|-- cBioPortal (RNA Seq Data)/
|   |-- data_mrna_seq_v2_rsem.txt
|   `-- data_clinical_patient.txt
|-- Gtex (normal samples)/
|   |-- gene_tpm_v11_lung.gct.gz
|   `-- GTEx_Analysis_v11_Annotations_SampleAttributesDS - LUAD.txt
`-- LCGene (Labeled LUAD Data)/
    `-- LCGene_human_LUAD_filtered.tsv
```

## Installation

For package users, install GLOOM from Bioconda:

```bash
conda install -c bioconda -c conda-forge gloom
```

For direct script usage from the repository root:

```bash
pip install -r src/gloom/pipeline/requirements.txt
```

Optional dependencies:

- `xgboost` for the XGBoost model in step 11
- `plotly` for step 17 interactive HTML outputs
- `gseapy` for step 19 KEGG enrichment
- `inmoose` or `rpy2 + R sva` for step 1b batch correction

## Usage

Run from inside `src/gloom/pipeline/`, or use the full paths from the repository root.

### Full run

```bash
python run_pipeline.py
```

### Full run without optional steps

```bash
python run_pipeline.py --skip-optional
```

### Resume from a step

```bash
python run_pipeline.py --from 4
```

### Run a single step

```bash
python run_pipeline.py --only 11
python run_pipeline.py --only 11b
python run_pipeline.py --only 19
```

### Run step scripts directly

```bash
python step11_model_training.py
python step11b_pu_bagging.py
python step19_kegg_enrichment.py
```

Equivalent commands from the repository root:

```bash
python src/gloom/pipeline/run_pipeline.py
python src/gloom/pipeline/run_pipeline.py --skip-optional
python src/gloom/pipeline/run_pipeline.py --from 4
python src/gloom/pipeline/run_pipeline.py --only 11b
python src/gloom/pipeline/step11_model_training.py
```

## Step Summary

- `1b`: optional batch-effect correction before downstream expression analysis
- `6b`: optional normal-tissue co-expression network
- `7b`: optional tumor-vs-normal differential network features
- `11`: supervised model training and best-model selection
- `11b`: PU bagging over the full gene universe
- `19`: KEGG enrichment for ranked candidate genes

## Outputs

All generated outputs live outside the code folder under `outputs/` at the `trial/` root:

```text
outputs/
|-- figures/
|-- logs/
|-- models/
`-- results/
    |-- enrichment/
    |-- network/
    `-- reports/
```

Key files:

- `outputs/results/gene_rankings.csv`
- `outputs/results/novel_candidates.csv`
- `outputs/results/network/annotated_network.graphml`
- `outputs/results/enrichment/kegg_all_candidates.csv`
- `outputs/results/enrichment/kegg_upregulated.csv`
- `outputs/results/enrichment/kegg_downregulated.csv`
- `outputs/results/enrichment/kegg_summary.csv`
- `outputs/results/reports/pipeline_summary_table.csv`
- `outputs/results/reports/pipeline_report.txt`
- `outputs/results/reports/pipeline_summary_figure.png`
- `outputs/figures/interactive_dashboard.html`
- `outputs/models/best_model.joblib`
- `outputs/models/cv_results.csv`

Additional optional outputs:

- `outputs/results/pu_bagging_scores.csv`
- `outputs/results/pu_bagging_metrics.csv`
- `outputs/figures/pu_bagging_score_distribution.png`
- `outputs/figures/interactive_kegg.html`

## Models Used

### Step 11: supervised training

- Random Forest
- Gradient Boosting
- Logistic Regression
- Extra Trees
- XGBoost, if `xgboost` is installed

Hyperparameter tuning uses `RandomizedSearchCV` when enabled in `config.py`.
The current tuning grid includes:

- Random Forest
- Gradient Boosting
- Logistic Regression
- Extra Trees
- XGBoost

### Step 11b: PU bagging

- Random Forest base learners inside a positive-unlabeled bagging ensemble

## Current Modeling Notes

- `SVM` is not part of the current implemented model list.
- `step14` prefers `pu_bagging_scores.csv` from `step11b` when available.
- The primary CV selection metric is currently `AUPRC`, not `AUROC`.
- Tuning is currently configured in `config.py` with:
  - `USE_HYPERPARAMETER_TUNING = True`
  - `HP_N_ITER = 20`

## Configuration Highlights

Main settings live in `config.py`.

Commonly edited options:

- `COEXPR_CORRELATION_CUTOFF`
- `CV_FOLDS`
- `DE_LOG2FC_THRESHOLD`
- `DE_PVALUE_THRESHOLD`
- `USE_BATCH_CORRECTION`
- `USE_HYPERPARAMETER_TUNING`
- `HP_N_ITER`
- `USE_SMOTE`
- `USE_UNDERSAMPLING`
- `PU_FRAMING`
- `PU_N_ESTIMATORS`
- `PU_SUBSAMPLE_RATIO`
- `PU_BASE_N_TREES`

## Notes

- On Windows, the pipeline now configures console output to avoid `UnicodeEncodeError` from characters that the active code page cannot encode.
- Log files are written with UTF-8 encoding.
- `outputs/results/enrichment/kegg_summary.csv` records subset-specific status and explanatory notes when a KEGG subset returns empty after thresholds or pathway filtering.
- If you skip optional steps, remove stale optional output files if you want a strictly optional-free run.
