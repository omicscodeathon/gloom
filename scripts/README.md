# LUAD ML Pipeline Scripts

This directory contains the executable LUAD pipeline for data loading, preprocessing, differential expression, network construction, machine learning, reporting, enrichment, and optional external validation.

The old README described an earlier 18-step version. The current code exposes 20 logical steps through `run_pipeline.py`, including optional modules that were added later:

- `step1b_batch_correction.py`
- `step6b_normal_network.py`
- `step7b_differential_network_features.py`
- `step11b_pu_bagging.py`
- `step19_kegg_enrichment.py`
- `step15_independent_validation.py` (used as step 20 by the runner)

## Support files

- `config.py` - central paths, thresholds, flags, and output locations
- `run_pipeline.py` - master runner for full runs or selected step ranges
- `requirements.txt` - base Python dependencies for this folder
- `__init__.py` - package marker

## Pipeline map

| Key | Script | Purpose | Optional |
|---|---|---|---|
| `0` | `config.py` | Create output directories and load central settings | No |
| `1` | `step1_data_loading.py` | Load cBioPortal, GTEx, and LCGene raw inputs into validated CSVs | No |
| `1b` | `step1b_batch_correction.py` | Apply optional batch correction before QC and DE | Yes |
| `2` | `step2_preprocessing.py` | Log-transform, filter, QC, and align expression matrices and metadata | No |
| `3` | `step3_harmonization.py` | Harmonize gene symbols and align tumor and normal matrices | No |
| `4` | `step4_differential_expression.py` | Run differential expression with Welch t-test, BH FDR, and effect size | No |
| `5` | `step5_expression_features.py` | Build per-gene expression and contrast features | No |
| `6` | `step6_coexpression_network.py` | Build the tumor co-expression network | No |
| `6b` | `step6b_normal_network.py` | Build the normal-tissue co-expression network | Yes |
| `7` | `step7_network_features.py` | Extract tumor network topology features | No |
| `7b` | `step7b_differential_network_features.py` | Compute tumor-vs-normal differential network features | Yes |
| `8` | `step8_feature_integration.py` | Merge, clean, and scale expression and network features | No |
| `9` | `step9_label_construction.py` | Build LCGene positive vs unlabeled labels and gene annotations | No |
| `10` | `step10_train_val_split.py` | Create the train/validation split and CV fold assignments | No |
| `11` | `step11_model_training.py` | Train core ML models and save the calibrated best model | No |
| `11b` | `step11b_pu_bagging.py` | Run Mordelet-Vert PU bagging for ranking genes in a PU setting | Yes |
| `12` | `step12_model_evaluation.py` | Evaluate trained models on the held-out validation set | No |
| `13` | `step13_feature_importance.py` | Compute native and permutation feature importance | No |
| `14` | `step14_gene_ranking.py` | Score all genes and export ranked novel candidates | No |
| `15` | `step15_network_annotation.py` | Annotate network nodes and edges with ML and DE results | No |
| `16` | `step16_network_export.py` | Export annotated network in multiple graph and tabular formats | No |
| `17` | `step17_interactive_visualization.py` | Build interactive HTML visualizations and a combined dashboard | No |
| `18` | `step18_final_report.py` | Write summary tables, text report, and summary figures | No |
| `19` | `step19_kegg_enrichment.py` | Run KEGG enrichment on high-confidence candidate sets | Yes |
| `20` | `step15_independent_validation.py` | Score an external cohort and report independent validation metrics | Yes |

Note: step 20 is implemented in `step15_independent_validation.py`. The filename is older than the current runner numbering.

## Required input data

The current `config.py` expects these raw inputs:

```text
data/raw/
  cBioPortal (RNA Seq Data)/
    data_mrna_seq_v2_rsem.txt
    data_clinical_patient.txt
  Gtex (normal samples)/
    gene_tpm_v11_lung.gct.gz
    GTEx_Analysis_v11_Annotations_SampleAttributesDS - LUAD.txt
  LCGene (Labeled LUAD Data)/
    LCGene_human_LUAD_filtered.tsv
```

Optional external validation inputs are configured in `config.py`:

- `EXTERNAL_EXPR_FILE`
- `EXTERNAL_LABELS_FILE`

## Installation

From the repository root:

```bash
pip install -r scripts/requirements.txt
```

Base requirements include `pandas`, `numpy`, `scipy`, `scikit-learn`, `networkx`, `matplotlib`, `statsmodels`, `plotly`, `joblib`, `openpyxl`, and `gseapy`.

Optional extras used by some paths:

- `pip install xgboost` for the optional XGBoost model in step 11
- `pip install imbalanced-learn` if `USE_SMOTE` or undersampling is enabled
- `pip install inmoose` for the Python ComBat-seq path in step 1b
- `pip install rpy2` plus an R installation with `sva` if you want the R ComBat-seq path instead

## Running the pipeline

From the repository root:

```bash
python scripts/run_pipeline.py
python scripts/run_pipeline.py --from 4
python scripts/run_pipeline.py --from 6 --to 8
python scripts/run_pipeline.py --only 12
python scripts/run_pipeline.py --only 11b
python scripts/run_pipeline.py --only 20
```

From inside `scripts/`:

```bash
cd scripts
python run_pipeline.py
```

All step files also expose a direct entrypoint, so ad hoc reruns are supported:

```bash
python scripts/step11_model_training.py
python scripts/step11b_pu_bagging.py
python scripts/step19_kegg_enrichment.py
```

## Current defaults worth knowing

These defaults are set in the current `config.py`:

- `DE_LOG2FC_THRESHOLD = 2.0`
- `DE_PVALUE_THRESHOLD = 0.001`
- `COEXPR_CORRELATION_CUTOFF = 0.60`
- `CV_FOLDS = 5`
- `NOVEL_PROB_THRESHOLD = 0.85`
- `NOVEL_PROB_THRESHOLD_SENS = 0.50`
- `USE_BATCH_CORRECTION = False`
- `EXTERNAL_EXPR_FILE = None`
- `EXTERNAL_LABELS_FILE = None`

Labeling is based on the LCGene LUAD file, not the older Cancer Gene Census input referenced by the previous README.

## Models currently trained

The implemented model set in `step11_model_training.py` is:

- Random Forest
- Gradient Boosting
- Logistic Regression
- Extra Trees
- XGBoost, only if the package is installed

`step11_model_training.py` still mentions SVM in its top docstring, but SVM is not present in the active model definitions anymore.

## Key outputs

Exact base directories are controlled by `config.py`, but the current scripts write the following logical outputs.

Processed and intermediate tables:

- `tumor_expression_raw.csv`
- `normal_expression_raw.csv`
- `tumor_expression_processed.csv`
- `normal_expression_processed.csv`
- `expression_features.csv`
- `network_features.csv`
- `differential_network_features.csv` when step `7b` is available
- `integrated_features.csv`
- `integrated_features_scaled.csv`
- `gene_labels.csv`
- `gene_annotation_table.csv`

Modeling and ranking outputs:

- `models/cv_results.csv`
- `results/model_metrics.csv`
- `results/feature_importance.csv`
- `results/gene_rankings.csv`
- `results/novel_candidates.csv`
- `results/novel_candidates_sensitivity.csv`
- `results/query_gene_rankings.csv`
- `results/pu_bagging_scores.csv` when step `11b` is run
- `results/pu_bagging_metrics.csv` when step `11b` is run

Network outputs:

- `results/network/coexpression_network.graphml`
- `results/network/normal_coexpression_network.graphml` when step `6b` is run
- `results/network/annotated_network.graphml`
- `results/network/annotated_nodes.csv`
- `results/network/annotated_edges.csv`
- export files from step 16 including GraphML, GML, TSV, and Cytoscape JSON

Interactive outputs from step 17:

- `figures/interactive_volcano.html`
- `figures/interactive_ranking.html`
- `figures/interactive_network.html`
- `figures/interactive_feature_importance.html`
- `figures/interactive_novel_candidates.html`
- `figures/interactive_kegg.html` when KEGG data is present
- `figures/interactive_dashboard.html`

Reporting and validation outputs:

- `results/reports/pipeline_summary_table.csv`
- `results/reports/pipeline_report.txt`
- `results/reports/pipeline_summary_figure.png`
- `results/reports/pipeline_summary_figure.pdf`
- `results/enrichment/kegg_all_candidates.csv`
- `results/enrichment/kegg_upregulated.csv`
- `results/enrichment/kegg_downregulated.csv`
- `results/enrichment/kegg_summary.csv`
- `results/independent_validation_scores.csv`
- `results/independent_validation_metrics.csv`

## Notes on optional steps

- Step `1b` only does work when `USE_BATCH_CORRECTION=True`.
- After step `1b`, `run_pipeline.py` checks whether the harmonized expression paths still point at the uncorrected files and warns if you need to reroute them before resuming from step `4`.
- Step `8` auto-includes differential network features if the `step7b` output file exists.
- Step `11b` is the main PU-learning path for whole-universe gene ranking.
- Step `19` depends on `novel_candidates.csv`, so step `14` must run first.
- Step `20` requires `EXTERNAL_EXPR_FILE`, and uses `EXTERNAL_LABELS_FILE` when provided.
