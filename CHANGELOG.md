# Changelog

All notable changes to **gloom** are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Version numbers follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Planned
- Additional disease models beyond LUAD
- Support for single-cell RNA-seq input
- Docker / Singularity container recipes

---

## [0.1.0] - 2026-04-25

### Added

#### CLI commands
- `gloom prioritize` — run the full 20-step ML pipeline for LUAD gene prioritization
  - Options: `--genes`, `--disease`, `--output`, `--data-dir`, `--from-step`, `--to-step`,
    `--top-k`, `--format`, `--fdr`, `--log2fc`, `--prob-threshold`, `--no-cache`,
    `--labels`, `--dry-run`, `--verbose`
- `gloom run` — run the pipeline with user-supplied expression data (bypasses Step 1)
  - Options: `--tumor-expr`, `--normal-expr`, `--tumor-meta`, `--normal-meta`,
    `--output`, `--genes`, `--labels`, `--from-step`, `--to-step`,
    `--fdr`, `--log2fc`, `--prob-threshold`, `--top-k`, `--format`, `--verbose`
- `gloom validate` — check that all required raw input files are present
- `gloom info` — display version, reference file status, and default thresholds
- `gloom diseases` — list supported disease contexts and their data sources
- `gloom cache clear` — delete cached intermediate files to force a full re-run

#### Pipeline (20 steps)
- **Step 1** — Data loading: TCGA/cBioPortal tumor expression, GTEx normal expression,
  clinical metadata, LCGene known LUAD gene database
- **Step 2** — Preprocessing: QC filtering, log2 transformation, low-expression removal
- **Step 3** — Harmonization: gene symbol unification, sample overlap validation
- **Step 4** — Differential expression: Welch t-test, Benjamini-Hochberg FDR correction,
  log2 fold-change thresholding, volcano plot
- **Step 5** — Expression feature construction: 15 per-gene features (mean, variance,
  log2FC rank, DE significance scores, etc.)
- **Step 6** — Co-expression network: Pearson correlation, configurable cutoff,
  GraphML export
- **Step 7** — Network feature extraction: degree, betweenness, closeness, clustering
  coefficient, eigenvector centrality, PageRank
- **Step 8** — Feature integration: expression + network features, RobustScaler,
  PCA visualization, high-correlation flagging
- **Step 9** — Label construction: positive-unlabeled (PU) framing, LCGene positives,
  class balance visualization
- **Step 10** — Train/validation split: stratified split, gene-level holdout
- **Step 11** — Model training: RandomForest, GradientBoosting, LogisticRegression,
  ExtraTrees; stratified k-fold CV; optional SMOTE and undersampling
- **Step 12** — Model evaluation: AUROC, AUPRC, F1, MCC, Brier score; ROC and PR curves
- **Step 13** — Feature importance: permutation importance, SHAP values, grouped ranking
- **Step 14** — Gene ranking: full genome scoring, novel candidate identification,
  configurable probability threshold
- **Step 15** — Network annotation: overlay ML scores onto co-expression graph
- **Step 16** — Network export: GraphML + Cytoscape-compatible XML
- **Step 17** — Interactive visualization: Plotly-based tabbed dashboard
  (Volcano, Gene Ranking, Network, Feature Importance, Novel Candidates, KEGG Pathways)
- **Step 18** — Final report: pipeline summary CSV
- **Step 19** — KEGG pathway enrichment: Enrichr REST API query, lung/cancer pathway
  filter, bar plots, dot plots, combined overview

#### Output structure
- `candidates/` — ranked_candidates, novel_candidates (CSV / Excel / JSON)
- `tables/` — DE results, network features
- `models/` — best_model.joblib, model_card.json
- `plots/` — interactive HTML: volcano, feature importance, ROC/PR curves, KEGG pathways
- `network/` — coexpression.graphml, coexpression.cytoscape.xml
- `report.html` — self-contained interactive dashboard (offline-ready)

#### Packaging
- `src/` layout with `pyproject.toml` (PEP 517/518)
- Optional dependency groups: `interactive`, `kegg`, `excel`, `resampling`, `full`, `dev`
- `conda.recipe/meta.yaml` for conda-forge / Anaconda distribution
- `environment.yml` for conda environment setup
- MIT License

---

[Unreleased]: https://github.com/omicscodeathon/gloom/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/omicscodeathon/gloom/releases/tag/v0.1.0
