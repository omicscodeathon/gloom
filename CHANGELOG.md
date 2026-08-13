# Changelog

All notable changes to `gloom` are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Version numbers follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Added
- `gloom` is now available on Bioconda.

### Changed
- Updated installation documentation to recommend `conda install -c bioconda -c conda-forge gloom`.

### Planned
- Additional disease models beyond LUAD
- Support for single-cell RNA-seq input
- Docker / Singularity container recipes

---

## [0.1.0] - 2026-04-25

### Added

#### CLI commands
- `gloom prioritize` runs the full core LUAD pipeline using bundled reference data.
  Options: `--genes`, `--disease`, `--output`, `--data-dir`, `--from-step`, `--to-step`, `--skip-optional`, `--skip-step`, `--top-k`, `--format`, `--fdr`, `--log2fc`, `--prob-threshold`, `--no-cache`, `--labels`, `--dry-run`, `--verbose`
- `gloom run` runs the pipeline with user-supplied expression data and bypasses Step 1 data loading.
  Options: `--tumor-expr`, `--normal-expr`, `--tumor-meta`, `--normal-meta`, `--output`, `--genes`, `--labels`, `--from-step`, `--to-step`, `--skip-optional`, `--skip-step`, `--fdr`, `--log2fc`, `--prob-threshold`, `--top-k`, `--format`, `--verbose`
- `gloom validate` checks that all required raw input files are present.
- `gloom info` displays version, reference file status, and default thresholds.
- `gloom diseases` lists supported disease contexts and their data sources.
- `gloom cache clear` deletes cached intermediate files to force a full re-run.

#### Pipeline
- 20 core stages from Step 0 through Step 19.
- Optional refinement stages `1b`, `6b`, `7b`, and `11b` can be included to improve robustness and, when appropriate, ranking quality.
- Step 1 loads TCGA/cBioPortal tumor expression, GTEx normal expression, metadata, and LCGene labels.
- Step 1b optionally performs batch correction before downstream analysis.
- Step 6b optionally builds a normal/control co-expression reference network.
- Step 7b optionally adds tumor-vs-normal differential network features.
- Step 11 trains core models with cross-validation, calibration, optional resampling, and XGBoost when installed.
- Step 11b optionally performs PU bagging as a supplemental ranking refinement.
- Step 13 computes native model importance and permutation importance.
- Step 18 generates reporting artifacts including a summary table, text report, summary figure, and the dashboard-backed `report.html`.
- Step 19 performs KEGG pathway enrichment.

#### Output structure
- `candidates/` contains ranked candidates, novel candidates, and full gene rankings.
- `tables/` contains DE results, expression features, network features, integrated features, feature importance, and model metrics.
- `models/` contains the best model, supporting model artifacts, and model card metadata.
- `plots/` contains dashboard and visualization HTML outputs.
- `network/` contains `annotated_coexpression.graphml`, `annotated_coexpression.cytoscape.xml`, and network exports.
- `reports/` contains `pipeline_report.txt`, `pipeline_summary_table.csv`, and `pipeline_summary_figure.png`.
- `report.html` provides a self-contained interactive dashboard.

#### Packaging
- `src/` layout with `pyproject.toml` (PEP 517/518)
- Optional dependency groups: `interactive`, `kegg`, `excel`, `resampling`, `full`, `dev`
- `conda.recipe/meta.yaml` for conda packaging
- `environment.yml` for conda environment setup
- MIT License

---

[Unreleased]: https://github.com/omicscodeathon/gloom/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/omicscodeathon/gloom/releases/tag/v0.1.0
