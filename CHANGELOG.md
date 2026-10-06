# Changelog

All notable changes to `gloom` are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Version numbers follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Planned
- Quantitative head-to-head benchmark against Endeavour, ToppGene, GeneMANIA and BioRank
- Publication of a Bioconda recipe (not yet available)
- Additional disease models beyond LUAD
- Support for single-cell RNA-seq input
- Docker / Singularity container recipes

---

## [0.2.0] - 2026-10-06 (revision release)

Major revision in response to peer review. **Results produced by 0.1.x are not comparable** with
0.2.0 (zero handling, gene universe, evaluation protocol and the default ranking all changed) and
the LUAD case study must be re-run. Release date to be confirmed when the version is tagged.

### Fixed (analysis correctness)
- **Zero handling (step 2).** True zeros are kept (`log2(0 + 1) = 0`); NaN stays missing and is counted;
  negative / infinite values are marked invalid, set to NaN and reported with a warning. The v0.1.x call
  `df.where(df > 0, NaN)` turned genuine zeros into missing values and biased means, DE statistics and correlations.
  The low-expression filter now divides by the number of *observed* samples; the correlation z-score (step 6) is NaN-safe.
- **Label-dependent gene universe (step 9).** New `LABEL_INDEPENDENT_UNIVERSE = True` (default): the universe is defined
  before labelling with identical criteria for all genes and no |log2FC| filter is applied to unlabeled genes.
  The old behaviour (remove unlabeled genes with |log2FC| < 2, keep all positives) is reachable with
  `LABEL_INDEPENDENT_UNIVERSE = False` and is DEPRECATED (label-dependent, not recommended).
- **Training positives in the ranking (steps 11b, 14).** `pu_bagging_metrics.csv` now reports top-K counts for held-out
  genes (`lcgene_top*_heldout`, precision/recall on held-out genes only) separately from training positives
  (`lcgene_top*_training`, resubstitution). `gene_rankings.csv` gains `full_model_score` and `is_training_positive`.
- **Enrichment (step 19).** `kegg_all_candidates.csv` is now the complete, truly unfiltered result table (explicit background =
  analysis universe, hypergeometric test, Benjamini-Hochberg over all tested terms; columns `term`, `n_query`,
  `n_term_in_background`, `n_overlap`, `p_raw`, `p_adj`). The lung/cancer keyword filter was removed from the primary output;
  it only produces the interpretive `kegg_lung_cancer_subset.csv` (`ENRICHMENT_THEMATIC_SUBSET`).
- Step 18 text report no longer fails on the model-metrics column names written by step 12.

### Added
- **Cross-fitting (step 11c, `USE_CROSSFIT`).** Repeated stratified K-fold (default 5 x 5) PU cross-fitting; the out-of-fold
  scores (`results/oof_scores.csv`) are the PRIMARY ranking used by step 14.
- **Out-of-fold metrics (step 12b).** AUROC, AUPRC, average precision, Precision/Recall/enrichment factor @10/50/100 with
  1000-resample stratified bootstrap 95% CIs (`results/oof_metrics.csv`).
- **Feature-set ablation (step 13b).** |log2FC|, adjusted p-value, expression-only, network-only, expression + network,
  combined without DE statistics and the full model under the same cross-fitting, paired-bootstrap comparison and a
  plain-text verdict (`results/ablation_metrics.csv`, `ablation_paired_comparison.csv`, `ablation_verdict.txt`).
  Feature groups are defined in `config.py` (`FEATURE_GROUP_*`).
- **Network stability (step 7c).** Threshold sweep |r| 0.50-0.70 (edges, density, mean degree, largest component, clustering,
  isolated nodes, degree-centrality Spearman between thresholds), optional equal-N sub-sampling (`NETWORK_EQUAL_N`) and a sample
  bootstrap (edge retention, hub Jaccard) -> `results/network_stability_*.csv`.
- **Cohort design checks.** `COHORT_DESIGN` check in step 1 (perfect collinearity between cohort and group -> loud warning and
  `results/qc_cohort_warning.txt`); new step 2b with per-sample median/IQR tables, sample PCA, global median log2FC,
  fraction of DE genes, up:down ratio and automatic warnings (|median log2FC| > 1, > 80 % DE, up:down > 10 or < 0.1).
- **Uniformly processed TCGA design.** `DATA_SOURCE = "tcga_gdc"`, `scripts/fetch_gdc_tcga_luad.py` (GDC STAR counts of
  TCGA-LUAD primary tumors and adjacent normals), `docs/tcga_paired_design.md`.
- **Differential-expression options (step 4).** `DE_METHOD = "welch" | "paired" | "limma_voom"` (limma-voom through optional
  rpy2; skipped with a message when unavailable).
- **External evidence for candidates (step 20, optional).** Open Targets GraphQL + Europe PMC check of the top non-LCGene
  genes with retries and a local cache -> `results/candidate_evidence.csv`; classes `already_LUAD`, `lung_cancer_unspecified`,
  `other_cancer`, `indirect_only`, `no_association_found` (+ `not_assessed`).
- **Reproducibility.** Central `SEED` in `config.py` (`RANDOM_STATE` kept as alias) used by every estimator/split/bootstrap/layout;
  `config.set_global_seed()`; `scripts/benchmark_runtime.py` (run time and peak memory vs genes x samples, machine info).
- `gloom_utils.py` with pure helper functions and pytest unit tests in `test/` (zero handling, cross-fitting, metrics and
  bootstrap CIs, enrichment BH + background, QC verdicts, pairing, network summary, evidence classes).
- `CITATION.cff`, `docs/RELEASE_CHECKLIST.md` (tag, Zenodo archive, DOI), `docs/revision/` (response-to-reviewers draft,
  manuscript replacement text, re-run guide).
- `run_pipeline.py` accepts unpadded step keys (`--from 4`) and the new keys 2b, 7c, 11c, 12b, 13b, 20; the CLI accepts them too
  (`--to-step 20` for the optional evidence step).

### Changed
- **Terminology: "novel candidates" -> "non-LCGene candidates"** (absence from LCGene does not establish novelty). Renamed
  `novel_candidates*.csv` -> `non_lcgene_candidates*.csv`, columns `novel_candidate*` -> `non_lcgene_candidate*`, config
  `NOVEL_PROB_THRESHOLD*` -> `CANDIDATE_PROB_THRESHOLD*` (old names kept as deprecated aliases), subnetwork `novel` -> `non_lcgene`,
  and all labels in the dashboards, step 18 report, CLI and README.
- Step 14 ranks by out-of-fold scores (`predicted_prob`); the full-data model score is in `full_model_score`; top-K metrics are
  computed on out-of-fold scores (or on held-out genes only when cross-fitting is disabled).
- `USE_BATCH_CORRECTION` documented: batch correction cannot separate cohort from disease when every tumor comes from one cohort
  and every normal from another; step 1b is skipped for `DATA_SOURCE = "tcga_gdc"`.
- Step 6 edge extraction is vectorised (same edges, much faster).
- Version bumped to 0.2.0; pytest `testpaths` now `test/`.

### Dependencies / installation
- **xgboost** is now declared in `pyproject.toml`, `environment.yml`, `requirements.txt`, `src/gloom/pipeline/requirements.txt`
  and `conda.recipe/meta.yaml` (it was imported by step 11 but never declared). New optional extra `limma` (rpy2).
- Git LFS documented (`git lfs install && git lfs pull`); `git-lfs` added to `environment.yml`.
- README no longer claims `conda install -c bioconda gloom` works: the Bioconda recipe is not published yet; the
  `pip install -e .` / `mamba env create -f environment.yml` route is documented instead.
- README: removed unsupported superiority/scalability wording; v0.1.x results are flagged as superseded.

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

[Unreleased]: https://github.com/omicscodeathon/gloom/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/omicscodeathon/gloom/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/omicscodeathon/gloom/releases/tag/v0.1.0
