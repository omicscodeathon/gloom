# Re-run guide for the major revision (GLOOM 0.2.0)

Goal: regenerate every number, table and figure of the LUAD case study with the revised pipeline
on a machine that has **Python 3.12** (or in **Google Colab**), then send the result files listed in
section 7 back so that the manuscript and the response letter can be completed.

All commands are for a POSIX shell (Linux, macOS, WSL, Colab with a leading `!`).
Run everything from the repository root. **Run the pipeline from `src/gloom/pipeline/`** (the
packaged copy): its `config.py` resolves `data/` and `outputs/` relative to the repository root.
`scripts/` holds the identical step files plus the helper tools (`fetch_gdc_tcga_luad.py`, `benchmark_runtime.py`).

Estimated cost (indicative; depends on the machine): the cross-fitting (step 11c) and the ablation (step 13b)
dominate the run time (tens of minutes to a few hours with the default settings, see section 5 for a quick smoke test).

---

## 0. Which analyses are needed

| # | Analysis | Manuscript / reviewer point |
|---|---|---|
| A | Legacy TCGA(cBioPortal) vs GTEx design: **QC and diagnostics only** (steps 1-4) | Reviewer comment 1 (PCA, DE asymmetry, why the old comparison is confounded) |
| B | **Primary analysis**: TCGA-LUAD primary tumor vs TCGA-LUAD adjacent normal (GDC STAR counts), full pipeline | Comments 1-7 |
| C | Optional external evidence for candidates (step 20, internet) | Comment 7 |
| D | Runtime benchmark | "scalable" claim |

---

## 1. Setup

```bash
git lfs install
git clone https://github.com/omicscodeathon/gloom.git
cd gloom
git lfs pull                                   # fetches the bundled data (otherwise only LFS pointer files!)
git checkout <revision-branch-or-commit>       # the commit that contains version 0.2.0

# environment (Python 3.12)
mamba env create -f environment.yml            # or: conda env create -f environment.yml
mamba activate gloom
pip install -e ".[full]"                       # xgboost is a core dependency; [full] adds plotly, gseapy, openpyxl, imbalanced-learn
python --version && git rev-parse HEAD         # keep these two lines for the record
```

Colab:

```python
!apt-get -qq install git-lfs && git lfs install
!git clone https://github.com/omicscodeathon/gloom.git && cd gloom && git lfs pull
%cd gloom
!pip install -q -e ".[full]"
```

Check the installation and the unit tests:

```bash
pytest -q                                      # synthetic-data tests in test/ (should all pass)
python -c "import gloom, xgboost; print(gloom.__version__)"    # expect 0.2.0
```

Optional, only for `DE_METHOD = "limma_voom"`: R with the Bioconductor packages `limma` and `edgeR`, plus `pip install rpy2`.

---

## 2. Analysis A — diagnostics of the legacy TCGA(cBioPortal)-vs-GTEx design (steps 1-4)

This documents *why* the original contrast is confounded (loud cohort warning, sample PCA, global DE statistics).
It does not need the whole pipeline.

```bash
# config.py must contain DATA_SOURCE = "cbioportal_gtex" (the default) and DE_METHOD = "welch"
python src/gloom/pipeline/run_pipeline.py --from 1 --to 4
mv outputs outputs_legacy_qc                   # keep these results; the next run starts from a clean outputs/
```

Expected: `outputs_legacy_qc/results/qc_cohort_warning.txt` (cohort perfectly collinear with group, and the
global-artefact warning), `qc_sample_pca.csv`, `qc_group_summary.csv`, `qc_global_de_summary.csv`,
`figures/qc_sample_pca.png`, `figures/de_volcano_plot.png`. Because step 2 now keeps true zeros, the numbers differ from the
original manuscript (9,715 / 10,986 DE genes, 9,700 up, 15 down); report the new values.

---

## 3. Analysis B — primary analysis: TCGA-LUAD tumor vs adjacent normal

### 3.1 Download the uniformly processed data from the GDC

```bash
python scripts/fetch_gdc_tcga_luad.py --paired-only --workers 6
#   about 59 adjacent normals + their matched tumors (about 0.5 GB). Without --paired-only: all primary tumors
#   and all normals (about 2.5 GB). Add --query-only first to list the files without downloading.
ls data/raw/tcga_gdc/     # counts_matrix.csv  tpm_matrix.csv  sample_sheet.csv  gdc_files_manifest.csv  files/
```

(The GDC API is open access; no token is needed. If the download is interrupted, re-run the same command: files already downloaded are cached.)

### 3.2 Configure the design

Edit **`src/gloom/pipeline/config.py`** (do not commit these edits; keep a copy with `git diff > config_used.diff`):

```bash
CFG=src/gloom/pipeline/config.py
sed -i 's/^DATA_SOURCE = .*/DATA_SOURCE = "tcga_gdc"/' $CFG
sed -i 's/^DE_METHOD = .*/DE_METHOD = "paired"/' $CFG        # paired t-test on matched patients ("welch" if you used all tumors unpaired)
grep -n '^DATA_SOURCE\|^DE_METHOD\|^SEED\|^USE_CROSSFIT\|^LABEL_INDEPENDENT_UNIVERSE\|^USE_BATCH_CORRECTION' $CFG
```

Leave `LABEL_INDEPENDENT_UNIVERSE = True`, `USE_CROSSFIT = True`, `USE_BATCH_CORRECTION = False`, `SEED = 42`.

Optional second DE analysis with limma-voom on raw counts (needs rpy2 + R): run the whole pipeline once with
`DE_METHOD = "paired"` and, to report the sensitivity of the DE results, re-run only step 4 with `DE_METHOD = "limma_voom"` into a copy of the
results folder (`cp -r outputs outputs_paired && sed -i 's/^DE_METHOD = .*/DE_METHOD = "limma_voom"/' $CFG && python src/gloom/pipeline/run_pipeline.py --only 4 && cp outputs/results/differential_expression_results.csv outputs_paired/results/differential_expression_results_limma_voom.csv`, then restore `DE_METHOD = "paired"`).

### 3.3 Run the full pipeline

```bash
python src/gloom/pipeline/run_pipeline.py 2>&1 | tee run_full.log
```

Step order (new steps in bold): 1 data loading + **COHORT_DESIGN check**, 2 preprocessing (true zeros kept),
**2b sample QC**, 3 harmonization, 4 differential expression (`DE_METHOD`), 5 expression features, 6 tumor network,
6b normal network, 7 / 7b network features, **7c network stability**, 8 integration, 9 labels (label-independent universe),
10 split, 11 / 11b models, **11c cross-fitting**, 12 evaluation, **12b out-of-fold metrics**, 13 importance,
**13b ablation**, 14 ranking (out-of-fold primary), 15-18 network/visualization/report, 19 KEGG (complete table), **20 candidate evidence**.

Optional steps are marked `[optional]` in the log: a failure there prints `SKIPPED` and does not stop the run.
Mandatory steps stop the run on failure: send the last 80 lines of `outputs/logs/pipeline.log`.

To run pieces separately:

```bash
python src/gloom/pipeline/run_pipeline.py --from 11 --to 14     # models + cross-fitting + ranking
python src/gloom/pipeline/run_pipeline.py --only 13b            # ablation only
python src/gloom/pipeline/run_pipeline.py --only 7c             # network stability only
python src/gloom/pipeline/run_pipeline.py --only 20             # external evidence (needs internet)
```

### 3.4 Sanity checks before trusting the results

```bash
test -f outputs/results/qc_cohort_warning.txt && cat outputs/results/qc_cohort_warning.txt || echo "no cohort warning (expected for the GDC design)"
cat outputs/results/qc_global_de_summary.csv        # median log2FC near 0? plausible fraction DE? up:down not extreme?
cat outputs/results/oof_metrics.csv | head -20
cat outputs/results/ablation_verdict.txt
```

If a cohort warning is raised for the GDC design, do not continue: send the file back first.

---

## 4. Analysis C — external evidence for the non-LCGene candidates (internet)

Step 20 runs at the end of the full pipeline when internet is available. To (re-)run it alone:

```bash
python src/gloom/pipeline/run_pipeline.py --only 20
# or: python src/gloom/pipeline/step20_candidate_evidence.py
```

It queries Open Targets (GraphQL) and Europe PMC for the top 200 non-LCGene genes, with retries and a local cache
(`outputs/results/candidate_evidence_cache.json`; delete it to force fresh queries).
If the Open Targets schema changed and every GraphQL variant is rejected, genes get `ot_status = error: ...` and class `not_assessed`
(never `no_association_found`): send `candidate_evidence.csv` back so the query can be adapted.

---

## 5. Quick smoke test (about 10-15 minutes)

To verify the installation end to end before the long run, temporarily lower the cost in `config.py`:

```bash
sed -i 's/^CROSSFIT_REPEATS *=.*/CROSSFIT_REPEATS      = 1/' $CFG
sed -i 's/^CROSSFIT_PU_N_ESTIMATORS *=.*/CROSSFIT_PU_N_ESTIMATORS = 30/' $CFG
sed -i 's/^BOOTSTRAP_N *=.*/BOOTSTRAP_N           = 100/' $CFG
sed -i 's/^ABLATION_REPEATS *=.*/ABLATION_REPEATS         = 1/' $CFG
sed -i 's/^ABLATION_PU_N_ESTIMATORS *=.*/ABLATION_PU_N_ESTIMATORS = 30/' $CFG
sed -i 's/^NETWORK_BOOTSTRAP_N *=.*/NETWORK_BOOTSTRAP_N     = 3/' $CFG
python src/gloom/pipeline/run_pipeline.py 2>&1 | tail -40
```

**Restore the defaults (5 repeats, 300 estimators, 1000 bootstrap resamples,
3 ablation repeats, 20 network resamples) for the final run** — the smoke-test numbers must not be used in the manuscript.
Remove the smoke-test outputs with `rm -rf outputs` before the final run.

---

## 6. Analysis D — runtime benchmark

```bash
python scripts/benchmark_runtime.py --genes 2000 5000 10000 --samples 100 300 500 --repeats 3
# writes outputs/results/runtime_benchmark.csv (time, peak memory, machine information)
```

Run it on the machine that will be named in the manuscript; it takes a while for the largest sizes
(use `--genes 2000 5000 --samples 100 300` for a short version). Install `psutil` (`pip install psutil`) to record RAM and, on Windows, peak RSS.

---

## 7. What to send back

Create one archive and send it together with the environment record:

```bash
python -m pip freeze > pip_freeze.txt
git rev-parse HEAD > git_commit.txt
git diff -- src/gloom/pipeline/config.py > config_used.diff
mkdir -p return_package/results return_package/figures return_package/data
cp outputs/logs/pipeline.log return_package/
cp pip_freeze.txt git_commit.txt config_used.diff return_package/
# core results
cp outputs/results/{differential_expression_results.csv,gene_rankings.csv,non_lcgene_candidates.csv,non_lcgene_candidates_sensitivity.csv} return_package/results/
cp outputs/results/{ranking_metrics.csv,pu_bagging_metrics.csv,pu_bagging_scores.csv,model_metrics.csv,feature_importance.csv} return_package/results/
cp outputs/results/{oof_scores.csv,oof_metrics.csv,oof_fold_assignments.csv} return_package/results/
cp outputs/results/{ablation_scores.csv,ablation_metrics.csv,ablation_paired_comparison.csv,ablation_verdict.txt} return_package/results/
cp outputs/results/network_stability_*.csv return_package/results/
cp outputs/results/qc_*.csv return_package/results/ ; cp outputs/results/qc_cohort_warning.txt return_package/results/ 2>/dev/null || true
cp outputs/results/candidate_evidence.csv outputs/results/runtime_benchmark.csv return_package/results/ 2>/dev/null || true
cp -r outputs/results/enrichment return_package/results/
cp outputs/results/reports/pipeline_report.txt return_package/results/ 2>/dev/null || true
# figures
cp outputs/figures/{qc_sample_pca.png,qc_sample_medians.png,de_volcano_plot.png,de_log2fc_distribution.png,ranking_score_distribution.png,pu_bagging_score_distribution.png,kegg_overview.png} return_package/figures/ 2>/dev/null || true
# data provenance (NOT the big matrices)
cp data/raw/tcga_gdc/{sample_sheet.csv,gdc_files_manifest.csv} return_package/data/
# the legacy diagnostics (analysis A)
cp -r outputs_legacy_qc return_package/legacy_qc 2>/dev/null || true
zip -r gloom_rerun_package.zip return_package
```

Checklist of what must be in the package:

- [ ] `pipeline.log` (complete), `git_commit.txt`, `pip_freeze.txt`, `config_used.diff`
- [ ] `qc_cohort_warning.txt` if it exists (for the GDC design it should NOT exist) and `qc_*.csv` for **both** designs (legacy in `legacy_qc/`)
- [ ] `differential_expression_results.csv`, `qc_global_de_summary.csv`, `qc_sample_pca.csv`, `qc_pca_explained_variance.csv`, `qc_sample_pca.png`
- [ ] `oof_metrics.csv`, `oof_scores.csv`, `ranking_metrics.csv`, `pu_bagging_metrics.csv`
- [ ] `ablation_metrics.csv`, `ablation_paired_comparison.csv`, `ablation_verdict.txt`
- [ ] `network_stability_thresholds.csv`, `network_stability_density_ratio.csv`, `network_stability_bootstrap_summary.csv`, `network_stability_hubs.csv`
- [ ] `enrichment/kegg_all_candidates.csv` (complete table), `kegg_lung_cancer_subset.csv`, `kegg_summary.csv`
- [ ] `non_lcgene_candidates.csv`, `candidate_evidence.csv`
- [ ] `runtime_benchmark.csv`
- [ ] Number of samples per group, number of matched pairs (printed by step 1 and in `sample_sheet.csv`)

Do **not** send the GDC matrices or the `files/` folder (large, re-downloadable); the manifest identifies them.

---

## 8. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| Step 1 fails with "file not found" or the CSV looks like `version https://git-lfs...` | Git LFS files were not pulled: `git lfs install && git lfs pull` |
| `ModuleNotFoundError: xgboost` in step 11 | Install it: `pip install xgboost` (it is now a declared dependency; re-run `pip install -e .`) |
| Step 4 prints "limma_voom SKIPPED" | rpy2 / R packages not available; the Welch test was used instead (check the `de_method` column of the DE table) |
| Step 4 prints "only 0 tumor/normal pair(s)" with `DE_METHOD = "paired"` | `patient_id` missing: check `sample_sheet.csv`, or the TCGA barcodes of the sample ids |
| Step 11c / 13b are too slow | lower `CROSSFIT_REPEATS`, `CROSSFIT_PU_N_ESTIMATORS`, `ABLATION_*` for a pilot; use the defaults for the final run |
| Step 19 fails offline | gseapy needs internet to download the KEGG library; set `ENRICHMENT_GENESETS_FILE` to a local GMT file |
| Step 20 prints "No internet access" | expected offline; re-run `--only 20` on a connected machine |
| Memory error in step 6 / 7c | lower the number of genes or raise `COEXPR_CORRELATION_CUTOFF`; step 7c at |r| >= 0.50 needs the most memory |
