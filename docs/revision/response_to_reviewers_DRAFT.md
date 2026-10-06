# Response to Reviewer 2 — DRAFT

> **Status: DRAFT.** Every quantitative result is a placeholder of the form `[TODO: fill after rerun]`
> and must be filled in from the files produced by the re-run (`docs/revision/RERUN_GUIDE.md`). No number in this
> draft has been computed. Section / page / line references to the manuscript are also placeholders.
> File names refer to GLOOM 0.2.0 (`results/...` = `outputs/results/...` of the pipeline).

Manuscript: *GLOOM: [TODO: title as submitted]* — Reviewer 2 (second review).

---

## General response

We thank the reviewer for a careful and constructive report. The central concerns are correct and, in several
cases, require new analyses rather than new wording: (i) in the original LUAD case study, tumor status was perfectly
confounded with cohort and processing pipeline; (ii) genuine zero expression values had been converted to missing
values; (iii) the construction of the training universe depended on the labels; (iv) the reported top-ranked list included
the positives used for training; (v) the contribution of network features had not been tested by ablation and the network
density difference had not been checked for stability; (vi) the enrichment output had been filtered for cancer relevance before
being interpreted; (vii) "novel" was used for genes that are merely absent from LCGene; and (viii) statements of methodological
advantage over existing tools were not supported by a benchmark.

We therefore revised the software (GLOOM 0.2.0), repeated the LUAD analysis on uniformly processed data, and rewrote the affected
sections. **The original TCGA-GTEx analysis is no longer part of the primary results.** Where the new results differ from
those previously reported, the manuscript now reports the new ones. The principal changes are summarised below.

| # | Reviewer point | What we did | Software / output |
|---|---|---|---|
| 1 | TCGA vs GTEx confounding | Repeated the analysis on TCGA-LUAD primary tumors vs TCGA-LUAD adjacent normals (GDC STAR counts); added sample-level QC; cohort-collinearity check | `fetch_gdc_tcga_luad.py`, `step1` COHORT_DESIGN, `step2b`, `DE_METHOD` |
| 2 | Zeros treated as missing | True zeros kept; NaN/invalid handled separately; everything re-run | `step2`, `gloom_utils.sanitize_expression` |
| 3 | Label-dependent selection | Universe defined before labelling with identical criteria; no \|log2FC\| filter | `LABEL_INDEPENDENT_UNIVERSE`, `step9` |
| 4 | Training positives in the ranking | Repeated stratified cross-fitting; ranking and metrics from out-of-fold scores; bootstrap CIs | `step11c`, `step12b`, `step14` |
| 5 | Network value and stability | Ablation (7 scorers) + paired bootstrap; threshold sweep and resampling stability | `step13b`, `step7c` |
| 6 | Circular enrichment | Complete unfiltered table, explicit background, BH; thematic subset separate | `step19` |
| 7 | "Novel candidates" | Renamed; systematic external evidence classification | `step14` rename, `step20` |
| 8 | Comparison with other tools | Superiority claims removed; architectural differences described | Discussion |
| + | Installation, XGBoost, archive, "scalable", batch correction | See "Additional points" | `README`, manifests, `CITATION.cff`, benchmark |

---

## Point-by-point responses

### Comment 1 — TCGA tumors versus GTEx normals; asymmetry of differential expression

**Reviewer comment:**
> The LUAD tumors are obtained from TCGA/cBioPortal, whereas normal lung samples are from GTEx. These datasets originate from different projects and processing pipelines, yet 9,715 of 10,986 genes are reported as differentially expressed, with 9,700 upregulated and only 15 downregulated. This extreme asymmetry makes it difficult to distinguish biological signal from cohort/platform effects. The LUAD analysis should be repeated using uniformly processed tumor and normal data, such as TCGA tumor versus TCGA adjacent-normal samples or a uniformly processed TCGA-GTEx resource. PCA or comparable sample-level diagnostics should also be provided.

**Response:**
We thank the reviewer for identifying this important limitation. We agree that the original comparison between cBioPortal TCGA-LUAD RSEM tumor data and GTEx TPM normal lung data confounded biological status with data source, sample preparation, quantification pipeline and normalization. Because all tumors originated from TCGA and all controls from GTEx, cohort and phenotype were perfectly collinear, and this could not be resolved by including cohort as a covariate or by batch correction (ComBat / ComBat-seq), which cannot separate a batch effect from a biological effect when the batch variable equals the group variable. We therefore did not attempt to defend the original comparison with a batch correction.

We repeated the LUAD case study using uniformly processed TCGA-LUAD primary tumor and TCGA-LUAD adjacent-normal ("Solid Tissue Normal") samples from the GDC, quantified with the same STAR pipeline ([TODO: n tumor, n normal, n matched patients]). Differential expression was run with [TODO: paired t-test on matched patients / limma-voom on raw counts / Welch]; network construction, feature engineering, model fitting and ranking were then re-run from the harmonized analysis universe. We added sample-level quality control (per-sample median/IQR, sample PCA, global median log2FC, fraction of differentially expressed genes and up:down ratio) and an automatic cohort-collinearity check that warns when tumor and normal samples come from disjoint cohorts. For the original design, the diagnostics give a median log2FC of [TODO] and PC1 separates tumor from normal with AUROC [TODO]; for the new design they give [TODO: median log2FC], [TODO: n up / n down] and [TODO: PC1 explained variance and group separation].

The revised manuscript no longer interprets the original extreme upregulation pattern as a LUAD-specific biological finding. The original TCGA-GTEx analysis is removed from the primary results and discussed only as an example of why matched processing is essential when GLOOM is applied to external cohorts. [TODO (optional): a secondary comparison against a uniformly reprocessed TCGA-GTEx resource, if performed.]

**Changes in the manuscript:**
[TODO: sections/pages/lines] — *Data sources and preprocessing*, *Differential expression*, *Results*, *Discussion* and *Limitations* were rewritten (replacement text: `manuscript_replacement_text.md`, sections 1, 2, 9). New figure [TODO: Figure X] (sample PCA, per-sample median distributions) and Supplementary Figure [TODO: SX]; revised differential-expression results in Supplementary Table [TODO: SX].

**New analysis or supplementary material:**
`results/qc_sample_pca.csv`, `qc_pca_explained_variance.csv`, `qc_sample_summary.csv`, `qc_group_summary.csv`, `qc_global_de_summary.csv`, `qc_cohort_warning.txt` (legacy design), `differential_expression_results.csv` (new design); `data/raw/tcga_gdc/sample_sheet.csv` and the GDC file manifest (provenance); `docs/tcga_paired_design.md`.

---

### Comment 2 — Treatment of zero expression values (Supplementary Method S3)

**Reviewer comment:**
> Supplementary Method S3 states that non-positive values are converted to missing values before log transformation. A true zero expression value is not equivalent to a missing observation, and omitting these values can alter group means, differential-expression statistics, and correlation estimates. Genuine zeros should remain zero, missing/invalid values should be handled separately, and the affected downstream analyses should be rerun.

**Response:**
The reviewer is right. In the original implementation every value ≤ 0 was replaced by a missing value before the log2(x + 1) transformation, so genuine zeros were silently removed from group means, differential-expression statistics and correlations. We corrected the preprocessing: a true zero remains zero (log2(0 + 1) = 0); a missing value remains missing and is counted; negative or infinite values, which are impossible for count/TPM/RSEM data, are flagged as invalid, set to missing and reported with a warning. The low-expression filter now divides by the number of observed values, and the correlation step is robust to missing cells. The change is not presented as a methodological refinement: **all downstream analyses were re-run** (means and medians, log2FC, differential-expression tests, correlations, networks, expression features, models and rankings). In the data used here, [TODO: n true zeros kept, n NaN, n invalid values (from the step-2 log)]. The effect on the results is [TODO: describe, e.g. change in number of retained genes / DE genes / network edges relative to the zero-as-missing version].

**Changes in the manuscript:**
[TODO: sections/pages/lines] — Supplementary Method S3 rewritten (replacement text section 2 of `manuscript_replacement_text.md`); all numbers in the Results updated; the Methods now state explicitly how zeros, missing and invalid values are handled.

**New analysis or supplementary material:**
Code: `step2_preprocessing.py` / `gloom_utils.sanitize_expression`; unit tests `test/test_expression_zero_handling.py`; QC counts in `pipeline.log` and `qc_summary.csv`. [TODO (optional): Supplementary Table SX comparing key counts before/after the correction.]

---

### Comment 3 — Label-dependent selection of the training universe

**Reviewer comment:**
> Unlabeled genes with |log2FC| < 2 are removed, whereas LCGene-positive genes are retained regardless of fold change. Fold-change-derived variables are then used again as model predictors. This creates label-dependent selection and may inflate class separation. The analysis universe should be defined using identical criteria for all genes before assigning LCGene labels, followed by repeated model training and evaluation.

**Response:**
We agree that this procedure creates a label-dependent selection bias: because positives were kept regardless of fold change while unlabeled genes with a small fold change were removed, and fold-change-derived variables were then model features, part of the class separation was created by the selection itself. In GLOOM 0.2.0 the analysis universe is defined **before** labelling, with criteria identical for every gene (detectable expression and variance filters applied uniformly, and the intersection of genes available in both groups); **no |log2FC| filter is applied to unlabeled genes**, and labels are assigned afterwards. The previous behaviour remains available only as a deprecated option (`LABEL_INDEPENDENT_UNIVERSE = False`) and is not used for any reported result. Model training and evaluation are repeated (see Comment 4). The final universe contains [TODO: n genes], of which [TODO: n] are LCGene positives ([TODO: %]); the old trimming would have removed [TODO: n] unlabeled genes.

**Changes in the manuscript:**
[TODO: sections/pages/lines] — *Label construction* and *Data sources and preprocessing* rewritten (replacement text section 3); the sentence describing removal of unlabeled genes with |log2FC| < 2 was deleted; Results updated.

**New analysis or supplementary material:**
`step9_label_construction.py`, `config.LABEL_INDEPENDENT_UNIVERSE`, `gene_labels.csv`, `gene_annotation_table.csv`. [TODO (optional): sensitivity table with the deprecated label-dependent universe to quantify how much the bias inflated the previous performance.]

---

### Comment 4 — Evaluation includes training positives in the ranked list

**Reviewer comment:**
> The manuscript reports that 98 of the top 100 ranked genes are known LCGene positives. However, LCGene-positive genes used to train the PU models are also included in the final ranked list. This demonstrates recovery of training positives rather than fully independent ranking performance. Ranking should be evaluated using out-of-sample predictions, for example through cross-fitting or evaluation restricted to held-out genes.

**Response:**
We agree. The "98 of the top 100" statement measured how well the model re-scores genes it was trained on, and we have removed it. The primary ranking is now obtained by **repeated stratified cross-fitting**: all genes are split into five stratified folds; the same PU-bagging model is trained on the positives and unlabeled genes of four folds and predicts only the held-out fold; the five sets of held-out predictions are concatenated, and the procedure is repeated [TODO: 5] times with different seeds and the out-of-fold scores are averaged. Every gene is therefore scored by a model that never saw it. All reported metrics — AUROC, AUPRC, average precision, Precision@10/50/100, Recall@10/50/100 and the enrichment factor — are computed from these out-of-fold scores, with 95% confidence intervals from a stratified bootstrap over genes (1,000 resamples, fixed seed). The results are: AUROC [TODO] (95% CI [TODO]), AUPRC [TODO] ([TODO]), Precision@100 [TODO] ([TODO]), Recall@100 [TODO] ([TODO]), enrichment factor@100 [TODO] ([TODO]). The number of LCGene positives among the top 10/50/100 out-of-fold-ranked genes is [TODO]/[TODO]/[TODO]. The score of the model fitted on all data is retained in a separate column (`full_model_score`) and the training positives are marked (`is_training_positive`); the top-K statistics of the single 80/20 hold-out are reported separately for held-out and training positives, and we note that a single split is not sufficient on its own. Because "unlabeled" genes contain undiscovered positives, precision and AUROC computed against LCGene are conservative.

**Changes in the manuscript:**
[TODO: sections/pages/lines] — *PU bagging*, *Evaluation*, *Results* and *Abstract* rewritten (replacement text section 4); the "98 of 100" statement removed from the Abstract, Results and Discussion; Table [TODO: X] now reports out-of-fold metrics with confidence intervals.

**New analysis or supplementary material:**
`results/oof_scores.csv`, `oof_metrics.csv`, `oof_fold_assignments.csv`, `gene_rankings.csv` (columns `predicted_prob`, `oof_score`, `full_model_score`, `is_training_positive`), `pu_bagging_metrics.csv` (held-out vs training top-K); unit tests `test/test_crossfit.py`, `test/test_metrics.py`.

---

### Comment 5 — Contribution of network features and network density difference

**Reviewer comment:**
> The feature-importance results show that most of the predictive signal comes from tumor-expression, normal-expression, and differential-expression features, while network variables contribute relatively little. An ablation analysis comparing expression-only, network-only, combined expression+network, and a simple differential-expression baseline would establish whether network integration improves performance. The large difference in network density also requires additional examination. The tumor network contains 110,508 edges compared with 1,276,368 edges in the normal network using the same correlation threshold. Network stability across reasonable correlation thresholds or resampling should be shown before this difference is interpreted as disease-specific rewiring.

**Response:**
We thank the reviewer; the Discussion previously stated that network integration added value without a direct comparison. We now provide two analyses.

*Ablation.* Under the same cross-fitting protocol (same folds, same PU-bagging model, same bootstrap) we compared: (a) ranking by |log2FC| alone, (b) ranking by adjusted p-value alone, (c) expression-derived features only, (d) network-derived features only, (e) expression + network features, (f) the combined set without explicit log2FC-derived features, and (g) the full model. Out-of-fold AUROC/AUPRC (95% CI) were: [TODO: table]. The paired-bootstrap difference between the full model and the expression-only model was ΔAUPRC = [TODO] (95% CI [TODO]; p = [TODO]) and ΔAUROC = [TODO] ([TODO]). [TODO: choose one] — *If the CI excludes zero:* network features add a measurable gain of [TODO]. — *Otherwise:* we conclude that network-derived features primarily support contextual interpretation of prioritized genes, whereas predictive performance in the LUAD case study was driven predominantly by expression-derived variables, and we revised the Abstract and Discussion accordingly. (We also note that tumor and normal mean expression still encode fold change implicitly in variant (f).)

*Network stability.* We rebuilt the networks at |r| ≥ 0.50, 0.55, 0.60, 0.65 and 0.70 and report edges, density, mean degree, size of the largest component, clustering coefficient, isolated nodes and the Spearman correlation of degree centrality between consecutive thresholds. Because the number of samples alone changes the amount of spurious correlation, the sweep was repeated after sub-sampling the larger group to the size of the smaller group. The normal-to-tumor edge ratio at the working threshold was [TODO] (equal-N: [TODO]) and [TODO: range] across thresholds. In addition, a sample bootstrap ([TODO: 20] resamples) gave a mean edge retention frequency of [TODO] ([TODO]% of edges retained in ≥ 80% of resamples) and a hub Jaccard stability of [TODO] ± [TODO] (hubs = top 5% by degree). [TODO: state whether the density difference is robust; if it is not, remove the "disease-specific rewiring" interpretation and describe the differential-network features as exploratory.]

**Changes in the manuscript:**
[TODO: sections/pages/lines] — new *Ablation* and *Network stability* subsections (replacement text sections 5–6); Feature-contribution section, Discussion and Abstract revised; Figure [TODO: X] (ablation) and Supplementary Figures/Tables [TODO: SX] (stability).

**New analysis or supplementary material:**
`results/ablation_metrics.csv`, `ablation_paired_comparison.csv`, `ablation_verdict.txt`, `ablation_scores.csv`; `network_stability_thresholds.csv`, `network_stability_density_ratio.csv`, `network_stability_bootstrap_summary.csv`, `network_stability_bootstrap_replicates.csv`, `network_stability_hubs.csv`.

---

### Comment 6 — Circular interpretation of enrichment (Supplementary Method S11)

**Reviewer comment:**
> Supplementary Method S11 states that enrichment results are filtered to retain lung- and cancer-relevant pathways. Reporting cancer-related pathways after specifically filtering for cancer relevance introduces circular interpretation. The complete enrichment results should be reported first, with appropriate multiple-testing correction and gene background. A cancer-focused subset can then be discussed separately.

**Response:**
We agree, and we also found that the file named as the "all candidates" result was itself already filtered. Enrichment is now reported in two separate outputs. The **primary output** is the complete, unfiltered table for all KEGG terms with at least one gene in the background, using an **explicit background consisting of the [TODO: n] genes actually analysed by GLOOM** (not all KEGG genes), a one-sided hypergeometric test and Benjamini–Hochberg correction across all tested terms; the columns are term, number of query genes, number of term genes in the background, number of overlapping genes, raw p-value and adjusted p-value. Of [TODO: n] terms tested, [TODO: n] have p_adj < 0.05; the top terms are [TODO]. The **secondary output** is a cancer/lung keyword subset, clearly labelled as an aid to reading: *"For interpretive purposes, we additionally extracted pathways containing predefined lung- or cancer-related terms. This secondary display was not used to determine statistical significance or establish biological validity."* The adjusted p-values shown in that subset are those of the complete correction, not re-computed on the subset.

**Changes in the manuscript:**
[TODO: sections/pages/lines] — Supplementary Method S11 rewritten (replacement text section 7); Results and Figure [TODO: X] now show the complete enrichment results; the discussion of cancer-related pathways is based on the secondary subset and is worded accordingly.

**New analysis or supplementary material:**
`results/enrichment/kegg_all_candidates.csv`, `kegg_upregulated.csv`, `kegg_downregulated.csv` (complete tables), `kegg_lung_cancer_subset.csv` (interpretive subset), `kegg_summary.csv`; unit tests `test/test_enrichment.py`.

---

### Comment 7 — "Novel candidates" and absence from LCGene

**Reviewer comment:**
> Absence from LCGene does not establish biological novelty. Several highlighted genes, including CDC20, MYBL2, and MARCO, have already been associated with lung cancer or LUAD in published studies. The term "novel candidates" should therefore be reconsidered unless novelty is established through a systematic literature assessment.

**Response:**
We agree: absence from a curated reference set does not establish novelty. We replaced "novel candidates" throughout the manuscript, figures, tables and software with **"non-LCGene candidate genes"** (genes absent from the LCGene reference set that are prioritized by GLOOM). In addition, every non-LCGene candidate in the top [TODO: 200] was checked against Open Targets (association scores for lung adenocarcinoma EFO_0000571, lung carcinoma EFO_0001071 and cancer) and Europe PMC (co-mention counts with "lung adenocarcinoma", "lung cancer" and "cancer") and classified as already associated with LUAD ([TODO: n]), associated with lung cancer without LUAD specificity ([TODO: n]), associated with other cancers ([TODO: n]), indirect/weak evidence only ([TODO: n]) or no association found ([TODO: n]). Only genes in the last category are described, with caution, as *potentially novel*, and only with respect to these two sources. CDC20, MYBL2 and MARCO are classified as [TODO] and are no longer presented as novel. [TODO (optional): DisGeNET, IntOGen and CIViC were also searched / are listed as future work; a systematic PubMed search was / was not performed.]

**Changes in the manuscript:**
[TODO: sections/pages/lines] — terminology changed in the Abstract, Results, Discussion, tables and figure legends; new subsection on external evidence (replacement text section 8); Table [TODO: X] lists the top candidates with their evidence class.

**New analysis or supplementary material:**
`results/candidate_evidence.csv` (with cached raw query results), `non_lcgene_candidates.csv`; code `step20_candidate_evidence.py` (requires internet; never interrupts the main pipeline).

---

### Comment 8 — Comparison with Endeavour, ToppGene, GeneMANIA and BioRank

**Reviewer comment:**
> The Discussion states that GLOOM offers a methodological advantage over tools such as Endeavour, ToppGene, GeneMANIA, and BioRank, but no head-to-head comparison is presented. These statements should either be supported by quantitative benchmarking or revised to describe differences in workflow architecture rather than superior performance.

**Response:**
We agree that no quantitative head-to-head comparison was performed and that the original statements of methodological or practical advantage were not supported. We removed them and now describe the differences in workflow architecture: GLOOM integrates user-supplied expression matrices, differential-expression analysis, disease-conditioned co-expression networks, feature engineering, positive-unlabeled learning and reporting in one reproducible workflow, whereas the cited tools [TODO: verify against each tool's current documentation] operate on pre-computed knowledge bases or networks. We state explicitly that quantitative benchmarking (same gene universe, same training set, same held-out positives, Precision@K, Recall@K, AUPRC) is an important direction for future work.

**Changes in the manuscript:**
[TODO: sections/pages/lines] — Discussion rewritten (replacement text section 10); "methodological advantage" and "practical advantage" removed; Limitations and Future Work updated.

**New analysis or supplementary material:**
None (no benchmark was performed).

---

## Additional points raised in the accompanying review report

### A. Installation: Bioconda command and Git LFS

**Point:** the report states that `conda install -c bioconda gloom` does not work, that no functional Bioconda package or recipe is available, and that a fresh clone can fail without Git LFS, which is not documented.

**Response:** The report is correct: the Bioconda recipe has not been published. We removed the Bioconda installation command from the README and the manuscript and now give the routes that work (`mamba env create -f environment.yml` followed by `pip install -e .`, or `pip install -e .` in a Python 3.12 environment). Git LFS is documented (`git lfs install && git lfs pull`, with an explanation of what happens without it) and `git-lfs` is included in `environment.yml`.
**Changes:** README (Installation), manuscript *Availability* section [TODO: lines]; `CHANGELOG.md`.

### B. XGBoost dependency

**Point:** XGBoost appears in the benchmark but is not declared in any dependency file.
**Response:** Correct. `xgboost>=1.7.0` is now declared in `pyproject.toml`, `environment.yml`, `requirements.txt`, `src/gloom/pipeline/requirements.txt` and the conda recipe, and the README states it. (The model step skipped XGBoost with only a log warning when it was not installed; the benchmark reported in the manuscript was produced with xgboost [TODO: version from `pip_freeze.txt`].)
**Changes:** manifests, README, Methods [TODO: lines].

### C. Archived, citable version

**Point:** a stable release with a tag, an archive and a DOI cited in the manuscript.
**Response:** We prepared version 0.2.0, a `CITATION.cff` and a release checklist (`docs/RELEASE_CHECKLIST.md`). The release [TODO: v0.2.0 tag] was archived on Zenodo with DOI [TODO: 10.5281/zenodo.XXXXXXX], which is cited in the Data and Code Availability section.
**Changes:** manuscript *Code availability* [TODO: lines]; reference list.

### D. "Scalable"

**Point:** the word "scalable" is not defensible without performance benchmarks.
**Response:** We replaced "scalable" by "modular" and "designed for reproducible execution" throughout and removed any implied claim of scalability. [TODO — choose one: *We additionally report measured run time and peak memory for [TODO: n] gene × sample sizes on [TODO: machine] (Supplementary Table [TODO: SX], `runtime_benchmark.csv`).* / *No scalability claim is made.*]
**Changes:** Abstract, Introduction, Discussion [TODO: lines]; `scripts/benchmark_runtime.py`.

### E. Batch correction (ComBat / ComBat-seq) for the TCGA–GTEx design

**Point (implicit in comment 1):** can batch correction solve the cohort effect?
**Response:** No. In the original design every tumor came from TCGA and every normal from GTEx, so batch and biological group are perfectly collinear: a ComBat-type model cannot estimate both, and correcting for batch without the group covariate removes the disease signal together with the batch signal. We therefore did not use batch correction to rescue the original comparison. GLOOM now documents this limitation in the optional batch-correction step and the README, detects perfect collinearity between cohort and group at data loading (loud warning, `qc_cohort_warning.txt`), and the primary analysis uses a design in which tumor and normal samples come from the same project and pipeline.
**Changes:** Methods and Limitations [TODO: lines]; README; `step1b_batch_correction.py` documentation.

### F. Reproducibility (seeds)

**Response:** A single `SEED` (default 42) in `config.py` now governs every random component (splits, folds, bagging, bootstrap, network layout), and Python's and NumPy's global generators are also seeded at pipeline start; the seed is stated in the Methods. Unit tests with synthetic data cover the new functions.
**Changes:** Methods [TODO: lines]; `config.py`; `test/`.

---

## Items to complete before submission

- [ ] Fill every `[TODO: ...]` from the re-run output; delete the alternatives that do not apply.
- [ ] Insert the exact section / page / line references of the revised manuscript.
- [ ] Verify each statement about Endeavour, ToppGene, GeneMANIA and BioRank against the tools' current documentation.
- [ ] Confirm the Zenodo DOI and the tag in `CITATION.cff`, the README and the manuscript.
- [ ] Confirm that no sentence of the manuscript still uses "novel candidates", "scalable", "methodological advantage", the "98 of the top 100" claim, or the Bioconda installation line.
