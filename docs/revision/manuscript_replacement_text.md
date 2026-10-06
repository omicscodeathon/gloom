# Manuscript replacement text (ready to paste) — GLOOM 0.2.0 revision

> Paragraphs are written in the manuscript's English register. Every number is a placeholder
> `[TODO: ...]` to be filled from the re-run output (`RERUN_GUIDE.md`); **no value here has been computed**.
> Where two variants are given, keep the one that matches the results and delete the other.
> Section numbering below matches the references in `response_to_reviewers_DRAFT.md`.

---

## 1. Data sources and preprocessing — data

**Replace the description of the TCGA/cBioPortal and GTEx data by:**

> **Data sources.** Gene expression of lung adenocarcinoma (LUAD) and adjacent non-malignant lung tissue was obtained from
> the Genomic Data Commons (GDC) for the TCGA-LUAD project, using the "Gene Expression Quantification" files produced by the
> GDC "STAR - Counts" workflow. We retained primary tumor samples (sample type "Primary Tumor") and adjacent normal samples
> ("Solid Tissue Normal"), which were quantified with the same pipeline and gene model. The analysis included
> [TODO: n] primary tumors and [TODO: n] adjacent normal samples, of which [TODO: n] patients contributed both a tumor and
> an adjacent-normal sample. Expression was analysed as [TODO: TPM (tpm_unstranded) / raw unstranded counts for limma-voom];
> gene symbols were taken from the GDC gene annotation, restricted to [TODO: protein-coding] genes, and duplicated symbols were
> collapsed by retaining the entry with the highest mean expression. Known LUAD-associated genes (positive labels) were taken from
> the LCGene database (LUAD-filtered; [TODO: 517] gene symbols, of which [TODO: n] were present in the analysis universe).
> The data files and the sample sheet used are listed in Supplementary Table [TODO: SX] and can be re-created with
> `scripts/fetch_gdc_tcga_luad.py`.

**Add (why the original design was dropped):**

> An earlier version of this work compared TCGA-LUAD tumors obtained from cBioPortal (RSEM) with normal lung samples from GTEx
> (TPM). Because all tumors originated from one resource and all controls from another, tumor status was perfectly collinear with
> cohort, sample processing and quantification pipeline, and batch-correction methods cannot separate such effects (see Limitations).
> That comparison was therefore removed from the primary analysis; it is retained only as a documented example of a confounded
> design, for which GLOOM now raises an explicit warning.

---

## 2. Preprocessing (Supplementary Method S3) and differential expression

**Replace Supplementary Method S3 by:**

> **Value handling and normalization.** Each expression value was classified before transformation. True zeros were kept as zeros,
> so that log2(0 + 1) = 0, because a zero is an observed measurement and not a missing value. Missing values (NaN) were kept as
> missing and counted. Negative or infinite values, which are not possible in count, TPM or RSEM data, were flagged as invalid,
> set to missing and reported (in this dataset: [TODO: n] zeros, [TODO: n] missing and [TODO: n] invalid values out of [TODO: n]
> values). Expression was transformed as log2(x + 1). Genes were retained if at least [TODO: 10%] of the observed samples had
> expression above [TODO: log2(1 + 1)] (zeros counting as not expressed) and if their interquartile range exceeded the
> [TODO: 10th] percentile of the cohort; both filters were applied uniformly to all genes, independently of their labels.
> Only genes present in both groups were retained ([TODO: n] genes).

**Replace the Differential expression paragraph by (choose the variant matching `DE_METHOD`):**

> **Differential expression.** *[Variant paired]* Differential expression between primary tumors and adjacent normal tissue was
> assessed with a paired t-test on the [TODO: n] patients with both samples (log2(TPM + 1) values; log2 fold change = mean paired
> difference), with Benjamini–Hochberg correction. *[Variant limma-voom]* Differential expression was assessed on raw counts with
> limma-voom (edgeR normalization factors; patient as a blocking factor for the [TODO: n] paired patients), with
> Benjamini–Hochberg correction. *[Variant Welch]* Differential expression was assessed with Welch's t-test on log2(x + 1) values
> with Benjamini–Hochberg correction. Genes were called differentially expressed at FDR ≤ [TODO: 0.001] and |log2 fold change|
> ≥ [TODO: 2]; [TODO: n] genes met this criterion ([TODO: n] up, [TODO: n] down).

**Add (sample-level diagnostics):**

> **Sample-level quality control.** Before differential expression we computed per-sample median and interquartile range of
> log2 expression, a principal-component analysis of the samples (most variable [TODO: 5,000] genes), the global median log2
> fold change, the fraction of differentially expressed genes and the up:down ratio. A warning is raised automatically if
> |median log2 fold change| > 1, if more than 80% of genes are differentially expressed, or if the up:down ratio exceeds 10 or
> falls below 0.1. In the present design the median log2 fold change was [TODO], [TODO: %] of genes were differentially
> expressed (up:down = [TODO]), and PC1 ([TODO: %] of the variance) [TODO: did / did not] separate tumor from normal samples
> (Figure [TODO: X]).

---

## 3. Label construction and analysis universe

**Replace the Label construction paragraph by:**

> **Label construction.** The analysis universe was defined before labels were assigned, using criteria that are identical for
> every gene: detectable expression, a variance filter (both described above) and presence in both groups. Genes were not
> selected or removed according to their fold change, because fold-change-derived variables are used as model features and a
> label-dependent filter would inflate class separation. LCGene-listed genes in the universe were labelled positive
> (n = [TODO]; [TODO: %] of [TODO: n] genes) and all other genes were labelled unlabeled. Unlabeled genes are not confirmed
> negatives: they include true negatives and undiscovered positives, which is why a positive-unlabeled (PU) learning framework was used.

---

## 4. PU bagging, cross-fitting and evaluation

**Replace the PU bagging / evaluation / "98 of the top 100" text by:**

> **PU bagging.** We used Mordelet–Vert bagging: in each of [TODO: 300] iterations a random subset of unlabeled genes of the same size as
> the positive set was combined with all positives, a random forest ([TODO: 100] trees, balanced class weights) was trained, and the final
> score of a gene is the mean predicted probability over iterations.
>
> **Out-of-fold evaluation (cross-fitting).** To obtain scores that are independent of the training labels, all genes were split into five
> stratified folds. For each fold, the PU-bagging model was trained on the positives and unlabeled genes of the other four folds and used to
> predict only the held-out fold. The held-out predictions of the five folds were concatenated into one out-of-fold score per gene, and the
> whole procedure was repeated [TODO: 5] times with different random seeds; the out-of-fold score of a gene is the mean across repeats.
> Each gene is thus scored by models that never used it for training, and the final ranking and all top-K statistics are computed from
> these out-of-fold scores only. We report the area under the ROC curve (AUROC), the area under the precision–recall curve (AUPRC),
> average precision, Precision@K and Recall@K for K = 10, 50 and 100, and the enrichment factor (Precision@K divided by the proportion of
> positives); 95% confidence intervals were obtained by a stratified bootstrap over genes (positives and unlabeled genes resampled
> separately, 1,000 resamples, fixed seed). Because unlabeled genes contain undiscovered positives, precision and AUROC computed against
> LCGene are conservative estimates.

**Results sentence (replace "98 of the 100 top-ranked genes are known LCGene positives"):**

> Using out-of-fold scores, the model achieved an AUROC of [TODO] (95% CI [TODO]) and an AUPRC of [TODO] ([TODO]) against a positive
> prevalence of [TODO]; Precision@100 was [TODO] ([TODO]), Recall@100 was [TODO] ([TODO]) and the enrichment factor at 100 was [TODO]
> ([TODO]). [TODO: n] of the 100 highest-ranked genes were LCGene positives that had not been seen by the model that scored them.

**Abstract sentence (replace any performance claim based on training positives):**

> In cross-fitted (out-of-fold) evaluation, GLOOM ranked held-out LCGene genes with an AUROC of [TODO] (95% CI [TODO]) and an AUPRC of [TODO].

---

## 5. Ablation (feature contribution)

**New subsection (Methods):**

> **Feature-set ablation.** To test whether network integration improves prioritization, we compared, under the same cross-fitting
> protocol (identical folds, model and bootstrap), seven scorers: (a) |log2 fold change| alone; (b) adjusted p-value alone;
> (c) expression-derived features only; (d) network-derived features only; (e) expression and network features combined; (f) the
> combined set without the explicit differential-expression statistics (log2 fold change, |log2 fold change|, adjusted p-value
> and effect-size derived columns); and (g) the full model. Differences between models were assessed with a paired bootstrap over genes
> (the same resampled genes for both models; 1,000 resamples) that gives a 95% confidence interval and a two-sided p-value for the
> difference in AUPRC and AUROC. Variant (f) removes explicit fold-change columns only; per-group mean expression is still
> present and implicitly encodes fold change.

**Results — choose ONE variant:**

> *[Variant: no measurable gain]* The full model did not outperform the expression-only model: ΔAUPRC = [TODO] (95% CI [TODO]; p = [TODO])
> and ΔAUROC = [TODO] ([TODO]). Network-only features reached an AUROC of [TODO] ([TODO]), and the differential-expression baselines
> [TODO: |log2FC|, adjusted p] reached [TODO] and [TODO]. Network-derived features therefore primarily support contextual interpretation
> of prioritized genes, whereas predictive performance in the LUAD case study was driven predominantly by expression-derived variables.
>
> *[Variant: measurable gain]* The full model outperformed the expression-only model: ΔAUPRC = [TODO] (95% CI [TODO]; p = [TODO]) and ΔAUROC =
> [TODO] ([TODO]), indicating that network-derived features add predictive information beyond expression and differential-expression
> statistics. The gain was [TODO: small / moderate] in absolute terms; network-only features reached an AUROC of [TODO] ([TODO]).

**Discussion sentence (replace "network integration confers added value"):**

> The value of the network layer in this study lies [TODO: in the interpretation of prioritized genes (hubs, neighbourhoods, differential
> connectivity) / in a measurable improvement in AUPRC of X]; we do not claim that it is the main source of predictive performance.

---

## 6. Network stability

**New subsection (Methods):**

> **Network stability.** Co-expression networks were built from Pearson correlations of the log2-transformed expression and rebuilt at
> |r| ≥ 0.50, 0.55, 0.60, 0.65 and 0.70. For each threshold we report the number of edges, density, mean degree, size of the largest
> connected component, mean clustering coefficient, number of isolated nodes, and the Spearman correlation of degree centrality
> between consecutive thresholds. Because the number of samples influences the amount of spurious correlation, the analysis was repeated
> after randomly sub-sampling the larger group to the size of the smaller group (equal-N). Stability to sampling was assessed with
> [TODO: 20] bootstrap resamples of the samples: for each resample the network was rebuilt at |r| ≥ [TODO: 0.60], and we recorded the
> retention frequency of each reference edge and the Jaccard index between the hubs (top 5% of genes by degree) of the reference network and of each resample.

**Results:**

> The tumor and normal networks contained [TODO: n] and [TODO: n] edges at |r| ≥ [TODO: 0.60] (ratio [TODO]); the ratio was [TODO: range]
> across thresholds and [TODO] after equal-N sub-sampling (Table [TODO: SX]). Degree centralities of consecutive thresholds were correlated
> (Spearman ρ = [TODO: range]). In the bootstrap, reference edges were retained in [TODO: %] of resamples on average ([TODO: %] of edges in ≥ 80%
> of resamples) and hub sets had a mean Jaccard index of [TODO] ± [TODO]. [TODO — choose: The difference in network density between tumor and
> normal samples was robust to the threshold and to sample size and is reported as differential co-expression. / The difference in density
> depended on the threshold and/or the sample size and is therefore not interpreted as disease-specific rewiring.]

---

## 7. Enrichment analysis (Supplementary Method S11)

**Replace Supplementary Method S11 by:**

> **Pathway enrichment.** Over-representation of the [TODO: n] high-confidence non-LCGene candidates (and of their up- and down-regulated
> subsets) was tested for every KEGG pathway of the KEGG_2021_Human library that contained at least one gene of the analysis universe, with a
> one-sided hypergeometric test. The background was the [TODO: n] genes actually analysed by GLOOM, and pathway members outside the
> background were ignored. P-values were adjusted across all tested pathways with the Benjamini–Hochberg procedure. The complete
> unfiltered results (term, number of query genes, number of pathway genes in the background, overlap, raw and adjusted p-value) are provided in
> Supplementary Table [TODO: SX]. For interpretive purposes, we additionally extracted pathways containing predefined lung- or
> cancer-related terms. This secondary display was not used to determine statistical significance or establish biological validity.

**Results:**

> Of [TODO: n] pathways tested, [TODO: n] reached an adjusted p < 0.05 with at least three overlapping genes; the most enriched were [TODO].
> [TODO: Among these, pathways containing lung- or cancer-related terms were X (secondary display, Supplementary Table SX).]

---

## 8. Candidate terminology and external evidence

**Terminology — apply globally (find/replace, then re-read each sentence):**

| Replace | By |
|---|---|
| novel candidates / novel candidate genes | non-LCGene candidate genes (genes absent from the LCGene reference set that are prioritized by GLOOM) |
| novel genes, newly identified genes | non-LCGene candidate genes |
| identification of novel LUAD genes | prioritization of non-LCGene candidate genes |
| "novelty" | (delete; or "absence from the LCGene reference set") |

**New Methods paragraph:**

> **External evidence for non-LCGene candidates.** For each of the top [TODO: 200] non-LCGene candidates we queried the Open Targets
> Platform (association scores of the target with lung adenocarcinoma [EFO_0000571], lung carcinoma [EFO_0001071] and cancer; indirect
> associations included) and Europe PMC (number of records that mention the gene together with "lung adenocarcinoma", "lung cancer" or
> "cancer"). Genes were classified as *already associated with LUAD*, *associated with lung cancer without LUAD specificity*,
> *associated with other cancers*, *indirect/weak evidence only* or *no association found*, using an Open Targets score ≥ [TODO: 0.05]
> or ≥ [TODO: 3] Europe PMC records as the evidence threshold. Europe PMC counts are text co-mentions, not curated associations, and may be
> inflated for gene symbols that are common words. Only genes without any association were considered *potentially novel*, with respect to these sources only.

**Results:**

> Of the [TODO: 200] candidates, [TODO: n] were already associated with LUAD (including [TODO: e.g. CDC20, MYBL2, MARCO]), [TODO: n] with lung
> cancer without LUAD specificity, [TODO: n] with other cancers, [TODO: n] had indirect or weak evidence only, and [TODO: n] had no association
> in the searched resources (Table [TODO: X]); only the latter are proposed as potentially novel hypotheses for experimental follow-up.

---

## 9. Limitations — cohort confounding and related points

**Add to Limitations (and to the Discussion where the original design was discussed):**

> **Cohort confounding in public-resource comparisons.** In an earlier version of the analysis, tumors from TCGA (via cBioPortal, RSEM) were compared
> with normal lung from GTEx (TPM). Because every tumor came from one resource and every control from another, tumor status was perfectly
> confounded with cohort, sample preparation, quantification pipeline and normalization. The resulting contrast showed extreme asymmetry
> (9,715 of 10,986 genes differentially expressed, 9,700 up and 15 down), which we regard as largely technical. Batch-correction methods such as
> ComBat and ComBat-seq cannot resolve this situation, since the batch variable equals the biological group; they would remove disease signal
> together with the technical signal. The analysis was therefore repeated with TCGA-LUAD primary tumors and adjacent normal tissue processed by the same
> pipeline. Adjacent normal tissue is not entirely normal (field effects) and the number of adjacent normals is limited ([TODO: n]); GLOOM
> reports sample-level diagnostics and a cohort-collinearity warning so that users can detect this problem in their own data.

**Further Limitations sentences:**

> Labels were taken from LCGene, which captures expression-based LUAD biomarkers and is incomplete; unlabeled genes include undiscovered positives, so
> reported precision and AUROC are conservative. The contribution of network-derived features to prediction was [TODO: not demonstrated / modest] (see Ablation).
> No quantitative benchmark against existing prioritization tools was performed. Candidate genes require experimental validation.

---

## 10. Comparison with existing tools (Option B — replace the Discussion paragraph claiming advantage)

**Use exactly:**

> GLOOM differs from established prioritization tools in workflow architecture rather than being demonstrated here to outperform them. Its distinguishing
> feature is the integration of user-supplied expression matrices, differential-expression analysis, disease-conditioned co-expression networks, feature
> engineering, PU learning, and reporting within a single reproducible workflow. Quantitative head-to-head benchmarking remains an important direction for future work.

**Also delete / reword:** every occurrence of "methodological advantage", "practical advantage", "outperforms", "superior to" with respect to Endeavour, ToppGene,
GeneMANIA or BioRank (check the Abstract, Introduction, Discussion and Conclusion).

---

## 11. "Scalable" → "modular"

| Replace | By |
|---|---|
| a scalable framework / scalable pipeline | a modular framework / a modular pipeline |
| scalable and reproducible | modular and designed for reproducible execution |
| scalability | extensibility (new data types and disease contexts can be added as steps) |

> *[Only if the runtime benchmark is reported]* Run time and peak memory were measured on [TODO: machine: CPU, cores, RAM, OS] for
> [TODO: 2,000 / 5,000 / 10,000] genes and [TODO: 100 / 300 / 500] samples (Supplementary Table [TODO: SX]); the full LUAD analysis ([TODO: n] genes × [TODO: n]
> samples) required [TODO: min] and [TODO: GB] of memory. These measurements describe the sizes tested and are not extrapolated.

---

## 12. Availability, installation and reproducibility

**Replace the installation sentence (remove the Bioconda command):**

> GLOOM is implemented in Python (≥ 3.12) and is available at https://github.com/omicscodeathon/gloom under the MIT licence. Version [TODO: 0.2.0] used in this
> study is archived on Zenodo (DOI: [TODO: 10.5281/zenodo.XXXXXXX]). The bundled data files are stored with Git LFS (`git lfs install && git lfs pull`);
> installation is performed with `mamba env create -f environment.yml` followed by `pip install -e .`. Dependencies, including XGBoost, are declared in
> `pyproject.toml` and `environment.yml`.

**Add to Methods (reproducibility):**

> All random components (data splits, cross-validation folds, PU-bagging subsamples, random forests, bootstrap resampling, network layouts) are controlled by a single
> seed (42), which can be changed in the configuration file. The software includes unit tests on synthetic data for the zero-handling, cross-fitting, metric,
> bootstrap and enrichment functions. Software versions used: Python [TODO], scikit-learn [TODO], xgboost [TODO], pandas [TODO], numpy [TODO], scipy [TODO],
> networkx [TODO] (`pip freeze` in Supplementary File [TODO]).

---

## 13. Abstract and Conclusion — what to change

- Remove: the "98 of the top 100" statement, "novel candidates", "scalable", any claim of advantage over other tools, the statement that network integration adds value (unless the ablation supports it).
- Add: the cohort design (TCGA-LUAD tumor vs adjacent normal, GDC STAR), out-of-fold performance with 95% CIs, the ablation conclusion, and the external-evidence classification ("[TODO: n] of the top [TODO: 200] non-LCGene candidates had no association in Open Targets / Europe PMC").
- Suggested positioning sentence (Conclusion): *GLOOM is an integrated, modular and reproducible workflow whose main contributions are transparent preprocessing, out-of-sample evaluation of positive-unlabeled gene prioritization, and structured interpretation through network and pathway analyses; its predictive performance in the LUAD case study was [TODO: driven predominantly by expression-derived variables / improved by network features by X].*
