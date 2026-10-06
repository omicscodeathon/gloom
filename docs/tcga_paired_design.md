# TCGA-LUAD tumor vs adjacent-normal design (uniformly processed)

## Why

The original GLOOM LUAD case study compared **tumors from TCGA (via cBioPortal, RSEM)** with
**normal lung from GTEx (TPM)**. Tumor status was therefore perfectly collinear with cohort,
sample preparation, quantification pipeline and normalization. Batch correction (ComBat /
ComBat-seq) cannot separate such a cohort effect from the disease effect, and the resulting
contrast showed the typical signature of an artefact (9,700 genes up vs 15 down, global median
log2FC of about +4.4).

The cleanest remedy is to compare **TCGA-LUAD primary tumors** with **TCGA-LUAD adjacent
"Solid Tissue Normal" samples**, both quantified by the **same GDC STAR - Counts pipeline**.
The adjacent normals of TCGA-LUAD are paired with tumors from the same patients (59 normals at
the time of writing), which also allows a **paired** analysis.

A secondary comparison against a *uniformly reprocessed* TCGA-GTEx resource (for example UCSC Xena
TOIL "TcgaTargetGtex") can be added, but the pipeline, version and normalization of that resource
must be documented precisely.

## 1. Obtain the data from the GDC

### Option A — the helper script (recommended)

```bash
python scripts/fetch_gdc_tcga_luad.py --paired-only      # normals + their matched tumors (about 0.5 GB)
# or: python scripts/fetch_gdc_tcga_luad.py              # all primary tumors + all normals (about 2.5 GB)
```

The script queries the public GDC REST API (open-access files, no token needed):

| GDC filter | Value |
|---|---|
| `cases.project.project_id` | `TCGA-LUAD` |
| `data_type` | `Gene Expression Quantification` |
| `analysis.workflow_type` | `STAR - Counts` |
| `cases.samples.sample_type` | `Primary Tumor`, `Solid Tissue Normal` |

and downloads every file with `GET https://api.gdc.cancer.gov/data/<file_id>` (thread pool,
cached in `data/raw/tcga_gdc/files/`). Each STAR counts file is a ~4.2 MB tab-separated table with a comment line,
the columns `gene_id gene_name gene_type unstranded stranded_first stranded_second tpm_unstranded
fpkm_unstranded fpkm_uq_unstranded` and four `N_*` summary rows, which are skipped.

Outputs in `data/raw/tcga_gdc/`:

| File | Content |
|---|---|
| `counts_matrix.csv` | genes x samples, raw unstranded STAR counts (input for limma-voom / DESeq2 / edgeR) |
| `tpm_matrix.csv` | genes x samples, `tpm_unstranded` (default input of the GLOOM pipeline) |
| `sample_sheet.csv` | `sample_id`, `patient_id`, `group` (`tumor` / `normal`), `sample_type`, `file_id`, `file_name` |
| `gdc_files_manifest.csv` | the GDC file list (provenance) |

### Option B — manual download from the GDC Data Portal

1. Go to <https://portal.gdc.cancer.gov/> -> Repository.
2. Filters: *Project* = TCGA-LUAD; *Data Category* = Transcriptome Profiling; *Data Type* = Gene
   Expression Quantification; *Workflow Type* = STAR - Counts; *Sample Type* = Primary Tumor and
   Solid Tissue Normal.
3. Add all files to the cart, download the manifest and the sample sheet, and fetch the files with
   the `gdc-client` (`gdc-client download -m manifest.txt`).
4. Build `counts_matrix.csv`, `tpm_matrix.csv` and `sample_sheet.csv` in the format above
   (genes as symbols in the first column, one column per `sample_id`).

## 2. Configure GLOOM

In `config.py` (or the packaged `src/gloom/pipeline/config.py`):

```python
DATA_SOURCE           = "tcga_gdc"     # uses data/raw/tcga_gdc/{tpm_matrix.csv, sample_sheet.csv}
GDC_EXPRESSION_MATRIX = "tpm"          # "tpm" (default) or "counts"
DE_METHOD             = "paired"       # "welch" | "paired" | "limma_voom"
```

* `welch` — unpaired Welch t-test on log2(x + 1) values (default).
* `paired` — paired t-test on the patients that have both a tumor and an adjacent-normal sample
  (patient ids come from `sample_sheet.csv`, or from the TCGA barcode `TCGA-XX-XXXX`). If fewer than
  three pairs are found the step warns and falls back to Welch.
* `limma_voom` — limma-voom on the **raw counts** (`counts_matrix.csv`), with a patient blocking
  factor when pairs exist. It uses `rpy2` and R (packages `limma`, `edgeR`), which are **optional**:
  install with `pip install rpy2` and `BiocManager::install(c("limma", "edgeR"))`. When they are
  missing the step prints a clear message and falls back to Welch.

Leave `USE_BATCH_CORRECTION = False`: tumor and adjacent normal share one cohort and one pipeline.

## 3. Run and check

```bash
python src/gloom/pipeline/run_pipeline.py --from 1
```

Check `results/qc_cohort_warning.txt` (it should not exist), `results/qc_global_de_summary.csv`
(median log2FC close to 0, a plausible fraction of DE genes and up:down ratio) and
`figures/qc_sample_pca.png` (tumor and normal should not be separated by a purely technical axis).

## Notes and caveats

* Adjacent "normal" tissue is not perfectly normal (field cancerization); expect smaller fold changes than against GTEx.
* With only 59 normals, a model fitted on the paired subset has fewer samples than the original design. Report this limitation.
* Use the same design for the tumor and the normal co-expression networks; differences in sample size alter network density
  (see `step7c_network_stability.py`, `NETWORK_EQUAL_N`).
