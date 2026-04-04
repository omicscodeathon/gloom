
# Raw Data Directory


This directory contains the original, unmodified source data downloaded from public repositories. These files serve as the **input** to the preprocessing pipeline.


## Sources and Files


### 1. Cancer Gene Census (Labeled Data)

**Source**: COSMIC (Catalogue of Somatic Mutations in Cancer) Cancer Gene Census

**File**: `Census_allWed Mar 18 05_48_31 2026.csv`

**Purpose**: Provides a curated list of genes with documented roles in cancer. Used as a gold-standard set for labeling and feature prioritization.

**Note**: The filename includes a timestamp indicating the download date.


### 2. cBioPortal (RNA Seq Data) – **Tumor Samples**

**Source**: cBioPortal for Cancer Genomics (TCGA Lung Adenocarcinoma - LUAD)

**Files**:

-`data_mrna_seq_v2_rsem.txt`: RNA-Seq expression data (log2-transformed RSEM values) for tumor samples.

-`data_clinical_patient.txt`: Clinical metadata (patient age, stage, survival, etc.) for the tumor cohort.

**Purpose**: Primary source of lung adenocarcinoma tumor gene expression and associated clinical data.


### 3. GTEx (normal samples) – **Normal/Control Samples**

**Source**: GTEx (Genotype-Tissue Expression) Consortium, V11 release

**Files**:

-`gene_tpm_v11_lung.gct.gz`: Gene expression data (TPM values) for normal lung tissue samples.

-`GTEx_Analysis_v11_Annotations_SampleAttributesDD.xlsx`: Sample-level metadata (tissue details, processing batch).

-`GTEx_Analysis_v11_Annotations_SubjectPhenotypesDD.xlsx`: Donor/subject-level metadata (demographics).

**Purpose**: Source of healthy lung tissue expression data for comparison/control.


## Important Notes

-**Do Not Modify**: These files should never be edited manually. They are the reproducible starting point.

-**Format**: Files are in their native formats (.txt, .csv, .xlsx, .gct.gz).

-**License & Citation**: When using results from this pipeline, ensure proper citation of COSMIC, cBioPortal/TCGA, and GTEx according to their data usage policies.

-**Storage**: The `gene_tpm_v11_lung.gct.gz` is compressed. The pipeline script decompresses it during processing.


## Pipeline Input

The preprocessing scripts (`scripts/preprocessing/`) are configured to read directly from these files and paths.

