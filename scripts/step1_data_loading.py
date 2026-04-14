"""
step1_data_loading.py
---------------------
Loads all raw input files into clean pandas DataFrames and saves
validated intermediate copies to the processed directory.

Handles:
  - cBioPortal RSEM tumor expression (tab-separated, Hugo_Symbol + Entrez cols)
  - cBioPortal clinical metadata (has comment lines starting with '#')
  - GTEx GCT v1.2 gzipped normal expression (3-line header)
  - GTEx sample attributes metadata (TSV: GTEx_Analysis_v11_Annotations_SampleAttributesDS_-_LUAD.txt)
  - LCGene LUAD labeled gene list (TSV, GeneSymbol column)

GTEx metadata filtering strategy (based on actual file contents):
  - Total rows           : 1892
  - RNA:Total RNA rows   : 1467
  - SMAFRZE == RNASEQ    :  604  ← the 604 bulk RNA-seq lung samples
  - SMLRNA               :  546  ← small RNA, excluded
  - EXCLUDE / WGS / WES  :  excluded
  Applying BOTH filters (SMAFRZE==RNASEQ AND ANALYTE_TYPE==RNA:Total RNA)
  yields exactly 604 clean RNA-seq samples whose SAMPIDs match the GTEx
  expression matrix columns.

Outputs (saved to data/processed/):
  - tumor_expression_raw.csv
  - tumor_metadata_raw.csv
  - normal_expression_raw.csv
  - normal_metadata_raw.csv    ← 604 RNA-seq samples, key QC columns retained
  - cancer_genes_raw.csv
  - lcgene_luad_full.csv       ← full LCGene table with Regulation + log2FC
"""

import gzip
import logging
import sys
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config

config.create_output_dirs()

logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(config.LOG_FILE),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)

# ------------------------------------------------------------------
# Columns to retain from the GTEx metadata after filtering.
# Keeps sample ID + tissue confirmation + the most useful RNA-seq
# QC metrics; drops the ~100 low-level alignment/sequencing columns
# that are not needed downstream.
# ------------------------------------------------------------------
NORMAL_META_KEEP_COLS = [
    "SAMPID",       # sample identifier (becomes index)
    "SMTS",         # broad tissue type  (should all be "Lung")
    "SMTSD",        # detailed tissue type
    "SMRIN",        # RNA Integrity Number  (higher = better quality)
    "SMATSSCR",     # autolysis score       (0=none … 3=severe; lower = better)
    "SMRDTTL",      # total reads
    "SMMAPRT",      # overall mapping rate
    "SMRRNART",     # ribosomal RNA rate    (lower = better)
    "SMRDLGTH",     # read length
    "SMNABTCHT",    # batch type (for batch-effect awareness)
    "SMGEBTCHT",    # genotyping batch type
]


def load_tumor_expression(path: Path) -> pd.DataFrame:
    log.info(f"Loading tumor expression from: {path}")
    if not path.exists():
        raise FileNotFoundError(f"Tumor expression file not found: {path}")
    df = pd.read_csv(path, sep="\t", low_memory=False)
    log.info(f"  Raw shape: {df.shape}")
    if config.TUMOR_GENE_ID_COL not in df.columns:
        raise ValueError(f"Expected gene ID column '{config.TUMOR_GENE_ID_COL}' not found.")
    # Drop Entrez ID column - not needed downstream
    entrez_col = "Entrez_Gene_Id"
    if entrez_col in df.columns:
        df = df.drop(columns=[entrez_col])
    df = df.set_index(config.TUMOR_GENE_ID_COL)
    df.index.name = "gene"
    # Drop rows with blank Hugo_Symbol (Entrez-only rows that cannot be
    # matched to GTEx gene symbols)
    blank_mask = df.index.astype(str).str.strip() == ""
    if blank_mask.sum() > 0:
        log.warning(
            f"  Dropping {blank_mask.sum()} genes with blank Hugo_Symbol "
            f"(Entrez-only rows). They cannot be matched to GTEx gene symbols."
        )
        df = df.loc[~blank_mask]
    df = df.apply(pd.to_numeric, errors="coerce")
    nan_count = df.isna().sum().sum()
    if nan_count > 0:
        log.warning(f"  {nan_count} NaN values detected after numeric coercion.")
    if df.index.duplicated().any():
        n_dup = df.index.duplicated().sum()
        log.warning(f"  {n_dup} duplicate gene symbols found - keeping highest-mean row.")
        df["_mean"] = df.mean(axis=1)
        df = df.sort_values("_mean", ascending=False)
        df = df[~df.index.duplicated(keep="first")]
        df = df.drop(columns=["_mean"])
    log.info(f"  Final shape: {df.shape}")
    return df


def load_tumor_metadata(path: Path) -> pd.DataFrame:
    log.info(f"Loading tumor metadata from: {path}")
    if not path.exists():
        raise FileNotFoundError(f"Tumor metadata file not found: {path}")
    df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    log.info(f"  Raw shape: {df.shape}")
    sample_col = config.TUMOR_META_SAMPLE_COL
    if sample_col not in df.columns:
        col_map = {c.upper(): c for c in df.columns}
        if sample_col.upper() in col_map:
            df = df.rename(columns={col_map[sample_col.upper()]: sample_col})
        else:
            raise ValueError(f"Sample ID column '{sample_col}' not found.")
    df = df.set_index(sample_col)
    df.index.name = "sample_id"
    df = df.dropna(how="all").dropna(axis=1, how="all")
    log.info(f"  Final shape: {df.shape}")
    return df


def load_normal_expression(path: Path) -> pd.DataFrame:
    log.info(f"Loading normal expression from: {path}")
    if not path.exists():
        raise FileNotFoundError(f"Normal expression file not found: {path}")
    open_fn = gzip.open if str(path).endswith(".gz") else open
    with open_fn(path, "rt") as fh:
        version_line = fh.readline().strip()
        dim_line     = fh.readline().strip()
        log.info(f"  GCT version line : {version_line}")
        log.info(f"  GCT dimension line: {dim_line}")
    df = pd.read_csv(
        path, sep="\t", skiprows=2, low_memory=False,
        compression="gzip" if str(path).endswith(".gz") else None,
    )
    log.info(f"  Raw shape (after preamble skip): {df.shape}")
    gene_col = config.NORMAL_GENE_ID_COL   # "Description"
    name_col = "Name"                       # Ensembl ID column - drop it
    if gene_col not in df.columns:
        raise ValueError(f"Expected gene symbol column '{gene_col}' not found.")
    if name_col in df.columns:
        df = df.drop(columns=[name_col])
    df = df.set_index(gene_col)
    df.index.name = "gene"
    df = df.apply(pd.to_numeric, errors="coerce")
    if df.index.duplicated().any():
        n_dup = df.index.duplicated().sum()
        log.warning(f"  {n_dup} duplicate gene symbols - keeping highest-mean row.")
        df["_mean"] = df.mean(axis=1)
        df = df.sort_values("_mean", ascending=False)
        df = df[~df.index.duplicated(keep="first")]
        df = df.drop(columns=["_mean"])
    log.info(f"  Final shape: {df.shape}")
    return df


def load_normal_metadata(path: Path) -> pd.DataFrame:
    """
    Load GTEx_Analysis_v11_Annotations_SampleAttributesDS_-_LUAD.txt.

    File facts (verified from actual data):
      - 1892 total rows (samples across all assay types)
      - 1467 rows with ANALYTE_TYPE == 'RNA:Total RNA'
      - 604  rows with SMAFRZE == 'RNASEQ'  (bulk RNA-seq only)
      - Applying both filters gives exactly 604 RNA-seq samples whose
        SAMPIDs match the column names of the GTEx expression matrix.

    Filtering steps applied here:
      1. SMAFRZE == 'RNASEQ'         → keeps bulk RNA-seq, drops small-RNA,
                                        WGS, WES, DEEPWGS, OMNI, EXCLUDE
      2. ANALYTE_TYPE == 'RNA:Total RNA' → secondary guard against any
                                        non-RNA rows with RNASEQ freeze label

    Only a curated subset of columns (NORMAL_META_KEEP_COLS) is retained
    to keep the saved file compact and readable.  Key QC fields kept:
      SMRIN      - RNA Integrity Number (quality indicator; higher = better)
      SMATSSCR   - autolysis score (0=none, 3=severe; lower = better)
      SMMAPRT    - mapping rate
      SMRRNART   - ribosomal RNA fraction (lower = better)
      SMRDTTL    - total read count
    """
    log.info(f"Loading normal metadata from: {path}")
    if not path.exists():
        raise FileNotFoundError(f"Normal metadata file not found: {path}")

    # Plain TSV - NOT Excel
    df = pd.read_csv(path, sep="\t", low_memory=False)
    log.info(f"  Raw shape: {df.shape}  (all assay types)")

    sample_col = config.NORMAL_META_SAMPLE_COL  # "SAMPID"
    if sample_col not in df.columns:
        col_map = {c.upper(): c for c in df.columns}
        if sample_col.upper() in col_map:
            df = df.rename(columns={col_map[sample_col.upper()]: sample_col})
        else:
            raise ValueError(f"Sample ID column '{sample_col}' not found.")

    # ---- Filter 1: bulk RNA-seq freeze only ----
    if "SMAFRZE" in df.columns:
        n_before = len(df)
        df = df[df["SMAFRZE"] == "RNASEQ"].copy()
        log.info(f"  SMAFRZE=='RNASEQ'         : {n_before:>4} → {len(df):>4} samples "
                 f"(dropped small-RNA / WGS / WES / OMNI / EXCLUDE)")
    else:
        log.warning("  Column 'SMAFRZE' not found - skipping freeze filter.")

    # ---- Filter 2: RNA analyte type guard ----
    if "ANALYTE_TYPE" in df.columns:
        n_before = len(df)
        df = df[df["ANALYTE_TYPE"] == "RNA:Total RNA"].copy()
        log.info(f"  ANALYTE_TYPE=='RNA:Total RNA': {n_before:>4} → {len(df):>4} samples")
    else:
        log.warning("  Column 'ANALYTE_TYPE' not found - skipping analyte filter.")

    # ---- Retain only the useful columns ----
    keep = [c for c in NORMAL_META_KEEP_COLS if c in df.columns]
    dropped = set(NORMAL_META_KEEP_COLS) - set(keep)
    if dropped:
        log.warning(f"  Expected metadata columns not found (skipped): {sorted(dropped)}")
    # Keep all columns if none of the curated list is present (safety fallback)
    if keep:
        df = df[keep].copy()
        log.info(f"  Retaining {len(keep)} QC columns: {keep}")

    # ---- Set index ----
    if sample_col in df.columns:
        df = df.set_index(sample_col)
    df.index.name = "sample_id"

    # ---- Report QC column summaries ----
    if "SMRIN" in df.columns:
        log.info(f"  SMRIN  (RNA integrity) : mean={df['SMRIN'].mean():.2f}  "
                 f"min={df['SMRIN'].min():.1f}  max={df['SMRIN'].max():.1f}")
    if "SMATSSCR" in df.columns:
        log.info(f"  SMATSSCR (autolysis)   : {df['SMATSSCR'].value_counts().sort_index().to_dict()}")
    if "SMMAPRT" in df.columns:
        log.info(f"  SMMAPRT  (map rate)    : mean={df['SMMAPRT'].mean():.4f}  "
                 f"min={df['SMMAPRT'].min():.4f}")

    log.info(f"  Final shape: {df.shape}  ({len(df)} RNA-seq samples × {df.shape[1]} QC columns)")
    return df


def load_cancer_genes(path: Path) -> pd.Series:
    """
    Load the labeled LUAD gene list from LCGene_human_LUAD_filtered.tsv.

    Reads a tab-separated LCGene file and extracts gene symbols from the
    'GeneSymbol' column.  Also saves the full table (with Regulation
    direction and log2FoldChange) to processed/lcgene_luad_full.csv so
    downstream label-construction steps can use directionality if needed.
    """
    log.info(f"Loading labeled LUAD gene list from: {path}")
    if not path.exists():
        raise FileNotFoundError(f"Labeled gene file not found: {path}")

    # LCGene file is tab-separated
    df = pd.read_csv(path, sep="\t", low_memory=False)
    log.info(f"  Raw shape: {df.shape}")

    # Primary gene-symbol column in LCGene TSV
    gene_col = "GeneSymbol"
    if gene_col not in df.columns:
        col_map = {c.strip().upper(): c for c in df.columns}
        alt = "GENESYMBOL"
        if alt in col_map:
            gene_col = col_map[alt]
        else:
            raise ValueError(
                f"Expected column 'GeneSymbol' not found in {path}. "
                f"Available columns: {list(df.columns)}"
            )

    genes = (
        df[gene_col]
        .dropna()
        .str.strip()
        .str.upper()
        .drop_duplicates()
        .reset_index(drop=True)
    )
    genes.name = "gene"
    log.info(f"  Unique labeled LUAD genes loaded: {len(genes)}")

    # Save full LCGene table with upper-cased symbols for downstream use
    lcgene_out = config.PROCESSED_DIR / "lcgene_luad_full.csv"
    df_save = df.copy()
    df_save[gene_col] = df_save[gene_col].astype(str).str.strip().str.upper()
    df_save.to_csv(lcgene_out, index=False)
    log.info(f"  Full LCGene table saved to: {lcgene_out}")

    return genes


def run_data_loading() -> dict:
    log.info("=" * 60)
    log.info("STEP 1 - DATA LOADING")
    log.info("=" * 60)
    config.validate_input_files()
    tumor_expr   = load_tumor_expression(config.TUMOR_EXPR_FILE)
    tumor_meta   = load_tumor_metadata(config.TUMOR_META_FILE)
    normal_expr  = load_normal_expression(config.NORMAL_EXPR_FILE)
    normal_meta  = load_normal_metadata(config.NORMAL_META_FILE)
    cancer_genes = load_cancer_genes(config.CANCER_GENE_FILE)
    log.info("\nSaving raw-loaded files to processed directory …")
    tumor_expr.to_csv(config.PROCESSED_DIR / "tumor_expression_raw.csv")
    tumor_meta.to_csv(config.PROCESSED_DIR / "tumor_metadata_raw.csv")
    normal_expr.to_csv(config.PROCESSED_DIR / "normal_expression_raw.csv")
    normal_meta.to_csv(config.PROCESSED_DIR / "normal_metadata_raw.csv")
    cancer_genes.to_csv(config.PROCESSED_DIR / "cancer_genes_raw.csv", index=False, header=True)
    log.info("\n" + "=" * 60)
    log.info("STEP 1 SUMMARY")
    log.info("=" * 60)
    log.info(f"  Tumor expression   : {tumor_expr.shape[0]:>6} genes  x {tumor_expr.shape[1]:>4} samples")
    log.info(f"  Normal expression  : {normal_expr.shape[0]:>6} genes  x {normal_expr.shape[1]:>4} samples")
    log.info(f"  Normal metadata    : {normal_meta.shape[0]:>6} RNA-seq samples  x {normal_meta.shape[1]:>3} QC columns")
    log.info(f"  Labeled LUAD genes : {len(cancer_genes):>6} unique gene symbols (LCGene)")
    log.info("STEP 1 COMPLETE")
    return {
        "tumor_expr"   : tumor_expr,
        "tumor_meta"   : tumor_meta,
        "normal_expr"  : normal_expr,
        "normal_meta"  : normal_meta,
        "cancer_genes" : cancer_genes,
    }


if __name__ == "__main__":
    data = run_data_loading()
    print("\n--- Tumor expression (first 3 genes, first 4 samples) ---")
    print(data["tumor_expr"].iloc[:3, :4])
    print("\n--- Normal metadata (first 3 rows) ---")
    print(data["normal_meta"].iloc[:3])
    print(f"\n  Total RNA-seq samples: {data['normal_meta'].shape[0]}")
    print("\n--- Labeled LUAD genes (first 10) ---")
    print(data["cancer_genes"].head(10).tolist())
