"""
step1_data_loading.py
---------------------
Loads all raw input files into clean pandas DataFrames and saves
validated intermediate copies to the processed directory.

Handles:
  - cBioPortal RSEM tumor expression (tab-separated, Hugo_Symbol + Entrez cols)
  - cBioPortal clinical metadata (has comment lines starting with '#')
  - GTEx GCT v1.2 gzipped normal expression (3-line header)
  - GTEx sample attributes metadata (Excel)
  - Cancer Gene Census CSV

Outputs (saved to data/processed/):
  - tumor_expression_raw.csv
  - tumor_metadata_raw.csv
  - normal_expression_raw.csv
  - normal_metadata_raw.csv
  - cancer_genes_raw.csv
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


def load_tumor_expression(path: Path) -> pd.DataFrame:
    log.info(f"Loading tumor expression from: {path}")
    if not path.exists():
        raise FileNotFoundError(f"Tumor expression file not found: {path}")
    df = pd.read_csv(path, sep="\t", low_memory=False)
    log.info(f"  Raw shape: {df.shape}")
    if config.TUMOR_GENE_ID_COL not in df.columns:
        raise ValueError(f"Expected gene ID column '{config.TUMOR_GENE_ID_COL}' not found.")
    entrez_col = "Entrez_Gene_Id"
    if entrez_col in df.columns:
        df = df.drop(columns=[entrez_col])
    df = df.set_index(config.TUMOR_GENE_ID_COL)
    df.index.name = "gene"
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
    df = pd.read_csv(path, sep="\t", skiprows=2, low_memory=False,
                     compression="gzip" if str(path).endswith(".gz") else None)
    log.info(f"  Raw shape (after preamble skip): {df.shape}")
    gene_col = config.NORMAL_GENE_ID_COL
    name_col = "Name"
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
    log.info(f"Loading normal metadata from: {path}")
    if not path.exists():
        raise FileNotFoundError(f"Normal metadata file not found: {path}")
    df = pd.read_excel(path, sheet_name=0, engine="openpyxl")
    log.info(f"  Raw shape: {df.shape}")
    sample_col = config.NORMAL_META_SAMPLE_COL
    if sample_col not in df.columns:
        col_map = {c.upper(): c for c in df.columns}
        if sample_col.upper() in col_map:
            df = df.rename(columns={col_map[sample_col.upper()]: sample_col})
        else:
            log.warning(f"  Column '{sample_col}' not found - returning as-is.")
            return df
    df = df.set_index(sample_col)
    df.index.name = "sample_id"
    df = df.dropna(how="all").dropna(axis=1, how="all")
    log.info(f"  Final shape: {df.shape}")
    return df


def load_cancer_genes(path: Path) -> pd.Series:
    log.info(f"Loading cancer gene list from: {path}")
    if not path.exists():
        raise FileNotFoundError(f"Cancer gene file not found: {path}")
    df = pd.read_csv(path, low_memory=False)
    log.info(f"  Raw shape: {df.shape}")
    gene_col = config.CGC_GENE_COL
    if gene_col not in df.columns:
        col_map = {c.strip().upper(): c for c in df.columns}
        if gene_col.upper() in col_map:
            gene_col = col_map[gene_col.upper()]
        else:
            raise ValueError(f"Gene symbol column '{config.CGC_GENE_COL}' not found.")
    genes = (
        df[gene_col]
        .dropna()
        .str.strip()
        .str.upper()
        .drop_duplicates()
        .reset_index(drop=True)
    )
    genes.name = "gene"
    log.info(f"  Unique cancer genes loaded: {len(genes)}")
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
    log.info(f"  Cancer genes       : {len(cancer_genes):>6} unique symbols")
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
    print("\n--- Cancer genes (first 10) ---")
    print(data["cancer_genes"].head(10).tolist())
