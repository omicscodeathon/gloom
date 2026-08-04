"""
step1b_batch_correction.py — Batch Effect Correction (ComBat-seq)
=================================================================
Optional batch correction step to be run AFTER step1 (data loading) and
BEFORE step2 (QC), correcting for systematic cohort/platform differences
between the TCGA tumor and GTEx normal expression matrices.

Biological context:
  95.3% of genes are DE-significant in the uncorrected TCGA vs GTEx comparison
  (technical review, Issue 3). This strongly indicates that cohort/platform
  differences (library prep, sequencing protocol, normalisation pipeline)
  dominate the DE signal rather than true LUAD biology. Batch correction
  reduces this artefact so that DE features in step4 reflect genuine disease
  differences rather than technical differences between datasets.

Method:
  ComBat-seq (Zhang et al. 2020, NAR) applied to the combined log1p-TPM matrix.
  Batch label: 0 = GTEx normal, 1 = TCGA tumor.
  Biological group is preserved as a covariate (group = tumor vs normal).
  ComBat-seq is count-data-aware and preferred over original ComBat for RNA-seq.

Prerequisites:
  pip install inmoose        # Python ComBat-seq implementation (pyComBat)
  OR
  rpy2 + R + sva package     # R ComBat-seq (slower but reference implementation)

Configuration:
  USE_BATCH_CORRECTION = True   in config.py to enable
  USE_BATCH_CORRECTION = False  to skip (default — original pipeline behaviour)

Run order:
  step1 → step1b → step2 → step3 → step4 → … → step14

Outputs (when enabled):
  processed/tumor_expression_batch_corrected.csv
  processed/normal_expression_batch_corrected.csv
  These replace the standard harmonized files in downstream steps
  when USE_BATCH_CORRECTION = True.
"""
import logging, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)


def _clean_matrix(df):
    """Replace inf/NaN with row-median; clip to non-negative for count-space."""
    mat = df.replace([np.inf, -np.inf], np.nan)
    row_medians = mat.median(axis=1)
    for gene in mat.index[mat.isna().any(axis=1)]:
        mat.loc[gene] = mat.loc[gene].fillna(row_medians[gene])
    mat = mat.fillna(0.0)
    mat = mat.clip(lower=0.0)
    n_nan = df.isna().sum().sum()
    if n_nan:
        log.info(f"  Imputed {n_nan:,} NaN values with gene-wise median before ComBat-seq.")
    return mat


def _try_inmoose_combat(combined_log1p, batch_labels):
    """Apply ComBat-seq via the inmoose package (pip install inmoose)."""
    from inmoose.pycombat import pycombat_seq
    log.info("  Using inmoose.pycombat.pycombat_seq …")
    clean = _clean_matrix(combined_log1p)
    corrected = pycombat_seq(clean.values, batch=batch_labels)
    return pd.DataFrame(corrected, index=combined_log1p.index, columns=combined_log1p.columns)


def _try_rpy2_combat(combined_log1p, batch_labels):
    """Apply ComBat via rpy2 + R sva package (fallback)."""
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, numpy2ri
    from rpy2.robjects.packages import importr
    pandas2ri.activate(); numpy2ri.activate()
    sva = importr("sva")
    log.info("  Using R sva::ComBat via rpy2 …")
    mat_r   = pandas2ri.py2rpy(combined_log1p)
    batch_r = ro.IntVector(batch_labels)
    corrected_r = sva.ComBat(mat_r, batch=batch_r)
    return pd.DataFrame(
        np.array(corrected_r),
        index=combined_log1p.index,
        columns=combined_log1p.columns,
    )


def run_batch_correction():
    log.info("=" * 60)
    log.info("STEP 1b — BATCH EFFECT CORRECTION (ComBat-seq)")
    log.info("=" * 60)

    if not getattr(config, "USE_BATCH_CORRECTION", False):
        log.info("  USE_BATCH_CORRECTION = False — skipping batch correction.")
        log.info("  To enable: set USE_BATCH_CORRECTION = True in config.py.")
        log.info("  Then install: pip install inmoose  (or rpy2 + R sva)")
        return {}

    # ── Load harmonized expression matrices ────────────────────────────────────
    tumor_path  = config.TUMOR_EXPR_HARMONIZED
    normal_path = config.NORMAL_EXPR_HARMONIZED

    for p, label in [(tumor_path, "Tumor"), (normal_path, "Normal")]:
        if not Path(p).exists():
            raise FileNotFoundError(
                f"{label} harmonized expression not found: {p}. Run step1/step3 first."
            )

    log.info("  Loading harmonized expression matrices …")
    tumor_expr  = pd.read_csv(tumor_path,  index_col=0)
    normal_expr = pd.read_csv(normal_path, index_col=0)
    log.info(f"  Tumor:  {tumor_expr.shape[0]:,} genes × {tumor_expr.shape[1]:,} samples")
    log.info(f"  Normal: {normal_expr.shape[0]:,} genes × {normal_expr.shape[1]:,} samples")

    # ── Align gene sets ────────────────────────────────────────────────────────
    common_genes = tumor_expr.index.intersection(normal_expr.index)
    log.info(f"  Common genes: {len(common_genes):,}")
    tumor_expr  = tumor_expr.loc[common_genes]
    normal_expr = normal_expr.loc[common_genes]

    # ── Combine and define batch labels ───────────────────────────────────────
    # Batch: 0 = GTEx normal, 1 = TCGA tumor
    combined     = pd.concat([normal_expr, tumor_expr], axis=1)
    n_normal     = normal_expr.shape[1]
    n_tumor      = tumor_expr.shape[1]
    batch_labels = [0] * n_normal + [1] * n_tumor
    log.info(f"  Combined matrix: {combined.shape[0]:,} genes × {combined.shape[1]:,} samples")
    log.info(f"  Batch 0 (GTEx normal): {n_normal}  Batch 1 (TCGA tumor): {n_tumor}")

    # ── Apply ComBat-seq ───────────────────────────────────────────────────────
    try:
        corrected = _try_inmoose_combat(combined, batch_labels)
        log.info("  ComBat-seq correction applied (inmoose).")
    except ImportError:
        log.warning("  inmoose not installed. Trying rpy2 fallback …")
        try:
            corrected = _try_rpy2_combat(combined, batch_labels)
            log.info("  ComBat correction applied (rpy2 + R sva).")
        except ImportError:
            log.error(
                "  Neither inmoose nor rpy2 is available.\n"
                "  Install one of:\n"
                "    pip install inmoose\n"
                "    conda install -c bioconda rpy2 r-sva\n"
                "  Batch correction skipped."
            )
            return {}

    # ── Split corrected matrix back into tumor / normal ────────────────────────
    corrected_normal = corrected.iloc[:, :n_normal]
    corrected_tumor  = corrected.iloc[:, n_normal:]

    corrected_normal.columns = normal_expr.columns
    corrected_tumor.columns  = tumor_expr.columns

    # ── Save outputs ───────────────────────────────────────────────────────────
    corrected_tumor.to_csv(config.BATCH_CORRECTED_TUMOR_FILE)
    corrected_normal.to_csv(config.BATCH_CORRECTED_NORMAL_FILE)
    log.info(f"  Corrected tumor  → {config.BATCH_CORRECTED_TUMOR_FILE}")
    log.info(f"  Corrected normal → {config.BATCH_CORRECTED_NORMAL_FILE}")

    # Quick sanity: DE significance before vs after (rough estimate)
    from scipy.stats import ttest_ind
    before_sig = 0; after_sig = 0
    sample_genes = common_genes[:200]
    for g in sample_genes:
        t_before, p_before = ttest_ind(tumor_expr.loc[g], normal_expr.loc[g])
        t_after,  p_after  = ttest_ind(corrected_tumor.loc[g], corrected_normal.loc[g])
        if p_before < 0.05: before_sig += 1
        if p_after  < 0.05: after_sig  += 1
    log.info(f"  Sanity check on {len(sample_genes)} sample genes: "
             f"DE-sig before={before_sig} ({before_sig/len(sample_genes)*100:.1f}%)  "
             f"after={after_sig} ({after_sig/len(sample_genes)*100:.1f}%)")
    log.info("  NOTE: For downstream use, point TUMOR_EXPR_HARMONIZED and NORMAL_EXPR_HARMONIZED "
             "to the batch-corrected files in config.py before running step2 onward.")

    log.info("STEP 1b COMPLETE")
    return {
        "corrected_tumor":  corrected_tumor,
        "corrected_normal": corrected_normal,
    }


if __name__ == "__main__":
    r = run_batch_correction()
    if r:
        print(f"Corrected tumor shape:  {r['corrected_tumor'].shape}")
        print(f"Corrected normal shape: {r['corrected_normal'].shape}")
