"""
gloom/cli.py
------------
Entry point for the gloom CLI.

Usage:
    gloom info
    gloom diseases
    gloom validate  --data-dir /path/to/data
    gloom cache clear [--output DIR]
    gloom prioritize --genes genes.txt --disease luad --output results/
    gloom prioritize --genes genes.txt --disease luad --output results/ --from-step 11
    gloom prioritize --genes genes.txt --output results/ --fdr 0.05 --log2fc 1.0 --prob-threshold 0.5
    gloom prioritize --genes genes.txt --output results/ --no-cache
    gloom prioritize --genes genes.txt --output results/ --labels my_known_genes.csv
    gloom prioritize --genes genes.txt --output results/ --dry-run
    gloom run --tumor-expr tumor.csv --normal-expr normal.csv --output results/
    gloom run --tumor-expr tumor.csv --normal-expr normal.csv \\
              --tumor-meta tumor_meta.csv --normal-meta normal_meta.csv \\
              --genes candidates.txt --labels known_genes.tsv --output results/
"""
import shutil
import sys
import logging
import time
from pathlib import Path

# Force UTF-8 on Windows consoles (cp1252/cp1256) so Unicode characters
# in pipeline log messages (→ ≥ ± ├── …) don't raise UnicodeEncodeError.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

import click

# LUAD_ML_Pipeline sits one level above this package
_PIPELINE_DIR = Path(__file__).resolve().parent / "pipeline"

# ── Step registry (number, module, function) ─────────────────────────────────
_STEPS = [
    (0,  "config",                            "create_output_dirs"),
    (1,  "step1_data_loading",                "run_data_loading"),
    (2,  "step2_preprocessing",               "run_preprocessing"),
    (3,  "step3_harmonization",               "run_harmonization"),
    (4,  "step4_differential_expression",     "run_differential_expression"),
    (5,  "step5_expression_features",         "run_expression_feature_construction"),
    (6,  "step6_coexpression_network",        "run_coexpression_network"),
    (7,  "step7_network_features",            "run_network_feature_extraction"),
    (8,  "step8_feature_integration",         "run_feature_integration"),
    (9,  "step9_label_construction",          "run_label_construction"),
    (10, "step10_train_val_split",            "run_train_val_split"),
    (11, "step11_model_training",             "run_model_training"),
    (12, "step12_model_evaluation",           "run_model_evaluation"),
    (13, "step13_feature_importance",         "run_feature_importance"),
    (14, "step14_gene_ranking",               "run_gene_ranking"),
    (15, "step15_network_annotation",         "run_network_annotation"),
    (16, "step16_network_export",             "run_network_export"),
    (17, "step17_interactive_visualization",  "run_interactive_visualization"),
    (18, "step18_final_report",               "run_final_report"),
    (19, "step19_kegg_enrichment",            "run_kegg_enrichment"),
]

# ── Disease registry ──────────────────────────────────────────────────────────
_DISEASES = {
    "luad": "Lung Adenocarcinoma (LUAD) — TCGA/cBioPortal tumour + GTEx normal",
}


# ── CLI root ──────────────────────────────────────────────────────────────────

@click.group()
@click.version_option()
def main():
    """gloom — Gene prioritization pipeline for lung adenocarcinoma (LUAD)."""
    pass


# ── gloom prioritize ──────────────────────────────────────────────────────────

@main.command()
@click.option(
    "--genes", required=True, type=click.Path(exists=True), metavar="FILE",
    help="Gene list: one symbol per line (.txt) or LCGene-format TSV.",
)
@click.option(
    "--disease", default="luad", show_default=True,
    type=click.Choice(["luad"], case_sensitive=False),
    help="Disease context (only 'luad' is currently supported).",
)
@click.option(
    "--output", "output_dir", required=True, type=click.Path(), metavar="DIR",
    help="Output directory (created if it does not exist).",
)
@click.option(
    "--data-dir", default=None, type=click.Path(exists=True), metavar="DIR",
    help="Root data directory containing a raw/ subfolder. "
         "Defaults to ../data relative to the pipeline folder.",
)
@click.option(
    "--from-step", default=0, type=click.IntRange(0, 19), show_default=True,
    help="Resume pipeline from this step number (0–19).",
)
@click.option(
    "--to-step", default=19, type=click.IntRange(0, 19), show_default=True,
    help="Stop after this step number (0–19).",
)
@click.option("--verbose", is_flag=True, help="Enable DEBUG-level logging.")
@click.option(
    "--top-k", default=None, type=click.IntRange(1), metavar="N",
    help="Keep only the top N candidates in ranked_candidates output (by predicted probability).",
)
@click.option(
    "--format", "fmt", default="csv", show_default=True,
    type=click.Choice(["csv", "excel", "json"], case_sensitive=False),
    help="Output file format for candidate tables.",
)
@click.option(
    "--fdr", default=None, type=click.FloatRange(0.0, 1.0), metavar="FLOAT",
    help="Adjusted p-value (FDR) cutoff for differential expression (default: 0.05).",
)
@click.option(
    "--log2fc", "log2fc_threshold", default=None, type=click.FloatRange(0.0), metavar="FLOAT",
    help="Log2 fold-change threshold for differential expression (default: 1.0).",
)
@click.option(
    "--prob-threshold", default=None, type=click.FloatRange(0.0, 1.0), metavar="FLOAT",
    help="Minimum ML probability to flag a gene as a novel candidate (default: 0.5).",
)
@click.option(
    "--no-cache", is_flag=True,
    help="Delete any existing intermediate cache before running (forces full re-run from --from-step).",
)
@click.option(
    "--labels", "labels_file", default=None, type=click.Path(exists=True), metavar="FILE",
    help="Custom gene label file (CSV/TSV with a GeneSymbol column) replacing the built-in LCGene database.",
)
@click.option(
    "--dry-run", is_flag=True,
    help="Validate inputs and print the execution plan without running the pipeline.",
)
def prioritize(genes, disease, output_dir, data_dir, from_step, to_step, verbose, top_k, fmt,
               fdr, log2fc_threshold, prob_threshold, no_cache, labels_file, dry_run):
    """Prioritize genes for LUAD using ML and co-expression networks.

    \b
    Examples:
        gloom prioritize --genes genes.txt --disease luad --output results/
        gloom prioritize --genes genes.txt --output results/ --from-step 11
        gloom prioritize --genes genes.txt --output results/ --top-k 200
        gloom prioritize --genes genes.txt --output results/ --format excel
        gloom prioritize --genes genes.txt --output results/ --fdr 0.05 --log2fc 1.0
        gloom prioritize --genes genes.txt --output results/ --prob-threshold 0.6
        gloom prioritize --genes genes.txt --output results/ --no-cache
        gloom prioritize --genes genes.txt --output results/ --labels my_known_genes.csv
        gloom prioritize --genes genes.txt --output results/ --dry-run
    """
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    _setup_logging(verbose, output_path)
    log = logging.getLogger("gloom")
    log.info(f"gloom {__import__('gloom').__version__}  disease={disease}")
    log.info(f"Output directory : {output_path}")

    # 1. Inject pipeline directory into sys.path and patch config
    sys.path.insert(0, str(_PIPELINE_DIR))
    import config  # noqa: PLC0415

    # 2. Handle --no-cache before patching
    if no_cache:
        cache_dir = output_path / ".pipeline_cache"
        if cache_dir.exists():
            shutil.rmtree(cache_dir)
            log.info("Cache cleared (--no-cache).")
        else:
            log.info("--no-cache: no existing cache found.")

    _patch_config(
        config=config,
        data_dir=Path(data_dir).resolve() if data_dir else None,
        output_path=output_path,
        genes_file=Path(genes).resolve(),
        fdr=fdr,
        log2fc_threshold=log2fc_threshold,
        prob_threshold=prob_threshold,
        labels_file=Path(labels_file).resolve() if labels_file else None,
    )

    # 3. Dry-run: validate + print plan, then exit
    if dry_run:
        _print_dry_run(config, from_step, to_step, output_path)
        return

    # 4. Run pipeline steps
    completed, failed = _run_steps(from_step, to_step, log)

    # 5. Print execution summary
    _print_summary(completed, log)

    if failed:
        log.error(f"Pipeline stopped at step {failed}.")
        sys.exit(1)

    # 6. Reorganize outputs into clean user-facing structure
    log.info("Organizing outputs …")
    from gloom.output import organize_outputs  # noqa: PLC0415
    organize_outputs(config, output_path, disease=disease, top_k=top_k, fmt=fmt)

    # 7. Final status
    ext = "xlsx" if fmt == "excel" else fmt
    click.echo("")
    click.echo(f"  Results ready : {output_path}")
    click.echo(f"  Main result   : {output_path / 'candidates' / f'ranked_candidates.{ext}'}")
    click.echo(f"  Report        : {output_path / 'report.html'}")


# ── gloom validate ────────────────────────────────────────────────────────────

@main.command()
@click.option(
    "--data-dir", required=True, type=click.Path(exists=True), metavar="DIR",
    help="Root data directory to validate (must contain raw/ subfolder).",
)
def validate(data_dir):
    """Check that all required input files are present before running."""
    sys.path.insert(0, str(_PIPELINE_DIR))
    import config  # noqa: PLC0415

    data_path = Path(data_dir).resolve()
    config.DATA_ROOT = data_path
    config.RAW_DIR   = data_path / "raw"
    config.TUMOR_EXPR_FILE  = config.RAW_DIR / "cBioPortal (RNA Seq Data)" / "data_mrna_seq_v2_rsem.txt"
    config.TUMOR_META_FILE  = config.RAW_DIR / "cBioPortal (RNA Seq Data)" / "data_clinical_patient.txt"
    config.NORMAL_EXPR_FILE = config.RAW_DIR / "Gtex (normal samples)" / "gene_tpm_v11_lung.gct.gz"
    config.NORMAL_META_FILE = config.RAW_DIR / "Gtex (normal samples)" / \
                              "GTEx_Analysis_v11_Annotations_SampleAttributesDS - LUAD.txt"
    config.CANCER_GENE_FILE = config.RAW_DIR / "LCGene (Labeled LUAD Data)" / \
                              "LCGene_human_LUAD_filtered.tsv"
    try:
        config.validate_input_files()
        click.echo("All required input files found.")
    except FileNotFoundError as exc:
        click.echo(str(exc), err=True)
        sys.exit(1)


# ── gloom info ────────────────────────────────────────────────────────────────

@main.command()
def info():
    """Show version, bundled reference paths, and any cached pipeline runs."""
    import gloom as _gloom

    click.echo(f"gloom version  : {_gloom.__version__}")
    click.echo(f"Python         : {sys.version.split()[0]}")
    click.echo(f"Platform       : {sys.platform}")
    click.echo("")

    # Reference data status
    sys.path.insert(0, str(_PIPELINE_DIR))
    try:
        import config  # noqa: PLC0415
        click.echo("Bundled reference files:")
        refs = [
            ("Tumor expression",  config.TUMOR_EXPR_FILE),
            ("Tumor metadata",    config.TUMOR_META_FILE),
            ("Normal expression", config.NORMAL_EXPR_FILE),
            ("Normal metadata",   config.NORMAL_META_FILE),
            ("LCGene labels",     config.CANCER_GENE_FILE),
        ]
        for label, path in refs:
            status = "FOUND  " if Path(path).exists() else "MISSING"
            click.echo(f"  [{status}] {label:<20} {path}")
        click.echo("")
        click.echo("Default thresholds:")
        click.echo(f"  FDR (adj p-value)  : {config.DE_PVALUE_THRESHOLD}  (override: --fdr)")
        click.echo(f"  Log2 fold-change   : {config.DE_LOG2FC_THRESHOLD}  (override: --log2fc)")
        click.echo(f"  Correlation cutoff : {config.COEXPR_CORRELATION_CUTOFF}")
        click.echo(f"  Novel prob cutoff  : 0.50  (override: --prob-threshold)")
        click.echo(f"  CV folds           : {config.CV_FOLDS}")
        click.echo(f"  CV primary metric  : {config.CV_METRIC_PRIMARY.upper()}")
    except Exception as exc:
        click.echo(f"  (Could not load config: {exc})")

    click.echo("")
    # Cached pipeline runs under cwd
    cache_dirs = sorted(Path(".").glob("*/.pipeline_cache"))
    if cache_dirs:
        click.echo("Cached pipeline runs (current directory):")
        for c in cache_dirs:
            click.echo(f"  {c.parent.resolve()}")
    else:
        click.echo("No cached pipeline runs found in the current directory.")


# ── gloom diseases ────────────────────────────────────────────────────────────

@main.command()
def diseases():
    """List supported --disease values and their data sources."""
    click.echo("Supported diseases:")
    click.echo("")
    for key, desc in _DISEASES.items():
        click.echo(f"  {key:<10} {desc}")
    click.echo("")
    click.echo("Pass the key to --disease when running 'gloom prioritize'.")
    click.echo("Additional disease models are planned for future releases.")


# ── gloom cache ───────────────────────────────────────────────────────────────

@main.group()
def cache():
    """Manage cached intermediate pipeline files."""
    pass


@cache.command("clear")
@click.option(
    "--output", "output_dir", default=None, type=click.Path(), metavar="DIR",
    help="Output directory whose cache to clear. "
         "Defaults to searching all subdirectories of the current directory.",
)
@click.confirmation_option(prompt="Delete all cached intermediate files?")
def cache_clear(output_dir):
    """Delete cached intermediate files to force a full pipeline re-run."""
    if output_dir:
        targets = [Path(output_dir).resolve() / ".pipeline_cache"]
    else:
        targets = sorted(Path(".").glob("*/.pipeline_cache"))

    if not targets:
        click.echo("No cache directories found.")
        return

    deleted = 0
    for cache_path in targets:
        if cache_path.exists():
            shutil.rmtree(cache_path)
            click.echo(f"Deleted: {cache_path}")
            deleted += 1
        else:
            click.echo(f"Not found (skipped): {cache_path}")

    click.echo(f"\n{deleted} cache director{'y' if deleted == 1 else 'ies'} deleted.")


# ── gloom run ────────────────────────────────────────────────────────────────

@main.command("run")
@click.option(
    "--tumor-expr", required=True, type=click.Path(exists=True), metavar="FILE",
    help="Tumor expression matrix: CSV with gene symbols as row index, sample IDs as columns.",
)
@click.option(
    "--normal-expr", required=True, type=click.Path(exists=True), metavar="FILE",
    help="Normal expression matrix: CSV with gene symbols as row index, sample IDs as columns.",
)
@click.option(
    "--tumor-meta", default=None, type=click.Path(exists=True), metavar="FILE",
    help="Tumor sample metadata: CSV with sample IDs as row index. "
         "If omitted, all samples in --tumor-expr are used.",
)
@click.option(
    "--normal-meta", default=None, type=click.Path(exists=True), metavar="FILE",
    help="Normal sample metadata: CSV with sample IDs as row index. "
         "If omitted, all samples in --normal-expr are used.",
)
@click.option(
    "--output", "output_dir", required=True, type=click.Path(), metavar="DIR",
    help="Output directory (created if it does not exist).",
)
@click.option(
    "--genes", default=None, type=click.Path(exists=True), metavar="FILE",
    help="Optional query gene list: one symbol per line (.txt) or TSV with GeneSymbol column. "
         "These genes are scored as candidates. If omitted, all genes are ranked.",
)
@click.option(
    "--labels", "labels_file", default=None, type=click.Path(exists=True), metavar="FILE",
    help="Known positive gene list for training labels: CSV/TSV with a GeneSymbol column. "
         "Falls back to the bundled LCGene database if omitted.",
)
@click.option(
    "--from-step", default=2, type=click.IntRange(2, 19), show_default=True,
    help="Resume from this step (minimum 2 — step 1 is always bypassed in gloom run).",
)
@click.option(
    "--to-step", default=19, type=click.IntRange(2, 19), show_default=True,
    help="Stop after this step number (2–19).",
)
@click.option(
    "--fdr", default=None, type=click.FloatRange(0.0, 1.0), metavar="FLOAT",
    help="Adjusted p-value (FDR) cutoff for differential expression (default: 0.05).",
)
@click.option(
    "--log2fc", "log2fc_threshold", default=None, type=click.FloatRange(0.0), metavar="FLOAT",
    help="Log2 fold-change threshold for differential expression (default: 1.0).",
)
@click.option(
    "--prob-threshold", default=None, type=click.FloatRange(0.0, 1.0), metavar="FLOAT",
    help="Minimum ML probability to flag a gene as a novel candidate (default: 0.5).",
)
@click.option(
    "--top-k", default=None, type=click.IntRange(1), metavar="N",
    help="Keep only the top N candidates in ranked_candidates output.",
)
@click.option(
    "--format", "fmt", default="csv", show_default=True,
    type=click.Choice(["csv", "excel", "json"], case_sensitive=False),
    help="Output file format for candidate tables.",
)
@click.option("--verbose", is_flag=True, help="Enable DEBUG-level logging.")
def run(tumor_expr, normal_expr, tumor_meta, normal_meta, output_dir,
        genes, labels_file, from_step, to_step,
        fdr, log2fc_threshold, prob_threshold, top_k, fmt, verbose):
    """Run the full pipeline with your own expression data files.

    Step 1 (data loading) is bypassed — your files are staged directly
    into the pipeline cache and preprocessing starts at Step 2.
    Steps 2–19 (QC, harmonization, DE, network, ML, ranking) run unchanged.

    \b
    Input file format:
      Expression matrices : CSV, gene symbols as row index, sample IDs as columns
      Metadata files      : CSV, sample IDs as row index (any extra columns OK)
      Labels file         : CSV/TSV with a GeneSymbol column (known positive genes)

    \b
    Examples:
        gloom run --tumor-expr tumor.csv --normal-expr normal.csv --output results/

        gloom run \\
          --tumor-expr  tumor_expr.csv \\
          --normal-expr normal_expr.csv \\
          --tumor-meta  tumor_meta.csv \\
          --normal-meta normal_meta.csv \\
          --output results/

        gloom run \\
          --tumor-expr  tumor.csv --normal-expr normal.csv \\
          --genes       my_candidates.txt \\
          --labels      known_genes.tsv \\
          --fdr 0.05 --log2fc 1.0 --prob-threshold 0.5 \\
          --output results/
    """
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    _setup_logging(verbose, output_path)
    log = logging.getLogger("gloom")
    log.info(f"gloom {__import__('gloom').__version__}  mode=run (custom data)")
    log.info(f"Output directory : {output_path}")

    # 1. Validate expression file shapes before doing anything else
    _validate_expression_file(Path(tumor_expr),  "tumor-expr")
    _validate_expression_file(Path(normal_expr), "normal-expr")

    # 2. Inject pipeline dir and patch config (no data-dir, no step-1 paths needed)
    sys.path.insert(0, str(_PIPELINE_DIR))
    import config  # noqa: PLC0415

    _patch_config(
        config=config,
        data_dir=None,
        output_path=output_path,
        genes_file=Path(genes).resolve() if genes else None,
        fdr=fdr,
        log2fc_threshold=log2fc_threshold,
        prob_threshold=prob_threshold,
        labels_file=Path(labels_file).resolve() if labels_file else None,
    )

    # 3. Stage user files into the cache as Step 1 outputs
    log.info("Staging user-supplied data files …")
    _stage_user_files(
        config=config,
        tumor_expr_path=Path(tumor_expr).resolve(),
        normal_expr_path=Path(normal_expr).resolve(),
        tumor_meta_path=Path(tumor_meta).resolve() if tumor_meta else None,
        normal_meta_path=Path(normal_meta).resolve() if normal_meta else None,
        labels_file=Path(labels_file).resolve() if labels_file else None,
        log=log,
    )

    # 4. Run pipeline from step 2 onward (step 1 is bypassed)
    completed, failed = _run_steps(from_step, to_step, log)

    # 5. Print execution summary
    _print_summary(completed, log)

    if failed:
        log.error(f"Pipeline stopped at step {failed}.")
        sys.exit(1)

    # 6. Organize outputs
    log.info("Organizing outputs …")
    from gloom.output import organize_outputs  # noqa: PLC0415
    organize_outputs(config, output_path, disease="luad", top_k=top_k, fmt=fmt)

    ext = "xlsx" if fmt == "excel" else fmt
    click.echo("")
    click.echo(f"  Results ready : {output_path}")
    click.echo(f"  Main result   : {output_path / 'candidates' / f'ranked_candidates.{ext}'}")
    click.echo(f"  Report        : {output_path / 'report.html'}")


# ── Helpers ───────────────────────────────────────────────────────────────────

def _setup_logging(verbose: bool, output_path: Path) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    log_file = output_path / ".pipeline_cache" / "logs" / "pipeline.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )


def _patch_config(config, data_dir, output_path: Path,
                  genes_file: Path | None = None,
                  fdr=None, log2fc_threshold=None, prob_threshold=None,
                  labels_file=None) -> None:
    """Redirect all config paths into a hidden cache inside output_path."""
    cache = output_path / ".pipeline_cache"

    # ── Data source paths ────────────────────────────────────────────────────
    if data_dir:
        config.DATA_ROOT = data_dir
        config.RAW_DIR   = data_dir / "raw"
        # Rebuild the raw file paths from the new RAW_DIR
        config.TUMOR_EXPR_FILE  = config.RAW_DIR / "cBioPortal (RNA Seq Data)" / "data_mrna_seq_v2_rsem.txt"
        config.TUMOR_META_FILE  = config.RAW_DIR / "cBioPortal (RNA Seq Data)" / "data_clinical_patient.txt"
        config.NORMAL_EXPR_FILE = config.RAW_DIR / "Gtex (normal samples)" / "gene_tpm_v11_lung.gct.gz"
        config.NORMAL_META_FILE = config.RAW_DIR / "Gtex (normal samples)" / \
                                  "GTEx_Analysis_v11_Annotations_SampleAttributesDS - LUAD.txt"

    # ── Custom label file (overrides built-in LCGene) ────────────────────────
    if labels_file is not None:
        config.CANCER_GENE_FILE = labels_file

    # ── DE / scoring thresholds ──────────────────────────────────────────────
    if fdr is not None:
        config.DE_PVALUE_THRESHOLD = fdr
    if log2fc_threshold is not None:
        config.DE_LOG2FC_THRESHOLD = log2fc_threshold
    if prob_threshold is not None:
        config.NOVEL_PROB_THRESHOLD = prob_threshold

    # ── Genes file: these are QUERY genes (unknown LUAD relationship) ────────
    # The LCGene database (config.CANCER_GENE_FILE) is always used for
    # training labels. The user's gene list is stored separately as
    # QUERY_GENES_FILE so the model can score them as candidates.
    # genes_file is None when the user does not pass --genes (e.g. gloom run).
    if genes_file is not None:
        if genes_file.suffix.lower() == ".txt":
            import pandas as pd  # noqa: PLC0415
            symbols = [
                ln.strip() for ln in genes_file.read_text().splitlines()
                if ln.strip() and not ln.startswith("#")
            ]
            tmp_tsv = cache / "user_genes.tsv"
            tmp_tsv.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame({"GeneSymbol": symbols}).to_csv(tmp_tsv, sep="\t", index=False)
            config.QUERY_GENES_FILE = tmp_tsv
        else:
            config.QUERY_GENES_FILE = genes_file
    # NOTE: config.CANCER_GENE_FILE is intentionally left unchanged unless
    # --labels was passed — LCGene always provides the training labels by default.

    # ── Redirect all output paths into cache ─────────────────────────────────
    p   = cache / "processed"
    r   = cache / "results"
    f   = cache / "figures"
    m   = cache / "models"
    n   = cache / "results" / "network"
    rp  = cache / "results" / "reports"
    e   = cache / "results" / "enrichment"
    lg  = cache / "logs"

    config.PROCESSED_DIR            = p
    config.OUTPUTS_ROOT             = cache
    config.RESULTS_DIR              = r
    config.FIGURES_DIR              = f
    config.MODELS_DIR               = m
    config.LOGS_DIR                 = lg
    config.NETWORK_DIR              = n
    config.REPORTS_DIR              = rp
    config.ENRICHMENT_DIR           = e

    config.TUMOR_EXPR_PROCESSED     = p / "tumor_expression_processed.csv"
    config.NORMAL_EXPR_PROCESSED    = p / "normal_expression_processed.csv"
    config.TUMOR_EXPR_HARMONIZED    = p / "tumor_expression_harmonized.csv"
    config.NORMAL_EXPR_HARMONIZED   = p / "normal_expression_harmonized.csv"
    config.DE_RESULTS_FILE          = r / "differential_expression_results.csv"
    config.EXPR_FEATURES_FILE       = p / "expression_features.csv"
    config.NETWORK_EDGES_FILE       = n / "coexpression_network_edges.csv"
    config.NETWORK_GRAPH_FILE       = n / "coexpression_network.graphml"
    config.NETWORK_FEATURES_FILE    = p / "network_features.csv"
    config.INTEGRATED_FEATURES_FILE = p / "integrated_features.csv"
    config.LABELS_FILE              = p / "gene_labels.csv"
    config.TRAIN_FEATURES_FILE      = p / "train_features.csv"
    config.VAL_FEATURES_FILE        = p / "val_features.csv"
    config.TRAIN_LABELS_FILE        = p / "train_labels.csv"
    config.VAL_LABELS_FILE          = p / "val_labels.csv"
    config.GENE_RANKINGS_FILE       = r / "gene_rankings.csv"
    config.FEATURE_IMPORTANCE_FILE  = r / "feature_importance.csv"
    config.MODEL_METRICS_FILE       = r / "model_metrics.csv"
    config.ANNOTATED_NETWORK_FILE   = n / "annotated_network.graphml"
    config.ANNOTATED_EDGES_FILE     = n / "annotated_edges.csv"
    config.ANNOTATED_NODES_FILE     = n / "annotated_nodes.csv"
    config.INTERACTIVE_HTML_FILE    = f / "interactive_network.html"
    config.FINAL_REPORT_FILE        = rp / "pipeline_summary_report.csv"
    config.KEGG_ALL_FILE            = e / "kegg_all_candidates.csv"
    config.KEGG_UP_FILE             = e / "kegg_upregulated.csv"
    config.KEGG_DOWN_FILE           = e / "kegg_downregulated.csv"
    config.KEGG_SUMMARY_FILE        = e / "kegg_summary.csv"
    config.LOG_FILE                 = lg / "pipeline.log"

    config.create_output_dirs()


def _print_dry_run(config, from_step: int, to_step: int, output_path: Path) -> None:
    """Print validation results and execution plan without running anything."""
    click.echo("")
    click.echo("DRY RUN — nothing will be executed")
    click.echo("=" * 60)

    # Input file validation
    click.echo("\nInput files:")
    inputs = [
        ("Tumor expression",  config.TUMOR_EXPR_FILE),
        ("Tumor metadata",    config.TUMOR_META_FILE),
        ("Normal expression", config.NORMAL_EXPR_FILE),
        ("Normal metadata",   config.NORMAL_META_FILE),
        ("Gene labels (ref)", config.CANCER_GENE_FILE),
        ("Query genes",       getattr(config, "QUERY_GENES_FILE", None)),
    ]
    all_ok = True
    for label, path in inputs:
        if path is None:
            click.echo(f"  [MISSING ] {label}")
            all_ok = False
            continue
        found = Path(path).exists()
        status = "FOUND  " if found else "MISSING"
        click.echo(f"  [{status}] {label:<22} {path}")
        if not found:
            all_ok = False

    # Effective thresholds
    click.echo("\nEffective thresholds:")
    click.echo(f"  FDR (adj p-value)   : {config.DE_PVALUE_THRESHOLD}")
    click.echo(f"  Log2 fold-change    : {config.DE_LOG2FC_THRESHOLD}")
    click.echo(f"  Novel prob cutoff   : {getattr(config, 'NOVEL_PROB_THRESHOLD', 0.50)}")

    # Steps that would run
    click.echo(f"\nSteps that would run (--from-step {from_step} to --to-step {to_step}):")
    for step_num, module_name, _ in _STEPS:
        if from_step <= step_num <= to_step:
            click.echo(f"  Step {step_num:>2}  {module_name}")

    # Output path
    click.echo(f"\nOutput directory : {output_path}")
    click.echo("")
    if all_ok:
        click.echo("All input files found. Ready to run.")
    else:
        click.echo("WARNING: Some input files are missing. Run 'gloom validate' for details.")
    click.echo("=" * 60)


def _run_steps(from_step: int, to_step: int, log) -> tuple[list, str | None]:
    completed = []
    failed_step = None

    for step_num, module_name, func_name in _STEPS:
        if step_num < from_step or step_num > to_step:
            continue

        label = f"Step {step_num:>2} — {module_name}"
        log.info(f">>> {label}")
        t0 = time.time()
        try:
            module = __import__(module_name)
            func   = getattr(module, func_name)
            func()
            elapsed = time.time() - t0
            log.info(f"<<< {label}  DONE ({elapsed:.1f}s)")
            completed.append((label, elapsed, "OK"))
        except Exception as exc:
            elapsed = time.time() - t0
            log.error(f"<<< {label}  FAILED: {type(exc).__name__}: {exc}")
            completed.append((label, elapsed, f"FAILED: {exc}"))
            failed_step = label
            break

    return completed, failed_step


def _print_summary(completed: list, log) -> None:
    log.info("=" * 60)
    log.info("EXECUTION SUMMARY")
    log.info("=" * 60)
    total = sum(e for _, e, _ in completed)
    for label, elapsed, status in completed:
        icon = "OK  " if status == "OK" else "FAIL"
        log.info(f"  [{icon}] {label:<50} {elapsed:>7.1f}s")
    log.info(f"  Total: {total:.1f}s  ({total / 60:.1f} min)")
    log.info("=" * 60)


def _validate_expression_file(path: Path, flag: str) -> None:
    """Quick sanity-check that a user-supplied expression CSV is readable and non-empty."""
    import pandas as pd  # noqa: PLC0415
    try:
        df = pd.read_csv(path, index_col=0, nrows=5)
    except Exception as exc:
        raise click.BadParameter(
            f"Cannot read file: {exc}", param_hint=f"--{flag}"
        )
    if df.shape[1] == 0:
        raise click.BadParameter(
            "File has no sample columns — expected a genes × samples CSV "
            "(gene symbols as row index, sample IDs as column headers).",
            param_hint=f"--{flag}",
        )
    if df.shape[0] == 0:
        raise click.BadParameter(
            "File has no gene rows.", param_hint=f"--{flag}"
        )


def _stage_user_files(config, tumor_expr_path: Path, normal_expr_path: Path,
                      tumor_meta_path: Path | None, normal_meta_path: Path | None,
                      labels_file: Path | None, log) -> None:
    """
    Write user-supplied files into the pipeline cache as Step 1 outputs so that
    Step 2 (preprocessing) can read them without any changes to the pipeline code.

    Expected input format
    ---------------------
    Expression CSVs : gene symbols as row index, sample IDs as column headers,
                      numeric expression values (any unit — Step 2 will log2-transform).
    Metadata CSVs   : sample IDs as row index, any additional columns are fine.
                      If omitted, a dummy metadata with all expression samples is created
                      so that Step 2's sample-alignment check passes cleanly.
    Labels file     : CSV/TSV with a GeneSymbol column.
                      Falls back to the bundled LCGene database when not provided.
    """
    import pandas as pd  # noqa: PLC0415

    proc = Path(config.PROCESSED_DIR)

    # ── Expression matrices ───────────────────────────────────────────────────
    log.info(f"  Reading tumor expression  : {tumor_expr_path.name}")
    tumor_expr = pd.read_csv(tumor_expr_path, index_col=0)
    tumor_expr.index.name = "gene"
    log.info(f"    Shape: {tumor_expr.shape[0]:,} genes × {tumor_expr.shape[1]:,} samples")
    tumor_expr.to_csv(proc / "tumor_expression_raw.csv")

    log.info(f"  Reading normal expression : {normal_expr_path.name}")
    normal_expr = pd.read_csv(normal_expr_path, index_col=0)
    normal_expr.index.name = "gene"
    log.info(f"    Shape: {normal_expr.shape[0]:,} genes × {normal_expr.shape[1]:,} samples")
    normal_expr.to_csv(proc / "normal_expression_raw.csv")

    # ── Metadata ─────────────────────────────────────────────────────────────
    if tumor_meta_path:
        log.info(f"  Reading tumor metadata    : {tumor_meta_path.name}")
        tumor_meta = pd.read_csv(tumor_meta_path, index_col=0)
        tumor_meta.index.name = "sample_id"
        log.info(f"    Shape: {tumor_meta.shape}")
    else:
        log.info("  No --tumor-meta provided — using all tumor expression samples.")
        tumor_meta = pd.DataFrame(index=pd.Index(tumor_expr.columns, name="sample_id"))
    tumor_meta.to_csv(proc / "tumor_metadata_raw.csv")

    if normal_meta_path:
        log.info(f"  Reading normal metadata   : {normal_meta_path.name}")
        normal_meta = pd.read_csv(normal_meta_path, index_col=0)
        normal_meta.index.name = "sample_id"
        log.info(f"    Shape: {normal_meta.shape}")
    else:
        log.info("  No --normal-meta provided — using all normal expression samples.")
        normal_meta = pd.DataFrame(index=pd.Index(normal_expr.columns, name="sample_id"))
    normal_meta.to_csv(proc / "normal_metadata_raw.csv")

    # ── Gene labels (training positives) ─────────────────────────────────────
    cancer_genes_dst = proc / "cancer_genes_raw.csv"

    if labels_file is not None:
        # User provided their own label file
        log.info(f"  Reading gene labels       : {labels_file.name}")
        sep = "\t" if labels_file.suffix.lower() in (".tsv", ".txt") else ","
        ldf = pd.read_csv(labels_file, sep=sep, low_memory=False)
        # Find the gene-symbol column (flexible naming)
        gene_col = next(
            (c for c in ldf.columns
             if c.strip().lower() in ("genesymbol", "gene_symbol", "gene", "symbol", "hgnc_symbol")),
            ldf.columns[0],
        )
        genes = (ldf[gene_col].dropna().astype(str)
                 .str.strip().str.upper().drop_duplicates().reset_index(drop=True))
        genes.name = "gene"
        genes.to_csv(cancer_genes_dst, index=False, header=True)
        log.info(f"    {len(genes):,} labeled genes loaded.")

    elif Path(config.CANCER_GENE_FILE).exists():
        # Fall back to the bundled LCGene database
        log.info(f"  Using bundled LCGene database: {Path(config.CANCER_GENE_FILE).name}")
        ldf = pd.read_csv(config.CANCER_GENE_FILE, sep="\t", low_memory=False)
        gene_col = getattr(config, "LCGENE_GENE_COL", "GeneSymbol")
        genes = (ldf[gene_col].dropna().astype(str)
                 .str.strip().str.upper().drop_duplicates().reset_index(drop=True))
        genes.name = "gene"
        genes.to_csv(cancer_genes_dst, index=False, header=True)
        # Also write the full LCGene table for downstream steps
        ldf_save = ldf.copy()
        ldf_save[gene_col] = ldf_save[gene_col].astype(str).str.strip().str.upper()
        ldf_save.to_csv(proc / "lcgene_luad_full.csv", index=False)
        log.info(f"    {len(genes):,} LCGene genes staged.")

    else:
        # No labels at all — write an empty list and warn loudly
        log.warning(
            "  No --labels file and no bundled LCGene database found.\n"
            "  All training labels will be 0 (no known positives).\n"
            "  The model will still run but gene ranking will not be meaningful.\n"
            "  Provide --labels <file> with a GeneSymbol column to fix this."
        )
        pd.Series([], name="gene", dtype=str).to_csv(
            cancer_genes_dst, index=False, header=True
        )

    log.info("  Staging complete — starting pipeline at Step 2.")
