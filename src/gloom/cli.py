"""
gloom/cli.py
------------
Entry point for the gloom CLI.

Usage:
    gloom info
    gloom diseases
    gloom validate --data-dir /path/to/data
    gloom cache clear [--output DIR]
    gloom prioritize --genes genes.txt --disease luad --output results/
    gloom prioritize --genes genes.txt --output results/ --from-step 11
    gloom prioritize --genes genes.txt --output results/ --skip-optional
    gloom prioritize --genes genes.txt --output results/ --skip-step 11b --skip-step 19
    gloom run --tumor-expr tumor.csv --normal-expr normal.csv --output results/
"""

import importlib
import logging
import re
import shutil
import sys
import time
from pathlib import Path

import click

# Force UTF-8 on Windows consoles (cp1252/cp1256) so Unicode characters
# in pipeline log messages do not raise UnicodeEncodeError.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")


_PIPELINE_DIR = Path(__file__).resolve().parent / "pipeline"
_DISEASES = {
    "luad": "Lung Adenocarcinoma (LUAD) - TCGA/cBioPortal tumour + GTEx normal",
}


def _ensure_pipeline_on_path() -> None:
    if str(_PIPELINE_DIR) not in sys.path:
        sys.path.insert(0, str(_PIPELINE_DIR))


def _load_pipeline_steps():
    """Use the packaged pipeline runner registry as the single step source."""
    _ensure_pipeline_on_path()
    runner = importlib.import_module("run_pipeline")
    return runner.STEPS


def _normalize_step_key(
    value: str | int | None,
    *,
    minimum: int | None = None,
    maximum: int | None = None,
) -> str | None:
    if value is None:
        return None
    raw = str(value).strip().lower()
    match = re.fullmatch(r"(\d+)(b?)", raw)
    if not match:
        raise click.BadParameter(
            f"Invalid step key '{value}'. Use keys like 4, 6b, 11, or 11b."
        )
    number = int(match.group(1))
    suffix = match.group(2)
    if minimum is not None and number < minimum:
        raise click.BadParameter(f"Step key '{value}' is below the minimum allowed step {minimum}.")
    if maximum is not None and number > maximum:
        raise click.BadParameter(f"Step key '{value}' is above the maximum allowed step {maximum}.")
    return f"{number:02d}{suffix}"


def _step_order_key(step_key: str) -> tuple[int, int]:
    match = re.fullmatch(r"(\d+)(b?)", step_key.lower())
    if not match:
        raise ValueError(f"Unsupported step key: {step_key}")
    return int(match.group(1)), 1 if match.group(2) else 0


def _select_steps(
    steps,
    *,
    from_key: str | None = None,
    to_key: str | None = None,
    skip_optional: bool = False,
    skip_keys: tuple[str, ...] = (),
):
    available = {step[0] for step in steps}
    requested = [key for key in (from_key, to_key, *skip_keys) if key is not None]
    invalid = [key for key in requested if key not in available]
    if invalid:
        raise click.BadParameter(f"Unknown step key(s): {', '.join(sorted(invalid))}")

    if from_key and to_key and _step_order_key(from_key) > _step_order_key(to_key):
        raise click.BadParameter("--from-step must be earlier than or equal to --to-step.")

    selected = list(steps)
    if from_key:
        start_idx = next(i for i, step in enumerate(selected) if step[0] == from_key)
        selected = selected[start_idx:]
    if to_key:
        end_idx = next(i for i, step in enumerate(selected) if step[0] == to_key)
        selected = selected[: end_idx + 1]
    if skip_optional:
        selected = [step for step in selected if not step[4]]
    if skip_keys:
        skip_set = set(skip_keys)
        selected = [step for step in selected if step[0] not in skip_set]
    if not selected:
        raise click.ClickException("No steps selected to run after applying the current filters.")
    return selected


@click.group()
@click.version_option()
def main():
    """gloom - Gene prioritization pipeline for lung adenocarcinoma (LUAD)."""


@main.command()
@click.option(
    "--genes",
    required=True,
    type=click.Path(exists=True),
    metavar="FILE",
    help="Gene list: one symbol per line (.txt) or LCGene-format TSV.",
)
@click.option(
    "--disease",
    default="luad",
    show_default=True,
    type=click.Choice(["luad"], case_sensitive=False),
    help="Disease context (only 'luad' is currently supported).",
)
@click.option(
    "--output",
    "output_dir",
    required=True,
    type=click.Path(),
    metavar="DIR",
    help="Output directory (created if it does not exist).",
)
@click.option(
    "--data-dir",
    default=None,
    type=click.Path(exists=True),
    metavar="DIR",
    help="Root data directory containing a raw/ subfolder. Defaults to ../data relative to the pipeline folder.",
)
@click.option(
    "--from-step",
    default="0",
    show_default=True,
    metavar="STEP",
    help="Resume pipeline from this step key (examples: 4, 6b, 11b).",
)
@click.option(
    "--to-step",
    default="19",
    show_default=True,
    metavar="STEP",
    help="Stop after this step key (examples: 8, 11, 19).",
)
@click.option(
    "--skip-optional",
    is_flag=True,
    help="Exclude all steps marked optional from the selected run.",
)
@click.option(
    "--skip-step",
    "skip_steps",
    multiple=True,
    metavar="STEP",
    help="Skip one specific step key. May be passed multiple times.",
)
@click.option("--verbose", is_flag=True, help="Enable DEBUG-level logging.")
@click.option(
    "--top-k",
    default=None,
    type=click.IntRange(1),
    metavar="N",
    help="Keep only the top N candidates in ranked_candidates output (by predicted probability).",
)
@click.option(
    "--format",
    "fmt",
    default="csv",
    show_default=True,
    type=click.Choice(["csv", "excel", "json"], case_sensitive=False),
    help="Output file format for candidate tables.",
)
@click.option(
    "--fdr",
    default=None,
    type=click.FloatRange(0.0, 1.0),
    metavar="FLOAT",
    help="Adjusted p-value (FDR) cutoff for differential expression (default: 0.05).",
)
@click.option(
    "--log2fc",
    "log2fc_threshold",
    default=None,
    type=click.FloatRange(0.0),
    metavar="FLOAT",
    help="Log2 fold-change threshold for differential expression (default: 1.0).",
)
@click.option(
    "--prob-threshold",
    default=None,
    type=click.FloatRange(0.0, 1.0),
    metavar="FLOAT",
    help="Minimum ML probability to flag a gene as a novel candidate (default: 0.5).",
)
@click.option(
    "--no-cache",
    is_flag=True,
    help="Delete any existing intermediate cache before running (forces a full re-run).",
)
@click.option(
    "--labels",
    "labels_file",
    default=None,
    type=click.Path(exists=True),
    metavar="FILE",
    help="Custom gene label file (CSV/TSV with a GeneSymbol column) replacing the built-in LCGene database.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Validate inputs and print the execution plan without running the pipeline.",
)
def prioritize(
    genes,
    disease,
    output_dir,
    data_dir,
    from_step,
    to_step,
    skip_optional,
    skip_steps,
    verbose,
    top_k,
    fmt,
    fdr,
    log2fc_threshold,
    prob_threshold,
    no_cache,
    labels_file,
    dry_run,
):
    """Prioritize genes for LUAD using ML and co-expression networks."""
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    cache_dir = output_path / ".pipeline_cache"
    cache_cleared = False
    if no_cache and cache_dir.exists():
        shutil.rmtree(cache_dir)
        cache_cleared = True

    _setup_logging(verbose, output_path)
    log = logging.getLogger("gloom")
    log.info(f"gloom {__import__('gloom').__version__}  disease={disease}")
    log.info(f"Output directory : {output_path}")
    if no_cache:
        if cache_cleared:
            log.info("Cache cleared (--no-cache).")
        else:
            log.info("--no-cache: no existing cache found.")

    _ensure_pipeline_on_path()
    import config  # noqa: PLC0415

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

    steps = _select_steps(
        _load_pipeline_steps(),
        from_key=_normalize_step_key(from_step, minimum=0, maximum=19),
        to_key=_normalize_step_key(to_step, minimum=0, maximum=19),
        skip_optional=skip_optional,
        skip_keys=tuple(
            _normalize_step_key(step, minimum=0, maximum=19) for step in skip_steps
        ),
    )

    if dry_run:
        _print_dry_run(
            config,
            steps,
            output_path,
            skip_optional=skip_optional,
            skip_step_keys=tuple(
                _normalize_step_key(step, minimum=0, maximum=19) for step in skip_steps
            ),
        )
        return

    completed, failed = _run_steps(steps, log)
    _print_summary(completed, log)

    if failed:
        log.error(f"Pipeline stopped at {failed}.")
        sys.exit(1)

    log.info("Organizing outputs ...")
    from gloom.output import organize_outputs  # noqa: PLC0415

    organize_outputs(
        config,
        output_path,
        disease=disease,
        top_k=top_k,
        fmt=fmt,
        include_models=False,
    )

    ext = "xlsx" if fmt == "excel" else fmt
    click.echo("")
    click.echo(f"  Results ready : {output_path}")
    click.echo(f"  Main result   : {output_path / 'candidates' / f'ranked_candidates.{ext}'}")
    click.echo(f"  Report        : {output_path / 'report.html'}")


@main.command()
@click.option(
    "--data-dir",
    required=True,
    type=click.Path(exists=True),
    metavar="DIR",
    help="Root data directory to validate (must contain raw/ subfolder).",
)
def validate(data_dir):
    """Check that all required input files are present before running."""
    _ensure_pipeline_on_path()
    import config  # noqa: PLC0415

    data_path = Path(data_dir).resolve()
    config.DATA_ROOT = data_path
    config.RAW_DIR = data_path / "raw"
    config.TUMOR_EXPR_FILE = config.RAW_DIR / "cBioPortal (RNA Seq Data)" / "data_mrna_seq_v2_rsem.txt"
    config.TUMOR_META_FILE = config.RAW_DIR / "cBioPortal (RNA Seq Data)" / "data_clinical_patient.txt"
    config.NORMAL_EXPR_FILE = config.RAW_DIR / "Gtex (normal samples)" / "gene_tpm_v11_lung.gct.gz"
    config.NORMAL_META_FILE = (
        config.RAW_DIR / "Gtex (normal samples)" / "GTEx_Analysis_v11_Annotations_SampleAttributesDS - LUAD.txt"
    )
    config.CANCER_GENE_FILE = (
        config.RAW_DIR / "LCGene (Labeled LUAD Data)" / "LCGene_human_LUAD_filtered.tsv"
    )

    try:
        config.validate_input_files()
        click.echo("All required input files found.")
    except FileNotFoundError as exc:
        click.echo(str(exc), err=True)
        sys.exit(1)


@main.command()
def info():
    """Show version, bundled reference paths, and any cached pipeline runs."""
    import gloom as _gloom

    click.echo(f"gloom version  : {_gloom.__version__}")
    click.echo(f"Python         : {sys.version.split()[0]}")
    click.echo(f"Platform       : {sys.platform}")
    click.echo("")

    _ensure_pipeline_on_path()
    try:
        import config  # noqa: PLC0415

        click.echo("Bundled reference files:")
        refs = [
            ("Tumor expression", config.TUMOR_EXPR_FILE),
            ("Tumor metadata", config.TUMOR_META_FILE),
            ("Normal expression", config.NORMAL_EXPR_FILE),
            ("Normal metadata", config.NORMAL_META_FILE),
            ("LCGene labels", config.CANCER_GENE_FILE),
        ]
        for label, path in refs:
            status = "FOUND  " if Path(path).exists() else "MISSING"
            click.echo(f"  [{status}] {label:<20} {path}")
        click.echo("")
        click.echo("Default thresholds:")
        click.echo(f"  FDR (adj p-value)  : {config.DE_PVALUE_THRESHOLD}  (override: --fdr)")
        click.echo(f"  Log2 fold-change   : {config.DE_LOG2FC_THRESHOLD}  (override: --log2fc)")
        click.echo(f"  Correlation cutoff : {config.COEXPR_CORRELATION_CUTOFF}")
        click.echo("  Novel prob cutoff  : 0.50  (override: --prob-threshold)")
        click.echo(f"  CV folds           : {config.CV_FOLDS}")
        click.echo(f"  CV primary metric  : {config.CV_METRIC_PRIMARY.upper()}")
    except Exception as exc:
        click.echo(f"  (Could not load config: {exc})")

    click.echo("")
    cache_dirs = sorted(Path(".").glob("*/.pipeline_cache"))
    if cache_dirs:
        click.echo("Cached pipeline runs (current directory):")
        for cache_dir in cache_dirs:
            click.echo(f"  {cache_dir.parent.resolve()}")
    else:
        click.echo("No cached pipeline runs found in the current directory.")


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


@main.group()
def cache():
    """Manage cached intermediate pipeline files."""


@cache.command("clear")
@click.option(
    "--output",
    "output_dir",
    default=None,
    type=click.Path(),
    metavar="DIR",
    help="Output directory whose cache to clear. Defaults to searching all subdirectories of the current directory.",
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


@main.command("run")
@click.option(
    "--tumor-expr",
    required=True,
    type=click.Path(exists=True),
    metavar="FILE",
    help="Tumor expression matrix: CSV with gene symbols as row index, sample IDs as columns.",
)
@click.option(
    "--normal-expr",
    required=True,
    type=click.Path(exists=True),
    metavar="FILE",
    help="Normal expression matrix: CSV with gene symbols as row index, sample IDs as columns.",
)
@click.option(
    "--tumor-meta",
    default=None,
    type=click.Path(exists=True),
    metavar="FILE",
    help="Tumor sample metadata: CSV with sample IDs as row index. If omitted, all samples in --tumor-expr are used.",
)
@click.option(
    "--normal-meta",
    default=None,
    type=click.Path(exists=True),
    metavar="FILE",
    help="Normal sample metadata: CSV with sample IDs as row index. If omitted, all samples in --normal-expr are used.",
)
@click.option(
    "--output",
    "output_dir",
    required=True,
    type=click.Path(),
    metavar="DIR",
    help="Output directory (created if it does not exist).",
)
@click.option(
    "--genes",
    default=None,
    type=click.Path(exists=True),
    metavar="FILE",
    help="Optional query gene list: one symbol per line (.txt) or TSV with GeneSymbol column. These genes are scored as candidates. If omitted, all genes are ranked.",
)
@click.option(
    "--labels",
    "labels_file",
    default=None,
    type=click.Path(exists=True),
    metavar="FILE",
    help="Known positive gene list for training labels: CSV/TSV with a GeneSymbol column. Falls back to the bundled LCGene database if omitted.",
)
@click.option(
    "--from-step",
    default="2",
    show_default=True,
    metavar="STEP",
    help="Resume from this step key (minimum 2 because step 1 is bypassed here).",
)
@click.option(
    "--to-step",
    default="19",
    show_default=True,
    metavar="STEP",
    help="Stop after this step key.",
)
@click.option(
    "--skip-optional",
    is_flag=True,
    help="Exclude all steps marked optional from the selected run.",
)
@click.option(
    "--skip-step",
    "skip_steps",
    multiple=True,
    metavar="STEP",
    help="Skip one specific step key. May be passed multiple times.",
)
@click.option(
    "--fdr",
    default=None,
    type=click.FloatRange(0.0, 1.0),
    metavar="FLOAT",
    help="Adjusted p-value (FDR) cutoff for differential expression (default: 0.05).",
)
@click.option(
    "--log2fc",
    "log2fc_threshold",
    default=None,
    type=click.FloatRange(0.0),
    metavar="FLOAT",
    help="Log2 fold-change threshold for differential expression (default: 1.0).",
)
@click.option(
    "--prob-threshold",
    default=None,
    type=click.FloatRange(0.0, 1.0),
    metavar="FLOAT",
    help="Minimum ML probability to flag a gene as a novel candidate (default: 0.5).",
)
@click.option(
    "--top-k",
    default=None,
    type=click.IntRange(1),
    metavar="N",
    help="Keep only the top N candidates in ranked_candidates output.",
)
@click.option(
    "--format",
    "fmt",
    default="csv",
    show_default=True,
    type=click.Choice(["csv", "excel", "json"], case_sensitive=False),
    help="Output file format for candidate tables.",
)
@click.option("--verbose", is_flag=True, help="Enable DEBUG-level logging.")
def run(
    tumor_expr,
    normal_expr,
    tumor_meta,
    normal_meta,
    output_dir,
    genes,
    labels_file,
    from_step,
    to_step,
    skip_optional,
    skip_steps,
    fdr,
    log2fc_threshold,
    prob_threshold,
    top_k,
    fmt,
    verbose,
):
    """Run the pipeline with your own expression data files."""
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    _setup_logging(verbose, output_path)
    log = logging.getLogger("gloom")
    log.info(f"gloom {__import__('gloom').__version__}  mode=run (custom data)")
    log.info(f"Output directory : {output_path}")

    _validate_expression_file(Path(tumor_expr), "tumor-expr")
    _validate_expression_file(Path(normal_expr), "normal-expr")

    _ensure_pipeline_on_path()
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

    log.info("Staging user-supplied data files ...")
    _stage_user_files(
        config=config,
        tumor_expr_path=Path(tumor_expr).resolve(),
        normal_expr_path=Path(normal_expr).resolve(),
        tumor_meta_path=Path(tumor_meta).resolve() if tumor_meta else None,
        normal_meta_path=Path(normal_meta).resolve() if normal_meta else None,
        labels_file=Path(labels_file).resolve() if labels_file else None,
        log=log,
    )

    steps = _select_steps(
        _load_pipeline_steps(),
        from_key=_normalize_step_key(from_step, minimum=2, maximum=19),
        to_key=_normalize_step_key(to_step, minimum=2, maximum=19),
        skip_optional=skip_optional,
        skip_keys=tuple(
            _normalize_step_key(step, minimum=2, maximum=19) for step in skip_steps
        ),
    )

    completed, failed = _run_steps(steps, log)
    _print_summary(completed, log)

    if failed:
        log.error(f"Pipeline stopped at {failed}.")
        sys.exit(1)

    log.info("Organizing outputs ...")
    from gloom.output import organize_outputs  # noqa: PLC0415

    organize_outputs(
        config,
        output_path,
        disease="luad",
        top_k=top_k,
        fmt=fmt,
        include_models=True,
    )

    ext = "xlsx" if fmt == "excel" else fmt
    click.echo("")
    click.echo(f"  Results ready : {output_path}")
    click.echo(f"  Main result   : {output_path / 'candidates' / f'ranked_candidates.{ext}'}")
    click.echo(f"  Report        : {output_path / 'report.html'}")


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
        force=True,
    )


def _patch_config(
    config,
    data_dir,
    output_path: Path,
    genes_file: Path | None = None,
    fdr=None,
    log2fc_threshold=None,
    prob_threshold=None,
    labels_file=None,
) -> None:
    """Redirect all config paths into a hidden cache inside output_path."""
    cache = output_path / ".pipeline_cache"

    if data_dir:
        config.DATA_ROOT = data_dir
        config.RAW_DIR = data_dir / "raw"
        config.TUMOR_EXPR_FILE = config.RAW_DIR / "cBioPortal (RNA Seq Data)" / "data_mrna_seq_v2_rsem.txt"
        config.TUMOR_META_FILE = config.RAW_DIR / "cBioPortal (RNA Seq Data)" / "data_clinical_patient.txt"
        config.NORMAL_EXPR_FILE = config.RAW_DIR / "Gtex (normal samples)" / "gene_tpm_v11_lung.gct.gz"
        config.NORMAL_META_FILE = (
            config.RAW_DIR / "Gtex (normal samples)" / "GTEx_Analysis_v11_Annotations_SampleAttributesDS - LUAD.txt"
        )

    if labels_file is not None:
        config.CANCER_GENE_FILE = labels_file

    if fdr is not None:
        config.DE_PVALUE_THRESHOLD = fdr
    if log2fc_threshold is not None:
        config.DE_LOG2FC_THRESHOLD = log2fc_threshold
    if prob_threshold is not None:
        config.NOVEL_PROB_THRESHOLD = prob_threshold

    if genes_file is not None:
        if genes_file.suffix.lower() == ".txt":
            import pandas as pd  # noqa: PLC0415

            symbols = [
                line.strip()
                for line in genes_file.read_text().splitlines()
                if line.strip() and not line.startswith("#")
            ]
            tmp_tsv = cache / "user_genes.tsv"
            tmp_tsv.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame({"GeneSymbol": symbols}).to_csv(tmp_tsv, sep="\t", index=False)
            config.QUERY_GENES_FILE = tmp_tsv
        else:
            config.QUERY_GENES_FILE = genes_file
    else:
        config.QUERY_GENES_FILE = None

    p = cache / "processed"
    r = cache / "results"
    f = cache / "figures"
    m = cache / "models"
    n = cache / "results" / "network"
    rp = cache / "results" / "reports"
    e = cache / "results" / "enrichment"
    lg = cache / "logs"

    config.PROCESSED_DIR = p
    config.OUTPUTS_ROOT = cache
    config.RESULTS_DIR = r
    config.FIGURES_DIR = f
    config.MODELS_DIR = m
    config.LOGS_DIR = lg
    config.NETWORK_DIR = n
    config.REPORTS_DIR = rp
    config.ENRICHMENT_DIR = e

    config.TUMOR_EXPR_PROCESSED = p / "tumor_expression_processed.csv"
    config.NORMAL_EXPR_PROCESSED = p / "normal_expression_processed.csv"
    config.TUMOR_EXPR_HARMONIZED = p / "tumor_expression_batch_corrected.csv"
    config.NORMAL_EXPR_HARMONIZED = p / "normal_expression_batch_corrected.csv"
    config.BATCH_CORRECTED_TUMOR_FILE = p / "tumor_expression_batch_corrected.csv"
    config.BATCH_CORRECTED_NORMAL_FILE = p / "normal_expression_batch_corrected.csv"
    config.DE_RESULTS_FILE = r / "differential_expression_results.csv"
    config.EXPR_FEATURES_FILE = p / "expression_features.csv"
    config.NETWORK_EDGES_FILE = n / "coexpression_network_edges.csv"
    config.NETWORK_GRAPH_FILE = n / "coexpression_network.graphml"
    config.NORMAL_NETWORK_EDGES_FILE = n / "normal_coexpression_network_edges.csv"
    config.NORMAL_NETWORK_GRAPH_FILE = n / "normal_coexpression_network.graphml"
    config.NETWORK_FEATURES_FILE = p / "network_features.csv"
    config.DIFFERENTIAL_NETWORK_FEATURES_FILE = p / "differential_network_features.csv"
    config.INTEGRATED_FEATURES_FILE = p / "integrated_features.csv"
    config.LABELS_FILE = p / "gene_labels.csv"
    config.TRAIN_FEATURES_FILE = p / "train_features.csv"
    config.VAL_FEATURES_FILE = p / "val_features.csv"
    config.TRAIN_LABELS_FILE = p / "train_labels.csv"
    config.VAL_LABELS_FILE = p / "val_labels.csv"
    config.GENE_RANKINGS_FILE = r / "gene_rankings.csv"
    config.FEATURE_IMPORTANCE_FILE = r / "feature_importance.csv"
    config.MODEL_METRICS_FILE = r / "model_metrics.csv"
    config.ANNOTATED_NETWORK_FILE = n / "annotated_network.graphml"
    config.ANNOTATED_EDGES_FILE = n / "annotated_edges.csv"
    config.ANNOTATED_NODES_FILE = n / "annotated_nodes.csv"
    config.INTERACTIVE_HTML_FILE = f / "interactive_network.html"
    config.FINAL_REPORT_FILE = rp / "pipeline_summary_report.csv"
    config.KEGG_ALL_FILE = e / "kegg_all_candidates.csv"
    config.KEGG_UP_FILE = e / "kegg_upregulated.csv"
    config.KEGG_DOWN_FILE = e / "kegg_downregulated.csv"
    config.KEGG_SUMMARY_FILE = e / "kegg_summary.csv"
    config.LOG_FILE = lg / "pipeline.log"

    config.create_output_dirs()


def _print_dry_run(
    config,
    steps,
    output_path: Path,
    *,
    skip_optional: bool = False,
    skip_step_keys: tuple[str, ...] = (),
) -> None:
    """Print validation results and execution plan without running anything."""
    click.echo("")
    click.echo("DRY RUN - nothing will be executed")
    click.echo("=" * 60)

    click.echo("\nInput files:")
    inputs = [
        ("Tumor expression", config.TUMOR_EXPR_FILE),
        ("Tumor metadata", config.TUMOR_META_FILE),
        ("Normal expression", config.NORMAL_EXPR_FILE),
        ("Normal metadata", config.NORMAL_META_FILE),
        ("Gene labels (ref)", config.CANCER_GENE_FILE),
        ("Query genes", getattr(config, "QUERY_GENES_FILE", None)),
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

    click.echo("\nEffective thresholds:")
    click.echo(f"  FDR (adj p-value)   : {config.DE_PVALUE_THRESHOLD}")
    click.echo(f"  Log2 fold-change    : {config.DE_LOG2FC_THRESHOLD}")
    click.echo(f"  Novel prob cutoff   : {getattr(config, 'NOVEL_PROB_THRESHOLD', 0.50)}")

    click.echo("\nSelected steps:")
    if skip_optional:
        click.echo("  Optional steps will be excluded.")
    if skip_step_keys:
        click.echo(f"  Explicitly skipped step keys: {', '.join(skip_step_keys)}")
    for _, step_label, _, _, _ in steps:
        click.echo(f"  {step_label}")

    click.echo(f"\nOutput directory : {output_path}")
    click.echo("")
    if all_ok:
        click.echo("All input files found. Ready to run.")
    else:
        click.echo("WARNING: Some input files are missing. Run 'gloom validate' for details.")
    click.echo("=" * 60)


def _run_steps(steps, log) -> tuple[list, str | None]:
    completed = []
    failed_step = None

    for step_key, step_label, module_name, func_name, optional in steps:
        log.info(f">>> {step_label}")
        t0 = time.time()
        try:
            module = importlib.import_module(module_name)
            importlib.reload(module)
            func = getattr(module, func_name)
            func()
            elapsed = time.time() - t0
            log.info(f"<<< {step_label} DONE ({elapsed:.1f}s)")
            completed.append((step_label, elapsed, "OK", optional))
        except Exception as exc:
            elapsed = time.time() - t0
            if optional:
                log.warning(f"<<< {step_label} SKIPPED ({type(exc).__name__}: {exc})")
                completed.append((step_label, elapsed, f"SKIPPED: {exc}", optional))
            else:
                log.error(f"<<< {step_label} FAILED: {type(exc).__name__}: {exc}")
                completed.append((step_label, elapsed, f"FAILED: {exc}", optional))
                failed_step = step_label
                break

    return completed, failed_step


def _print_summary(completed: list, log) -> None:
    log.info("=" * 60)
    log.info("EXECUTION SUMMARY")
    log.info("=" * 60)
    total = sum(elapsed for _, elapsed, _, _ in completed)
    for label, elapsed, status, _ in completed:
        if status == "OK":
            icon = "OK  "
        elif status.startswith("SKIPPED"):
            icon = "SKIP"
        else:
            icon = "FAIL"
        log.info(f"  [{icon}] {label:<50} {elapsed:>7.1f}s")
    log.info(f"  Total: {total:.1f}s  ({total / 60:.1f} min)")
    log.info("=" * 60)


def _validate_expression_file(path: Path, flag: str) -> None:
    """Quick sanity-check that a user-supplied expression CSV is readable and non-empty."""
    import pandas as pd  # noqa: PLC0415

    try:
        df = pd.read_csv(path, index_col=0, nrows=5)
    except Exception as exc:
        raise click.BadParameter(f"Cannot read file: {exc}", param_hint=f"--{flag}")
    if df.shape[1] == 0:
        raise click.BadParameter(
            "File has no sample columns - expected a genes x samples CSV (gene symbols as row index, sample IDs as column headers).",
            param_hint=f"--{flag}",
        )
    if df.shape[0] == 0:
        raise click.BadParameter("File has no gene rows.", param_hint=f"--{flag}")


def _stage_user_files(
    config,
    tumor_expr_path: Path,
    normal_expr_path: Path,
    tumor_meta_path: Path | None,
    normal_meta_path: Path | None,
    labels_file: Path | None,
    log,
) -> None:
    """
    Write user-supplied files into the pipeline cache as Step 1 outputs so that
    Step 2 (preprocessing) can read them without any changes to the pipeline code.
    """
    import pandas as pd  # noqa: PLC0415

    proc = Path(config.PROCESSED_DIR)

    log.info(f"  Reading tumor expression  : {tumor_expr_path.name}")
    tumor_expr = pd.read_csv(tumor_expr_path, index_col=0)
    tumor_expr.index.name = "gene"
    log.info(f"    Shape: {tumor_expr.shape[0]:,} genes x {tumor_expr.shape[1]:,} samples")
    tumor_expr.to_csv(proc / "tumor_expression_raw.csv")

    log.info(f"  Reading normal expression : {normal_expr_path.name}")
    normal_expr = pd.read_csv(normal_expr_path, index_col=0)
    normal_expr.index.name = "gene"
    log.info(f"    Shape: {normal_expr.shape[0]:,} genes x {normal_expr.shape[1]:,} samples")
    normal_expr.to_csv(proc / "normal_expression_raw.csv")

    if tumor_meta_path:
        log.info(f"  Reading tumor metadata    : {tumor_meta_path.name}")
        tumor_meta = pd.read_csv(tumor_meta_path, index_col=0)
        tumor_meta.index.name = "sample_id"
        log.info(f"    Shape: {tumor_meta.shape}")
    else:
        log.info("  No --tumor-meta provided - using all tumor expression samples.")
        tumor_meta = pd.DataFrame(index=pd.Index(tumor_expr.columns, name="sample_id"))
    tumor_meta.to_csv(proc / "tumor_metadata_raw.csv")

    if normal_meta_path:
        log.info(f"  Reading normal metadata   : {normal_meta_path.name}")
        normal_meta = pd.read_csv(normal_meta_path, index_col=0)
        normal_meta.index.name = "sample_id"
        log.info(f"    Shape: {normal_meta.shape}")
    else:
        log.info("  No --normal-meta provided - using all normal expression samples.")
        normal_meta = pd.DataFrame(index=pd.Index(normal_expr.columns, name="sample_id"))
    normal_meta.to_csv(proc / "normal_metadata_raw.csv")

    cancer_genes_dst = proc / "cancer_genes_raw.csv"

    if labels_file is not None:
        log.info(f"  Reading gene labels       : {labels_file.name}")
        sep = "\t" if labels_file.suffix.lower() in (".tsv", ".txt") else ","
        ldf = pd.read_csv(labels_file, sep=sep, low_memory=False)
        gene_col = next(
            (
                column
                for column in ldf.columns
                if column.strip().lower() in ("genesymbol", "gene_symbol", "gene", "symbol", "hgnc_symbol")
            ),
            ldf.columns[0],
        )
        genes = (
            ldf[gene_col]
            .dropna()
            .astype(str)
            .str.strip()
            .str.upper()
            .drop_duplicates()
            .reset_index(drop=True)
        )
        genes.name = "gene"
        genes.to_csv(cancer_genes_dst, index=False, header=True)
        log.info(f"    {len(genes):,} labeled genes loaded.")
    elif Path(config.CANCER_GENE_FILE).exists():
        log.info(f"  Using bundled LCGene database: {Path(config.CANCER_GENE_FILE).name}")
        ldf = pd.read_csv(config.CANCER_GENE_FILE, sep="\t", low_memory=False)
        gene_col = getattr(config, "LCGENE_GENE_COL", "GeneSymbol")
        genes = (
            ldf[gene_col]
            .dropna()
            .astype(str)
            .str.strip()
            .str.upper()
            .drop_duplicates()
            .reset_index(drop=True)
        )
        genes.name = "gene"
        genes.to_csv(cancer_genes_dst, index=False, header=True)
        ldf_save = ldf.copy()
        ldf_save[gene_col] = ldf_save[gene_col].astype(str).str.strip().str.upper()
        ldf_save.to_csv(proc / "lcgene_luad_full.csv", index=False)
        log.info(f"    {len(genes):,} LCGene genes staged.")
    else:
        log.warning(
            "  No --labels file and no bundled LCGene database found.\n"
            "  All training labels will be 0 (no known positives).\n"
            "  The model will still run but gene ranking will not be meaningful.\n"
            "  Provide --labels <file> with a GeneSymbol column to fix this."
        )
        pd.Series([], name="gene", dtype=str).to_csv(cancer_genes_dst, index=False, header=True)

    log.info("  Staging complete - starting pipeline at Step 2.")
