"""
gloom/output.py
---------------
Collects pipeline outputs from the hidden .pipeline_cache/ directory and
reorganizes them into a cleaner user-facing structure.
"""

import base64
import json
import logging
import shutil
from datetime import date
from pathlib import Path

log = logging.getLogger("gloom.output")


def organize_outputs(
    config,
    output_path: Path,
    disease: str = "luad",
    top_k: int | None = None,
    fmt: str = "csv",
    include_models: bool = True,
) -> None:
    """Map pipeline cache outputs to the public package output directory."""
    _make_dirs(output_path, include_models=include_models)

    _candidates(config, output_path, top_k=top_k, fmt=fmt)
    _tables(config, output_path)
    if include_models:
        _models(config, output_path, disease)
    else:
        _remove_if_exists(output_path / "models")
    _plots(config, output_path)
    _network(config, output_path)
    _report(config, output_path)
    _reports(config, output_path)

    log.info("Output structure ready.")
    _print_tree(output_path)


def _candidates(config, out: Path, top_k: int | None = None, fmt: str = "csv") -> None:
    import pandas as pd  # noqa: PLC0415

    d = out / "candidates"
    query_rankings = Path(config.RESULTS_DIR) / "query_gene_rankings.csv"
    primary_rankings = (
        query_rankings
        if query_rankings.exists() and query_rankings.stat().st_size > 0
        else Path(config.GENE_RANKINGS_FILE)
    )

    files = [
        (primary_rankings, "ranked_candidates"),
        (Path(config.RESULTS_DIR) / "novel_candidates.csv", "novel_candidates"),
        (Path(config.RESULTS_DIR) / "novel_candidates_sensitivity.csv", "novel_candidates_sensitivity"),
        (Path(config.GENE_RANKINGS_FILE), "full_gene_rankings"),
    ]

    for src_path, dst_stem in files:
        if not src_path.exists():
            log.warning(f"  MISSING (skipped): {src_path}")
            continue
        df = pd.read_csv(src_path, index_col=0)
        if top_k is not None and dst_stem in {
            "ranked_candidates",
            "novel_candidates",
            "novel_candidates_sensitivity",
        }:
            df = df.head(top_k)
        _save_table(df, d / dst_stem, fmt)
        log.info(f"  {src_path.name:50s} -> {dst_stem}.{_ext(fmt)}  ({len(df)} rows)")


def _tables(config, out: Path) -> None:
    d = out / "tables"
    table_map = [
        (config.DE_RESULTS_FILE, d / "de_results.csv"),
        (config.EXPR_FEATURES_FILE, d / "expression_features.csv"),
        (config.NETWORK_FEATURES_FILE, d / "network_features.csv"),
        (getattr(config, "DIFFERENTIAL_NETWORK_FEATURES_FILE", None), d / "differential_network_features.csv"),
        (config.INTEGRATED_FEATURES_FILE, d / "integrated_features.csv"),
        (config.LABELS_FILE, d / "gene_labels.csv"),
        (config.FEATURE_IMPORTANCE_FILE, d / "feature_importance.csv"),
        (config.MODEL_METRICS_FILE, d / "model_metrics.csv"),
        (Path(config.RESULTS_DIR) / "ranking_metrics.csv", d / "ranking_metrics.csv"),
        (Path(config.RESULTS_DIR) / "pu_bagging_scores.csv", d / "pu_bagging_scores.csv"),
        (Path(config.RESULTS_DIR) / "pu_bagging_metrics.csv", d / "pu_bagging_metrics.csv"),
        (getattr(config, "KEGG_SUMMARY_FILE", None), d / "kegg_summary.csv"),
    ]
    for src, dst in table_map:
        _copy(src, dst)


def _models(config, out: Path, disease: str) -> None:
    d = out / "models"
    static_files = [
        (config.MODELS_DIR / "best_model.joblib", d / "best_model.joblib"),
        (config.MODELS_DIR / "best_model_calibrated.joblib", d / "best_model_calibrated.joblib"),
        (config.MODELS_DIR / "robust_scaler.joblib", d / "robust_scaler.joblib"),
        (config.MODELS_DIR / "best_model_name.txt", d / "best_model_name.txt"),
        (config.MODELS_DIR / "cv_results.csv", d / "cv_results.csv"),
    ]
    for src, dst in static_files:
        _copy(src, dst)
    for model_file in sorted(Path(config.MODELS_DIR).glob("model_*.joblib")):
        _copy(model_file, d / model_file.name)
    _build_model_card(config, d / "model_card.json", disease)


def _plots(config, out: Path) -> None:
    d = out / "plots"
    html_map = [
        (config.FIGURES_DIR / "interactive_volcano.html", d / "volcano_plot.html"),
        (config.FIGURES_DIR / "interactive_feature_importance.html", d / "feature_importance.html"),
        (config.FIGURES_DIR / "interactive_network.html", d / "network_view.html"),
        (config.FIGURES_DIR / "interactive_ranking.html", d / "ranking_view.html"),
        (config.FIGURES_DIR / "interactive_dashboard.html", d / "dashboard.html"),
    ]
    for src, dst in html_map:
        _copy(src, dst)

    kegg_html = Path(config.FIGURES_DIR) / "interactive_kegg.html"
    if not kegg_html.exists():
        _build_kegg_html(config, kegg_html)
    _copy(kegg_html, d / "kegg_pathways.html")
    _build_roc_pr_html(config, d / "roc_pr_curves.html")

    png_map = [
        ("de_volcano_plot.png", "de_volcano_plot.png"),
        ("feature_importance_grouped.png", "feature_importance_grouped.png"),
        ("ranking_score_distribution.png", "ranking_score_distribution.png"),
        ("pu_bagging_score_distribution.png", "pu_bagging_score_distribution.png"),
        ("cv_auroc_auprc_comparison.png", "cv_auroc_auprc_comparison.png"),
        ("eval_roc_curves.png", "eval_roc_curves.png"),
        ("eval_pr_curves.png", "eval_pr_curves.png"),
        ("eval_confusion_matrices.png", "eval_confusion_matrices.png"),
        ("eval_confusion_multithreshold.png", "eval_confusion_multithreshold.png"),
        ("split_feature_pca.png", "split_feature_pca.png"),
        ("integrated_feature_pca.png", "integrated_feature_pca.png"),
        ("annotated_degree_vs_rank.png", "annotated_degree_vs_rank.png"),
    ]
    for src_name, dst_name in png_map:
        _copy(Path(config.FIGURES_DIR) / src_name, d / dst_name)


def _network(config, out: Path) -> None:
    d = out / "network"
    annotated_graph = Path(config.ANNOTATED_NETWORK_FILE)
    _copy(annotated_graph, d / "annotated_coexpression.graphml")
    _copy(annotated_graph, d / "annotated_coexpression.cytoscape.xml")
    _copy(config.ANNOTATED_NODES_FILE, d / "annotated_nodes.csv")
    _copy(config.ANNOTATED_EDGES_FILE, d / "annotated_edges.csv")
    _copy(config.NETWORK_GRAPH_FILE, d / "tumor_coexpression.graphml")
    _copy(config.NETWORK_EDGES_FILE, d / "tumor_coexpression_edges.csv")
    _copy(getattr(config, "NORMAL_NETWORK_GRAPH_FILE", None), d / "normal_coexpression.graphml")
    _copy(getattr(config, "NORMAL_NETWORK_EDGES_FILE", None), d / "normal_coexpression_edges.csv")
    _copy_tree(Path(config.NETWORK_DIR) / "exports", d / "exports")


def _report(config, out: Path) -> None:
    _copy(config.FIGURES_DIR / "interactive_dashboard.html", out / "report.html")
    report = out / "report.html"
    if report.exists():
        _inject_kegg_tab(report, config)


def _reports(config, out: Path) -> None:
    d = out / "reports"
    report_map = [
        (Path(config.REPORTS_DIR) / "pipeline_report.txt", d / "pipeline_report.txt"),
        (Path(config.REPORTS_DIR) / "pipeline_summary_table.csv", d / "pipeline_summary_table.csv"),
        (Path(config.REPORTS_DIR) / "pipeline_summary_figure.png", d / "pipeline_summary_figure.png"),
    ]
    for src, dst in report_map:
        _copy(src, dst)


def _build_model_card(config, dst: Path, disease: str) -> None:
    """Generate a lightweight model card for the best pipeline model."""
    import pandas as pd  # noqa: PLC0415

    card: dict = {
        "name": "gloom-luad-prioritizer",
        "version": __import__("gloom").__version__,
        "disease": disease,
        "generated": str(date.today()),
    }

    name_file = Path(config.MODELS_DIR) / "best_model_name.txt"
    if name_file.exists():
        card["best_model"] = name_file.read_text().strip()

    metrics_file = Path(config.MODEL_METRICS_FILE)
    if metrics_file.exists():
        df = pd.read_csv(metrics_file)
        if not df.empty:
            row = df.iloc[0]
            scalar_cols = [
                "auroc",
                "auprc",
                "accuracy",
                "f1_optimal",
                "mcc",
                "brier_score",
                "threshold_optimal",
                "recall_at_50",
                "recall_at_100",
                "recall_at_200",
            ]
            card["metrics"] = {
                col: round(float(row[col]), 4)
                for col in scalar_cols
                if col in df.columns and not pd.isna(row[col])
            }

    card["training_params"] = {
        "cv_folds": config.CV_FOLDS,
        "cv_metric_primary": config.CV_METRIC_PRIMARY,
        "threshold_strategy": config.THRESHOLD_STRATEGY,
        "use_smote": config.USE_SMOTE,
        "use_undersampling": config.USE_UNDERSAMPLING,
        "random_state": config.RANDOM_STATE,
    }

    dst.write_text(json.dumps(card, indent=2))
    log.info(f"  model_card.json -> {dst}")


def _build_roc_pr_html(config, dst: Path) -> None:
    """Build a simple HTML page embedding ROC and PR PNG outputs."""
    roc_png = Path(config.FIGURES_DIR) / "eval_roc_curves.png"
    pr_png = Path(config.FIGURES_DIR) / "eval_pr_curves.png"

    sections = []
    for title, png_path in [("ROC Curves", roc_png), ("Precision-Recall Curves", pr_png)]:
        if not png_path.exists():
            log.warning(f"  {png_path.name} not found - skipping from roc_pr_curves.html")
            continue
        b64 = base64.b64encode(png_path.read_bytes()).decode()
        sections.append(
            f"<h2>{title}</h2>"
            f'<img src="data:image/png;base64,{b64}" '
            f'style="max-width:100%;border:1px solid #e0e0e0;border-radius:6px;">'
        )

    if not sections:
        log.warning("  No ROC/PR PNGs found - roc_pr_curves.html not created.")
        return

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>gloom - ROC &amp; PR Curves</title>
  <style>
    body  {{ font-family: system-ui, sans-serif; max-width: 1100px;
             margin: 2rem auto; padding: 0 1.5rem; color: #222; }}
    h1    {{ font-size: 1.6rem; border-bottom: 2px solid #4a90d9; padding-bottom: .4rem; }}
    h2    {{ font-size: 1.2rem; color: #4a90d9; margin-top: 2.5rem; }}
    img   {{ display: block; margin-bottom: 2rem; }}
  </style>
</head>
<body>
  <h1>Model Evaluation - ROC &amp; PR Curves</h1>
  {"".join(sections)}
</body>
</html>"""
    dst.write_text(html, encoding="utf-8")
    log.info(f"  roc_pr_curves.html -> {dst}")


def _build_kegg_figure(config):
    """Build a Plotly KEGG enrichment figure from step 19 CSV outputs."""
    try:
        import numpy as np  # noqa: PLC0415
        import pandas as pd  # noqa: PLC0415
        import plotly.graph_objects as go  # noqa: PLC0415
    except ImportError:
        log.warning("  Plotly not installed - KEGG figure not built.")
        return None

    padj_cutoff = 0.05
    top_n = 15
    subset_colors = {
        "All candidates": "#3A7DBF",
        "Upregulated": "#E8524A",
        "Downregulated": "#2CA07C",
    }

    kegg_files = [
        ("All candidates", getattr(config, "KEGG_ALL_FILE", None)),
        ("Upregulated", getattr(config, "KEGG_UP_FILE", None)),
        ("Downregulated", getattr(config, "KEGG_DOWN_FILE", None)),
    ]

    subsets = {}
    for label, fpath in kegg_files:
        if fpath and Path(fpath).exists():
            df = pd.read_csv(fpath)
            if not df.empty:
                if "neg_log10_padj" not in df.columns and "padj" in df.columns:
                    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
                subsets[label] = df

    if not subsets:
        log.info("  No KEGG enrichment data found - skipping KEGG figure.")
        return None

    fig = go.Figure()
    subset_names = list(subsets.keys())
    has_traces = False

    for i, (label, df) in enumerate(subsets.items()):
        sig = df[df["padj"] < padj_cutoff].head(top_n).sort_values("neg_log10_padj")
        if sig.empty:
            continue
        has_traces = True
        color = subset_colors.get(label, "#888888")

        def _hover(row, subset_label=label):
            gene_str = str(row.get("genes", ""))
            gene_preview = gene_str[:80] + "..." if len(gene_str) > 80 else gene_str
            return (
                f"<b>{row['pathway']}</b><br>Subset: {subset_label}<br>"
                f"adj P: {row['padj']:.2e}<br>"
                f"-log10(padj): {row['neg_log10_padj']:.2f}<br>"
                f"Overlap: {row.get('overlap_genes', '?')}/{row.get('pathway_size', '?')}<br>"
                f"Odds Ratio: {float(row.get('odds_ratio', 0)):.2f}<br>"
                f"Genes: {gene_preview}"
            )

        fig.add_trace(
            go.Bar(
                x=sig["neg_log10_padj"],
                y=sig["pathway"],
                orientation="h",
                name=label,
                visible=(i == 0),
                marker=dict(
                    color=sig["neg_log10_padj"],
                    colorscale=[[0, "lightyellow"], [0.5, color], [1.0, color]],
                    cmin=0,
                    cmax=sig["neg_log10_padj"].max() or 1,
                    showscale=(i == 0),
                    colorbar=dict(title="-log10(padj)", thickness=14, len=0.7),
                ),
                text=sig["neg_log10_padj"].round(2),
                textposition="outside",
                customdata=sig.apply(_hover, axis=1).tolist(),
                hovertemplate="%{customdata}<extra></extra>",
            )
        )

    if not has_traces:
        log.info("  No significant KEGG pathways found - skipping KEGG figure.")
        return None

    buttons = []
    for i, name in enumerate(subset_names):
        buttons.append(
            dict(
                label=name,
                method="update",
                args=[
                    {"visible": [j == i for j in range(len(subset_names))]},
                    {"title.text": f"<b>KEGG Pathway Enrichment - {name}</b>"},
                ],
            )
        )
    buttons.append(
        dict(
            label="All subsets",
            method="update",
            args=[
                {"visible": [True] * len(subset_names)},
                {"title.text": "<b>KEGG Pathway Enrichment - All Subsets</b>"},
            ],
        )
    )

    fig.add_vline(
        x=-np.log10(padj_cutoff),
        line_dash="dash",
        line_color="grey",
        line_width=1.5,
        annotation_text="padj=0.05",
        annotation_position="top right",
        annotation_font_size=11,
    )

    fig.update_layout(
        title=f"<b>KEGG Pathway Enrichment - {subset_names[0]}</b>",
        xaxis_title="-log10(adjusted p-value)",
        yaxis_title="KEGG Pathway",
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=max(520, top_n * 36 + 160),
        width=1000,
        margin=dict(l=300, r=80, t=110, b=60),
        updatemenus=[
            dict(
                type="buttons",
                buttons=buttons,
                direction="right",
                showactive=True,
                x=0.0,
                xanchor="left",
                y=1.14,
                yanchor="top",
                bgcolor="#EEF3FA",
                bordercolor="#2D6A9F",
                font=dict(size=12),
            )
        ],
        legend=dict(orientation="h", y=-0.18),
    )
    return fig


def _build_kegg_html(config, out_path: Path) -> None:
    """Write an interactive_kegg.html file from KEGG CSV outputs."""
    fig = _build_kegg_figure(config)
    if fig is None:
        return
    fig.write_html(str(out_path), include_plotlyjs="cdn")
    log.info(f"  interactive_kegg.html built from KEGG CSVs -> {out_path.name}")


def _inject_kegg_tab(report: Path, config) -> None:
    """Embed the KEGG chart inline into report.html if the dashboard is present."""
    import plotly.io as pio  # noqa: PLC0415

    html = report.read_text(encoding="utf-8")
    if "KEGG Pathways" in html:
        return

    fig = _build_kegg_figure(config)
    if fig is None:
        log.info("  KEGG tab not added to report.html - no enrichment data.")
        return

    kegg_div = pio.to_html(
        fig,
        full_html=False,
        include_plotlyjs=False,
        div_id="kegg_div",
        config={"responsive": True},
    )

    tab_button = "  <div class=\"tab\" onclick=\"showTab('kegg',this)\">KEGG Pathways</div>\n"
    old_bar_end = (
        "  <div class=\"tab\"         onclick=\"showTab('novel',this)\">Novel Candidates</div>\n"
        "</div>"
    )
    new_bar_end = (
        "  <div class=\"tab\"         onclick=\"showTab('novel',this)\">Novel Candidates</div>\n"
        + tab_button
        + "</div>"
    )
    if old_bar_end not in html:
        log.warning("  Could not inject KEGG tab - tab-bar marker not found in report.html")
        return

    html = html.replace(old_bar_end, new_bar_end, 1)
    kegg_tab_content = (
        '\n<div id="tab-kegg" class="tab-content">\n'
        '  <div class="plot-card">\n'
        '    <div class="plot-desc"><b>KEGG Pathway Enrichment</b> - Novel LUAD candidates tested against KEGG_2021_Human via Enrichr, filtered to lung/cancer-relevant pathways. Use the buttons above the chart to switch between subsets. Hover a bar for pathway details and overlapping gene list.</div>\n'
        f"    {kegg_div}\n"
        "  </div>\n"
        "</div>"
    )
    html = html.replace("<footer>LUAD ML Pipeline", kegg_tab_content + "\n<footer>LUAD ML Pipeline", 1)
    report.write_text(html, encoding="utf-8")
    log.info("  KEGG tab embedded inline into report.html")


def _make_dirs(out: Path, *, include_models: bool = True) -> None:
    subdirs = ["candidates", "tables", "plots", "network", "reports"]
    if include_models:
        subdirs.append("models")
    for sub in subdirs:
        (out / sub).mkdir(parents=True, exist_ok=True)


def _ext(fmt: str) -> str:
    return {"csv": "csv", "excel": "xlsx", "json": "json"}[fmt]


def _save_table(df, dst_stem: Path, fmt: str) -> None:
    ext = _ext(fmt)
    dst = dst_stem.with_suffix(f".{ext}")
    if fmt == "csv":
        df.to_csv(dst)
    elif fmt == "excel":
        df.to_excel(dst, engine="openpyxl")
    elif fmt == "json":
        df.to_json(dst, orient="records", indent=2)


def _copy(src, dst: Path) -> None:
    if src is None:
        return
    src = Path(src)
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        log.info(f"  {src.name:50s} -> {dst.relative_to(dst.parent.parent)}")
    else:
        log.warning(f"  MISSING (skipped): {src}")


def _copy_tree(src_dir: Path, dst_dir: Path) -> None:
    if not src_dir.exists():
        log.warning(f"  MISSING (skipped): {src_dir}")
        return
    for src_file in src_dir.rglob("*"):
        if src_file.is_file():
            rel = src_file.relative_to(src_dir)
            dst_file = dst_dir / rel
            dst_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src_file, dst_file)
            log.info(f"  {src_file.name:50s} -> {dst_file.relative_to(dst_dir.parent)}")


def _remove_if_exists(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
        log.info(f"  removed stale directory -> {path.name}/")


def _print_tree(out: Path) -> None:
    log.info("")
    log.info(f"{out.name}/")
    for file_path in sorted(out.rglob("*")):
        if file_path.is_file() and ".pipeline_cache" not in str(file_path):
            rel = file_path.relative_to(out)
            depth = len(rel.parts) - 1
            indent = "    " * depth
            log.info(f"  {indent}|-- {file_path.name}")
    log.info("")
