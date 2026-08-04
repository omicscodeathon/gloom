## `output/`

This directory is the human-readable summary layer for the project outputs.

It currently contains:

- `report.txt`
  A plain-text overview of the latest generated results, figures, ranking outputs, network statistics, and enrichment findings.
- `readme.md`
  This file, which explains the purpose of the directory and where the full pipeline artifacts live.

### Purpose

Use `output/` when you want a quick summary of what the code produced without opening every CSV, GraphML, PNG, or HTML file individually.

### Relationship To Other Folders

- `results/`
  Stores the structured outputs of the pipeline, including model metrics, gene rankings, candidate tables, enrichment results, classification reports, and network exports.
- `figures/`
  Stores the generated static figures and interactive HTML visualizations.
- `results/reports/`
  Stores the larger report assets and presentation-ready summary figures.

### Suggested Reading Order

1. Open `report.txt` for the summary.
2. Open `../results/model_metrics.csv` for model performance details.
3. Open `../results/gene_rankings.csv` and `../results/novel_candidates.csv` for the main prioritization outputs.
4. Open `../figures/interactive_dashboard.html` or `../results/reports/pipeline_summary_figure.png` for visual review.

### Notes

- The files in this folder are documentation-style outputs.
- The source data and full generated artifacts remain in `results/` and `figures/`.
- This folder is a good place to keep future exported summaries or user-facing deliverables.
