# How to use GLOOM

## A. Run the full pipeline

Use this when you want everything from raw data to final results.

```bash
cd script
python run_pipeline.py
```

### Best for

- first full analysis
- final production run
- reproducible end-to-end results

---

## B. Run step by step

Use this when you want control or debugging.

```bash
python config.py
python step1_data_loading.py
python step2_preprocessing.py
python step3_harmonization.py
python step4_differential_expression.py
python step5_expression_features.py
python step6_coexpression_network.py
python step7_network_features.py
python step8_feature_integration.py
python step9_label_construction.py
python step10_train_val_split.py
python step11_model_training.py
python step12_model_evaluation.py
python step13_feature_importance.py
python step14_gene_ranking.py
python step15_network_annotation.py
python step16_network_export.py
python step17_interactive_visualization.py
python step18_final_report.py
```

### Best for

- debugging
- changing one parameter and rerunning only later steps
- teaching/demo use

---

## C. Use GLOOM only for gene ranking

If your real question is: **“Which genes should I study first?”**

Use:

- `step14_gene_ranking.py`
- `results/gene_rankings.csv`
- `results/novel_candidates.csv`

### Best for

- candidate gene discovery
- shortlist generation
- downstream manual review

---

## D. Use GLOOM only for model comparison

If your real question is: **“Which machine-learning model worked best?”**

Use:

- `step12_model_evaluation.py`
- `results/model_metrics.csv`
- ROC / PR figures

### Best for

- ML benchmarking
- method comparison
- choosing the model for deployment

---

## E. Use GLOOM only for interpretation

If your question is: **“Why did the model choose these genes?”**

Use:

- `step13_feature_importance.py`
- `results/feature_importance.csv`
- grouped importance figure

### Best for

- explainability
- reports and manuscripts
- understanding model behavior

---

## F. Use GLOOM only for network exploration

If your question is: **“How are these genes connected?”**

Use:

- `step15_network_annotation.py`
- `step16_network_export.py`
- `results/network/annotated_network.graphml`
- `results/network/exports/`

### Best for

- Cytoscape/Gephi work
- hub gene analysis
- subnetwork discovery

---

## G. Use GLOOM only for interactive browsing

If you want an easy visual interface without touching Python tables.

Open:

- `figures/interactive_dashboard.html`

### Best for

- collaborators
- presentations
- fast exploration

---

## H. Use GLOOM only for final reporting

If you want a compact end product.

Use:

- `results/reports/pipeline_report.txt`
- `results/reports/pipeline_summary_table.csv`
- `results/reports/pipeline_summary_figure.png`

### Best for

- posters
- slides
- reporting to supervisors or collaborators
