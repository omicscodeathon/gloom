# What GLOOM is

**GLOOM** is a Python package for finding important LUAD genes from transcriptomic data.

In simple words, GLOOM does this:

1. Reads raw tumor and normal gene-expression files.
2. Cleans the data.
3. Keeps genes that can be fairly compared.
4. Builds many clues (features) about each gene.
5. Trains machine-learning models on known cancer genes.
6. Scores every gene.
7. Builds and annotates a gene network.
8. Exports tables, graphs, reports, and interactive dashboards.

Short version:

**raw data → cleaned data → features → machine learning → gene ranking → network annotation → interactive visuals → final report**

---

## Real run summary (reference execution)

The reference pipeline run completed all 18 steps successfully in about **480.2 seconds (~8.0 minutes)**.

### Main numbers from the successful run

- Raw tumor expression after duplicate cleanup: **20,512 genes × 510 samples**
- Raw normal expression after duplicate cleanup: **73,321 genes × 604 samples**
- After preprocessing:
  - tumor: **16,570 genes × 510 samples**
  - normal: **21,159 genes × 604 samples**
- Shared genes after harmonization: **10,986**
- Positive labeled genes retained from Cancer Gene Census: **488**
- Integrated features: **10,986 genes × 42 features**
- Final network size: **10,986 nodes, 36,943 edges**
- Novel candidates found by ranking: **3,754**
- Best model selected during training: **logistic_regression**

Why these numbers matter:

- they prove the pipeline ran correctly,
- they tell you the scale of the data,
- and they show what the package really produced.

---

# Folder-by-folder explanation

## 1) `data/`

This folder stores the data used by the package.

### Purpose

This is where GLOOM gets its input and where it saves cleaned data that later steps need.

### Subfolders

- `data/raw/` = original files from outside sources
- `data/processed/` = cleaned and transformed files created by GLOOM

### Very simple idea

- **raw** = ingredients from the market
- **processed** = ingredients already washed and cut

---

## 2) `data/raw/`

This folder stores the original input files.

### Purpose

These are the files the package starts from. They should not be edited by hand.

### Expected inputs

- cBioPortal LUAD tumor RNA-seq matrix
- cBioPortal tumor metadata
- GTEx lung normal expression matrix
- GTEx normal metadata
- Cancer Gene Census file for known cancer genes

### Why this matters

If these files are wrong or missing, the whole pipeline cannot start.

### When to use this folder

Use this folder when:

- you want to run the full pipeline from scratch,
- you want to replace the data with a new cohort,
- you want to check data provenance.

---

## 3) `data/processed/`

This folder stores cleaned and prepared files.

### Purpose

This folder makes the rest of the pipeline easier and safer. Instead of using messy raw files again and again, later steps use standardized processed files.

### Main file types here

- cleaned tumor expression
- cleaned normal expression
- harmonized tumor and normal matrices
- expression feature tables
- network feature tables
- integrated feature matrix
- labels
- train/validation split files

### Why this matters

This folder is the bridge between raw data and machine learning.

### When to use this folder

Use this folder when:

- you want to inspect how the data changed after cleaning,
- you want to rerun only later steps,
- you want to debug feature creation or splitting.

---

## 4) `results/`

This folder stores the main answers from GLOOM.

### Purpose

This is where the package gives you the final useful outputs.

### Key files

- `gene_rankings.csv`
- `novel_candidates.csv`
- `feature_importance.csv`
- `model_metrics.csv`

### Why this matters

If someone asks, “What did the pipeline discover?”, the answer is usually in this folder.

### When to use this folder

Use this folder when:

- you want the final ranked genes,
- you want model performance,
- you want interpretation of important features.

---

## 5) `results/network/`

This folder stores network outputs.

### Purpose

GLOOM does not only rank genes. It also shows how genes are connected.

### Key files

- `annotated_network.graphml`
- `annotated_nodes.csv`
- `annotated_edges.csv`

### Why this matters

Some genes may be important not only because of their expression, but also because of where they sit inside a gene network.

### When to use this folder

Use this folder when:

- you want to view the network in Cytoscape or Gephi,
- you want to study hub genes,
- you want to inspect which top genes connect together.

---

## 6) `results/network/exports/`

This folder stores the network in many file formats.

### Purpose

Different tools like different formats. This folder lets you reuse the network anywhere.

### Key files

- `network_full.graphml`
- `network_full.gml`
- `network_full_edgelist.tsv`
- `network_nodes.tsv`
- `network_cytoscape.json`
- small subnetworks such as top100, cgc, novel, candidates
- network statistics report files

### Why this matters

It makes sharing and reuse much easier.

### When to use this folder

Use this folder when:

- you need a format for another program,
- you want a smaller subnetwork,
- you want a statistics summary of the graph.

---

## 7) `results/reports/`

This folder stores final summaries.

### Purpose

This folder gives you the easy-to-read final outputs.

### Key files

- `pipeline_summary_table.csv`
- `pipeline_report.txt`
- `pipeline_summary_figure.png`
- `pipeline_summary_figure.pdf`

### Why this matters

Not everybody wants raw tables. Some people want one summary they can read, share, or present.

### When to use this folder

Use this folder when:

- you need a quick report,
- you are preparing a poster or slide deck,
- you want a compact summary of the whole analysis.

---

## 8) `figures/`

This folder stores visual outputs.

### Purpose

It helps humans understand results faster.

### Key files

- static plots
- interactive HTML pages
- dashboard files

### Examples

- `interactive_volcano.html`
- `interactive_ranking.html`
- `interactive_network.html`
- `interactive_feature_importance.html`
- `interactive_dashboard.html`

### When to use this folder

Use this folder when:

- you want to inspect results visually,
- you want to present findings,
- you want to browse results without looking at raw tables.

---

## 9) `models/`

This folder stores trained machine-learning models.

### Purpose

Training can take time. This folder keeps the trained models so you can reuse them.

### Typical contents

- `.joblib` model files
- `best_model.joblib`
- `best_model_name.txt`
- cross-validation results

### When to use this folder

Use this folder when:

- you want to score new data,
- you want to reload a trained model,
- you want to compare models or deploy one later.

---

## 10) `logs/`

This folder stores execution logs.

### Purpose

This folder is the diary of the pipeline.

### Main file

- `pipeline.log`

### When to use this folder

Use this folder when:

- something failed,
- you want runtimes,
- you want proof that the full run completed,
- you want reference values for documentation.

---

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

---

# Important result files explained simply

## `model_metrics.csv`

This is the **report card of the models**.

Use it to compare models.

Simple idea:

- higher AUROC / AUPRC / MCC / F1 = better
- lower Brier score = better
- `optimal_threshold` = decision cut-off used later in ranking

---

## `feature_importance.csv`

This tells you **which clues the AI trusted most**.

Use it to explain the model.

Simple idea:

- `rank = 1` means strongest clue
- `mean_importance` = average usefulness of that clue
- `feature_group` tells whether the clue came from tumor stats, normal stats, differential expression, network, etc.

---

## `gene_rankings.csv`

This is the **big leaderboard of genes**.

Use it to pick genes to study first.

Simple idea:

- each row = one gene
- higher `predicted_prob` = more suspicious / more interesting
- smaller `rank` = closer to the top

---

## `novel_candidates.csv`

This is the **shortlist of new interesting genes**.

Use it for discovery.

Simple idea:

- genes here are strongly scored by the model
- but they are not already known cancer genes in the label list

---

# What do I do with the `.graphml` file?

## Simple answer

A `.graphml` file is a **network file**.
It stores genes as **nodes** and links between genes as **edges**, plus extra attributes.

## What you can do with it

### 1. Open it in Cytoscape

Use this when you want a desktop biological network viewer.

Typical workflow:

1. Open Cytoscape.
2. Go to **File → Import → Network → File**.
3. Select the `.graphml` file.
4. Map node color to `node_color`.
5. Map node size to `node_size` or `predicted_prob`.
6. Use `rank`, `novel_candidate`, or `is_cgc_gene` to filter/highlight nodes.

### 2. Open it in Gephi

Use this when you want fast graph layout and visual exploration.

Typical workflow:

1. Open Gephi.
2. Import the `.graphml` file.
3. Run a layout such as ForceAtlas2.
4. Color nodes by category.
5. Resize nodes by degree or predicted probability.

### 3. Load it in Python with NetworkX

Use this when you want programmatic work.

```python
import networkx as nx
G = nx.read_graphml("annotated_network.graphml")
print(G.number_of_nodes(), G.number_of_edges())
```

### 4. Use small subnetworks instead of the full graph

If the full graph is too big, use files such as:

- `subnetwork_top100.graphml`
- `subnetwork_cgc.graphml`
- `subnetwork_novel.graphml`
- `subnetwork_candidates.graphml`

## Who uses `.graphml` most?

- network biologists
- Cytoscape users
- graph analysts
- anyone building visual biology stories

---

# What do I do with the `.joblib` model file?

## Simple answer

A `.joblib` file stores a **trained machine-learning model**.
It saves time because you do not need to train again.

## What you can do with it

### 1. Load the model again in Python

```python
import joblib
model = joblib.load("best_model.joblib")
```

### 2. Use it to score new genes or new feature tables

If you build features for another compatible dataset, you can use the saved model to predict on that new table.

### 3. Reuse it in later pipeline steps

The package itself reloads the model for evaluation, feature importance, and gene ranking.

### 4. Keep it for reproducibility

It is a frozen copy of the trained model used in the reference analysis.

## Important warning

Only load `.joblib` files from a trusted source.

## Who uses `.joblib` most?

- ML developers
- bioinformaticians doing repeated scoring
- people deploying or reusing the best model

---

# Comparative study: GLOOM vs alternatives

## Why compare?

No single tool does everything. GLOOM is strongest when you need **one reproducible workflow** that goes from data to ranking to network to report.

---

## 1. Cytoscape

### What it is good at

- strong interactive network visualization
- large app ecosystem
- very popular in biology
- good for manual visual exploration

### Limits compared with GLOOM

- not mainly an end-to-end machine-learning ranking workflow
- you still need other scripts/tools for preprocessing, model training, and ranking

### Best use case

Use Cytoscape **with GLOOM** after GLOOM exports the `.graphml` file.

---

## 2. Gephi

### What it is good at

- fast graph visualization
- good layouts for large networks
- good for visual storytelling and publications

### Limits compared with GLOOM

- strong on graph display, weak on genomics-specific preprocessing and ML ranking
- not a complete cancer-gene prioritization pipeline

### Best use case

Use Gephi when your main goal is **graph exploration and graph design**.

---

## 3. PINTA

### What it is good at

- network-based gene prioritization from expression data
- good historical example of web-based prioritization

### Limits compared with GLOOM

- web-server style workflow
- less flexible for custom scripting and reproducibility
- less integrated with local custom visualization and model interpretation

### Best use case

Use PINTA when you want a **quick gene-prioritization web workflow** and do not need a local custom Python project.

---

## 4. WGCNA

### What it is good at

- strong co-expression network analysis
- module detection
- module–trait relationships
- widely used in transcriptomics

### Limits compared with GLOOM

- mainly an R network-analysis framework, not a full ML ranking + reporting package
- usually needs extra scripting for classification, ranking, and final dashboards

### Best use case

Use WGCNA when your main goal is **module-based co-expression biology**, not full end-to-end ML-driven ranking.

---

## 5. NetworkX alone

### What it is good at

- flexible Python graph operations
- reading/writing many graph formats
- custom graph algorithms

### Limits compared with GLOOM

- no ready-made omics workflow by itself
- you must build preprocessing, differential expression, feature creation, ranking, and reporting by hand

### Best use case

Use NetworkX when you want **custom Python graph programming**.

---

## Best position of GLOOM

GLOOM is strongest when you need:

- one pipeline instead of many disconnected tools,
- reproducibility,
- gene ranking + network annotation together,
- final tables + final visuals + reports in one place.

---

# How to do this analysis without GLOOM

You *can* do it manually, but it takes more work.

## Manual path without GLOOM

### Step A: load raw tumor and normal expression data

Tools:

- Python (pandas)
- or R

### Step B: clean the data

You would need to:

- remove duplicates
- replace invalid values
- log-transform
- filter low-expression genes
- filter low-variance genes

### Step C: harmonize gene IDs

You must make sure tumor and normal use the same gene names.

### Step D: run differential expression

You can do this with:

- Python statistics
- statsmodels
- or R packages like limma / DESeq2 / edgeR (if you change workflow)

### Step E: build expression features

You would manually calculate means, medians, variability, ranks, contrasts, etc.

### Step F: build co-expression network

You would manually calculate correlations and keep strong edges.

### Step G: extract network features

You would use NetworkX or igraph to compute:

- degree
- betweenness
- clustering
- eigenvector centrality
- component size

### Step H: create labels

You would manually match known cancer genes from CGC to your gene list.

### Step I: split train and validation data

You would use scikit-learn by hand.

### Step J: train many ML models

You would manually define models, parameters, cross-validation, and model persistence.

### Step K: evaluate models

You would manually compute AUROC, AUPRC, MCC, F1, threshold selection, plots.

### Step L: rank all genes

You would load the best model, score the full feature matrix, and build the candidate tables yourself.

### Step M: annotate and export the network

You would manually attach node attributes and save GraphML/JSON/TSV.

### Step N: build reports and figures

You would manually write reporting code and interactive HTML pages.

## Simple truth

Without GLOOM, the work is possible, but it is more fragmented, more error-prone, and harder to reproduce.

---

# Who benefits most from GLOOM?

## 1. Bioinformaticians

Because it reduces many steps into one reproducible workflow.

## 2. Computational biologists

Because it links transcriptomics, networks, and machine learning in one place.

## 3. Cancer researchers

Because it helps find candidate genes and see how they connect in the network.

## 4. Students and trainees

Because it gives a full example of how an omics ML pipeline is structured.

## 5. Wet-lab collaborators

Because they can use final outputs like `gene_rankings.csv`, `novel_candidates.csv`, and the dashboard without understanding every technical step.

## 6. Presentation/report users

Because the package produces final figures, dashboards, and summary reports automatically.

---

# Which tools can help *with* GLOOM?

## For raw analysis and reruns

- Python
- pandas
- scikit-learn
- matplotlib
- plotly
- NetworkX

## For network visualization after export

- Cytoscape
- Gephi

## For spreadsheet review

- Excel
- LibreOffice Calc
- Google Sheets (for smaller files)

## For notebooks

- Jupyter
- VS Code notebooks

## For report/presentation work

- PowerPoint
- Word
- PDF viewers
- browser for HTML dashboards

---

# Recommended workflow by user type

## If you are a researcher and want biological candidates

1. Run the full pipeline.
2. Open `gene_rankings.csv`.
3. Inspect `novel_candidates.csv`.
4. Open `interactive_dashboard.html`.
5. Use Cytoscape on `annotated_network.graphml`.

## If you are an ML person

1. Check `model_metrics.csv`.
2. Check `feature_importance.csv`.
3. Load the `.joblib` best model.
4. Reuse the model on compatible new feature matrices.

## If you are a graph/network person

1. Use `annotated_network.graphml`.
2. Use `network_full.graphml` or `subnetwork_*.graphml`.
3. Open in Cytoscape/Gephi.
4. Explore hubs and subnetworks.

## If you are preparing a poster or report

1. Use `pipeline_report.txt`.
2. Use `pipeline_summary_figure.png`.
3. Use `interactive_dashboard.html` for screenshots.
4. Use `novel_candidates.csv` to highlight discoveries.

---

# Final advice

If you want the easiest way to understand and use GLOOM, start with these files:

1. `results/gene_rankings.csv`
2. `results/novel_candidates.csv`
3. `figures/interactive_dashboard.html`
4. `results/network/annotated_network.graphml`
5. `results/reports/pipeline_report.txt`

These five files usually give the fastest understanding of what the package found.
