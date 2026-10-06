# GLOOM

## Gene Network some are truly important for disease biology.

## Gene Network Learning and Organization through Optimized Machine Intelligence

The idea behind **GLOOM** is simple:

- take **tumor and normal gene expression data**
- compare them
- learn patterns from **known cancer genes**
- rank all genes by how likely they are to matter
- place the results inside a **gene network**
- make everything easy to inspect through tables, figures, and interactive dashboards

In short:

**raw data → cleaned data → features → machine learning → ranked genes → annotated network → report**

---

# 2. Why GLOOM Was Built

Many existing workflows are fragmented:

- one tool for preprocessing
- another for differential expression
- another for networks
- another for machine learning
- another for visualization

This makes analysis harder to reproduce and harder to explain.

**GLOOM was built to unify these steps in one reproducible Python workflow.**

---

# 3. Biological Problem

The current GLOOM case study focuses on **LUAD**
(**Lung Adenocarcinoma**), a major form of lung cancer.

Main challenge:

> Among thousands of genes, which ones are most likely to be biologically important?

Traditional differential expression alone may miss genes that are not the most strongly changed, but still sit in important parts of the network.

So GLOOM combines:

- **expression evidence**
- **network evidence**
- **machine learning prioritization**

---

# 4. What Makes GLOOM Different

GLOOM is not just a model.

It is a **full analysis package** that includes:

- raw data loading
- preprocessing and quality control
- gene harmonization
- differential expression
- feature construction
- co-expression network building
- machine learning training
- model evaluation
- gene ranking
- network annotation
- graph export
- interactive visualization
- final reporting

That means GLOOM is both:

- a **scientific workflow**
- and a **usable software package**

---

# 5. Main Inputs

GLOOM starts from three main information sources:

### Tumor data

- LUAD tumor RNA-seq expression matrix
- tumor metadata

### Normal data

- GTEx lung normal expression
- normal metadata

### Known labels

- Cancer Gene Census gene list

These inputs let the package compare tumor vs normal and learn from known cancer genes.

---

# 6. Core Workflow

## GLOOM Pipeline

```
graph TD
1. Load raw tumor and normal datasets
2. Clean and preprocess the data
3. Keep genes that can be compared in both datasets
4. Run differential expression analysis
5. Build expression features
6. Build co-expression network
7. Extract network features
8. Merge all features
9. Construct labels from known cancer genes
10. Split into training and validation sets
11. Train machine learning models
12. Evaluate model performance
13. Compute feature importance
14. Rank all genes
15. Annotate the network
16. Export the network in reusable formats
17. Create interactive visualizations
18. Write final report and summary outputs
```

1. Load raw tumor and normal datasets
2. Clean and preprocess the data
3. Keep genes that can be compared in both datasets
4. Run differential expression analysis
5. Build expression features
6. Build co-expression network
7. Extract network features
8. Merge all features
9. Construct labels from known cancer genes
10. Split into training and validation sets
11. Train machine learning models
12. Evaluate model performance
13. Compute feature importance
14. Rank all genes
15. Annotate the network
16. Export the network in reusable formats
17. Create interactive visualizations
18. Write final report and summary outputs

---

# 7. What the Model Learns From

GLOOM uses two big families of clues:

## A. Expression-based clues

Examples:

- tumor mean
- normal mean
- tumor variability
- log2 fold change
- adjusted p-value
- statistical contrast features

## B. Network-based clues

Examples:

- degree
- weighted degree
- betweenness centrality
- clustering coefficient
- component size
- edge-weight summaries

This lets the package see both:

- how a gene behaves
- and where it sits in the network

---

# 8. Machine Learning Strategy

GLOOM trains several models and compares them.

In the reference run, the trained models included:

- Random Forest
- Gradient Boosting
- Logistic Regression
- SVM
- Extra Trees

The workflow then selects the best model and uses it to score all genes genome-wide.

---

# 9. Key Outputs

GLOOM produces outputs for different kinds of users.

## Main result tables

- `gene_rankings.csv`
- `novel_candidates.csv`
- `feature_importance.csv`
- `model_metrics.csv`

## Network outputs

- `annotated_network.graphml`
- `annotated_nodes.csv`
- `annotated_edges.csv`

## Visual outputs

- static figures
- interactive HTML visualizations
- combined dashboard

## Final reports

- summary table
- text report
- summary figure
- master runner script

---

# 10. What the Gene Ranking Means

The ranking step gives every gene a model score.

Each gene gets:

- a **predicted probability**
- a **rank**
- a **predicted label**
- optional biological annotations
- a flag showing whether it is a **novel candidate**

This helps answer:

> Which genes should we study first?

It also helps separate:

- genes already known in cancer biology
- genes that may represent new discoveries

---

# 11. What the Network Part Adds

The network is important because genes do not act alone.

GLOOM builds a co-expression network and then annotates each node with:

- model score
- rank
- known cancer status
- novel candidate status
- differential expression information
- visualization attributes like color and size

This means the final network is not just a graph.

It becomes a **biological map of ranked evidence**.

---

# 12. Interactivity and Usability

GLOOM also produces interactive browser-based outputs:

- interactive volcano plot
- interactive gene ranking plot
- interactive co-expression network
- interactive feature importance view
- combined dashboard

So the package is not only for programmers.

It can also be used by:

- collaborators
- supervisors
- biologists
- students
- presentation audiences

---

# 13. What Happened in the Reference Run

In the successful reference run:

- all **18 steps completed**
- total runtime was about **480.2 seconds (~8 minutes)**
- harmonized gene set contained **10,986 genes**
- integrated feature table had **42 features**
- final network had **10,986 nodes** and **36,943 edges**
- the ranking step found **3,754 novel candidates**

This shows that GLOOM is not only an idea — it is already a working end-to-end pipeline.

---

# 14. Why This Matters

GLOOM matters because it helps turn raw omics data into usable scientific insight.

It helps:

- prioritize candidate genes
- interpret machine learning decisions
- connect expression changes with network structure
- create reusable visual and report outputs
- improve reproducibility

So the package is useful not only for analysis, but also for:

- communication
- collaboration
- teaching
- and downstream biological discovery

---

# 15. Who GLOOM Is For

GLOOM is especially useful for:

- **bioinformaticians** who want one workflow instead of many disconnected scripts
- **computational biologists** who want both ranking and network context
- **cancer researchers** who want candidate genes
- **students** learning how an omics ML pipeline works
- **wet-lab collaborators** who mainly need the results
- **presentation/report users** who want ready-made figures and dashboards

---

# 16. One-Sentence Summary

**GLOOM is a reproducible Python workflow that combines transcriptomics, machine learning, and network biology to rank LUAD genes and make the results easy to explore, visualize, and report.**

---

# 17. Thank You

## GLOOM

### From raw expression data to ranked genes, annotated networks, and interactive biological insight

**Next step:**Refine this draft into a polished presentation with:

- stronger visuals
- fewer words
- slide-specific titles
- speaker notes
- and a poster/pitch version

### An interpretable machine learning workflow for LUAD gene prioritization and co-expression network visualization

**Presenter:**
Rahma Yasser and collaborators

**Project theme:**
Transcriptomics • Machine Learning • Network Biology • Reproducible Bioinformatics

*GLOOM integrates raw expression data, feature engineering, model training, gene ranking, network annotation, and reporting into one workflow.*
