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
