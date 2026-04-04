# `/results/` Directory

This directory contains all analytical outputs, statistical results, network analyses, and comprehensive reports generated throughout the project workflow.

### **Top-Level Results Files**

| File Name                               | Description                                                                                                                                      | Key Content                                                                                                                        | Usage                                                                                           |
| --------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------- |
| `differential_expression_results.csv` | **Differential Expression Analysis** - Statistical comparison of gene expression between experimental conditions (e.g., case vs. control). | Gene IDs, log2 fold change, p-values, adjusted p-values (FDR), significance flags, expression means per group.                     | Identifies significantly dysregulated genes for downstream analysis and validation.             |
| `feature_importance.csv`              | **Model Feature Importance** - Quantitative ranking of features (genes) by their contribution to model predictions.                        | Feature names, importance scores (e.g., Gini, permutation, SHAP), ranks, standard deviations (if applicable).                      | Highlights top predictive biomarkers; used for feature selection and biological interpretation. |
| `gene_rankings.csv`                   | **Comprehensive Gene Prioritization** - Aggregated ranking of genes across multiple criteria (DE, network centrality, ML importance).      | Gene IDs, DE rank, network rank, ML rank, combined score, final rank, annotation columns.                                          | Provides unified gene prioritization for candidate selection and experimental follow-up.        |
| `model_metrics.csv`                   | **Model Performance Evaluation** - Quantitative assessment of trained ML models using multiple metrics.                                    | Model names, AUROC, accuracy, precision, recall, F1-score, training time, hyperparameters.                                         | Compares model performance; informs model selection and optimization decisions.                 |
| `novel_candidates.csv`                | **Novel Candidate Biomarkers** - List of prioritized genes not previously associated with the disease/condition.                           | Candidate gene IDs, supporting evidence (DE p-value, ML importance, network metrics), novelty score, literature validation status. | Identifies new therapeutic targets or biomarkers for experimental validation.                   |

---

### **`/network/` Subdirectory**

Contains network construction, analysis, and export files for gene co-expression and protein-protein interaction networks.

#### **Core Network Files**

| File Name                          | Description                                                                                                               | Format & Content                                                                            |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| `coexpression_network.graphml`   | **Raw Co-expression Network** - Undirected graph of gene-gene correlations above threshold.                         | GraphML format: nodes=genes, edges=correlation weights, no annotations.                     |
| `coexpression_network_edges.csv` | **Edge List** - Tabular representation of network connections with correlation metrics.                             | CSV: Source, target, correlation coefficient, p-value, FDR-adjusted p-value.                |
| `annotated_network.graphml`      | **Annotated Network** - Network enriched with gene annotations, centrality scores, and community detection results. | GraphML: Nodes contain gene symbols, DE status, ML importance, betweenness centrality, etc. |
| `annotated_nodes.csv`            | **Node Attributes** - Complete annotation table for all genes in the network.                                       | CSV: All node metadata from annotated network, suitable for statistical analysis.           |
| `annotated_edges.csv`            | **Annotated Edge List** - Enhanced edge table with additional network metrics.                                      | CSV: Edge weights, interaction types, shared pathways, co-expression significance.          |

#### **`/exports/` Subdirectory**

Network files exported in multiple formats for different analysis tools and visualization platforms.

| File Name                         | Description                                                                                       | Target Tool/Use Case                                                                |
| --------------------------------- | ------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| `network_cytoscape.json`        | **Cytoscape-compatible JSON** - Network formatted for Cytoscape visualization and analysis. | Cytoscape desktop application for advanced network visualization.                   |
| `network_full.graphml`          | **Complete Network (GraphML)** - Full annotated network in standard GraphML format.         | Gephi, yEd, or any GraphML-compatible network analyzer.                             |
| `network_full.gml`              | **Complete Network (GML)** - Alternative format for network analysis tools.                 | NetworkX, igraph, or other GML-supporting libraries.                                |
| `network_full_edgelist.tsv`     | **Tab-separated Edge List** - Simple TSV format for quick loading into various tools.       | R, Python, or custom scripts for network analysis.                                  |
| `network_nodes.tsv`             | **Node Attribute Table** - Tab-separated node metadata for external analysis.               | Statistical analysis in R/Python or integration with other datasets.                |
| `network_statistics_report.csv` | **Network Statistics (CSV)** - Quantitative network metrics in machine-readable format.     | Automated reporting or integration with downstream analyses.                        |
| `network_statistics_report.txt` | **Network Statistics (Text)** - Human-readable summary of network properties.               | Quick reference: node count, edge count, density, diameter, clustering coefficient. |

#### **Subnetwork Exports**

| File Name                              | Description                                                                                | Selection Criteria                                    |
| -------------------------------------- | ------------------------------------------------------------------------------------------ | ----------------------------------------------------- |
| `subnetwork_top100.graphml`          | **Top 100 Genes Network** - Subnetwork containing only the top-ranked genes.         | Based on combined ranking from `gene_rankings.csv`. |
| `subnetwork_top100_edgelist.tsv`     | **Edge list for Top 100**                                                            | Corresponding edges for the top 100 gene subnetwork.  |
| `subnetwork_cgc.graphml`             | **Cancer Gene Census Network** - Subnetwork of genes from COSMIC Cancer Gene Census. | Genes annotated in COSMIC CGC database.               |
| `subnetwork_cgc_edgelist.tsv`        | **Edge list for CGC Network**                                                        | Corresponding edges for CGC gene subnetwork.          |
| `subnetwork_novel.graphml`           | **Novel Candidates Network** - Subnetwork of novel candidate biomarkers.             | Genes from `novel_candidates.csv`.                  |
| `subnetwork_novel_edgelist.tsv`      | **Edge list for Novel Network**                                                      | Corresponding edges for novel candidate subnetwork.   |
| `subnetwork_candidates.graphml`      | **All Candidate Network** - Subnetwork containing all candidate biomarkers.          | Union of top-ranked and novel candidates.             |
| `subnetwork_candidates_edgelist.tsv` | **Edge list for Candidates**                                                         | Corresponding edges for candidate subnetwork.         |

---

### **`/reports/` Subdirectory**

Contains comprehensive reports and summary visualizations of the complete analysis pipeline.

| File Name                       | Description                                                                                      | Key Contents                                                                                                                                                                                                                                                                              |
| ------------------------------- | ------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `pipeline_report.txt`         | **Complete Pipeline Report** - Detailed textual summary of all analysis steps and results. | 1. Executive summary`<br>`2. Data preprocessing statistics`<br>`3. DE analysis summary (# significant genes)`<br>`4. Model performance comparison`<br>`5. Top biomarker findings`<br>`6. Network statistics`<br>`7. Novel candidate list`<br>`8. Limitations and next steps |
| `pipeline_summary_table.csv`  | **Summary Statistics Table** - Aggregated metrics from all analysis stages.                | Sample counts, DE gene counts, model AUROCs, network metrics, candidate counts.                                                                                                                                                                                                           |
| `pipeline_summary_figure.png` | **Summary Figure (PNG)** - Visual overview of key findings for presentations.              | Multi-panel figure: DE volcano plot, model ROC curves, network visualization, candidate heatmap.                                                                                                                                                                                          |
| `pipeline_summary_figure.pdf` | **Summary Figure (PDF)** - High-resolution vector version for publications.                | Same content as PNG but in scalable vector format for print/publication.                                                                                                                                                                                                                  |

---

### **Key Analytical Workflows**

1. **Differential Expression → Feature Importance Integration**

   ```
   differential_expression_results.csv → gene_rankings.csv (DE rank)
   feature_importance.csv → gene_rankings.csv (ML rank)
   ```
2. **Network-Enhanced Candidate Discovery**

   ```
   coexpression_network.graphml → annotated_network.graphml → gene_rankings.csv (network rank)
   gene_rankings.csv → novel_candidates.csv (filter by literature)
   ```
3. **Multi-Format Network Analysis**

   ```
   annotated_network.graphml → /exports/ (various formats) → tool-specific analysis
   ```
4. **Comprehensive Reporting**

   ```
   All CSV results + network statistics → pipeline_report.txt + pipeline_summary_figure.png
   ```

---

### **Usage Guidelines**

- **For Quick Insights**: Start with `pipeline_summary_figure.png` and `pipeline_report.txt`
- **For Candidate Selection**: Use `gene_rankings.csv` for prioritized lists, `novel_candidates.csv` for new discoveries
- **For Network Analysis**: Use appropriate format from `/exports/` for your preferred tool
- **For Model Comparison**: Refer to `model_metrics.csv` for quantitative performance assessment
- **For Publication**: Use PDF figures and CSV tables for supplementary materials

---

### **File Relationships**

```
differential_expression_results.csv → Input for network construction and gene ranking
feature_importance.csv → Input for gene ranking and candidate selection
coexpression_network.graphml → Base for all annotated and subnetworks
gene_rankings.csv → Source for top100 and candidate subnetworks
All results → Compiled into pipeline_report.txt and summary_figure.png
```
