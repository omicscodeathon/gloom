## `/docs/` Directory

This directory contains supplementary documentation, reference materials, and external resources that support the project's scientific context, methodology, and dissemination.

### **Documentation Files**

| File Name          | Description                                                                                                                                                             | Format & Contents                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | Primary Use                                                                                    |
| ------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| `notes.txt`      | **Research Notes & Methodology Documentation** - Detailed technical notes, parameter settings, methodological decisions, and analytical rationale.                | Plain text with structured sections:`<br>`• **Experimental Design**: Sample sizes, conditions, replicates `<br>`• **Preprocessing Parameters**: Normalization methods, filtering thresholds, batch correction `<br>`• **Statistical Settings**: DE analysis parameters (FDR cutoff, fold change threshold)`<br>`• **Model Specifications**: ML algorithms used, hyperparameter ranges, validation strategy `<br>`• **Network Parameters**: Correlation thresholds, community detection algorithms `<br>`• **Interpretation Notes**: Biological context, assumptions, limitations | Reference for methodology replication, parameter justification, and technical troubleshooting. |
| `papers_url.txt` | **Reference Literature Database** - Curated collection of relevant scientific publications supporting the project's biological context and analytical approaches. | Text file with structured entries:`<br>`• **Format**: `[PMID/DOI]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | [Title]                                                                                        |
| `slides_url.txt` | **Presentation & Dissemination Resources** - Links to project presentations, conference materials, and educational resources.                                     | Text file with categorized entries:`<br>`• **Format**: `[Title]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | [Presenter/Author]                                                                             |

---

### **Detailed Contents**

#### **`notes.txt` - Technical Documentation**

```
# PROJECT NOTES - [Project Name]
# Last Updated: [Date]
# Maintainer: [Name/Team]

## 1. EXPERIMENTAL DESIGN
- Sample Collection: [Description]
- Conditions: [Case/Control definitions]
- Replicates: [Number per condition]
- Quality Control: [QC metrics and thresholds]

## 2. DATA PREPROCESSING
- Normalization: [Method: TPM/FPKM/RPKM with justification]
- Filtering: [Low-expression cutoff: >1 CPM in >50% samples]
- Batch Correction: [ComBat/limma with parameters]
- Missing Value Imputation: [kNN/mean with k=10]

## 3. DIFFERENTIAL EXPRESSION ANALYSIS
- Tool: [limma/edgeR/DESeq2]
- Parameters: [FDR cutoff: 0.05, logFC threshold: 1.0]
- Contrasts: [Specific comparisons made]
- Multiple Testing Correction: [Benjamini-Hochberg]

## 4. MACHINE LEARNING PIPELINE
- Feature Selection: [Top 5000 most variable genes]
- Algorithms: [Random Forest, SVM, XGBoost with hyperparameters]
- Validation: [5-fold CV repeated 10 times]
- Performance Metrics: [AUROC, accuracy, precision, recall]

## 5. NETWORK ANALYSIS
- Construction: [WGCNA/co-expression with r > 0.7]
- Annotation Sources: [GO, KEGG, STRING]
- Centrality Metrics: [Betweenness, degree, closeness]
- Community Detection: [Louvain method with resolution=1.0]

## 6. CANDIDATE PRIORITIZATION
- Ranking Criteria: [DE p-value (40%), ML importance (30%), network centrality (30%)]
- Novelty Assessment: [PubMed search, DisGeNET, OMIM]
- Validation Strategy: [Experimental plan]
```

#### **`papers_url.txt` - Literature Database**

```
# REFERENCE LITERATURE
# Format: ID | Title | Authors | Source | URL | Tags

## DISEASE BIOLOGY
31801896 | Molecular pathogenesis of Disease X | Chen et al. | Cell, 2020 | https://doi.org/xxx | [pathogenesis][mechanism]
24523456 | Genetic susceptibility factors | Johnson et al. | Nature Genetics, 2018 | https://doi.org/xxx | [genetics][risk_factors]

## METHODOLOGY
27899658 | limma powers differential expression analyses | Ritchie et al. | Nucleic Acids Research, 2015 | https://doi.org/xxx | [DE_method][R_package]
12345678 | Network-based biomarker discovery | Zhang et al. | Bioinformatics, 2021 | https://doi.org/xxx | [network_analysis][biomarkers]

## REVIEW ARTICLES
34567890 | Comprehensive review of ML in genomics | Wang et al. | Nature Reviews Genetics, 2022 | https://doi.org/xxx | [review][ML_genomics]
```

#### **`slides_url.txt` - Presentation Resources**

```
[1] H. Shadman, S. Gomrok, C. Litle, Q. Cheng, Y. Jiang, X. Huang, J. D. Ziebarth, and Y. Wang, “A machine learning-based investigation of integrin expression patterns in cancer and metastasis,” Scientific Reports, vol. 15, art. no. 5270, 2025.
[2] E. Glaab, J. Bacardit, J. M. Garibaldi, and N. Krasnogor, “Using rule-based machine learning for candidate disease gene prioritization and sample classification of cancer gene expression data,” PLoS ONE, vol. 7, no. 7, e39932, 2012.
[3] J. Figueroa-Martínez, D. M. Saz-Navarro, A. López-Fernández, D. S. Rodríguez-Baena, and F. A. Gómez-Vela, “Computational ensemble gene co-expression networks for the analysis of cancer biomarkers,” 2024.
[4] M. Khalsan, L. R. Machado, E. S. Al-Shamery, S. Ajit, K. Anthony, M. U. Mu, and M. O. Agyeman, “A survey of machine learning approaches applied to gene expression analysis for cancer prediction,” IEEE Access, vol. 10, pp. 27522–27534, 2022, doi: 10.1109/ACCESS.2022.3146312.
[5] G. Kallah-Dagadu, M. Mohammed, J. B. Nasejje, N. N. Mchunu, H. S. Twabi, J. M. Batidzirai, G. C. Singini, P. Nevhungoni, and I. Maposa, “Breast cancer prediction based on gene expression data using interpretable machine learning techniques,” Scientific Reports, vol. 15, art. no. 7594, 2025.

```

---

### **Usage Guidelines**

#### **For New Project Members**

1. **Start with `notes.txt`** to understand methodology and parameters
2. **Review `papers_url.txt`** for biological context and methodological references
3. **Check `slides_url.txt`** for previous presentations and training materials

#### **For Manuscript Preparation**

- **Methodology Section**: Reference relevant entries in `papers_url.txt`
- **Technical Details**: Cite specific parameters from `notes.txt`
- **Supplementary Materials**: Link to presentations in `slides_url.txt` for dissemination history

#### **For Analysis Replication**

- **Exact Parameters**: Use values documented in `notes.txt`
- **Tool Citations**: Reference papers from `papers_url.txt`
- **Validation**: Compare with methods from referenced presentations

#### **For Project Continuity**

- **Update `notes.txt`** when changing analysis parameters
- **Add to `papers_url.txt`** when discovering new relevant literature
- **Record in `slides_url.txt`** after each presentation or dissemination activity

---

### **Maintenance & Updates**

| File               | Update Frequency                 | Responsibility   | Version Control                    |
| ------------------ | -------------------------------- | ---------------- | ---------------------------------- |
| `notes.txt`      | After each major analysis change | Lead Analyst     | Git commit with change description |
| `papers_url.txt` | Monthly literature review        | All Team Members | Add new entries with date stamps   |
| `slides_url.txt` | After each presentation          | Presenter        | Add within 1 week of presentation  |

**Version Convention**: Append date to filenames for major updates (e.g., `notes_2026-04.txt`)
