## **Data Organization**

This project follows a structured data workflow that separates **raw input files**, **processed data**, **intermediate feature tables**, and **final results**. The goal is to keep the pipeline reproducible, easy to inspect, and easy to debug. All key file locations are defined in `config.py`, and the analysis is executed stepwise from raw transcriptomic input to final gene ranking and network visualization outputs.

### **1. Raw Data**

Raw data are the original external input files used by the pipeline. These files should not be modified manually. They are stored under `data/raw/` and include tumor expression, tumor metadata, normal lung expression, normal metadata, and the Cancer Gene Census labels.

**Expected raw data structure**

```text
data/raw/
├── cBioPortal (RNA Seq Data)/
│   ├── data_mrna_seq_v2_rsem.txt
│   └── data_clinical_patient.txt
├── Gtex (normal samples)/
│   ├── gene_tpm_v11_lung.gct.gz
│   └── GTEx_Analysis_v11_Annotations_SampleAttributesDD.xlsx
└── Cancer Gene Census (Labeled Data)/
    └── Census_allWed Mar 18 05_48_31 2026.csv
```

This layout matches both the project README and the paths declared in `config.py`.

**Raw data description**

* `data_mrna_seq_v2_rsem.txt`: LUAD tumor RNA-seq expression matrix from cBioPortal.
* `data_clinical_patient.txt`: clinical metadata for tumor samples.
* `gene_tpm_v11_lung.gct.gz`: GTEx lung normal expression matrix.
* `GTEx_Analysis_v11_Annotations_SampleAttributesDD.xlsx`: GTEx sample annotation metadata.
* `Census_allWed Mar 18 05_48_31 2026.csv`: Cancer Gene Census file used to build binary labels for supervised learning.

---

### **2. Processed Data**

Processed data are cleaned and standardized tables generated after quality control, normalization, and harmonization. These files are stored in `data/processed/` and are used as the main structured inputs for downstream feature engineering and machine learning.

**Processed data files**

```text
data/processed/
├── tumor_expression_processed.csv
├── normal_expression_processed.csv
├── tumor_expression_harmonized.csv
├── normal_expression_harmonized.csv
├── expression_features.csv
├── network_features.csv
├── integrated_features.csv
├── gene_labels.csv
├── train_features.csv
├── val_features.csv
├── train_labels.csv
└── val_labels.csv
```

These filenames correspond directly to the processed and intermediate path definitions in `config.py`.

**Processed data description**

* `tumor_expression_processed.csv`: tumor expression matrix after QC and preprocessing.
* `normal_expression_processed.csv`: normal expression matrix after QC and preprocessing.
* `tumor_expression_harmonized.csv`: tumor matrix restricted to the harmonized shared gene set.
* `normal_expression_harmonized.csv`: normal matrix restricted to the harmonized shared gene set.
* `expression_features.csv`: per-gene expression-derived features such as differential expression statistics.
* `network_features.csv`: per-gene network topology features extracted from the co-expression network.
* `integrated_features.csv`: merged feature matrix combining expression and network-derived features.
* `gene_labels.csv`: binary target labels constructed from the Cancer Gene Census.
* `train_features.csv`, `val_features.csv`, `train_labels.csv`, `val_labels.csv`: stratified machine learning train/validation split outputs.

---

### **3. Intermediate Data**

Intermediate data are files produced between major pipeline stages. They are not necessarily final outputs, but they are essential for traceability and debugging. In this project, intermediate artifacts include differential expression results, feature tables, train/validation split tables, and co-expression network files generated before final annotation and export.

**Examples of intermediate files**

* `results/differential_expression_results.csv`: output of the differential expression step.
* `results/network/coexpression_network_edges.csv`: edge list of the initial co-expression network.
* `results/network/coexpression_network.graphml`: graph representation of the co-expression network before final annotation.
* `results/feature_importance.csv`: model interpretation output generated after training and evaluation.
* `results/model_metrics.csv`: model-level evaluation metrics.

**Why keep intermediate files?**
Intermediate outputs make the pipeline more reproducible and easier to inspect. They allow users to validate preprocessing, troubleshoot feature integration, and re-run only selected steps without restarting the full workflow.

---

### **4. Final Results**

Final results are the main deliverables of the pipeline. These outputs summarize the machine learning prioritization, annotated network, figures, and reporting products generated in the final analysis stages. According to the pipeline README, the key deliverables include gene rankings, novel candidate genes, annotated networks, interactive visualization, and human-readable reports.

**Final result directories**

```text
results/
├── gene_rankings.csv
├── feature_importance.csv
├── model_metrics.csv
├── network/
│   ├── annotated_network.graphml
│   ├── annotated_edges.csv
│   └── annotated_nodes.csv
└── reports/
    └── pipeline_summary_report.csv

figures/
└── interactive_network.html
```

The file paths above follow the outputs declared in `config.py`.

**Final result description**

* `results/gene_rankings.csv`: all genes ranked by model-predicted probability.
* `results/feature_importance.csv`: feature-importance scores used for interpretation of model behavior.
* `results/model_metrics.csv`: evaluation metrics for the trained models.
* `results/network/annotated_network.graphml`: annotated network ready for graph tools such as Cytoscape or Gephi.
* `results/network/annotated_edges.csv` and `results/network/annotated_nodes.csv`: tabular exports of the annotated network.
* `figures/interactive_network.html`: browser-based interactive visualization of the network and ML results.
* `results/reports/pipeline_summary_report.csv`: structured summary report of the pipeline.

---

### **5. Models, Figures, and Logs**

The pipeline also organizes supporting outputs into dedicated folders for trained models, generated figures, and execution logs. These are defined in `config.py` as `models/`, `figures/`, and `logs/`.

**Support directories**

* `models/`: stores trained model objects and related serialized artifacts.
* `figures/`: stores visual outputs including the interactive HTML network view.
* `logs/`: stores pipeline execution logs such as `pipeline.log`.
