# EMGD-ML

**Feature engineering pipeline for predicting gene dependency from label-free proteomics**

---

## Overview

**EMGD-ML** is a framework for deriving proteomics-based features associated with gene dependency using label-free proteomic data.

The method combines cancer cell line (CCL) label-free proteomics with CRISPR–Cas9 gene dependency probabilities (DepMap) to identify proteins consistently associated with gene essentiality states.

This repository currently provides the **feature engineering backbone** of the EMGD-ML framework:
- Empirical Marker of Gene Dependency (EMGD) construction
- Marker identification from essential vs non-essential sample groups
- D-value computation for downstream machine learning

These components generate interpretable, robust features that capture pathway-level signals of gene dependency.

---

## Data Sources

This pipeline utilises gene dependency and label-free proteomics data from the following sources:

- **DepMap CRISPR–Cas9 gene dependency data**  
  Broad Institute DepMap Portal  
  https://depmap.org/portal/
  Dempster JM, Pacini C, Pantel S et al. Agreement between two large pan-cancer CRISPR-Cas9 gene dependency datasets. Nat Commun 2021;12:1255.

- **Cancer cell line proteomics**  
  Gonçalves E, Poulos RC, Cai Z et al. Pan-cancer proteomic map of 949 human cell lines. Cancer Cell 2022;40:835–49.

All datasets are publicly available and can be obtained seperately from source

---

## Example Data

This repository includes a reduced gene dependency training file containing the first 10 genes for demonstration and testing purposes. 
The accompanying proteomics data contains training samples required to run the pipeline. 
This example data is intended for reproducibility and code testing only, and does not represent the full dataset used in the study.

---

## Current Scope

This version of the repository includes:

### 1. EMGD Group Construction

For each gene, samples are stratified into a **essential (ESS)** or **non-essential (NES)** group using multiple strategies:
- Median 
- Quartile 
- Fixed high/low dependency thresholds
- Extremes (top/bottom fraction)
- K-means clustering

### 2. Marker Generation

Protein markers associated with ESS and NES states are identified via:
- Repeated statistical comparisons
- Directional consistency filtering
- Combined p-value ranking

### 3. D-value Computation

Marker sets are reduced into **D-values**, defined as:

\[
D = (Q2_{ESS} - Q2_{NES}) + (Q3_{ESS} - Q3_{NES})
\]

These dimensionality reduced D values represent proteomic patterns per gene and sample.

---
## Installation

```
git clone https://github.com/robert-nkwo/EMGD-ML.git
cd EMGD-ML

python -m venv venv
source venv/bin/activate   # Mac/Linux
# OR
venv\Scripts\activate      # Windows

pip install -r requirements.txt
```

---

## Input Data Format

| Data Type       | Format             |
| --------------- | ------------------ |
| Gene dependency | genes × samples    |
| Proteomics      | samples × proteins |


## Usage Example

Step 1 - Build EMGD Groups
```
python scripts/02_build_emgd_groups.py \
  --dependency data/processed/gene_dependency_train_example.csv \
  --proteomics data/processed/proteomics_train.csv \
  --output-dir results/emgd_groups
```
Step 2 - Generate Marker Database
```
python scripts/03_generate_emgd_markers.py \
  --dependency data/processed/gene_dependency_train_example.csv \
  --proteomics data/processed/proteomics_train.csv \
  --groups results/emgd_groups/emgd_median_groups.json \
  --output results/emgd_markers/emgd_median_marker_db.csv
```
Step 3 - Calculate D-values
```
python scripts/04_calculate_d_values.py \
  --markers results/emgd_markers/emgd_median_marker_db.csv \
  --proteomics data/processed/proteomics_train.csv \
  --output results/d_values/emgd_median_dvalues.csv \
  --missing-at-random
```
Run Full Pipeline:
```
bash run_emgd_pipeline.sh
```
---

## Current Repository Structure

```
EMGD-ML/
│
├── data/
│   └── processed/
│
├── scripts/
│   ├── 02_build_emgd_groups.py
│   ├── 03_generate_emgd_markers.py
│   └── 04_compute_d_values.py
│
├── results/
│   ├── emgd_groups/
│   ├── emgd_markers/
│   └── d_values/
│
├── run_emgd_pipeline.sh
├── requirements.txt
└── README.md
```

## Current Key Features

- Uses label-free proteomics only
- Captures gene dependency signals without genomic input
- Generates interpretable feature representations (D-values)
- Robust to missing protein measurements
- Modular and reproducible pipeline

---

## Roadmap

Future updates will include:

- Machine learning model training (RF, XGBoost, NN, SVR)
- Model evaluation and benchmarking
- Application to clinical HCC proteomics datasets
- Visualisation tools/scripts

---

## Citation


A formal citation will be provided at a later date

---

## License

MIT License

---
## Contact
rtnkwo@gmail.com