# MGP_HNSC32

**A H&E Morphology-Derived 32-Gene Prognostic Signature in Head and Neck Squamous Cell Carcinoma**

This repository provides code for pathology feature extraction, gene expression prediction, and the development and validation of MGP_HNSC32, a morphology-derived prognostic signature for head and neck squamous cell carcinoma (HNSCC).

The original pathology feature extraction and gene expression prediction modules are retained. Additional R scripts for MGP_HNSC32 are organized in the `R/` directory.

## Repository Structure

```text
MGP_HNSC32/
├── ctran.py
├── Ctranspath_h5.ipynb
├── Transformer_feature_extraction.ipynb
├── scoring.r
├── features.csv
├── R/
│   ├── 01_Morphology_Associated_Gene_Identification.r
│   ├── 02_RSF_Gene_Ranking.r
│   ├── 03_Cox_Model_Construction.r
│   ├── 04_Risk_Score_and_External_Validation.r
│   ├── 05_Model_Stability.r
│   └── 06_Propensity_Score_Matching.r
└── README.md
```

**Note:** `rf_models.rds` is described below but is not included in the repository.

# 1. Pathology Feature Extraction (WSI-level)

### `ctran.py`

Defines the **CTransPath** backbone used for pathology feature extraction.

A convolutional stem is implemented and integrated into a Swin Transformer architecture (via `timm`) to generate patch-level embeddings from histopathology images.

CTransPath pretrained weights (`ctranspath.phth`) can be downloaded here:

https://drive.google.com/file/d/1dhysqcv_Ct_A96qOF8i6COTK3jLb56vx/view

### `Ctranspath_h5.ipynb`

Uses the CTransPath model to extract patch-level features from WSIs and save them in HDF5 (`.h5`) format.

These features serve as inputs for downstream slide-level modeling.

### `Transformer_feature_extraction.ipynb`

Performs feature aggregation and preprocessing for **Transformer**, a transformer-based aggregation framework.

Extracted features are organized at the slide level for immune-related stratification.

### Module Overview

The scripts in this module enable:

- CTransPath-based feature extraction from histopathology patches.
- Transformer-ready representations at the slide level.

This module is designed to interface with scRNA-seq- and TCGA-based analyses, facilitating cross-modal immune- and pathology-informed studies.

### References

If you use the CTransPath and Transformer-based aggregation framework, please cite the following papers:

Wang X, Yang S, Zhang J, Wang M, Zhang J, Yang W, Huang J, Han X. Transformer-based unsupervised contrastive learning for histopathological image classification. *Medical Image Analysis*. 2022 Oct;81:102559. doi: [10.1016/j.media.2022.102559](https://doi.org/10.1016/j.media.2022.102559). Epub 2022 Jul 30. PMID: 35952419.

Wagner SJ, Reisenbüchler D, West NP, Niehues JM, Zhu J, Foersch S, Veldhuizen GP, Quirke P, Grabsch HI, van den Brandt PA, Hutchins GGA, Richman SD, Yuan T, Langer R, Jenniskens JCA, Offermans K, Mueller W, Gray R, Gruber SB, Greenson JK, Rennert G, Bonner JD, Schmolze D, Jonnagaddala J, Hawkins NJ, Ward RL, Morton D, Seymour M, Magill L, Nowak M, Hay J, Koelzer VH, Church DN; TransSCOT consortium; Matek C, Geppert C, Peng C, Zhi C, Ouyang X, James JA, Loughrey MB, Salto-Tellez M, Brenner H, Hoffmeister M, Truhn D, Schnabel JA, Boxberg M, Peng T, Kather JN. Transformer-based biomarker prediction from colorectal cancer histology: A large-scale multicentric study. *Cancer Cell*. 2023 Sep 11;41(9):1650-1661.e4. doi: [10.1016/j.ccell.2023.08.002](https://doi.org/10.1016/j.ccell.2023.08.002). Epub 2023 Aug 30. PMID: 37652006; PMCID: PMC10507381.

# 2. Gene Expression Prediction

### `scoring.r`

This script implements the gene scoring procedure.

It takes pathology-derived features as input and applies a set of pretrained models to estimate gene expression levels.

Specifically, it generates predicted expression values for **671 risk-associated genes** based on the input features.

### `features.csv`

This file contains the input feature matrix.

Features are extracted from H&E whole-slide images (WSIs) using a combination of CTransPath and Transformer-based models, resulting in high-dimensional pathology representations for each sample.

### `rf_models.rds`

This file stores the pretrained model parameters used for gene expression prediction.

Due to size and sharing constraints, it is not included in this repository. Researchers with reasonable requests may contact the authors to obtain access.

# 3. MGP_HNSC32 Statistical Analysis

The `R/` directory contains six scripts covering morphology-associated gene identification, prognostic model development, external validation, model stability and propensity score matching.

| Script | Description |
|---|---|
| `01_Morphology_Associated_Gene_Identification.r` | Identification of morphology-associated genes using pathology features and gene expression prediction. |
| `02_RSF_Gene_Ranking.r` | Random survival forest ranking of candidate genes and selection of the top 100 genes. |
| `03_Cox_Model_Construction.r` | Construction of the prognostic signature using stepwise multivariable Cox regression. |
| `04_Risk_Score_and_External_Validation.r` | Risk score calculation, TCGA-derived cutoff determination and external validation in HANCOCK. |
| `05_Model_Stability.r` | Bootstrap assessment of gene selection stability. |
| `06_Propensity_Score_Matching.r` | Propensity score matching and post-matching survival analysis. |

## 3.1. Morphology-Associated Gene Identification

Histopathological features derived from H&E WSIs are used to predict gene expression.

Candidate genes are screened using Spearman correlation. Gene-specific random forest models are subsequently evaluated in an independent validation subset.

Morphology-associated genes are identified using the following validation criteria:

- Spearman correlation coefficient ≥ 0.30.
- Benjamini–Hochberg-adjusted FDR ≤ 0.05.

## 3.2. Random Survival Forest Gene Ranking

Candidate morphology-associated genes are ranked according to their prognostic importance using random survival forest.

The top 100 genes are retained for subsequent Cox regression analysis.

## 3.3. Cox Model Construction

Stepwise multivariable Cox proportional hazards regression is used to construct the MGP_HNSC32 prognostic signature.

The final model contains 32 genes.

## 3.4. Risk Score Calculation and External Validation

Risk scores are calculated using the fitted Cox model.

The optimal cutoff is determined in the TCGA-HNSC cohort using the maximally selected rank statistic implemented in the `maxstat` package:

```r
smethod = "LogRank"
minprop = 0.2
maxprop = 0.8
pmethod = "exactGauss"
```

The TCGA-derived cutoff is applied unchanged to the independent HANCOCK cohort.

## 3.5. Model Stability

Gene selection stability is evaluated using 1,000 bootstrap iterations.

The selection frequency of each candidate gene is calculated across successfully completed iterations.

## 3.6. Propensity Score Matching

Propensity score matching is performed using the `MatchIt` package.

The matching variables are T stage, N stage, resection margin status and radiotherapy.

Matching parameters include:

- 1:1 nearest-neighbor matching.
- Logistic regression-based propensity score estimation.
- A caliper of 0.2 standard deviations.
- Matching without replacement.

Overall survival (OS), progression-free interval (PFI) and disease-specific survival (DSS) are evaluated after matching.

# 4. Data Availability

TCGA-HNSC data are available through the [NCI Genomic Data Commons](https://portal.gdc.cancer.gov/).

The HANCOCK cohort is used for independent external validation. Data access is subject to the original dataset's availability and usage conditions.

# 5. Reproducibility

The R scripts are numbered according to the analysis workflow.

To reproduce the full analysis, the corresponding pathology features, gene expression data, survival information, trained prediction models, Cox coefficients and TCGA-derived cutoff are required.

Input paths should be configured before running the scripts.

# 6. Citation

For the original pathology feature extraction and Transformer-based aggregation modules, please cite the corresponding publications listed in Section 1.
