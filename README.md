# Pancancer_LACTATE

Pathology feature extraction (WSI-level)
- **`ctran.py`**  
  Defines the **CTransPath** backbone used for pathology feature extraction.  
  A convolutional stem is implemented and integrated into a Swin Transformer architecture (via `timm`) to generate patch-level embeddings from histopathology images.
CTransPath pretrained weights (ctranspath.phth) can be downloaded here: https://drive.google.com/file/d/1dhysqcv_Ct_A96qOF8i6COTK3jLb56vx/view
- **`Ctranspath_h5.ipynb`**  
  Uses the CTransPath model to extract patch-level features from WSIs and save them in HDF5 (`.h5`) format.  
  These features serve as inputs for downstream slide-level modeling.

- **`Transformer_feature_extraction.ipynb`**  
  Performs feature aggregation and preprocessing for **Transformer**, a transformer-based aggregation framework.  
  Extracted features are organized at the slide level for immune-related stratification.


The scripts in this folder enable:
- CTransPath-based feature extraction from histopathology patches  
- Transformer–ready representations at the slide level  

This module is designed to interface with the scRNA-seq– and TCGA-based analyses, facilitating cross-modal immune and pathology-informed studies.
If you use the Ctranspath and Transformer-based aggregation framework, please cite the following papers:

Wang X, Yang S, Zhang J, Wang M, Zhang J, Yang W, Huang J, Han X. Transformer-based unsupervised contrastive learning for histopathological image classification. Med Image Anal. 2022 Oct;81:102559. doi: 10.1016/j.media.2022.102559. Epub 2022 Jul 30. PMID: 35952419.

Wagner SJ, Reisenbüchler D, West NP, Niehues JM, Zhu J, Foersch S, Veldhuizen GP, Quirke P, Grabsch HI, van den Brandt PA, Hutchins GGA, Richman SD, Yuan T, Langer R, Jenniskens JCA, Offermans K, Mueller W, Gray R, Gruber SB, Greenson JK, Rennert G, Bonner JD, Schmolze D, Jonnagaddala J, Hawkins NJ, Ward RL, Morton D, Seymour M, Magill L, Nowak M, Hay J, Koelzer VH, Church DN; TransSCOT consortium; Matek C, Geppert C, Peng C, Zhi C, Ouyang X, James JA, Loughrey MB, Salto-Tellez M, Brenner H, Hoffmeister M, Truhn D, Schnabel JA, Boxberg M, Peng T, Kather JN. Transformer-based biomarker prediction from colorectal cancer histology: A large-scale multicentric study. Cancer Cell. 2023 Sep 11;41(9):1650-1661.e4. doi: 10.1016/j.ccell.2023.08.002. Epub 2023 Aug 30. PMID: 37652006; PMCID: PMC10507381.

