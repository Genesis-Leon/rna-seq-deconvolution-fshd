# RNA-seq deconvolution analysis in FSHD muscle samples
RNA-seq deconvolution pipeline using MuSiC and CIBERSORT

## Objective

This project aims to model skeletal muscle differentiation and to investigate the contribution of different cell populations in muscle impairment in Facioscapulohumeral Muscular Dystrophy (FSHD type 1 and 2).

The goal is to compare cell-type deconvolution results using two computational approaches:
- MuSiC2
- CIBERSORT

Different single-cell reference datasets were tested to evaluate their impact on deconvolution results.

---

## Methods

Two deconvolution methods were applied:

### 1. MuSiC2
Implemented in R using multiple single-cell reference datasets.

### 2. CIBERSORT
Signature matrices were generated using the same single-cell reference datasets as in the MuSiC2 analysis.
The deconvolution step was performed using the official CIBERSORT web interface (https://cibersortx.stanford.edu/)

---

## Reference datasets used

Three independent single-cell RNA-seq references were used:

### Human skeletal muscle tissue atlas
De Micheli AJ et al., 2020  
*A reference single-cell transcriptomic atlas of human skeletal muscle tissue*  
Skeletal Muscle.  
doi: 10.1186/s13395-020-00236-3


### Human skeletal muscle development atlas
Filippis et al., 2020  
*A Human Skeletal Muscle Atlas Identifies the Trajectories of Stem and Progenitor Cells across Development*  
Cell Stem Cell.  
doi: 10.1016/j.stem.2020.04.017


### Neuromuscular organoid model
Martins JMF et al., 2020  
*Self-Organizing 3D Human Trunk Neuromuscular Organoids*  
Cell Stem Cell.  
doi: 10.1016/j.stem.2019.12.007

---

## Project structure

RNAseq_deconvolution_project/  
│  
├── music2/  
│ ├── Music_script_GSE147457_ADT.R  
│ ├── Music_script_GSE128357.R  
│ ├── Music_script_GSE143704.R  
│  
├── CIBERSORTx/  
│ ├── cibersort_script_signature_matrix_GSE147457_ADT.R  
│ ├── cibersort_script_signature_matrix_GSE128357.R  
│ ├── cibersort_script_signature_matrix_GSE143704.R  
│  
├── docs/  
│ ├── references.md  
│ ├── pipeline_schema.png  
│  
└── README.md  
