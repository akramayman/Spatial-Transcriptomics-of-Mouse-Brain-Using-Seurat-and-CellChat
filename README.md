<p align="left">
  <img align="left"
       src="https://uploads-ssl.webflow.com/5f02d7013319a34e214498a6/62d62e7969615a7d137014e9_Icon%20Design-08.png"
       width="120"
       alt="Logo">

# Spatial Transcriptomics of Mouse Brain Using Seurat and CellChat
</p>

<br clear="left"/>

<p align="left">
  <img src="https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white" alt="R">
  <img src="https://img.shields.io/badge/Seurat-Spatial%20Analysis-blueviolet?style=for-the-badge" alt="Seurat">
  <img src="https://img.shields.io/badge/CellChat-Cell--Cell%20Communication-orange?style=for-the-badge" alt="CellChat">
  <img src="https://img.shields.io/badge/Mouse-Brain-red?style=for-the-badge" alt="Mouse Brain">
</p>

## 🧬 Project Description

This repository contains an end-to-end analysis pipeline for **Spatial Transcriptomics** data generated using the **10x Genomics Visium** platform. The project investigates spatial gene expression patterns, tissue architecture, and cell–cell communication within mouse brain tissue sections using a combination of **Seurat**, **CellChat**, and **SCDC** in R.

The workflow provides a comprehensive framework for processing, visualizing, and interpreting spatial transcriptomics data, enabling the identification of biologically meaningful spatial domains and cellular interactions across brain regions.

---

## 🧠 Overview

Spatial transcriptomics enables transcriptome-wide profiling of gene expression while preserving the spatial context of cells within intact tissue. By integrating gene expression data with tissue morphology, researchers can uncover spatially organized biological processes and cellular interactions.

In this project, two mouse brain tissue sections are analyzed to:

- 🔍 Explore spatial gene expression patterns
- 📍 Visualize tissue-specific expression landscapes
- 🧩 Identify spatially distinct clusters
- 🔗 Integrate multiple Visium tissue sections
- 📈 Detect differentially expressed genes (DEGs)
- 🏷️ Annotate cell types using reference single-cell RNA-seq datasets
- 💬 Infer cell–cell communication networks
- 🧬 Perform cell-type deconvolution analysis

---

## 📊 Methods and Tools

### Main Packages

| Package | Purpose |
|----------|----------|
| **Seurat** | Spatial transcriptomics preprocessing, clustering, integration, and visualization |
| **CellChat** | Inference and analysis of cell–cell communication networks |
| **SCDC** | Cell-type deconvolution using reference scRNA-seq datasets |

## ⚙️ System Setup

Before running the analysis, ensure your R environment is properly configured with all required packages.

---

## 🖥️ 1. General Setup (All Users)

Open an R session and run the following commands:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install Seurat and core packages
install.packages(c("Seurat", "patchwork", "ggplot2", "dplyr", "Matrix"))

# Install Bioconductor manager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install CellChat from GitHub
devtools::install_github("sqjin/CellChat")

# Install SCDC for deconvolution analysis
devtools::install_github("meichendong/SCDC")
```

---

## 🖥️ 2. Additional Dependencies (Platform-Specific)

Some packages may require manual installation depending on your operating system.

### 🪟 Windows / Compatibility Fixes

```r
# Needed for SCTransform and variance stabilization
devtools::install_github("const-ae/glmGamPoi")

# Required for Seurat integration and nearest-neighbor search
BiocManager::install("BiocNeighbors")
```
---

## Citation 

If you use this repository, please cite:

Akram Abushmais.
